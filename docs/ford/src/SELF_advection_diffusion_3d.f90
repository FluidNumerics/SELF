module self_advection_diffusion_3d

    use self_model
    use self_model3d
    use self_mesh
    
    implicit none
    
      type, extends(model3d) :: advection_diffusion_3d
        real(prec) :: nu ! diffusion coefficient
        real(prec) :: u  ! constant x-component of velocity
        real(prec) :: v  ! constant y-component of velocity
        real(prec) :: w  ! constant z-component of velocity
    
        contains
        procedure :: pretendency => pretendency_advection_diffusion_3d
        procedure :: setboundarycondition => setboundarycondition_advection_diffusion_3d
        procedure :: riemannsolver => riemannsolver_advection_diffusion_3d
        procedure :: fluxmethod => fluxmethod_advection_diffusion_3d
    
      end type advection_diffusion_3d
    
    ! Remember, the order of operations for the tendency calculation is
    !
    !    solution % BoundaryInterp(this % gpuAccel)
    !    solution % SideExchange(this % mesh,this % decomp,this % gpuAccel)
    !    PreTendency()
    !    SetBoundaryCondition()
    !    SourceMethod()
    !    RiemannSolver()
    !    FluxMethod()
    !    CalculateFluxDivergence()
    
    contains
    
      subroutine pretendency_advection_diffusion_3d(this)
      ! Here, we use the pre-tendency method to calculate the
      ! derivative of the solution using a bassi-rebay method
      ! We then do a boundary interpolation and side exchange
      ! on the gradient field
      implicit none
      class(advection_diffusion_3d), intent(inout) :: this
      ! local
      integer :: i, j, ivar, iEl, iSide, e2    
    
        do iEl = 1,this % solution % nElem ! Loop over all elements
            do iSide = 1,6 ! Loop over all sides

                !bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
                e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID

                if (e2 == 0) then
                    do ivar = 1, this % solution % nvar
                        do j = 0,this % solution % interp % N ! Loop over quadrature point
                            do i = 0,this % solution % interp % N ! Loop over quadrature points
                                this % solution % extBoundary % hostData(i,j,ivar,iSide,iEl) = 0.0_prec
                            enddo
                        enddo
                    enddo
                  
                end if
    
            end do
        end do
    
        ! calculate the averages of the solutions on the element
        ! boundaries and store is this % solution % avgBoundary
        call this % solution % BassiRebaySides(this % gpuaccel)
    
        ! calculate the derivative using the bassi-rebay form
        call this % solution % Gradient(this % geometry, &
                this % solutionGradient, selfWeakBRForm, &
                this % gpuaccel)
    
        ! interpolate the solutiongradient to the element boundaries
        call this % solutionGradient % BoundaryInterp(this % gpuaccel)
      
        ! perform the side exchange to populate the 
        ! solutionGradient % extBoundary attribute
        call this % solutionGradient % SideExchange(this % mesh, &
               this % decomp, this % gpuaccel) 
      
      end subroutine pretendency_advection_diffusion_3d
    
      subroutine setboundarycondition_advection_diffusion_3d(this)
      ! Here, we set the boundary conditions for the 
      ! solution and the solution gradient at the left
      ! and right most boundaries.
      ! 
      ! Here, we use periodic boundary conditions
      implicit none
      class(advection_diffusion_3d), intent(inout) :: this
      ! local
      integer :: i, j, ivar, iEl, iSide, e2 

        do iEl = 1,this % solution % nElem ! Loop over all elements
            do iSide = 1,6 ! Loop over all sides

                !bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
                e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID

                if (e2 == 0) then

                    do ivar = 1, this % solution % nvar
                        do j = 0,this % solution % interp % N ! Loop over quadrature points
                            do i = 0,this % solution % interp % N ! Loop over quadrature points
                                this % solutionGradient % extBoundary % hostData(1:3,i,j,ivar,iSide,iEl) = 0.0_prec
                            enddo
                        enddo
                    enddo
                  
                end if
    
            end do
        end do

        call this % solutionGradient % BassiRebaySides(this % gpuaccel)

      end subroutine setboundarycondition_advection_diffusion_3d
    
      subroutine fluxmethod_advection_diffusion_3d(this)
      implicit none
      class(advection_diffusion_3d), intent(inout) :: this
      ! Local
      integer :: iel
      integer :: ivar
      integer :: i
      integer :: j
      integer :: k 
      real(prec) :: u, v, w, nu, f, dfdx, dfdy, dfdz
    
        u = this % u
        v = this % v
        nu = this % nu
        do iel = 1, this % mesh % nelem
          do ivar = 1, this % solution % nvar
            do k = 0, this % solution % interp % N
                do j = 0, this % solution % interp % N
                    do i = 0, this % solution % interp % N
            
                    f = this % solution % interior % hostdata(i,j,k,ivar,iel)
                    dfdx = this % solutionGradient % interior % hostdata(1,i,j,k,ivar,iel)
                    dfdy = this % solutionGradient % interior % hostdata(2,i,j,k,ivar,iel)
                    dfdz = this % solutionGradient % interior % hostdata(3,i,j,k,ivar,iel)
            
                    this % flux % interior % hostdata(1,i,j,k,ivar,iel) = u*f - nu*dfdx  ! advective flux + diffusive flux (x-component)
                    this % flux % interior % hostdata(2,i,j,k,ivar,iel) = v*f - nu*dfdy  ! advective flux + diffusive flux (y-component)
                    this % flux % interior % hostdata(3,i,j,k,ivar,iel) = w*f - nu*dfdz  ! advective flux + diffusive flux (z-component)
            
                    enddo
                enddo
            enddo
          enddo
        enddo
    
      end subroutine fluxmethod_advection_diffusion_3d
    
      subroutine riemannsolver_advection_diffusion_3d(this)
      ! this method uses an linear upwind solver for the
      ! advective flux and the bassi-rebay method for the 
      ! diffusive fluxes
      implicit none
      class(advection_diffusion_3d), intent(inout) :: this
      ! Local
      integer :: iel
      integer :: ivar
      integer :: iside
      integer :: i, j
      real(prec) :: fin, fout, dfdn, un
      real(prec) :: nhat(1:3), nmag
        
      do iEl = 1,this % solution % nElem
        do iSide = 1,6
          do ivar = 1, this % solution % nvar
            do j = 0,this % solution % interp % N
                do i = 0,this % solution % interp % N

                ! Get the boundary normals on cell edges from the mesh geometry
                nhat(1:3) = this % geometry % nHat % boundary % hostData(1:3,i,j,1,iSide,iEl)
                nmag = this % geometry % nScale % boundary % hostData(i,j,1,iSide,iEl)

                un = this % u*nhat(1) + this % v*nhat(2) + this % w*nhat(3)
                dfdn = this % solutionGradient % boundary % hostData(1,i,j,ivar,iSide,iEl)*nhat(1) +&
                    this % solutionGradient % boundary % hostData(2,i,j,ivar,iSide,iEl)*nhat(2) +&
                    this % solutionGradient % boundary % hostData(3,i,j,ivar,iSide,iEl)*nhat(3)

                fin = this % solution % boundary % hostdata(i,j,ivar,iside,iel) ! interior solution
                fout = this % solution % extboundary % hostdata(i,j,ivar,iside,iel) ! exterior solution

                this % flux % boundaryNormal % hostData(i,j,1,iSide,iEl) =  ( 0.5_prec*( &
                    un*(fin + fout) + abs(un)*(fin - fout) ) -& ! advective flux
                    this % nu*dfdn )*nmag
                enddo
            enddo
          enddo
        enddo
      enddo
    
      end subroutine riemannsolver_advection_diffusion_3d
    
      
    end module self_advection_diffusion_3d