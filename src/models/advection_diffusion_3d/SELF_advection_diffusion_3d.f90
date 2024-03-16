module self_advection_diffusion_3d

    use self_model
    use self_dgmodel3d
    use self_mesh
    
    implicit none
    
      type, extends(dgmodel3d) :: advection_diffusion_3d
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
      integer :: i, j, ivar, iEl, k, e2    

      do ivar = 1, this % solution % nvar
        do iEl = 1,this % solution % nElem ! Loop over all elements
          do k = 1,6 ! Loop over all sides

            !bcid = this % mesh % sideInfo(5,k,iEl) ! Boundary Condition ID
            e2 = this % mesh % sideInfo(3,k,iEl) ! Neighboring Element ID

            if (e2 == 0) then
              do j = 1,this % solution % interp % N+1 ! Loop over quadrature point
                do i = 1,this % solution % interp % N+1 ! Loop over quadrature points
                  this % solution % extBoundary(i,j,k,iEl,iVar) = 0.0_prec
                enddo
              enddo
            end if
          end do
        end do
      end do
    
        ! calculate the averages of the solutions on the element
        ! boundaries and store is this % solution % avgBoundary
        call this % solution % BassiRebaySides()
    
        ! calculate the derivative using the bassi-rebay form
        call this % solution % BRGradient(this % geometry,this % solutionGradient)
    
        ! interpolate the solutiongradient to the element boundaries
        call this % solutionGradient % BoundaryInterp()
      
        ! perform the side exchange to populate the solutionGradient % extBoundary attribute
        call this % solutionGradient % SideExchange(this % mesh, this % decomp) 
      
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
      integer :: i, j, ivar, iEl, k, e2 

      do ivar = 1, this % solution % nvar
        do iEl = 1,this % solution % nElem ! Loop over all elements
          do k = 1,6 ! Loop over all sides

            !bcid = this % mesh % sideInfo(5,k,iEl) ! Boundary Condition ID
            e2 = this % mesh % sideInfo(3,k,iEl) ! Neighboring Element ID

            if (e2 == 0) then
              do j = 1,this % solution % interp % N+1 ! Loop over quadrature point
                do i = 1,this % solution % interp % N+1 ! Loop over quadrature points
                  this % solutionGradient % extBoundary(i,j,k,iEl,iVar,1:3) = 0.0_prec
                enddo
              enddo
            end if
          end do
        end do
      end do

        call this % solutionGradient % BassiRebaySides()

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
        do ivar = 1, this % solution % nvar
          do iel = 1, this % mesh % nelem
            do k = 1, this % solution % interp % N+1
              do j = 1, this % solution % interp % N+1
                do i = 1, this % solution % interp % N+1
            
                  f = this % solution % interior(i,j,k,iel,ivar)
                  dfdx = this % solutionGradient % interior(i,j,k,iel,ivar,1)
                  dfdy = this % solutionGradient % interior(i,j,k,iel,ivar,2)
                  dfdz = this % solutionGradient % interior(i,j,k,iel,ivar,3)
          
                  this % flux % interior(i,j,k,iel,ivar,1) = u*f - nu*dfdx  ! advective flux + diffusive flux (x-component)
                  this % flux % interior(i,j,k,iel,ivar,2) = v*f - nu*dfdy  ! advective flux + diffusive flux (y-component)
                  this % flux % interior(i,j,k,iel,ivar,3) = w*f - nu*dfdz  ! advective flux + diffusive flux (z-component)
          
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
      integer :: k
      integer :: i, j
      real(prec) :: fin, fout, dfdn, un
      real(prec) :: nhat(1:3), nmag
        
      do ivar = 1, this % solution % nvar
        do iEl = 1,this % solution % nElem
          do k = 1,6
            do j = 1,this % solution % interp % N+1
              do i = 1,this % solution % interp % N+1

                ! Get the boundary normals on cell edges from the mesh geometry
                nhat(1:3) = this % geometry % nHat % boundary(i,j,k,iEl,1,1:3)

                un = this % u*nhat(1) + this % v*nhat(2) + this % w*nhat(3)
                dfdn = this % solutionGradient % avgBoundary(i,j,k,iEl,iVar,1)*nhat(1) +&
                    this % solutionGradient % avgBoundary(i,j,k,iEl,iVar,2)*nhat(2) +&
                    this % solutionGradient % avgBoundary(i,j,k,iEl,iVar,3)*nhat(3)

                fin = this % solution % boundary(i,j,k,iEl,iVar) ! interior solution
                fout = this % solution % extboundary(i,j,k,iEl,iVar) ! exterior solution

                nmag = this % geometry % nScale % boundary(i,j,k,iEl,1)

                this % flux % boundaryNormal(i,j,k,iEl,1) =  ( 0.5_prec*( &
                    un*(fin + fout) + abs(un)*(fin - fout) ) -& ! advective flux
                    this % nu*dfdn )*nmag
              enddo
            enddo
          enddo
        enddo
      enddo
    
      end subroutine riemannsolver_advection_diffusion_3d
    
      
    end module self_advection_diffusion_3d