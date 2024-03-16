module self_advection_diffusion_2d

    use self_model
    use self_dgmodel2d
    use self_mesh
    
    implicit none
    
      type, extends(dgmodel2d) :: advection_diffusion_2d
        real(prec) :: nu ! diffusion coefficient
        real(prec) :: u  ! constant x-component of velocity
        real(prec) :: v  ! constant y-component of velocity
    
        contains
        procedure :: pretendency => pretendency_advection_diffusion_2d
        procedure :: setboundarycondition => setboundarycondition_advection_diffusion_2d
        procedure :: riemannsolver => riemannsolver_advection_diffusion_2d
        procedure :: fluxmethod => fluxmethod_advection_diffusion_2d
    
      end type advection_diffusion_2d
    
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
    
      subroutine pretendency_advection_diffusion_2d(this)
      ! Here, we use the pre-tendency method to calculate the
      ! derivative of the solution using a bassi-rebay method
      ! We then do a boundary interpolation and side exchange
      ! on the gradient field
      implicit none
      class(advection_diffusion_2d), intent(inout) :: this
      ! local
      integer :: i, ivar, iEl, j, e2    
    
      do ivar = 1, this % solution % nvar
        do iEl = 1,this % solution % nElem ! Loop over all elements
          do j = 1,4 ! Loop over all sides

            !bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
            e2 = this % mesh % sideInfo(3,j,iEl) ! Neighboring Element ID

            if (e2 == 0) then
              do i = 1,this % solution % interp % N+1 ! Loop over quadrature points
                this % solution % extBoundary(i,j,iEl,ivar) = 0.0_prec
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
    
    
      end subroutine pretendency_advection_diffusion_2d
    
      subroutine setboundarycondition_advection_diffusion_2d(this)
      ! Here, we set the boundary conditions for the 
      ! solution and the solution gradient at the left
      ! and right most boundaries.
      ! 
      ! Here, we use periodic boundary conditions
      implicit none
      class(advection_diffusion_2d), intent(inout) :: this
      ! local
      integer :: i, ivar, iEl, j, e2 

        do ivar = 1, this % solution % nvar

          do iEl = 1,this % solution % nElem ! Loop over all elements
            do j = 1,4 ! Loop over all sides

              !bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
              e2 = this % mesh % sideInfo(3,j,iEl) ! Neighboring Element ID

              if (e2 == 0) then

                do i = 1,this % solution % interp % N+1 ! Loop over quadrature points
                    this % solutionGradient % extBoundary(i,j,iEl,ivar,1:2) = 0.0_prec
                end do
                
              end if
    
            end do
          end do
        end do

        call this % solutionGradient % BassiRebaySides()

      end subroutine setboundarycondition_advection_diffusion_2d
    
      subroutine fluxmethod_advection_diffusion_2d(this)
      implicit none
      class(advection_diffusion_2d), intent(inout) :: this
      ! Local
      integer :: iel
      integer :: ivar
      integer :: i
      integer :: j
      real(prec) :: u, v, nu, f, dfdx, dfdy
    
        u = this % u
        v = this % v
        nu = this % nu
        do ivar = 1, this % solution % nvar
          do iel = 1, this % mesh % nelem
            do j = 1, this % solution % interp % N+1
              do i = 1, this % solution % interp % N+1
        
                f = this % solution % interior(i,j,iel,ivar)
                dfdx = this % solutionGradient % interior(i,j,iel,ivar,1)
                dfdy = this % solutionGradient % interior(i,j,iel,ivar,2)
        
                this % flux % interior(i,j,iel,ivar,1) = u*f - nu*dfdx  ! advective flux + diffusive flux (x-component)
                this % flux % interior(i,j,iel,ivar,2) = v*f - nu*dfdy  ! advective flux + diffusive flux (y-component)
        
              enddo
            enddo
          enddo
        enddo
    
      end subroutine fluxmethod_advection_diffusion_2d
    
      subroutine riemannsolver_advection_diffusion_2d(this)
      ! this method uses an linear upwind solver for the
      ! advective flux and the bassi-rebay method for the 
      ! diffusive fluxes
      implicit none
      class(advection_diffusion_2d), intent(inout) :: this
      ! Local
      integer :: iel
      integer :: ivar
      integer :: j
      integer :: i
      real(prec) :: fin, fout, dfdn, un
      real(prec) :: nhat(1:2), nmag
        
      do ivar = 1, this % solution % nvar
        do iEl = 1,this % solution % nElem
          do j = 1,4
            do i = 1,this % solution % interp % N+1

              ! Get the boundary normals on cell edges from the mesh geometry
              nhat(1:2) = this % geometry % nHat % boundary(i,j,iEl,1,1:2)
      
              un = this % u*nhat(1) + this % v*nhat(2)
              dfdn = this % solutionGradient % avgBoundary(i,j,iEl,ivar,1)*nhat(1) +&
                this % solutionGradient % avgBoundary(i,j,iEl,ivar,2)*nhat(2)

              fin = this % solution % boundary(i,j,iel,ivar) ! interior solution
              fout = this % solution % extboundary(i,j,iel,ivar) ! exterior solution
      
              nmag = this % geometry % nScale % boundary(i,j,iEl,1)

              this % flux % boundaryNormal(i,j,iEl,1) =  ( 0.5_prec*( &
                  un*(fin + fout) + abs(un)*(fin - fout) ) -& ! advective flux
                  this % nu*dfdn )*nmag

            enddo
          enddo
        enddo
      enddo
    
      end subroutine riemannsolver_advection_diffusion_2d
    
      
    end module self_advection_diffusion_2d