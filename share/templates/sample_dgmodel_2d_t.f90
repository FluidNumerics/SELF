! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module sample_dgmodel_2d_t

    use self_model
    use self_dgmodel2d
    use self_mesh
  
    implicit none
  
    type,extends(dgmodel2d) :: sample_dgmodel2d_t
      ! Here you can declare any attributes you may need in addition to the 
      ! attributes provided by the dgmodel2d class. These could be things like
      ! physical constants, additional diagnostic fields, procedure pointers,
      ! whatever you need. 
  
    contains
      ! This is the list of procedures you will need to define at a minimum
      ! to define your conservation law solver.
      !
      procedure :: setboundarycondition => setboundarycondition_sample_dgmodel2d_t
      procedure :: setgradientboundarycondition => setgradientboundarycondition_sample_dgmodel2d_t
      procedure :: riemannsolver => riemannsolver_sample_dgmodel2d_t
      procedure :: fluxmethod => fluxmethod_sample_dgmodel2d_t
      procedure :: sourcemethod => sourcemethod_sample_dgmodel2d_t
      procedure :: calculateentropy => calculateentropy_sample_dgmodel2d_t
  
    endtype sample_dgmodel2d_t
  
  contains
  
    subroutine calculateentropy_sample_dgmodel2d_t(this)
      ! Here we define a measure of the model's entropy.
      ! The "entropy" is a mathematical entropy function.
      ! It is a sign definite convex function that depends 
      ! on the prognostic state variables of the conservation
      ! law you are solving. For physical systems like the shallow
      ! water equations, the total energy (kinetic plus potential)
      ! is one possible entropy function.
      ! In this routine, we specifically compute the integral of the
      ! entropy function over the model domain.
      implicit none
      class(sample_dgmodel2d_t),intent(inout) :: this
      ! Local
      integer :: iel,i,j,ivar
      real(prec) :: e,ei,jac
      real(prec) :: s(1:this%solution%nvar)
  
      e = 0.0_prec ! The entropy is initialized to zero
      do iel = 1,this%geometry%nelem
        do j = 1,this%solution%interp%N+1
          do i = 1,this%solution%interp%N+1 
            ! These tightly nested loops loop over all of the points
            ! in the model space.
            jac = this%geometry%J%interior(i,j,iel,1) ! This is the jacobian of the coordinate transformation and is needed for computing integrals
            s(1:this%solution%nvar) = this%solution%interior(i,j,iel,1:this%solution%nvar) ! Here, we are getting the prognostic variables at each point
            ! ei =  ! Uncomment this line and define the entropy at this point in terms of the prognostic variables.
            e = e + ei*jac
          enddo
        enddo
      enddo
  
      this%entropy = e
  
    endsubroutine calculateentropy_sample_dgmodel2d_t
  
    subroutine setboundarycondition_sample_dgmodel2d_t(this)
      !! Boundary conditions for the solution are set to 
      !! 0 for the external state to provide radiation type
      !! boundary conditions.
      implicit none
      class(sample_dgmodel2d_t),intent(inout) :: this
      ! local
      integer :: i,ivar,iEl,j,e2,bcid
      real(prec) :: internalstate(1:this%solution%nvar)
      real(prec) :: nx,ny
  

        do iEl = 1,this%solution%nElem ! Loop over all elements
          do j = 1,4 ! Loop over all sides of each element
  
            bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
            e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID
  
            if(e2 == 0) then ! When the neighboring element id is zero, then there is no neighboring element
              ! SELF has a few pre-defined integers that can be used for distinguishing different
              ! boundary conditions. You can use these conditionals as you need. Or, if you only want
              ! one boundary condition type, just remove these conditionals on the `bcid` value.

              if(bcid == SELF_BC_NONORMALFLOW)then

                  do i = 1,this%solution%interp%N+1 ! Loop over quadrature points on this boundary
                    ! To set the boundary condition, we need to set the external state (solution%extboundary)
                    ! using our knowledge of the internal state.
                    internalstate = this%solution%boundary(i,j,iel,1:this%solution%nvar)  ! Get the internal state

                    ! You may sometimes need to know the boundary normal direction to compute your
                    ! external state. These next two lines get the x and y components of the bounary normal
                    nx = this % geometry %nhat%boundary(i,j,iel,1,1) ! x-component of the boundary unit normal
                    ny = this % geometry %nhat%boundary(i,j,iel,1,2) ! y-component of the boundary unit normal

                    !this%solution%extBoundary(i,j,iEl,ivar) = 0.0_prec ! Uncomment this line and define the external state in terms of the internal state.
                  enddo
                enddo
                
            !   elseif(bcid == SELF_BC_PRESCRIBED)then
            !   elseif(bcid == SELF_BC_RADIATION)then
            !   else
              endif    
            endif
  
          enddo
        enddo
      enddo
  
    endsubroutine setboundarycondition_sample_dgmodel2d_t
  
    subroutine setgradientboundarycondition_sample_dgmodel2d_t(this)
      !! Boundary conditions on the solution gradient are set
      !! here. Gradient boundary conditions are often set for stress terms
      !! 
      implicit none
      class(sample_dgmodel2d_t),intent(inout) :: this
      integer :: i,ivar,iEl,j,e2,bcid
      real(prec) :: internalstate(1:this%solution%nvar)
      real(prec) :: nx,ny
  

        do iEl = 1,this%solution%nElem ! Loop over all elements
          do j = 1,4 ! Loop over all sides of each element
  
            bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
            e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID
  
            if(e2 == 0) then ! When the neighboring element id is zero, then there is no neighboring element
              ! SELF has a few pre-defined integers that can be used for distinguishing different
              ! boundary conditions. You can use these conditionals as you need. Or, if you only want
              ! one boundary condition type, just remove these conditionals on the `bcid` value.

              if(bcid == SELF_BC_NOSTRESS)then

                  do i = 1,this%solution%interp%N+1 ! Loop over quadrature points on this boundary
                    ! To set the boundary condition, we need to set the external state (solutiongradient%extboundary)
                    ! using our knowledge of the internal state.
                    internalstate = this%solutiongradient%boundary(i,j,iel,1:this%solution%nvar)  ! Get the internal state

                    ! You may sometimes need to know the boundary normal direction to compute your
                    ! external state. These next two lines get the x and y components of the bounary normal
                    nx = this % geometry %nhat%boundary(i,j,iel,1,1) ! x-component of the boundary unit normal
                    ny = this % geometry %nhat%boundary(i,j,iel,1,2) ! y-component of the boundary unit normal

                    !this%solutiongradient%extBoundary(i,j,iEl,ivar) = 0.0_prec ! Uncomment this line and define the external state in terms of the internal state.
                  enddo
                enddo
                
            !   elseif(bcid == SELF_BC_PRESCRIBED_STRESS)then
            !   else
              endif    
            endif
  
          enddo
        enddo
      enddo

    endsubroutine setgradientboundarycondition_sample_dgmodel2d_t
  
    subroutine fluxmethod_sample_dgmodel2d_t(this)
      implicit none
      class(sample_dgmodel2d_t),intent(inout) :: this
      ! Local
      integer :: iel
      integer :: ivar
      integer :: i
      integer :: j
      real(prec) :: s(1:this%solution%nvar)
  

        do iel = 1,this%mesh%nelem
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1
  
              f = this%solution%interior(i,j,iel,1:this%solution%nvarr)
              dfdx = this%solutionGradient%interior(i,j,iel,1:this%solution%nvar,1)
              dfdy = this%solutionGradient%interior(i,j,iel,1:this%solution%nvar,2)
  
              this%flux%interior(i,j,iel,1,1) = ! x-component of the flux for the first variable
              this%flux%interior(i,j,iel,1,2) = ! y-component of the flux for the first variable

              this%flux%interior(i,j,iel,2,1) = ! x-component of the flux for the second variable
              this%flux%interior(i,j,iel,2,2) = ! y-component of the flux for the second variable

              this%flux%interior(i,j,iel,3,1) = ! x-component of the flux for the third variable
              this%flux%interior(i,j,iel,3,2) = ! y-component of the flux for the third variable
  
            enddo
          enddo
        enddo

  
    endsubroutine fluxmethod_sample_dgmodel2d_t
  
    subroutine riemannsolver_sample_dgmodel2d_t(this)
      ! Set the riemann solver for the element edges
      implicit none
      class(sample_dgmodel2d_t),intent(inout) :: this
      ! Local
      integer :: iel
      integer :: ivar
      integer :: j
      integer :: i
      real(prec) :: fin,fout,dfdn,un
      real(prec) :: nx,ny,nmag
  
        do iEl = 1,this%solution%nElem
          do j = 1,4
            do i = 1,this%solution%interp%N+1
  
              ! Get the boundary normals on cell edges from the mesh geometry
              nx = this%geometry%nHat%boundary(i,j,iEl,1,1)
              ny = this%geometry%nHat%boundary(i,j,iEl,1,2)
  
              dfdn = this%solutionGradient%avgboundary(i,j,iEl,ivar,1)*nx+ &
                     this%solutionGradient%avgboundary(i,j,iEl,ivar,2)*ny
  
              fin = this%solution%boundary(i,j,iel,ivar) ! interior solution
              fout = this%solution%extboundary(i,j,iel,ivar) ! exterior solution
  
              nmag = this%geometry%nScale%boundary(i,j,iEl,1)
  
              this%flux%boundaryNormal(i,j,iEl,1) = ! Set the boundary normal component of the flux for the first variable 
              this%flux%boundaryNormal(i,j,iEl,2) = ! Set the boundary normal component of the flux for the second variable 
              this%flux%boundaryNormal(i,j,iEl,3) = ! Set the boundary normal component of the flux for the third variable 


  
            enddo
          enddo
        enddo
  
    endsubroutine riemannsolver_sample_dgmodel2d_t

    subroutine sourcemethod_sample_dgmodel2d_t(this)
      implicit none
      class(sample_dgmodel2d_t),intent(inout) :: this
      ! Local
      integer :: iel,i,j,ivar
      real(prec) :: s(1:this%solution%nvar)
  
      do iel = 1,this%geometry%nelem
        do j = 1,this%solution%interp%N+1
          do i = 1,this%solution%interp%N+1 
            ! These tightly nested loops loop over all of the points
            ! in the model space.
            s(1:this%solution%nvar) = this%solution%interior(i,j,iel,1:this%solution%nvar) ! Here, we are getting the prognostic variables at each point
            this%source%interior(i,j,iel,1) = ! Source term for first prognostic variable
            this%source%interior(i,j,iel,2) = ! Source term for second prognostic variable
            this%source%interior(i,j,iel,3) = ! Source term for third prognostic variable
          enddo
        enddo
      enddo
  
    endsubroutine sourcemethod_sample_dgmodel2d_t
  
  endmodule sample_dgmodel_2d_t
  