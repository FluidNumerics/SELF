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

module SELF_shallow_water_2d_t

    use self_model
    use self_dgmodel2d
    use self_mesh

    implicit none

    type,extends(dgmodel2d) :: shallow_water_2d_t
        real(prec) :: H = 1.0_prec ! uniform resting depth
        real(prec) :: g = 9.8_prec

    contains
        procedure :: setboundarycondition => setboundarycondition_shallow_water_2d_t
        ! procedure :: setgradientboundarycondition => setgradientboundarycondition_shallow_water_2d_t
        procedure :: riemannsolver => riemannsolver_shallow_water_2d_t
        procedure :: fluxmethod => fluxmethod_shallow_water_2d_t
        ! procedure :: sourcemethod => sourcemethod_shallow_water_2d_t
        procedure :: calculateentropy => calculateentropy_shallow_water_2d_t
        procedure :: Init => Init_shallow_water_2d_t 

    endtype shallow_water_2d_t

    contains
        subroutine setboundarycondition_shallow_water_2d_t(this)
            implicit none
            class(shallow_water_2d_t),intent(inout) :: this

            integer :: i,ivar,iEl,j,e2,bcid
            real(prec) :: s(1:this%solution%nvar)
            real(prec) :: nx, ny
            
            do iEl = 1,this%solution%nElem ! Loop over all elements
                do j = 1,4 ! Loop over all sides of each element
                    bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
                    e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID
                    
                    if (e2 == 0) then
                        if (bcid == SELF_BC_RADIATION) then
                            do i = 1,this%solution%interp%N+1
                                this%solution%extBoundary(i,j,iEl,1:3) = 0.0_prec
                            enddo
                        else if (bcid == SELF_BC_NONORMALFLOW) then
                            do i = 1,this%solution%interp%N+1
                                s = this%solution%boundary(i,j,iel,1:this%solution%nvar)
                                nx = this % geometry%nhat%boundary(i,j,iel,1,1)
                                ny = this % geometry%nhat%boundary(i,j,iel,1,2)

                                this%solution%extBoundary(i,j,iEl,1) = (ny * ny - nx * nx) * s(1) - 2 * nx * ny * s(2)
                                this%solution%extBoundary(i,j,iEl,2) = (nx * nx - ny * ny) * s(2) - 2 * nx * ny * s(1)
                                this%solution%extBoundary(i,j,iEl,3) = s(3)
                            enddo
                        end if
                    endif
                enddo
            enddo

        end subroutine setboundarycondition_shallow_water_2d_t
        ! subroutine setgradientboundarycondition_shallow_water_2d_t(this)
        !     implicit none
        !     class(shallow_water_2d_t),intent(inout) :: this

        !     integer :: i,ivar,iel,j,e2,bcid
        !     real(prec) :: internalstate(1:this%solution%nvar)
        !     ! real(prec) :: nx, ny

        !     do iEl = 1,this%solution%nElem ! Loop over all elements
        !         do j = 1,4 ! Loop over all sides of each element
        !             bcid = this % mesh % sideInfo(5,j,iEl) ! Boundary Condition ID
        !             e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID
                    
        !             if (e2 == 0) then
        !                 do i = 1,this%solution%interp%N+1
        !                     internalstate = this%solutiongradient%boundary(i,j,iel,1:this%solution%nvar)
        !                     ! nx = this % geometry%nhat%boundary(i,j,iel,1,1)
        !                     ! ny = this % geometry%nhat%boundary(i,j,iel,1,2)
        !                     this%solutiongradient%extBoundary(i,j,iEl,ivar) = -1 * this%solutiongradient%boundary(i,j,iel,1:this%solution%nvar)
        !                 enddo
        !             endif
        !         enddo
        !     enddo
        ! end subroutine
        subroutine fluxmethod_shallow_water_2d_t(this)
            implicit none
            class(shallow_water_2d_t),intent(inout) :: this
            ! Local
            integer :: iel, ivar, i, j
            real(prec) :: s(1:this%solution%nvar)

            do concurrent( &
                i = 1:this%solution%interp%N+1, &
                j = 1:this%solution%interp%N+1, &
                iel = 1:this%mesh%nelem &
            )
                s = this%solution%interior(i,j,iel,1:this%solution%nvar)


                this%flux%interior(i,j,iel,1,1) = this%g * s(3)
                this%flux%interior(i,j,iel,1,2) = 0.0_prec

                this%flux%interior(i,j,iel,2,1) = 0.0_prec
                this%flux%interior(i,j,iel,2,2) = this%g * s(3)

                this%flux%interior(i,j,iel,3,1) = this%H * s(1)
                this%flux%interior(i,j,iel,3,2) = this%H * s(2)
            end do
        end subroutine fluxmethod_shallow_water_2d_t
        subroutine riemannsolver_shallow_water_2d_t(this)
            implicit none
            class(shallow_water_2d_t),intent(inout) :: this
            ! Local
            integer :: iel, ivar, i, j
            real(prec) :: fin(1:3), fout(1:3), unl, unr
            real(prec) :: nx, ny, nmag
            real(prec) :: c

            c = sqrt(this%g * this%H)

            do concurrent( &
                i = 1:this%solution%interp%N+1, &
                j = 1:4, &
                iel = 1:this%mesh%nelem &
            )
                nx = this % geometry%nhat%boundary(i,j,iel,1,1)
                ny = this % geometry%nhat%boundary(i,j,iel,1,2)
                nmag = this%geometry%nScale%boundary(i,j,iEl,1)

                fin = this%solution%boundary(i,j,iel,1:this%solution%nvar)
                fout = this%solution%extboundary(i,j,iel,1:this%solution%nvar)
                
                unl = fin(1) * nx + fin(2) * ny
                unr = fout(1) * nx + fout(2) * ny

                this%flux%boundaryNormal(i,j,iEl,1) = 0.5_prec * (this%g * (fin(3) + fout(3)) + c * (unl - unr)) * nx * nmag
                this%flux%boundaryNormal(i,j,iEl,2) = 0.5_prec * (this%g * (fin(3) + fout(3)) + c * (unl - unr)) * ny * nmag
                this%flux%boundaryNormal(i,j,iEl,3) = 0.5_prec * (this%H * (unl + unr) + c * (fin(3) - fout(3))) * nmag

            end do

        end subroutine riemannsolver_shallow_water_2d_t

        subroutine calculateentropy_shallow_water_2d_t(this)
            implicit none
            class(shallow_water_2d_t),intent(inout) :: this
            ! Local
            integer :: iel, ivar, i, j
            real(prec) :: e, ei, jac
            real(prec) :: s(1:this%solution%nvar)

            e = 0.0_prec
            do iel = 1,this%geometry%nelem
                do j = 1,this%solution%interp%N+1
                    do i = 1,this%solution%interp%N+1
                        jac = this%geometry%J%interior(i,j,iel,1)
                        s(1:this%solution%nvar) = this%solution%interior(i,j,iel,1:this%solution%nvar)
                        ei =  0.5_prec * (this%H * s(1) * s(1) + this%H * s(2) * s(2)) + &
                              0.5_prec * this%g * s(3) * s(3)
                        e = e + ei * jac
                    enddo
                enddo
            enddo

            this%entropy = e

        end subroutine calculateentropy_shallow_water_2d_t

        subroutine Init_shallow_water_2d_t(this,nvar,mesh,geometry,decomp)
            implicit none
            class(shallow_water_2d_t),intent(out) :: this
            integer,intent(in) :: nvar
            type(Mesh2D),intent(in),target :: mesh
            type(SEMQuad),intent(in),target :: geometry
            type(MPILayer),intent(in),target :: decomp
            ! Local
            integer :: nvar_three = 3
            integer :: ivar
            character(LEN=3) :: ivarChar
            character(LEN=25) :: varname
        
            print*,__FILE__ // " init nvar = 3"

            this%decomp => decomp
            this%mesh => mesh
            this%geometry => geometry
        
            call this%solution%Init(geometry%x%interp,nvar_three,this%mesh%nElem)
            call this%workSol%Init(geometry%x%interp,nvar_three,this%mesh%nElem)
            call this%dSdt%Init(geometry%x%interp,nvar_three,this%mesh%nElem)
            call this%solutionGradient%Init(geometry%x%interp,nvar_three,this%mesh%nElem)
            call this%flux%Init(geometry%x%interp,nvar_three,this%mesh%nElem)
            call this%source%Init(geometry%x%interp,nvar_three,this%mesh%nElem)
            call this%fluxDivergence%Init(geometry%x%interp,nvar_three,this%mesh%nElem)
        
            ! set default metadata
            ! do ivar = 1,nvar_three
            !   write(ivarChar,'(I3.3)') ivar
            !   varname = "solution"//trim(ivarChar)
            !   call this%solution%SetName(ivar,varname)
            !   call this%solution%SetUnits(ivar,"[null]")
            ! enddo
            call this%solution%SetName(1,"u")
            call this%solution%SetUnits(1,"[null]")
            call this%solution%SetName(2,"v")
            call this%solution%SetUnits(2,"[null]")
            call this%solution%SetName(3,"eta")
            call this%solution%SetUnits(3,"[null]")

            call this%solution%AssociateGeometry(geometry)
            call this%solutionGradient%AssociateGeometry(geometry)
            call this%flux%AssociateGeometry(geometry)
            call this%fluxDivergence%AssociateGeometry(geometry)
        
          endsubroutine Init_shallow_water_2d_t

endmodule SELF_shallow_water_2d_t