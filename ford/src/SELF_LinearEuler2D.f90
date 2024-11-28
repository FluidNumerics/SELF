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

module self_LinearEuler2D

  use self_LinearEuler2D_t

  implicit none

  type,extends(LinearEuler2D_t) :: LinearEuler2D
  contains
    procedure :: setboundarycondition => setboundarycondition_LinearEuler2D
    procedure :: boundaryflux => boundaryflux_LinearEuler2D
    procedure :: fluxmethod => fluxmethod_LinearEuler2D
    procedure :: sourcemethod => sourcemethod_LinearEuler2D

  endtype LinearEuler2D

  interface
    subroutine setboundarycondition_LinearEuler2D_gpu(extboundary,boundary,sideinfo,nhat,N,nel,nvar) &
      bind(c,name="setboundarycondition_LinearEuler2D_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,sideinfo,nhat
      integer(c_int),value :: N,nel,nvar
    endsubroutine setboundarycondition_LinearEuler2D_gpu
  endinterface

  interface
    subroutine fluxmethod_LinearEuler2D_gpu(solution,flux,rho0,c,N,nel,nvar) &
      bind(c,name="fluxmethod_LinearEuler2D_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,flux
      real(c_prec),value :: rho0,c
      integer(c_int),value :: N,nel,nvar
    endsubroutine fluxmethod_LinearEuler2D_gpu
  endinterface

  interface
    subroutine boundaryflux_LinearEuler2D_gpu(fb,fextb,nhat,nscale,flux,rho0,c,N,nel,nvar) &
      bind(c,name="boundaryflux_LinearEuler2D_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,flux,nhat,nscale
      real(c_prec),value :: rho0,c
      integer(c_int),value :: N,nel,nvar
    endsubroutine boundaryflux_LinearEuler2D_gpu
  endinterface

contains

  subroutine sourcemethod_LinearEuler2D(this)
    implicit none
    class(LinearEuler2D),intent(inout) :: this

    return

  endsubroutine sourcemethod_LinearEuler2D

  subroutine boundaryflux_LinearEuler2D(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(LinearEuler2D),intent(inout) :: this

    call boundaryflux_LinearEuler2D_gpu(this%solution%boundary_gpu, &
                                        this%solution%extBoundary_gpu, &
                                        this%geometry%nhat%boundary_gpu, &
                                        this%geometry%nscale%boundary_gpu, &
                                        this%flux%boundarynormal_gpu, &
                                        this%rho0,this%c,this%solution%interp%N, &
                                        this%solution%nelem,this%solution%nvar)

  endsubroutine boundaryflux_LinearEuler2D

  subroutine fluxmethod_LinearEuler2D(this)
    implicit none
    class(LinearEuler2D),intent(inout) :: this

    call fluxmethod_LinearEuler2D_gpu(this%solution%interior_gpu, &
                                      this%flux%interior_gpu, &
                                      this%rho0,this%c,this%solution%interp%N,this%solution%nelem, &
                                      this%solution%nvar)

  endsubroutine fluxmethod_LinearEuler2D

  subroutine setboundarycondition_LinearEuler2D(this)
    !! Boundary conditions are set to periodic boundary conditions
    implicit none
    class(LinearEuler2D),intent(inout) :: this
    ! local
    integer :: i,iEl,j,e2,bcid
    real(prec) :: x(1:2)

    if(this%prescribed_bcs_enabled) then
      call gpuCheck(hipMemcpy(c_loc(this%solution%extboundary), &
                              this%solution%extboundary_gpu,sizeof(this%solution%extboundary), &
                              hipMemcpyDeviceToHost))

      ! Prescribed boundaries are still done on the GPU
      do iEl = 1,this%solution%nElem ! Loop over all elements
        do j = 1,4 ! Loop over all sides

          bcid = this%mesh%sideInfo(5,j,iEl) ! Boundary Condition ID
          e2 = this%mesh%sideInfo(3,j,iEl) ! Neighboring Element ID

          if(e2 == 0) then
            if(bcid == SELF_BC_PRESCRIBED) then

              do i = 1,this%solution%interp%N+1 ! Loop over quadrature points
                x = this%geometry%x%boundary(i,j,iEl,1,1:2)

                this%solution%extBoundary(i,j,iEl,1:this%nvar) = &
                  this%hbc2d_Prescribed(x,this%t)
              enddo

            endif
          endif

        enddo
      enddo

      call gpuCheck(hipMemcpy(this%solution%extBoundary_gpu, &
                              c_loc(this%solution%extBoundary), &
                              sizeof(this%solution%extBoundary), &
                              hipMemcpyHostToDevice))
    endif
    call setboundarycondition_LinearEuler2D_gpu(this%solution%extboundary_gpu, &
                                                this%solution%boundary_gpu, &
                                                this%mesh%sideInfo_gpu, &
                                                this%geometry%nhat%boundary_gpu, &
                                                this%solution%interp%N, &
                                                this%solution%nelem,this%solution%nvar)

  endsubroutine setboundarycondition_LinearEuler2D

endmodule self_LinearEuler2D
