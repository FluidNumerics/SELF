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

module self_shallow_water_2d

  use self_shallow_water_2d_t

  implicit none

  type,extends(shallow_water_2d_t) :: shallow_water_2d

  contains
    procedure :: setboundarycondition => setboundarycondition_shallow_water_2d
    procedure :: riemannsolver => riemannsolver_shallow_water_2d
    procedure :: fluxmethod => fluxmethod_shallow_water_2d
    procedure :: CalculateEntropy => CalculateEntropy_shallow_water_2d

  endtype shallow_water_2d

  interface
    subroutine setboundarycondition_shallow_water_2d_gpu(extboundary,boundary,sideinfo,nHat,N,nel,nvar) &
      bind(c,name="setboundarycondition_shallow_water_2d_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,sideinfo,nHat
      integer(c_int),value :: N,nel,nvar
    endsubroutine setboundarycondition_shallow_water_2d_gpu
  endinterface

  interface
    subroutine fluxmethod_shallow_water_2d_gpu(solution,flux,g,H,N,nel,nvar) &
      bind(c,name="fluxmethod_shallow_water_2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,flux
      real(c_prec),value :: g,H
      integer(c_int),value :: N,nel,nvar
    endsubroutine fluxmethod_shallow_water_2d_gpu
  endinterface

  interface
    subroutine riemannsolver_shallow_water_2d_gpu(fb,fextb,nhat,nscale,flux,g,H,N,nel,nvar) &
      bind(c,name="riemannsolver_shallow_water_2d_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fb,fextb,flux,nhat,nscale
      real(c_prec),value :: g,H
      integer(c_int),value :: N,nel,nvar
    endsubroutine riemannsolver_shallow_water_2d_gpu
  endinterface

contains
  subroutine CalculateEntropy_shallow_water_2d(this)
    implicit none
    class(shallow_water_2d),intent(inout) :: this
    ! Local
    integer :: iel,i,j,ivar
    real(prec) :: e,ei,jac
    real(prec) :: s(1:this%solution%nvar)

    call gpuCheck(hipMemcpy(c_loc(this%solution%interior), &
                            this%solution%interior_gpu,sizeof(this%solution%interior), &
                            hipMemcpyDeviceToHost))

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

  endsubroutine CalculateEntropy_shallow_water_2d

  subroutine setboundarycondition_shallow_water_2d(this)
    !! Boundary conditions are set to periodic boundary conditions
    implicit none
    class(shallow_water_2d),intent(inout) :: this

    call setboundarycondition_shallow_water_2d_gpu(this%solution%extboundary_gpu, &
                                                         this%solution%boundary_gpu,this%mesh%sideInfo_gpu,this%geometry%nHat%boundary_gpu,this%solution%interp%N, &
                                                         this%solution%nelem,this%solution%nvar)

  endsubroutine setboundarycondition_shallow_water_2d

  subroutine fluxmethod_shallow_water_2d(this)
    implicit none
    class(shallow_water_2d),intent(inout) :: this

    call fluxmethod_shallow_water_2d_gpu(this%solution%interior_gpu, &
                                               this%flux%interior_gpu, &
                                               this%g,this%H,this%solution%interp%N,this%solution%nelem, &
                                               this%solution%nvar)

  endsubroutine fluxmethod_shallow_water_2d

  subroutine riemannsolver_shallow_water_2d(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(shallow_water_2d),intent(inout) :: this

    call riemannsolver_shallow_water_2d_gpu(this%solution%boundary_gpu, &
                                                  this%solution%extBoundary_gpu, &
                                                  this%geometry%nhat%boundary_gpu,this%geometry%nscale%boundary_gpu, &
                                                  this%flux%boundarynormal_gpu,this%g,this%H,this%solution%interp%N, &
                                                  this%solution%nelem,this%solution%nvar)

  endsubroutine riemannsolver_shallow_water_2d

endmodule self_shallow_water_2d
