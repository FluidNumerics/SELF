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

module SELF_ColumnOperators_3D
!! Device backend of the column operators (see SELF_ColumnOperators_3D_t). Every
!! operator reads and writes the device copies (interior_gpu) of its fields and the
!! device copy of the covariant basis, and moves no data between host and device. The
!! caller uploads the inputs (UpdateDevice) and downloads the outputs (UpdateHost) when
!! it needs them on the host.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_2D
  use SELF_Scalar_3D
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_ColumnOperators_3D_t
  use SELF_GPU
  use iso_c_binding

  implicit none

  type,extends(ColumnOperators3D_t),public :: ColumnOperators3D
    type(c_ptr) :: upInverse_gpu = c_null_ptr
    type(c_ptr) :: downInverse_gpu = c_null_ptr

  contains

    procedure,public :: Init => Init_ColumnOperators3D
    procedure,public :: Free => Free_ColumnOperators3D
    procedure,public :: ColumnIntegral => ColumnIntegral_ColumnOperators3D
    procedure,public :: Broadcast => Broadcast_ColumnOperators3D
    procedure,public :: VerticalSweepUp => VerticalSweepUp_ColumnOperators3D
    procedure,public :: VerticalSweepDown => VerticalSweepDown_ColumnOperators3D

  endtype ColumnOperators3D

  interface
    subroutine ColumnIntegral_3D_gpu(q,dxds,qWeights,qbar,N,nLayers,nColumns,nvar) &
      bind(c,name="ColumnIntegral_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: q,dxds,qWeights,qbar
      integer(c_int),value :: N,nLayers,nColumns,nvar
    endsubroutine ColumnIntegral_3D_gpu
  endinterface

  interface
    subroutine ColumnBroadcast_3D_gpu(qbar,q,N,nLayers,nColumns,nvar) &
      bind(c,name="ColumnBroadcast_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: qbar,q
      integer(c_int),value :: N,nLayers,nColumns,nvar
    endsubroutine ColumnBroadcast_3D_gpu
  endinterface

  interface
    subroutine VerticalSweep_3D_gpu(aInv,r,inflow,q,N,nLayers,nColumns,nvar,upward) &
      bind(c,name="VerticalSweep_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: aInv,r,inflow,q
      integer(c_int),value :: N,nLayers,nColumns,nvar,upward
    endsubroutine VerticalSweep_3D_gpu
  endinterface

contains

  subroutine Init_ColumnOperators3D(this,interp,mesh,geometry)
    !! Init_ColumnOperators3D_t, followed by a copy of the two sweep inverses to the device
    implicit none
    class(ColumnOperators3D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    type(Mesh3D),intent(in) :: mesh
    type(SEMHex),target,intent(in) :: geometry

    call Init_ColumnOperators3D_t(this,interp,mesh,geometry)

    if(this%lobatto) then
      call gpuCheck(hipMalloc(this%upInverse_gpu,sizeof(this%upInverse)))
      call gpuCheck(hipMalloc(this%downInverse_gpu,sizeof(this%downInverse)))
      call gpuCheck(hipMemcpy(this%upInverse_gpu, &
                              c_loc(this%upInverse), &
                              sizeof(this%upInverse), &
                              hipMemcpyHostToDevice))
      call gpuCheck(hipMemcpy(this%downInverse_gpu, &
                              c_loc(this%downInverse), &
                              sizeof(this%downInverse), &
                              hipMemcpyHostToDevice))
    endif

  endsubroutine Init_ColumnOperators3D

  subroutine Free_ColumnOperators3D(this)
    implicit none
    class(ColumnOperators3D),intent(inout) :: this

    if(c_associated(this%upInverse_gpu)) call gpuCheck(hipFree(this%upInverse_gpu))
    if(c_associated(this%downInverse_gpu)) call gpuCheck(hipFree(this%downInverse_gpu))
    this%upInverse_gpu = c_null_ptr
    this%downInverse_gpu = c_null_ptr
    call Free_ColumnOperators3D_t(this)

  endsubroutine Free_ColumnOperators3D

  subroutine ColumnIntegral_ColumnOperators3D(this,q,qbar)
    !! Device version of ColumnIntegral_ColumnOperators3D_t
    implicit none
    class(ColumnOperators3D),intent(in) :: this
    class(Scalar3D),intent(in) :: q
    class(Scalar2D),intent(inout) :: qbar

    call CheckFields_ColumnOperators3D_t(this,q,qbar,"ColumnIntegral")
    call ColumnIntegral_3D_gpu(q%interior_gpu,this%geometry%dxds%interior_gpu, &
                               this%interp%qWeights_gpu,qbar%interior_gpu, &
                               this%N,this%nLayers,this%nColumns,q%nVar)

  endsubroutine ColumnIntegral_ColumnOperators3D

  subroutine Broadcast_ColumnOperators3D(this,qbar,q)
    !! Device version of Broadcast_ColumnOperators3D_t
    implicit none
    class(ColumnOperators3D),intent(in) :: this
    class(Scalar2D),intent(in) :: qbar
    class(Scalar3D),intent(inout) :: q

    call CheckFields_ColumnOperators3D_t(this,q,qbar,"Broadcast")
    call ColumnBroadcast_3D_gpu(qbar%interior_gpu,q%interior_gpu, &
                                this%N,this%nLayers,this%nColumns,q%nVar)

  endsubroutine Broadcast_ColumnOperators3D

  subroutine VerticalSweepUp_ColumnOperators3D(this,r,qBottom,q)
    !! Device version of VerticalSweepUp_ColumnOperators3D_t
    implicit none
    class(ColumnOperators3D),intent(in) :: this
    class(Scalar3D),intent(in) :: r
    class(Scalar2D),intent(in) :: qBottom
    class(Scalar3D),intent(inout) :: q

    call CheckSweep_ColumnOperators3D_t(this,r,qBottom,q,"VerticalSweepUp")
    call VerticalSweep_3D_gpu(this%upInverse_gpu,r%interior_gpu,qBottom%interior_gpu,q%interior_gpu, &
                              this%N,this%nLayers,this%nColumns,r%nVar,1)

  endsubroutine VerticalSweepUp_ColumnOperators3D

  subroutine VerticalSweepDown_ColumnOperators3D(this,r,qTop,q)
    !! Device version of VerticalSweepDown_ColumnOperators3D_t
    implicit none
    class(ColumnOperators3D),intent(in) :: this
    class(Scalar3D),intent(in) :: r
    class(Scalar2D),intent(in) :: qTop
    class(Scalar3D),intent(inout) :: q

    call CheckSweep_ColumnOperators3D_t(this,r,qTop,q,"VerticalSweepDown")
    call VerticalSweep_3D_gpu(this%downInverse_gpu,r%interior_gpu,qTop%interior_gpu,q%interior_gpu, &
                              this%N,this%nLayers,this%nColumns,r%nVar,0)

  endsubroutine VerticalSweepDown_ColumnOperators3D

endmodule SELF_ColumnOperators_3D
