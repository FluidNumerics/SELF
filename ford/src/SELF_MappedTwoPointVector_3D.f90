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
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_MappedTwoPointVector_3D

  use SELF_MappedTwoPointVector_3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(MappedTwoPointVector3D_t),public :: MappedTwoPointVector3D

  contains

    generic,public :: MappedDivergence => MappedDivergence_MappedTwoPointVector3D
    procedure,private :: MappedDivergence_MappedTwoPointVector3D

  endtype MappedTwoPointVector3D

contains

  subroutine MappedDivergence_MappedTwoPointVector3D(this,df)
    !! GPU implementation of the physical-space split-form divergence for
    !! a 3-D two-point vector on a curvilinear mesh.
    !!
    !! interior_gpu must already contain pre-projected SCALAR contravariant
    !! two-point fluxes.  Applies the reference-element split-form sum then
    !! divides by J.
    implicit none
    class(MappedTwoPointVector3D),intent(inout) :: this
    type(c_ptr),intent(out) :: df

    call TwoPointVectorDivergence_3D_gpu(this%interior_gpu,df, &
                                         this%interp%dSplitMatrix_gpu, &
                                         this%interp%N,this%nVar,this%nElem)

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu, &
                               this%interp%N,this%nVar,this%nElem)

  endsubroutine MappedDivergence_MappedTwoPointVector3D

endmodule SELF_MappedTwoPointVector_3D
