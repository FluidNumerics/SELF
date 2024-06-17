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

module SELF_Tensor_2D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5
  use SELF_Data

  use HDF5
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(SELF_DataObj),public :: Tensor2D

    real(prec),pointer,dimension(:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:) :: extBoundary

  contains

    procedure,public :: Init => Init_Tensor2D
    procedure,public :: Free => Free_Tensor2D

    procedure,public :: BoundaryInterp => BoundaryInterp_Tensor2D
    procedure,public :: Divergence => Divergence_Tensor2D
    procedure,public :: DGDivergence => DGDivergence_Tensor2D
    procedure,public :: Determinant => Determinant_Tensor2D

  endtype Tensor2D

contains

  subroutine Init_Tensor2D(this,interp,nVar,nElem)
    implicit none
    class(Tensor2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! local
    integer :: i

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:2,1:2), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2,1:2), &
             this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2,1:2))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:4*nVar))

    ! Initialize equation parser
    ! This is done to prevent segmentation faults that arise
    ! when building with amdflang that are traced back to
    ! feqparse_functions.f90 : finalize routine
    ! When the equation parser is not initialized, the
    ! functions are not allocated, which I think are the
    ! source of the segfault - joe@fluidnumerics.com
    do i = 1,4*nvar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)

  endsubroutine Init_Tensor2D

  subroutine Free_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%extBoundary)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)

    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_Tensor2D

  subroutine BoundaryInterp_Tensor2D(this)
    implicit none
    class(Tensor2D),intent(inout) :: this

    call this%interp%TensorBoundaryInterp_2D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Tensor2D

  subroutine Divergence_Tensor2D(this,df)
    implicit none
    class(Tensor2D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2)

    call this%interp%TensorDivergence_2D(this%interior, &
                                         df, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine Divergence_Tensor2D

  subroutine DGDivergence_Tensor2D(this,df)
    implicit none
    class(Tensor2D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:2)

    call this%interp%TensorDGDivergence_2D(this%interior, &
                                           this%boundary, &
                                           df, &
                                           this%nVar, &
                                           this%nElem)

  endsubroutine DGDivergence_Tensor2D

  subroutine Determinant_Tensor2D(this,det)
    implicit none
    class(Tensor2D),intent(in) :: this
    real(prec),intent(out) :: det(1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j

    do iVar = 1,this%nVar
      do iEl = 1,this%nElem
        do j = 1,this%interp%N+1
          do i = 1,this%interp%N+1

            det(i,j,iEl,iVar) = this%interior(i,j,iEl,iVar,1,1)* &
                                this%interior(i,j,iEl,iVar,2,2)- &
                                this%interior(i,j,iEl,iVar,1,2)* &
                                this%interior(i,j,iEl,iVar,2,1)

          enddo
        enddo
      enddo
    enddo

  endsubroutine Determinant_Tensor2D

endmodule SELF_Tensor_2D
