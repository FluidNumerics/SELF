! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Tensor_3D

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

  type,extends(SELF_DataObj),public :: Tensor3D

    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: interior
    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: boundary
    real(prec),pointer,dimension(:,:,:,:,:,:,:) :: extBoundary

  contains

    procedure,public :: Init => Init_Tensor3D
    procedure,public :: Free => Free_Tensor3D

    procedure,public :: BoundaryInterp => BoundaryInterp_Tensor3D
    procedure,public :: Divergence => Divergence_Tensor3D
    procedure,public :: DGDivergence => DGDivergence_Tensor3D
    procedure,public :: Determinant => Determinant_Tensor3D

  endtype Tensor3D

contains

  subroutine Init_Tensor3D(this,interp,nVar,nElem)
    implicit none
    class(Tensor3D),intent(out) :: this
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

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:3,1:3), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3,1:3), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar,1:3,1:3))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:9*nVar))


    ! Initialize equation parser
    ! This is done to prevent segmentation faults that arise
    ! when building with amdflang that are traced back to 
    ! feqparse_functions.f90 : finalize routine
    ! When the equation parser is not initialized, the 
    ! functions are not allocated, which I think are the 
    ! source of the segfault - joe@fluidnumerics.com
    do i = 1, 9*nvar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)

  endsubroutine Init_Tensor3D

  subroutine Free_Tensor3D(this)
    implicit none
    class(Tensor3D),intent(inout) :: this

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

  endsubroutine Free_Tensor3D

  subroutine BoundaryInterp_Tensor3D(this)
    implicit none
    class(Tensor3D),intent(inout) :: this

    call this%interp%TensorBoundaryInterp_3D(this%interior, &
                                             this%boundary, &
                                             this%nVar, &
                                             this%nElem)

  endsubroutine BoundaryInterp_Tensor3D

  subroutine Divergence_Tensor3D(this,df)
    implicit none
    class(Tensor3D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3)

    call this%interp%TensorDivergence_3D(this%interior, &
                                         df, &
                                         this%nVar, &
                                         this%nElem)

  endsubroutine Divergence_Tensor3D

  subroutine DGDivergence_Tensor3D(this,df)
    implicit none
    class(Tensor3D),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar,1:3)

    call this%interp%TensorDGDivergence_3D(this%interior, &
                                           this%boundary, &
                                           df, &
                                           this%nVar, &
                                           this%nElem)

  endsubroutine DGDivergence_Tensor3D

  subroutine Determinant_Tensor3D(this,det)
    implicit none
    class(Tensor3D),intent(in) :: this
    real(prec),intent(out) :: det(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nelem,1:this%nvar)
    ! Local
    integer :: iEl,iVar,i,j,k

    do iEl = 1,this%nElem
      do iVar = 1,this%nVar
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              det(i,j,k,iEl,iVar) = &
                this%interior(i,j,k,iEl,iVar,1,1)* &
                (this%interior(i,j,k,iEl,iVar,2,2)* &
                 this%interior(i,j,k,iEl,iVar,3,3)- &
                 this%interior(i,j,k,iEl,iVar,2,3)* &
                 this%interior(i,j,k,iEl,iVar,3,2))- &
                this%interior(i,j,k,iEl,iVar,2,1)* &
                (this%interior(i,j,k,iEl,iVar,1,2)* &
                 this%interior(i,j,k,iEl,iVar,3,3)- &
                 this%interior(i,j,k,iEl,iVar,1,3)* &
                 this%interior(i,j,k,iEl,iVar,3,2))+ &
                this%interior(i,j,k,iEl,iVar,3,1)* &
                (this%interior(i,j,k,iEl,iVar,1,2)* &
                 this%interior(i,j,k,iEl,iVar,2,3)- &
                 this%interior(i,j,k,iEl,iVar,1,3)* &
                 this%interior(i,j,k,iEl,iVar,2,2))

            enddo
          enddo
        enddo
      enddo
    enddo

  endsubroutine Determinant_Tensor3D

endmodule SELF_Tensor_3D
