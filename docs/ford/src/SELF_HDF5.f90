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

module SELF_HDF5

  use SELF_Constants
  use iso_fortran_env
  use HDF5
  use mpi

#ifdef DOUBLE_PRECISION
#define HDF5_IO_PREC H5T_IEEE_F64LE
#else
#define HDF5_IO_PREC H5T_IEEE_F32LE
#endif

  implicit none

  interface Open_HDF5
    module procedure :: Open_HDF5_serial
    module procedure :: Open_HDF5_parallel
  endinterface

  interface ReadAttribute_HDF5
    module procedure :: ReadAttribute_HDF5_int32
    module procedure :: ReadAttribute_HDF5_real
    module procedure :: ReadAttribute_HDF5_character
  endinterface

  interface WriteAttribute_HDF5
    module procedure :: WriteAttribute_HDF5_int32
  endinterface

  interface ReadArray_HDF5
    module procedure :: ReadArray_HDF5_real_r1_serial
    module procedure :: ReadArray_HDF5_real_r2_serial
    module procedure :: ReadArray_HDF5_real_r3_serial
    module procedure :: ReadArray_HDF5_real_r4_serial
    module procedure :: ReadArray_HDF5_real_r5_serial
    ! module procedure :: ReadArray_HDF5_real_r6_serial

    module procedure :: ReadArray_HDF5_int32_r1_serial
    module procedure :: ReadArray_HDF5_int32_r2_serial

    module procedure :: ReadArray_HDF5_real_r1_parallel
    module procedure :: ReadArray_HDF5_real_r2_parallel
    module procedure :: ReadArray_HDF5_real_r3_parallel
    module procedure :: ReadArray_HDF5_real_r4_parallel
    ! module procedure :: ReadArray_HDF5_real_r5_parallel
    ! module procedure :: ReadArray_HDF5_real_r6_parallel

    module procedure :: ReadArray_HDF5_int32_r1_parallel
    module procedure :: ReadArray_HDF5_int32_r2_parallel

  endinterface

  interface WriteCharacter_HDF5
    module procedure :: WriteCharacter_HDF5_serial
  endinterface WriteCharacter_HDF5

  interface WriteArray_HDF5
    module procedure :: WriteArray_HDF5_real_r1_serial
    module procedure :: WriteArray_HDF5_real_r2_serial
    module procedure :: WriteArray_HDF5_real_r3_serial
    module procedure :: WriteArray_HDF5_real_r4_serial
    module procedure :: WriteArray_HDF5_real_r5_serial
    ! module procedure :: WriteArray_HDF5_real_r6_serial

    module procedure :: WriteArray_HDF5_int32_r1_serial
    module procedure :: WriteArray_HDF5_int32_r2_serial
    module procedure :: WriteArray_HDF5_int32_r3_serial
    module procedure :: WriteArray_HDF5_int32_r4_serial

    module procedure :: WriteArray_HDF5_real_r3_parallel
    module procedure :: WriteArray_HDF5_real_r4_parallel

    !module procedure :: WriteArray_HDF5_int32_r3_parallel
    !module procedure :: WriteArray_HDF5_int32_r4_parallel

  endinterface

  private

  public :: Open_HDF5
  public :: Close_HDF5
  public :: CreateGroup_HDF5
  public :: ReadAttribute_HDF5
  public :: WriteAttribute_HDF5
  public :: ReadArray_HDF5
  public :: WriteArray_HDF5
  public :: WriteCharacter_HDF5

contains

  subroutine Open_HDF5_serial(fileName,accessFlag,fileId)
    implicit none
    character(*),intent(in) :: fileName
    integer,intent(in) :: accessFlag
    integer(HID_T),intent(inout) :: fileId
    ! Local
    integer :: error

    call h5open_f(error)

    if(accessFlag == H5F_ACC_TRUNC_F) then
      call h5fcreate_f(trim(fileName),accessFlag,fileId,error)
    else
      call h5fopen_f(trim(fileName),accessFlag,fileId,error)
    endif

    if(error == -1) then
      print*,'Failed to open '//trim(fileName)//'.'
      stop 1
    endif

  endsubroutine Open_HDF5_serial

  subroutine Open_HDF5_parallel(fileName,accessFlag,fileId,mpiComm)
    implicit none
    character(*),intent(in) :: fileName
    integer,intent(in) :: accessFlag
    integer(HID_T),intent(inout) :: fileId
    integer,intent(in) :: mpiComm
    ! Local
    integer(HID_T) :: plistId
    integer :: error

    call h5open_f(error)

    call h5pcreate_f(H5P_FILE_ACCESS_F,plistId,error)
    call h5pset_fapl_mpio_f(plistId,mpiComm,MPI_INFO_NULL,error)

    if(accessFlag == H5F_ACC_TRUNC_F) then
      call h5fcreate_f(trim(fileName),accessFlag,fileId,error,access_prp=plistId)
    else
      call h5fopen_f(trim(fileName),accessFlag,fileId,error,access_prp=plistId)
    endif
    call h5pclose_f(plistId,error)

    if(error == -1) then
      print*,'Failed to open '//trim(fileName)//'.'
      stop 1
    endif

  endsubroutine Open_HDF5_parallel

  subroutine Close_HDF5(fileId)
    implicit none
    integer(HID_T),intent(in) :: fileId
    ! Local
    integer :: error

    call h5fclose_f(fileId,error)
    call h5close_f(error)

  endsubroutine Close_HDF5

  subroutine CreateGroup_HDF5(fileId,groupName)
#undef __FUNC__
#define __FUNC__ "CreateGroup_HDF5"
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: groupName
    ! Local
    integer(HID_T) :: groupId
    logical :: groupExists
    integer :: error

    call h5lexists_f(fileId,trim(groupName),groupExists,error)
    if(error /= 0) then
      print*,__FILE__," : Link check failure for "//trim(groupName)
    else

      if(.not. groupExists) then
        ! Create groups
        call h5gcreate_f(fileId,trim(groupName),groupId,error)

        if(error /= 0) then
          print*,__FILE__," :Failed to create group "//trim(groupName)
        endif

        call h5gclose_f(groupId,error)

        if(error /= 0) then
          print*,__FILE__," :Failed to close group "//trim(groupName)
        endif

      endif

    endif

  endsubroutine CreateGroup_HDF5

  subroutine ReadAttribute_HDF5_int32(fileId,attributeName,attribute)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: attributeName
    integer,intent(out) :: attribute
    ! Local
    integer(HID_T) :: attrId
    integer(HID_T) :: typeId
    integer(HSIZE_T) :: dims(1:1)
    integer :: error

    dims(1) = 1
    call h5aopen_f(fileId,trim(attributeName),attrId,error)
    call h5aget_type_f(attrId,typeId,error)

    call h5aread_f(attrId,typeId,attribute,dims,error)

    call h5tclose_f(typeId,error)
    call h5aclose_f(attrId,error)

  endsubroutine ReadAttribute_HDF5_int32

  subroutine ReadAttribute_HDF5_real(fileId,attributeName,attribute)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: attributeName
    real(prec),intent(out) :: attribute
    ! Local
    integer(HID_T) :: attrId
    integer(HID_T) :: typeId
    integer(HSIZE_T) :: dims(1:1)
    integer :: error

    dims(1) = 1
    call h5aopen_f(fileId,trim(attributeName),attrId,error)
    call h5aget_type_f(attrId,typeId,error)

    call h5aread_f(attrId,typeId,attribute,dims,error)

    call h5tclose_f(typeId,error)
    call h5aclose_f(attrId,error)

  endsubroutine ReadAttribute_HDF5_real

  subroutine ReadAttribute_HDF5_character(fileId,attributeName,attribute)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: attributeName
    character(*),intent(out) :: attribute
    ! Local
    integer(HID_T) :: attrId
    integer(HID_T) :: typeId
    integer(HSIZE_T) :: dims(1:1)
    integer :: error

    dims(1) = 1
    call h5aopen_f(fileId,trim(attributeName),attrId,error)
    call h5aget_type_f(attrId,typeId,error)

    call h5aread_f(attrId,typeId,attribute,dims,error)

    call h5tclose_f(typeId,error)
    call h5aclose_f(attrId,error)

  endsubroutine ReadAttribute_HDF5_character

  subroutine WriteAttribute_HDF5_int32(fileId,attributeName,attribute)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: attributeName
    integer,intent(in) :: attribute
    ! Local
    integer(HID_T) :: aspaceId
    integer(HID_T) :: attrId
    integer(HSIZE_T) :: dims(1:1)
    integer :: error

    dims(1) = 1
    call h5screate_f(H5S_SCALAR_F,aspaceId,error)
    call h5acreate_f(fileId,trim(attributeName),H5T_STD_I32LE, &
                     aspaceId,attrId,error)
    call h5awrite_f(attrId,H5T_STD_I32LE,attribute,dims,error)
    call h5sclose_f(aspaceId,error)
    call h5aclose_f(attrId,error)

  endsubroutine WriteAttribute_HDF5_int32

  subroutine WriteArray_HDF5_real_r1_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(1,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    call h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_real_r1_serial

  subroutine WriteArray_HDF5_real_r2_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:,:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:2)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(2,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    call h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_real_r2_serial

  subroutine WriteArray_HDF5_real_r3_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:,:,:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:3)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(3,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    call h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_real_r3_serial

  subroutine WriteArray_HDF5_real_r4_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:,:,:,:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:4)
    integer :: error

    dims = shape(hfArray)

    call h5screate_simple_f(4,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    call h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_real_r4_serial

  subroutine WriteArray_HDF5_real_r5_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:,:,:,:,:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:5)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(5,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,memspace,dsetId,error)

    call h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_real_r5_serial

  ! subroutine WriteArray_HDF5_real_r6_serial(fileId,arrayName,hfArray)
  !   implicit none
  !   integer(HID_T),intent(in) :: fileId
  !   character(*),intent(in) :: arrayName
  !   real(prec),dimension(:,:,:,:,:,:),intent(in) :: hfArray
  !   ! Local
  !   integer(HID_T) :: dsetId
  !   integer(HID_T) :: memspace
  !   integer(HSIZE_T) :: dims(1:6)
  !   integer :: error

  !   dims = shape(hfArray)
  !   call h5screate_simple_f(6,dims,memspace,error)

  !   call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,memspace,dsetId,error)

  !   call h5dwrite_f(dsetId,HDF5_IO_PREC, &
  !                   hfArray,dims,error)

  !   call h5dclose_f(dSetId,error)
  !   call h5sclose_f(memspace,error)

  ! endsubroutine WriteArray_HDF5_real_r6_serial

  subroutine WriteArray_HDF5_int32_r1_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(int32),dimension(:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(1,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    call h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_int32_r1_serial

  subroutine WriteArray_HDF5_int32_r2_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(int32),dimension(:,:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:2)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(2,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    call h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_int32_r2_serial

  subroutine WriteArray_HDF5_int32_r3_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(int32),dimension(:,:,:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:3)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(3,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    call h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_int32_r3_serial

  subroutine WriteArray_HDF5_int32_r4_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(int32),dimension(:,:,:,:),intent(in) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:4)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(4,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    call h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray,dims,error)

    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_int32_r4_serial

  subroutine WriteArray_HDF5_real_r3_parallel(fileId,arrayName,hfArray,offset,globalDims)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(HID_T),intent(in) :: offset(1:3)
    real(prec),dimension(:,:,:),intent(in) :: hfArray
    integer(HID_T),intent(in) :: globalDims(1:3)
    ! Local
    integer(HID_T) :: plistId
    integer(HID_T) :: dsetId
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:3)
    integer :: error

    dims = shape(hfArray)

    call h5screate_simple_f(3,globalDims,filespace,error)
    call h5screate_simple_f(3,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    call h5sselect_hyperslab_f(filespace, &
                               H5S_SELECT_SET_F, &
                               offset, &
                               dims, &
                               error)

    call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    call h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error, &
                    mem_space_id=memspace,file_space_id=filespace,xfer_prp=plistId)

    if(error /= 0) then
      print*,'Failure to write dataset'
      stop 1
    endif

    call h5pclose_f(plistId,error)
    call h5sclose_f(filespace,error)
    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine WriteArray_HDF5_real_r3_parallel

  subroutine WriteArray_HDF5_real_r4_parallel(fileId,arrayName,hfArray,offset,globalDims)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(HID_T),intent(in) :: offset(1:4)
    real(prec),dimension(:,:,:,:),intent(in) :: hfArray
    integer(HID_T),intent(in) :: globalDims(1:4)

    ! Local
    integer(HID_T) :: plistId
    integer(HID_T) :: dsetId
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(1:4)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(4,globalDims,filespace,error)
    call h5screate_simple_f(4,dims,memspace,error)

    call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    call h5sselect_hyperslab_f(filespace, &
                               H5S_SELECT_SET_F, &
                               offset, &
                               dims, &
                               error)

    call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    call h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error, &
                    mem_space_id=memspace,file_space_id=filespace,xfer_prp=plistId)

    if(error /= 0) then
      print*,'Failure to write dataset'
      stop 1
    endif

    call h5pclose_f(plistId,error)
    call h5sclose_f(filespace,error)
    call h5dclose_f(dSetId,error)
    call h5sclose_f(memspace,error)
  endsubroutine WriteArray_HDF5_real_r4_parallel

  ! subroutine WriteArray_HDF5_int32_r3_parallel(fileId,arrayName,hfArray,offset,globalDims)
  !   implicit none
  !   integer(HID_T),intent(in) :: fileId
  !   character(*),intent(in) :: arrayName
  !   integer(HID_T),intent(in) :: offset(1:3)
  !   integer(int32),dimension(:,:,:),intent(in) :: hfArray
  !   integer(HID_T),intent(in) :: globalDims(1:3)
  !   ! Local
  !   integer(HID_T) :: plistId
  !   integer(HID_T) :: dsetId
  !   integer(HID_T) :: filespace
  !   integer(HID_T) :: memspace
  !   integer(HSIZE_T) :: dims(1:3)
  !   integer :: error

  !   dims = shape(hfArray)

  !   call h5screate_simple_f(3,globalDims,filespace,error)
  !   call h5screate_simple_f(3,dims,memspace,error)

  !   call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

  !   call h5sselect_hyperslab_f(filespace, &
  !                              H5S_SELECT_SET_F, &
  !                              offset, &
  !                              dims, &
  !                              error)

  !   call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
  !   call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
  !   call h5dwrite_f(dsetId,H5T_STD_I32LE,hfArray,dims,error, &
  !                   mem_space_id=memspace,file_space_id=filespace,xfer_prp=plistId)

  !   if(error /= 0) then
  !     print*,'Failure to write dataset'
  !     stop 1
  !   endif

  !   call h5pclose_f(plistId,error)
  !   call h5sclose_f(filespace,error)
  !   call h5dclose_f(dSetId,error)
  !   call h5sclose_f(memspace,error)

  ! endsubroutine WriteArray_HDF5_int32_r3_parallel

  ! subroutine WriteArray_HDF5_int32_r4_parallel(fileId,arrayName,hfArray,offset,globalDims)
  !   implicit none
  !   integer(HID_T),intent(in) :: fileId
  !   character(*),intent(in) :: arrayName
  !   integer(HID_T),intent(in) :: offset(1:4)
  !   integer(int32),dimension(:,:,:,:),intent(in) :: hfArray
  !   integer(HID_T),intent(in) :: globalDims(1:4)
  !   ! Local
  !   integer(HID_T) :: plistId
  !   integer(HID_T) :: dsetId
  !   integer(HID_T) :: filespace
  !   integer(HID_T) :: memspace
  !   integer(HSIZE_T) :: dims(1:4)
  !   integer :: error

  !   dims = shape(hfArray)

  !   call h5screate_simple_f(4,globalDims,filespace,error)
  !   call h5screate_simple_f(4,dims,memspace,error)

  !   call h5dcreate_f(fileId,trim(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

  !   call h5sselect_hyperslab_f(filespace, &
  !                              H5S_SELECT_SET_F, &
  !                              offset, &
  !                              dims, &
  !                              error)

  !   call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
  !   call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
  !   call h5dwrite_f(dsetId,H5T_STD_I32LE,hfArray,dims,error, &
  !                   mem_space_id=memspace,file_space_id=filespace,xfer_prp=plistId)

  !   if(error /= 0) then
  !     print*,'Failure to write dataset'
  !     stop 1
  !   endif

  !   call h5pclose_f(plistId,error)
  !   call h5sclose_f(filespace,error)
  !   call h5dclose_f(dSetId,error)
  !   call h5sclose_f(memspace,error)

  ! endsubroutine WriteArray_HDF5_int32_r4_parallel

  subroutine WriteCharacter_HDF5_serial(fileid,name,hfField)
    ! adapted from https://forum.hdfgroup.org/t/writing-a-string-array-as-attribute-in-fortran/8503/6
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(len=*),intent(in) :: name
    character(len=*),intent(in) :: hfField
    ! Local
    integer(HID_T) :: h5_strtype,h5_dspace,h5_dset
    integer(HSIZE_T),dimension(2) :: size
    character(len=len(hfField)+1),dimension(1) :: str_data
    integer(SIZE_T),dimension(1) :: str_len
    integer :: error

    ! string output requires to open a file local = non-parallel

    str_len(1) = len_trim(hfField)
    size(1) = str_len(1)
    size(2) = 1
    str_data(1) = hfField//char(0)

    ! create data space
    call H5Tcopy_f(H5T_STRING,h5_strtype,error)
    call H5Tset_strpad_f(h5_strtype,H5T_STR_NULLPAD_F,error)
    call h5screate_simple_f(1,size(2),h5_dspace,error)
    call h5dcreate_f(fileid,trim(name),h5_strtype,h5_dspace,h5_dset,error)
    call h5dwrite_vl_f(h5_dset,h5_strtype,str_data,size,str_len,error,h5_dspace)
    call h5dclose_f(h5_dset,error)
    call h5sclose_f(h5_dspace,error)

  endsubroutine WriteCharacter_HDF5_serial

  subroutine ReadArray_HDF5_real_r1_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: dims(1)
    integer :: error

    dims = shape(hfArray)

    call h5dopen_f(fileId,arrayName,dsetId,error)

    call h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    call h5dclose_f(dsetId,error)

  endsubroutine ReadArray_HDF5_real_r1_serial

  subroutine ReadArray_HDF5_real_r2_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: dims(1:2)
    integer :: error

    dims = shape(hfArray)

    call h5dopen_f(fileId,arrayName,dsetId,error)

    call h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    call h5dclose_f(dsetId,error)

  endsubroutine ReadArray_HDF5_real_r2_serial

  subroutine ReadArray_HDF5_real_r3_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:,:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: dims(1:3)
    integer :: error

    dims = shape(hfArray)

    call h5dopen_f(fileId,arrayName,dsetId,error)

    call h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    call h5dclose_f(dsetId,error)

  endsubroutine ReadArray_HDF5_real_r3_serial

  subroutine ReadArray_HDF5_real_r4_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:,:,:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: dims(1:4)
    integer :: error

    dims = shape(hfArray)

    call h5dopen_f(fileId,arrayName,dsetId,error)

    call h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    call h5dclose_f(dsetId,error)

  endsubroutine ReadArray_HDF5_real_r4_serial

  subroutine ReadArray_HDF5_real_r5_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    real(prec),dimension(:,:,:,:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: dims(1:5)
    integer :: error

    dims = shape(hfArray)

    call h5dopen_f(fileId,arrayName,dsetId,error)

    call h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    call h5dclose_f(dsetId,error)

  endsubroutine ReadArray_HDF5_real_r5_serial

  ! subroutine ReadArray_HDF5_real_r6_serial(fileId,arrayName,hfArray)
  !   implicit none
  !   integer(HID_T),intent(in) :: fileId
  !   character(*),intent(in) :: arrayName
  !   real(prec),dimension(:,:,:,:,:,:),intent(inout) :: hfArray
  !   ! Local
  !   integer(HID_T) :: dsetId
  !   integer(HID_T) :: dims(1:6)
  !   integer :: error

  !   dims = shape(hfArray)

  !   call h5dopen_f(fileId,arrayName,dsetId,error)

  !   call h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

  !   call h5dclose_f(dsetId,error)

  ! endsubroutine ReadArray_HDF5_real_r6_serial

  subroutine ReadArray_HDF5_int32_r1_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(int32),dimension(:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: dims(1)
    integer :: error

    dims = shape(hfArray)

    call h5dopen_f(fileId,arrayName,dsetId,error)

    call h5dread_f(dsetId,h5kind_to_type(int32,H5_INTEGER_KIND),hfArray,dims,error)

    call h5dclose_f(dsetId,error)

  endsubroutine ReadArray_HDF5_int32_r1_serial

  subroutine ReadArray_HDF5_int32_r2_serial(fileId,arrayName,hfArray)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(int32),dimension(:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: dsetId
    integer(HID_T) :: dims(1:2)
    integer :: error

    dims = shape(hfArray)

    call h5dopen_f(fileId,arrayName,dsetId,error)

    call h5dread_f(dsetId,h5kind_to_type(int32,H5_INTEGER_KIND),hfArray,dims,error)

    call h5dclose_f(dsetId,error)

  endsubroutine ReadArray_HDF5_int32_r2_serial

  subroutine ReadArray_HDF5_real_r1_parallel(fileId,arrayName,hfArray,offset)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(HID_T),intent(in) :: offset(1)
    real(prec),dimension(:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: plistId
    integer(HID_T) :: dsetId
    integer(HID_T) :: dtypeId
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HID_T) :: dims(1)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(1,dims,memspace,error)
    call h5dopen_f(fileId,arrayName,dsetId,error)
    call h5dget_space_f(dsetId,filespace,error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    call h5dget_type_f(dsetId,dtypeId,error)

    call h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    call h5tclose_f(dtypeId,error)
    call h5pclose_f(plistId,error)
    call h5sclose_f(filespace,error)
    call h5dclose_f(dsetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine ReadArray_HDF5_real_r1_parallel

  subroutine ReadArray_HDF5_real_r2_parallel(fileId,arrayName,hfArray,offset)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(HID_T),intent(in) :: offset(1:2)
    real(prec),dimension(:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: plistId
    integer(HID_T) :: dsetId
    integer(HID_T) :: dtypeId
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HID_T) :: dims(1:2)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(2,dims,memspace,error)
    call h5dopen_f(fileId,arrayName,dsetId,error)
    call h5dget_space_f(dsetId,filespace,error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    call h5dget_type_f(dsetId,dtypeId,error)

    call h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    call h5tclose_f(dtypeId,error)
    call h5pclose_f(plistId,error)
    call h5sclose_f(filespace,error)
    call h5dclose_f(dsetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine ReadArray_HDF5_real_r2_parallel

  subroutine ReadArray_HDF5_real_r3_parallel(fileId,arrayName,hfArray,offset)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(HID_T),intent(in) :: offset(1:3)
    real(prec),dimension(:,:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: plistId
    integer(HID_T) :: dsetId
    integer(HID_T) :: dtypeId
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HID_T) :: dims(1:3)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(3,dims,memspace,error)
    call h5dopen_f(fileId,arrayName,dsetId,error)
    call h5dget_space_f(dsetId,filespace,error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    call h5dget_type_f(dsetId,dtypeId,error)

    call h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    call h5tclose_f(dtypeId,error)
    call h5pclose_f(plistId,error)
    call h5sclose_f(filespace,error)
    call h5dclose_f(dsetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine ReadArray_HDF5_real_r3_parallel

  subroutine ReadArray_HDF5_real_r4_parallel(fileId,arrayName,hfArray,offset)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(HID_T),intent(in) :: offset(1:4)
    real(prec),dimension(:,:,:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: plistId
    integer(HID_T) :: dsetId
    integer(HID_T) :: dtypeId
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HID_T) :: dims(1:4)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(4,dims,memspace,error)
    call h5dopen_f(fileId,arrayName,dsetId,error)
    call h5dget_space_f(dsetId,filespace,error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    call h5dget_type_f(dsetId,dtypeId,error)

    call h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    call h5tclose_f(dtypeId,error)
    call h5pclose_f(plistId,error)
    call h5sclose_f(filespace,error)
    call h5dclose_f(dsetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine ReadArray_HDF5_real_r4_parallel

  ! subroutine ReadArray_HDF5_real_r5_parallel(fileId,arrayName,hfArray,offset)
  !   implicit none
  !   integer(HID_T),intent(in) :: fileId
  !   character(*),intent(in) :: arrayName
  !   integer(HID_T),intent(in) :: offset(1:5)
  !   real(prec),dimension(:,:,:,:,:),intent(inout) :: hfArray
  !   ! Local
  !   integer(HID_T) :: plistId
  !   integer(HID_T) :: dsetId
  !   integer(HID_T) :: dtypeId
  !   integer(HID_T) :: filespace
  !   integer(HID_T) :: memspace
  !   integer(HID_T) :: dims(1:5)
  !   integer :: error

  !   dims = shape(hfArray)
  !   call h5screate_simple_f(5,dims,memspace,error)
  !   call h5dopen_f(fileId,arrayName,dsetId,error)
  !   call h5dget_space_f(dsetId,filespace,error)
  !   call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
  !   call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
  !   call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
  !   call h5dget_type_f(dsetId,dtypeId,error)

  !   call h5dread_f(dsetId,dtypeId,hfArray,dims, &
  !                  error,memspace,filespace,plistId)

  !   call h5tclose_f(dtypeId,error)
  !   call h5pclose_f(plistId,error)
  !   call h5sclose_f(filespace,error)
  !   call h5dclose_f(dsetId,error)
  !   call h5sclose_f(memspace,error)

  ! endsubroutine ReadArray_HDF5_real_r5_parallel

  ! subroutine ReadArray_HDF5_real_r6_parallel(fileId,arrayName,hfArray,offset)
  !   implicit none
  !   integer(HID_T),intent(in) :: fileId
  !   character(*),intent(in) :: arrayName
  !   integer(HID_T),intent(in) :: offset(1:6)
  !   real(prec),dimension(:,:,:,:,:,:),intent(inout) :: hfArray
  !   ! Local
  !   integer(HID_T) :: plistId
  !   integer(HID_T) :: dsetId
  !   integer(HID_T) :: dtypeId
  !   integer(HID_T) :: filespace
  !   integer(HID_T) :: memspace
  !   integer(HID_T) :: dims(1:6)
  !   integer :: error

  !   dims = shape(hfArray)
  !   call h5screate_simple_f(6,dims,memspace,error)
  !   call h5dopen_f(fileId,arrayName,dsetId,error)
  !   call h5dget_space_f(dsetId,filespace,error)
  !   call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
  !   call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
  !   call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
  !   call h5dget_type_f(dsetId,dtypeId,error)

  !   call h5dread_f(dsetId,dtypeId,hfArray,dims, &
  !                  error,memspace,filespace,plistId)

  !   call h5tclose_f(dtypeId,error)
  !   call h5pclose_f(plistId,error)
  !   call h5sclose_f(filespace,error)
  !   call h5dclose_f(dsetId,error)
  !   call h5sclose_f(memspace,error)

  ! endsubroutine ReadArray_HDF5_real_r6_parallel

  subroutine ReadArray_HDF5_int32_r1_parallel(fileId,arrayName,hfArray,offset)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(HID_T),intent(in) :: offset(1)
    integer(int32),dimension(:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: plistId
    integer(HID_T) :: dsetId
    integer(HID_T) :: dtypeId
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HID_T) :: dims(1)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(1,dims,memspace,error)
    call h5dopen_f(fileId,arrayName,dsetId,error)
    call h5dget_space_f(dsetId,filespace,error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    call h5dget_type_f(dsetId,dtypeId,error)

    call h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    call h5tclose_f(dtypeId,error)
    call h5pclose_f(plistId,error)
    call h5sclose_f(filespace,error)
    call h5dclose_f(dsetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine ReadArray_HDF5_int32_r1_parallel

  subroutine ReadArray_HDF5_int32_r2_parallel(fileId,arrayName,hfArray,offset)
    implicit none
    integer(HID_T),intent(in) :: fileId
    character(*),intent(in) :: arrayName
    integer(HID_T),intent(in) :: offset(1:2)
    integer(int32),dimension(:,:),intent(inout) :: hfArray
    ! Local
    integer(HID_T) :: plistId
    integer(HID_T) :: dsetId
    integer(HID_T) :: dtypeId
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HID_T) :: dims(1:2)
    integer :: error

    dims = shape(hfArray)
    call h5screate_simple_f(2,dims,memspace,error)
    call h5dopen_f(fileId,arrayName,dsetId,error)
    call h5dget_space_f(dsetId,filespace,error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    call h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    call h5dget_type_f(dsetId,dtypeId,error)

    call h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    call h5tclose_f(dtypeId,error)
    call h5pclose_f(plistId,error)
    call h5sclose_f(filespace,error)
    call h5dclose_f(dsetId,error)
    call h5sclose_f(memspace,error)

  endsubroutine ReadArray_HDF5_int32_r2_parallel

endmodule SELF_HDF5
