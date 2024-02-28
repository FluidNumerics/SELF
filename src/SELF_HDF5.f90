!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_HDF5

  USE SELF_Constants
  USE ISO_FORTRAN_ENV
  USE HDF5

#ifdef DOUBLE_PRECISION
#define HDF5_IO_PREC H5T_IEEE_F64LE
#else
#define HDF5_IO_PREC H5T_IEEE_F32LE
#endif

IMPLICIT NONE

#include "SELF_Macros.h"

  INTERFACE Open_HDF5
    MODULE PROCEDURE :: Open_HDF5_serial
    MODULE PROCEDURE :: Open_HDF5_parallel
  END INTERFACE

  INTERFACE ReadAttribute_HDF5
    MODULE PROCEDURE :: ReadAttribute_HDF5_int32
    MODULE PROCEDURE :: ReadAttribute_HDF5_real
    MODULE PROCEDURE :: ReadAttribute_HDF5_character
  END INTERFACE

  INTERFACE WriteAttribute_HDF5
    MODULE PROCEDURE :: WriteAttribute_HDF5_int32
    ! MODULE PROCEDURE :: WriteAttribute_HDF5_real
    ! MODULE PROCEDURE :: WriteAttribute_HDF5_character
  END INTERFACE

  INTERFACE ReadArray_HDF5
    MODULE PROCEDURE :: ReadArray_HDF5_real_r1_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r2_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r3_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r4_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r5_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r6_serial
   ! MODULE PROCEDURE :: ReadArray_HDF5_real_r7_serial

    MODULE PROCEDURE :: ReadArray_HDF5_int32_r1_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r2_serial
  !  MODULE PROCEDURE :: ReadArray_HDF5_int32_r3_serial
  !  MODULE PROCEDURE :: ReadArray_HDF5_int32_r4_serial

    MODULE PROCEDURE :: ReadArray_HDF5_real_r1_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r2_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r3_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r4_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r5_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r6_parallel
   ! MODULE PROCEDURE :: ReadArray_HDF5_real_r7_parallel

    MODULE PROCEDURE :: ReadArray_HDF5_int32_r1_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r2_parallel
   ! MODULE PROCEDURE :: ReadArray_HDF5_int32_r3_parallel
   ! MODULE PROCEDURE :: ReadArray_HDF5_int32_r4_parallel

  END INTERFACE

  INTERFACE WriteCharacter_HDF5
    MODULE PROCEDURE :: WriteCharacter_HDF5_serial 
  END INTERFACE WriteCharacter_HDF5

  INTERFACE WriteArray_HDF5
    MODULE PROCEDURE :: WriteArray_HDF5_real_r1_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r2_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r3_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r4_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r5_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r6_serial
   ! MODULE PROCEDURE :: WriteArray_HDF5_real_r7_serial

    MODULE PROCEDURE :: WriteArray_HDF5_int32_r1_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r2_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r3_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r4_serial

    MODULE PROCEDURE :: WriteArray_HDF5_real_r1_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r2_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r3_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r4_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r5_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r6_parallel
   ! MODULE PROCEDURE :: WriteArray_HDF5_real_r7_parallel

    MODULE PROCEDURE :: WriteArray_HDF5_int32_r1_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r2_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r3_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r4_parallel

  END INTERFACE

  PRIVATE

  PUBLIC :: Open_HDF5
  PUBLIC :: Close_HDF5
  PUBLIC :: CreateGroup_HDF5
  PUBLIC :: ReadAttribute_HDF5
  PUBLIC :: WriteAttribute_HDF5
  PUBLIC :: ReadArray_HDF5
  PUBLIC :: WriteArray_HDF5
  PUBLIC :: WriteCharacter_HDF5

CONTAINS

  SUBROUTINE Open_HDF5_serial(fileName,accessFlag,fileId)
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fileName
    INTEGER,INTENT(in) :: accessFlag
    INTEGER(HID_T),INTENT(inout) :: fileId
    ! Local
    INTEGER :: error

    CALL h5open_f(error)

    IF (accessFlag == H5F_ACC_TRUNC_F) THEN
       CALL h5fcreate_f(TRIM(fileName),accessFlag,fileId,error)
    ELSE
      CALL h5fopen_f(TRIM(fileName),accessFlag,fileId,error)
    END IF

    IF (error == -1) THEN
      PRINT *, 'Failed to open '//TRIM(fileName)//'.'
      STOP -1
    END IF

  END SUBROUTINE Open_HDF5_serial

  SUBROUTINE Open_HDF5_parallel(fileName,accessFlag,fileId,mpiComm)
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fileName
    INTEGER,INTENT(in) :: accessFlag
    INTEGER(HID_T),INTENT(inout) :: fileId
    INTEGER,INTENT(in) :: mpiComm
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER :: error

    CALL h5open_f(error)

    CALL h5pcreate_f(H5P_FILE_ACCESS_F,plistId,error)
    CALL h5pset_fapl_mpio_f(plistId,mpiComm,MPI_INFO_NULL,error)

    IF (accessFlag == H5F_ACC_TRUNC_F) THEN
      CALL h5fcreate_f(TRIM(fileName),accessFlag,fileId,error,access_prp = plistId)
    ELSE
      CALL h5fopen_f(TRIM(fileName),accessFlag,fileId,error,access_prp = plistId)
    END IF
    CALL h5pclose_f(plistId,error)

    IF (error == -1) THEN
      PRINT *, 'Failed to open '//TRIM(fileName)//'.'
      STOP -1
    END IF

  END SUBROUTINE Open_HDF5_parallel

  SUBROUTINE Close_HDF5(fileId)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    ! Local
    INTEGER :: error

    CALL h5fclose_f(fileId,error)
    CALL h5close_f(error)

  END SUBROUTINE Close_HDF5

  SUBROUTINE CreateGroup_HDF5(fileId,groupName)
#undef __FUNC__
#define __FUNC__ "CreateGroup_HDF5"
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: groupName
    ! Local
    INTEGER(HID_T) :: groupId
    LOGICAL :: groupExists
    INTEGER :: error

    CALL h5lexists_f(fileId, TRIM(groupName), groupExists, error)
    IF( error /= 0 )THEN
      ERROR( "Link check failure for "//TRIM(groupName) )
    ELSE

      IF( .NOT. groupExists )THEN
        INFO("Creating group "//TRIM(groupName))
        ! Create groups
        CALL h5gcreate_f(fileId,TRIM(groupName),groupId,error)
    
        IF( error /= 0 )THEN
          ERROR( "Failed to create group "//TRIM(groupName) )
        ENDIF

        CALL h5gclose_f(groupId,error)

        IF( error /= 0 )THEN
          ERROR( "Failed to close group "//TRIM(groupName) )
        ENDIF

      ENDIF

    ENDIF

  END SUBROUTINE CreateGroup_HDF5

  SUBROUTINE ReadAttribute_HDF5_int32(fileId,attributeName,attribute)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: attributeName
    INTEGER,INTENT(out) :: attribute
    ! Local
    INTEGER(HID_T) :: attrId
    INTEGER(HID_T) :: typeId
    INTEGER(HSIZE_T) :: dims(1:1)
    INTEGER :: error

    dims(1) = 1
    CALL h5aopen_f(fileId,TRIM(attributeName),attrId,error)
    CALL h5aget_type_f(attrId,typeId,error)

    CALL h5aread_f(attrId,typeId,attribute,dims,error)

    CALL h5tclose_f(typeId,error)
    CALL h5aclose_f(attrId,error)

  END SUBROUTINE ReadAttribute_HDF5_int32

  SUBROUTINE ReadAttribute_HDF5_real(fileId,attributeName,attribute)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: attributeName
    REAL(prec),INTENT(out) :: attribute
    ! Local
    INTEGER(HID_T) :: attrId
    INTEGER(HID_T) :: typeId
    INTEGER(HSIZE_T) :: dims(1:1)
    INTEGER :: error

    dims(1) = 1
    CALL h5aopen_f(fileId,TRIM(attributeName),attrId,error)
    CALL h5aget_type_f(attrId,typeId,error)

    CALL h5aread_f(attrId,typeId,attribute,dims,error)

    CALL h5tclose_f(typeId,error)
    CALL h5aclose_f(attrId,error)

  END SUBROUTINE ReadAttribute_HDF5_real

  SUBROUTINE ReadAttribute_HDF5_character(fileId,attributeName,attribute)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: attributeName
    CHARACTER(*),INTENT(out) :: attribute
    ! Local
    INTEGER(HID_T) :: attrId
    INTEGER(HID_T) :: typeId
    INTEGER(HSIZE_T) :: dims(1:1)
    INTEGER :: error

    dims(1) = 1
    CALL h5aopen_f(fileId,TRIM(attributeName),attrId,error)
    CALL h5aget_type_f(attrId,typeId,error)

    CALL h5aread_f(attrId,typeId,attribute,dims,error)

    CALL h5tclose_f(typeId,error)
    CALL h5aclose_f(attrId,error)

  END SUBROUTINE ReadAttribute_HDF5_character

  SUBROUTINE WriteAttribute_HDF5_int32(fileId,attributeName,attribute)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: attributeName
    INTEGER,INTENT(in) :: attribute
    ! Local
    INTEGER(HID_T) :: aspaceId
    INTEGER(HID_T) :: attrId
    INTEGER(HSIZE_T) :: dims(1:1)
    INTEGER :: error

    dims(1) = 1
    CALL h5screate_f(H5S_SCALAR_F,aspaceId,error)
    CALL h5acreate_f(fileId,TRIM(attributeName),H5T_STD_I32LE, &
                     aspaceId,attrId,error)
    CALL h5awrite_f(attrId,H5T_STD_I32LE,attribute,dims,error)
    CALL h5sclose_f(aspaceId,error)
    CALL h5aclose_f(attrId,error)

  END SUBROUTINE WriteAttribute_HDF5_int32

  ! SUBROUTINE WriteAttribute_HDF5_real(fileId,attributeName,attribute)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: attributeName
  !   REAL(prec),INTENT(in) :: attribute
  !   ! Local
  !   INTEGER(HID_T) :: aspaceId
  !   INTEGER(HID_T) :: attrId
  !   INTEGER(HSIZE_T) :: dims(1:1)
  !   INTEGER :: error

  !   dims(1) = 1
  !   CALL h5screate_f(H5S_SCALAR_F,aspaceId,error)
  !   CALL h5acreate_f(fileId,TRIM(attributeName),HDF5_IO_PREC, &
  !                    aspaceId,attrId,error)
  !   CALL h5awrite_f(attrId,HDF5_IO_PREC,attribute,dims,error)
  !   CALL h5sclose_f(aspaceId,error)
  !   CALL h5aclose_f(attrId,error)

  ! END SUBROUTINE WriteAttribute_HDF5_real

  ! SUBROUTINE WriteAttribute_HDF5_character(fileId,attributeName,attribute)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: attributeName
  !   CHARACTER(*),INTENT(in) :: attribute
  !   ! Local
  !   INTEGER(HID_T) :: aspaceId
  !   INTEGER(HID_T) :: attrId
  !   INTEGER(HSIZE_T) :: dims(1:1)
  !   INTEGER :: error

  !   dims(1) = 1
  !   CALL h5screate_f(H5S_SCALAR_F,aspaceId,error)
  !   CALL h5acreate_f(fileId,TRIM(attributeName),H5T_STRING, &
  !                    aspaceId,attrId,error)
  !   CALL h5awrite_f(attrId,H5T_STRING,TRIM(attribute),dims,error)
  !   CALL h5sclose_f(aspaceId,error)
  !   CALL h5aclose_f(attrId,error)

  ! END SUBROUTINE WriteAttribute_HDF5_character

  SUBROUTINE WriteArray_HDF5_real_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r1_serial

  SUBROUTINE WriteArray_HDF5_real_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r2_serial

  SUBROUTINE WriteArray_HDF5_real_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:,:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r3_serial

  SUBROUTINE WriteArray_HDF5_real_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:,:,:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r4_serial

  SUBROUTINE WriteArray_HDF5_real_r5_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:,:,:,:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace,dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r5_serial

  SUBROUTINE WriteArray_HDF5_real_r6_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:,:,:,:,:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace,dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r6_serial

  ! SUBROUTINE WriteArray_HDF5_real_r7_serial(fileId,arrayName,hfArray)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: arrayName
  !   real(prec), dimension(:,:,:,:,:,:,:),INTENT(in) :: hfArray
  !   ! Local
  !   INTEGER(HID_T) :: dsetId
  !   INTEGER(HID_T) :: memspace
  !   INTEGER(HSIZE_T) :: dims(1:7)
  !   INTEGER :: error
  !   INTEGER :: aRank

  !   aRank = RANK(hfArray)

  !   dims = SHAPE(hfArray)
  !   CALL h5screate_simple_f(aRank,dims,memspace,error)

  !   CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
  !                    dsetId,error)

  !   CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
  !                   hfArray,dims,error)

  !   CALL h5dclose_f(dSetId,error)
  !   CALL h5sclose_f(memspace,error)

  ! END SUBROUTINE WriteArray_HDF5_real_r7_serial

  SUBROUTINE WriteArray_HDF5_int32_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    integer(int32), dimension(:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r1_serial

  SUBROUTINE WriteArray_HDF5_int32_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    integer(int32), dimension(:,:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r2_serial

  SUBROUTINE WriteArray_HDF5_int32_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    integer(int32), dimension(:,:,:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r3_serial

  SUBROUTINE WriteArray_HDF5_int32_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    integer(int32), dimension(:,:,:,:),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r4_serial

  SUBROUTINE WriteArray_HDF5_real_r1_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    real(prec), dimension(:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r1_parallel

  SUBROUTINE WriteArray_HDF5_real_r2_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:2)
    real(prec), dimension(:,:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:2)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r2_parallel

  SUBROUTINE WriteArray_HDF5_real_r3_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:3)
    real(prec), dimension(:,:,:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:3)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r3_parallel

  SUBROUTINE WriteArray_HDF5_real_r4_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:4)
    real(prec), dimension(:,:,:,:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:4)

    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r4_parallel

  SUBROUTINE WriteArray_HDF5_real_r5_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:5)
    real(prec), dimension(:,:,:,:,:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:5)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r5_parallel

  SUBROUTINE WriteArray_HDF5_real_r6_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:6)
    real(prec), dimension(:,:,:,:,:,:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:6)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r6_parallel

  ! SUBROUTINE WriteArray_HDF5_real_r7_parallel(fileId,arrayName,hfArray,offset,globalDims)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: arrayName
  !   INTEGER(HID_T),INTENT(in) :: offset(1:7)
  !   real(prec), dimension(:,:,:,:,:,:,:),INTENT(in) :: hfArray
  !   INTEGER(HID_T),INTENT(in) :: globalDims(1:7)
  !   ! Local
  !   INTEGER(HID_T) :: plistId
  !   INTEGER(HID_T) :: dsetId
  !   INTEGER(HID_T) :: filespace
  !   INTEGER(HID_T) :: memspace
  !   INTEGER(HSIZE_T) :: dims(1:7)
  !   INTEGER :: error
  !   INTEGER :: aRank

  !   aRank = RANK(hfArray)
  !   dims = SHAPE(hfArray)

  !   CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
  !   CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

  !   CALL h5screate_simple_f(aRank,globalDims,filespace,error)
  !   CALL h5screate_simple_f(aRank,dims,memspace,error)

  !   CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

  !   CALL h5sselect_hyperslab_f(filespace,&
  !                              H5S_SELECT_SET_F,&
  !                              offset,&
  !                              dims,&
  !                              error)
  !   CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

  !   IF( error /= 0 )THEN
  !     PRINT*, 'Failure to write dataset'
  !     STOP
  !   ENDIF

  !   CALL h5pclose_f(plistId,error)
  !   CALL h5sclose_f(filespace,error)
  !   CALL h5dclose_f(dSetId,error)
  !   CALL h5sclose_f(memspace,error)

  ! END SUBROUTINE WriteArray_HDF5_real_r7_parallel

  SUBROUTINE WriteArray_HDF5_int32_r1_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    integer(int32), dimension(:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r1_parallel

  SUBROUTINE WriteArray_HDF5_int32_r2_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:2)
    integer(int32), dimension(:,:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:2)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r2_parallel

  SUBROUTINE WriteArray_HDF5_int32_r3_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:3)
    integer(int32), dimension(:,:,:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:3)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r3_parallel

  SUBROUTINE WriteArray_HDF5_int32_r4_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:4)
    integer(int32), dimension(:,:,:,:),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:4)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)
    dims = SHAPE(hfArray)

    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_INDEPENDENT_F,error)

    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace,dsetId,error)

    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               dims,&
                               error)
    CALL h5dwrite_f(dsetId,HDF5_IO_PREC,hfArray,dims,error,memspace,filespace,plistId)

    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r4_parallel

  subroutine WriteCharacter_HDF5_serial(fileid, name, hfField)
    ! adapted from https://forum.hdfgroup.org/t/writing-a-string-array-as-attribute-in-fortran/8503/6
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    character (len=*), intent(in) :: name
    character (len=*), intent(in) :: hfField
    ! Local
    integer(HID_T) :: h5_strtype, h5_dspace, h5_dset
    integer(HSIZE_T), dimension(2) :: size
    character (len=len(hfField)+1), dimension(1) :: str_data
    integer(SIZE_T), dimension(1) :: str_len
    INTEGER :: error
  
    ! string output requires to open a file local = non-parallel
  
    str_len(1) = len_trim (hfField)
    size(1) = str_len(1)
    size(2) = 1
    str_data(1) = hfField//char(0)
  
    ! create data space
    call H5Tcopy_f (H5T_STRING, h5_strtype, error)
    call H5Tset_strpad_f (h5_strtype, H5T_STR_NULLPAD_F, error)
    call h5screate_simple_f (1, size(2), h5_dspace, error)
    call h5dcreate_f (fileid, trim (name), h5_strtype, h5_dspace, h5_dset, error)
    call h5dwrite_vl_f (h5_dset, h5_strtype, str_data, size, str_len, error, h5_dspace)
    call h5dclose_f (h5_dset, error)
    call h5sclose_f (h5_dspace, error)

  end subroutine WriteCharacter_HDF5_serial 

  SUBROUTINE ReadArray_HDF5_real_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error

    dims = SHAPE(hfArray)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r1_serial

  SUBROUTINE ReadArray_HDF5_real_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error

    dims = SHAPE(hfArray)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r2_serial

  SUBROUTINE ReadArray_HDF5_real_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:3)
    INTEGER :: error

    dims = SHAPE(hfArray)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r3_serial

  SUBROUTINE ReadArray_HDF5_real_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:,:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:4)
    INTEGER :: error

    dims = SHAPE(hfArray)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r4_serial

  SUBROUTINE ReadArray_HDF5_real_r5_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:,:,:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:5)
    INTEGER :: error

    dims = SHAPE(hfArray)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r5_serial

  SUBROUTINE ReadArray_HDF5_real_r6_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    real(prec), dimension(:,:,:,:,:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:6)
    INTEGER :: error

    dims = SHAPE(hfArray)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r6_serial

  ! SUBROUTINE ReadArray_HDF5_real_r7_serial(fileId,arrayName,hfArray)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: arrayName
  !   real(prec), dimension(:,:,:,:,:,:,:),INTENT(inout) :: hfArray
  !   ! Local
  !   INTEGER(HID_T) :: dsetId
  !   INTEGER(HID_T) :: dims(1:7)
  !   INTEGER :: error

  !   dims = SHAPE(hfArray)

  !   CALL h5dopen_f(fileId,arrayName,dsetId,error)

  !   CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray,dims,error)

  !   CALL h5dclose_f(dsetId,error)

  ! END SUBROUTINE ReadArray_HDF5_real_r7_serial

  SUBROUTINE ReadArray_HDF5_int32_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    integer(int32), dimension(:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error

    dims = SHAPE(hfArray)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r1_serial

  SUBROUTINE ReadArray_HDF5_int32_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    integer(int32), dimension(:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error

    dims = SHAPE(hfArray)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r2_serial

  ! SUBROUTINE ReadArray_HDF5_int32_r3_serial(fileId,arrayName,hfArray)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: arrayName
  !   integer(int32), dimension(:,:,:),INTENT(inout) :: hfArray
  !   ! Local
  !   INTEGER(HID_T) :: dsetId
  !   INTEGER(HID_T) :: dims(1:3)
  !   INTEGER :: error

  !   dims = SHAPE(hfArray)

  !   CALL h5dopen_f(fileId,arrayName,dsetId,error)

  !   CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray,dims,error)

  !   CALL h5dclose_f(dsetId,error)

  ! END SUBROUTINE ReadArray_HDF5_int32_r3_serial

  ! SUBROUTINE ReadArray_HDF5_int32_r4_serial(fileId,arrayName,hfArray)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: arrayName
  !   integer(int32), dimension(:,:,:,:),INTENT(inout) :: hfArray
  !   ! Local
  !   INTEGER(HID_T) :: dsetId
  !   INTEGER(HID_T) :: dims(1:4)
  !   INTEGER :: error

  !   dims = SHAPE(hfArray)

  !   CALL h5dopen_f(fileId,arrayName,dsetId,error)

  !   CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray,dims,error)

  !   CALL h5dclose_f(dsetId,error)

  ! END SUBROUTINE ReadArray_HDF5_int32_r4_serial

  SUBROUTINE ReadArray_HDF5_real_r1_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    real(prec), dimension(:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_real_r1_parallel

  SUBROUTINE ReadArray_HDF5_real_r2_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:2)
    real(prec), dimension(:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_real_r2_parallel

  SUBROUTINE ReadArray_HDF5_real_r3_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:3)
    real(prec), dimension(:,:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_real_r3_parallel

  SUBROUTINE ReadArray_HDF5_real_r4_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:4)
    real(prec), dimension(:,:,:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_real_r4_parallel

  SUBROUTINE ReadArray_HDF5_real_r5_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:5)
    real(prec), dimension(:,:,:,:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_real_r5_parallel

  SUBROUTINE ReadArray_HDF5_real_r6_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:6)
    real(prec), dimension(:,:,:,:,:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_real_r6_parallel

  ! SUBROUTINE ReadArray_HDF5_real_r7_parallel(fileId,arrayName,hfArray,offset)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: arrayName
  !   INTEGER(HID_T),INTENT(in) :: offset(1:7)
  !   real(prec), dimension(:,:,:,:,:,:,:),INTENT(inout) :: hfArray
  !   ! Local
  !   INTEGER(HID_T) :: plistId
  !   INTEGER(HID_T) :: dsetId
  !   INTEGER(HID_T) :: dtypeId
  !   INTEGER(HID_T) :: filespace
  !   INTEGER(HID_T) :: memspace
  !   INTEGER(HID_T) :: dims(1:7)
  !   INTEGER :: error
  !   INTEGER :: aRank

  !   aRank = RANK(hfArray)

  !   dims = SHAPE(hfArray)
  !   CALL h5screate_simple_f(aRank,dims,memspace,error)
  !   CALL h5dopen_f(fileId,arrayName,dsetId,error)
  !   CALL h5dget_space_f(dsetId,filespace,error)
  !   CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
  !   CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
  !   CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
  !   CALL h5dget_type_f(dsetId,dtypeId,error)

  !   CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
  !                  error,memspace,filespace,plistId)

  !   CALL h5tclose_f(dtypeId,error)
  !   CALL h5pclose_f(plistId,error)
  !   CALL h5sclose_f(filespace,error)
  !   CALL h5dclose_f(dsetId,error)
  !   CALL h5sclose_f(memspace,error)

  ! END SUBROUTINE ReadArray_HDF5_real_r7_parallel

  SUBROUTINE ReadArray_HDF5_int32_r1_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    integer(int32), dimension(:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(1,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int32_r1_parallel

  SUBROUTINE ReadArray_HDF5_int32_r2_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:2)
    integer(int32), dimension(:,:),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray)

    dims = SHAPE(hfArray)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int32_r2_parallel

  ! SUBROUTINE ReadArray_HDF5_int32_r3_parallel(fileId,arrayName,hfArray,offset)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: arrayName
  !   INTEGER(HID_T),INTENT(in) :: offset(1:3)
  !   integer(int32), dimension(:,:,:),INTENT(inout) :: hfArray
  !   ! Local
  !   INTEGER(HID_T) :: plistId
  !   INTEGER(HID_T) :: dsetId
  !   INTEGER(HID_T) :: dtypeId
  !   INTEGER(HID_T) :: filespace
  !   INTEGER(HID_T) :: memspace
  !   INTEGER(HID_T) :: dims(1:3)
  !   INTEGER :: error
  !   INTEGER :: aRank

  !   aRank = RANK(hfArray)

  !   dims = SHAPE(hfArray)
  !   CALL h5screate_simple_f(aRank,dims,memspace,error)
  !   CALL h5dopen_f(fileId,arrayName,dsetId,error)
  !   CALL h5dget_space_f(dsetId,filespace,error)
  !   CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
  !   CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
  !   CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
  !   CALL h5dget_type_f(dsetId,dtypeId,error)

  !   CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
  !                  error,memspace,filespace,plistId)

  !   CALL h5tclose_f(dtypeId,error)
  !   CALL h5pclose_f(plistId,error)
  !   CALL h5sclose_f(filespace,error)
  !   CALL h5dclose_f(dsetId,error)
  !   CALL h5sclose_f(memspace,error)

  ! END SUBROUTINE ReadArray_HDF5_int32_r3_parallel

  ! SUBROUTINE ReadArray_HDF5_int32_r4_parallel(fileId,arrayName,hfArray,offset)
  !   IMPLICIT NONE
  !   INTEGER(HID_T),INTENT(in) :: fileId
  !   CHARACTER(*),INTENT(in) :: arrayName
  !   INTEGER(HID_T),INTENT(in) :: offset(1:4)
  !   integer(int32), dimension(:,:,:,:),INTENT(inout) :: hfArray
  !   ! Local
  !   INTEGER(HID_T) :: plistId
  !   INTEGER(HID_T) :: dsetId
  !   INTEGER(HID_T) :: dtypeId
  !   INTEGER(HID_T) :: filespace
  !   INTEGER(HID_T) :: memspace
  !   INTEGER(HID_T) :: dims(1:4)
  !   INTEGER :: error
  !   INTEGER :: aRank

  !   aRank = RANK(hfArray)

  !   dims = SHAPE(hfArray)
  !   CALL h5screate_simple_f(aRank,dims,memspace,error)
  !   CALL h5dopen_f(fileId,arrayName,dsetId,error)
  !   CALL h5dget_space_f(dsetId,filespace,error)
  !   CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
  !   CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
  !   CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
  !   CALL h5dget_type_f(dsetId,dtypeId,error)

  !   CALL h5dread_f(dsetId,dtypeId,hfArray,dims, &
  !                  error,memspace,filespace,plistId)

  !   CALL h5tclose_f(dtypeId,error)
  !   CALL h5pclose_f(plistId,error)
  !   CALL h5sclose_f(filespace,error)
  !   CALL h5dclose_f(dsetId,error)
  !   CALL h5sclose_f(memspace,error)

  ! END SUBROUTINE ReadArray_HDF5_int32_r4_parallel

END MODULE SELF_HDF5
