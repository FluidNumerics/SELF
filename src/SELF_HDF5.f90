!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_HDF5

  USE SELF_Constants
  USE SELF_Memory
  USE ISO_FORTRAN_ENV
  USE HDF5

#ifdef DOUBLE_PRECISION
#define HDF5_IO_PREC H5T_IEEE_F64LE
#else
#define HDF5_IO_PREC H5T_IEEE_F32LE
#endif

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
    MODULE PROCEDURE :: WriteAttribute_HDF5_real
    MODULE PROCEDURE :: WriteAttribute_HDF5_character
  END INTERFACE

  INTERFACE ReadArray_HDF5
    MODULE PROCEDURE :: ReadArray_HDF5_real_r1_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r2_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r3_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r4_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r5_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r6_serial
    MODULE PROCEDURE :: ReadArray_HDF5_real_r7_serial

    MODULE PROCEDURE :: ReadArray_HDF5_int32_r1_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r2_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r3_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r4_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r5_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r6_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r7_serial

    MODULE PROCEDURE :: ReadArray_HDF5_int64_r1_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r2_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r3_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r4_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r5_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r6_serial
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r7_serial

    MODULE PROCEDURE :: ReadArray_HDF5_real_r1_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r2_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r3_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r4_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r5_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r6_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_real_r7_parallel

    MODULE PROCEDURE :: ReadArray_HDF5_int32_r1_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r2_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r3_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r4_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r5_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r6_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r7_parallel

    MODULE PROCEDURE :: ReadArray_HDF5_int64_r1_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r2_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r3_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r4_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r5_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r6_parallel
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r7_parallel
  END INTERFACE

  INTERFACE WriteArray_HDF5
    MODULE PROCEDURE :: WriteArray_HDF5_real_r1_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r2_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r3_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r4_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r5_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r6_serial
    MODULE PROCEDURE :: WriteArray_HDF5_real_r7_serial

    MODULE PROCEDURE :: WriteArray_HDF5_int32_r1_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r2_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r3_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r4_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r5_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r6_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r7_serial

    MODULE PROCEDURE :: WriteArray_HDF5_int64_r1_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r2_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r3_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r4_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r5_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r6_serial
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r7_serial

    MODULE PROCEDURE :: WriteArray_HDF5_real_r1_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r2_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r3_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r4_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r5_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r6_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_real_r7_parallel

    MODULE PROCEDURE :: WriteArray_HDF5_int32_r1_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r2_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r3_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r4_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r5_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r6_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int32_r7_parallel

    MODULE PROCEDURE :: WriteArray_HDF5_int64_r1_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r2_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r3_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r4_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r5_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r6_parallel
    MODULE PROCEDURE :: WriteArray_HDF5_int64_r7_parallel
  END INTERFACE

  PRIVATE

  PUBLIC :: Open_HDF5
  PUBLIC :: Close_HDF5
  PUBLIC :: CreateGroup_HDF5
  PUBLIC :: ReadAttribute_HDF5
  PUBLIC :: WriteAttribute_HDF5
  PUBLIC :: ReadArray_HDF5
  PUBLIC :: WriteArray_HDF5

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
      STOP - 1
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
      STOP - 1
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
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: groupName
    ! Local
    INTEGER(HID_T) :: groupId
    INTEGER :: error

    ! Create groups
    CALL h5gcreate_f(fileId,TRIM(groupName),groupId,error)
    CALL h5gclose_f(groupId,error)

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

  SUBROUTINE WriteAttribute_HDF5_real(fileId,attributeName,attribute)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: attributeName
    REAL(prec),INTENT(in) :: attribute
    ! Local
    INTEGER(HID_T) :: aspaceId
    INTEGER(HID_T) :: attrId
    INTEGER(HSIZE_T) :: dims(1:1)
    INTEGER :: error

    dims(1) = 1
    CALL h5screate_f(H5S_SCALAR_F,aspaceId,error)
    CALL h5acreate_f(fileId,TRIM(attributeName),HDF5_IO_PREC, &
                     aspaceId,attrId,error)
    CALL h5awrite_f(attrId,HDF5_IO_PREC,attribute,dims,error)
    CALL h5sclose_f(aspaceId,error)
    CALL h5aclose_f(attrId,error)

  END SUBROUTINE WriteAttribute_HDF5_real

  SUBROUTINE WriteAttribute_HDF5_character(fileId,attributeName,attribute)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: attributeName
    CHARACTER(*),INTENT(in) :: attribute
    ! Local
    INTEGER(HID_T) :: aspaceId
    INTEGER(HID_T) :: attrId
    INTEGER(HSIZE_T) :: dims(1:1)
    INTEGER :: error

    dims(1) = 1
    CALL h5screate_f(H5S_SCALAR_F,aspaceId,error)
    CALL h5acreate_f(fileId,TRIM(attributeName),H5T_STRING, &
                     aspaceId,attrId,error)
    CALL h5awrite_f(attrId,H5T_STRING,TRIM(attribute),dims,error)
    CALL h5sclose_f(aspaceId,error)
    CALL h5aclose_f(attrId,error)

  END SUBROUTINE WriteAttribute_HDF5_character

  SUBROUTINE WriteArray_HDF5_real_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r1),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r1_serial

  SUBROUTINE WriteArray_HDF5_real_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r2),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r2_serial

  SUBROUTINE WriteArray_HDF5_real_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r3),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r3_serial

  SUBROUTINE WriteArray_HDF5_real_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r4),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r4_serial

  SUBROUTINE WriteArray_HDF5_real_r5_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r5),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace,dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r5_serial

  SUBROUTINE WriteArray_HDF5_real_r6_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r6),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace,dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r6_serial

  SUBROUTINE WriteArray_HDF5_real_r7_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r7),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r7_serial

  SUBROUTINE WriteArray_HDF5_int32_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint32_r1),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r1_serial

  SUBROUTINE WriteArray_HDF5_int32_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint32_r2),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r2_serial

  SUBROUTINE WriteArray_HDF5_int32_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint32_r3),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r3_serial

  SUBROUTINE WriteArray_HDF5_int32_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint32_r4),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r4_serial

  SUBROUTINE WriteArray_HDF5_int32_r5_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint32_r5),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r5_serial

  SUBROUTINE WriteArray_HDF5_int32_r6_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint32_r6),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r6_serial

  SUBROUTINE WriteArray_HDF5_int32_r7_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint32_r7),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r7_serial

  SUBROUTINE WriteArray_HDF5_int64_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint64_r1),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r1_serial

  SUBROUTINE WriteArray_HDF5_int64_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint64_r2),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r2_serial

  SUBROUTINE WriteArray_HDF5_int64_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint64_r3),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r3_serial

  SUBROUTINE WriteArray_HDF5_int64_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint64_r4),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r4_serial

  SUBROUTINE WriteArray_HDF5_int64_r5_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint64_r5),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r5_serial

  SUBROUTINE WriteArray_HDF5_int64_r6_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint64_r6),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r6_serial

  SUBROUTINE WriteArray_HDF5_int64_r7_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfint64_r7),INTENT(in) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,memspace, &
                     dsetId,error)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error)

    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r7_serial

  SUBROUTINE WriteArray_HDF5_real_r1_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    TYPE(hfReal_r1),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER(HSIZE_T) :: strides(1)
    INTEGER(HSIZE_T) :: counts(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace, &
                     dsetId,error,plistId)

    strides = (/1/)
    counts = (/1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error,memspace,filespace)
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
    TYPE(hfReal_r2),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:2)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER(HSIZE_T) :: strides(1:2)
    INTEGER(HSIZE_T) :: counts(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1/)
    counts = (/1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error,memspace,filespace)
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
    TYPE(hfReal_r3),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:3)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER(HSIZE_T) :: strides(1:3)
    INTEGER(HSIZE_T) :: counts(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1/)
    counts = (/1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error,memspace,filespace)
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
    TYPE(hfReal_r4),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:4)

    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER(HSIZE_T) :: strides(1:4)
    INTEGER(HSIZE_T) :: counts(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1/)
    counts = (/1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r4_parallel

  SUBROUTINE WriteArray_HDF5_real_r5_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:5)
    TYPE(hfReal_r5),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:5)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:5)
    INTEGER(HSIZE_T) :: strides(1:5)
    INTEGER(HSIZE_T) :: counts(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1/)
    counts = (/1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace, &
                               H5S_SELECT_SET_F, &
                               offset, &
                               counts, &
                               error, &
                               strides, &
                               dims)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error,memspace,filespace)
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
    TYPE(hfReal_r6),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:6)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:6)
    INTEGER(HSIZE_T) :: strides(1:6)
    INTEGER(HSIZE_T) :: counts(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1,1/)
    counts = (/1,1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r6_parallel

  SUBROUTINE WriteArray_HDF5_real_r7_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:7)
    TYPE(hfReal_r7),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:7)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:7)
    INTEGER(HSIZE_T) :: strides(1:7)
    INTEGER(HSIZE_T) :: counts(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),HDF5_IO_PREC,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1,1,1/)
    counts = (/1,1,1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,HDF5_IO_PREC, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_real_r7_parallel

  SUBROUTINE WriteArray_HDF5_int32_r1_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    TYPE(hfint32_r1),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER(HSIZE_T) :: strides(1)
    INTEGER(HSIZE_T) :: counts(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1/)
    counts = (/1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
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
    TYPE(hfint32_r2),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:2)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER(HSIZE_T) :: strides(1:2)
    INTEGER(HSIZE_T) :: counts(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1/)
    counts = (/1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
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
    TYPE(hfint32_r3),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:3)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER(HSIZE_T) :: strides(1:3)
    INTEGER(HSIZE_T) :: counts(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1/)
    counts = (/1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
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
    TYPE(hfint32_r4),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:4)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER(HSIZE_T) :: strides(1:4)
    INTEGER(HSIZE_T) :: counts(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1/)
    counts = (/1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r4_parallel

  SUBROUTINE WriteArray_HDF5_int32_r5_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:5)
    TYPE(hfint32_r5),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:5)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:5)
    INTEGER(HSIZE_T) :: strides(1:5)
    INTEGER(HSIZE_T) :: counts(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1/)
    counts = (/1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r5_parallel

  SUBROUTINE WriteArray_HDF5_int32_r6_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:6)
    TYPE(hfint32_r6),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:6)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:6)
    INTEGER(HSIZE_T) :: strides(1:6)
    INTEGER(HSIZE_T) :: counts(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1,1/)
    counts = (/1,1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r6_parallel

  SUBROUTINE WriteArray_HDF5_int32_r7_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:7)
    TYPE(hfint32_r7),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:7)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:7)
    INTEGER(HSIZE_T) :: strides(1:7)
    INTEGER(HSIZE_T) :: counts(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I32LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1,1,1/)
    counts = (/1,1,1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I32LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int32_r7_parallel

  SUBROUTINE WriteArray_HDF5_int64_r1_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    TYPE(hfint64_r1),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER(HSIZE_T) :: strides(1)
    INTEGER(HSIZE_T) :: counts(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1/)
    counts = (/1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r1_parallel

  SUBROUTINE WriteArray_HDF5_int64_r2_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:2)
    TYPE(hfint64_r2),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:2)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:2)
    INTEGER(HSIZE_T) :: strides(1:2)
    INTEGER(HSIZE_T) :: counts(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1/)
    counts = (/1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r2_parallel

  SUBROUTINE WriteArray_HDF5_int64_r3_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:3)
    TYPE(hfint64_r3),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:3)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:3)
    INTEGER(HSIZE_T) :: strides(1:3)
    INTEGER(HSIZE_T) :: counts(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1/)
    counts = (/1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r3_parallel

  SUBROUTINE WriteArray_HDF5_int64_r4_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:4)
    TYPE(hfint64_r4),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:4)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:4)
    INTEGER(HSIZE_T) :: strides(1:4)
    INTEGER(HSIZE_T) :: counts(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1/)
    counts = (/1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r4_parallel

  SUBROUTINE WriteArray_HDF5_int64_r5_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:5)
    TYPE(hfint64_r5),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:5)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:5)
    INTEGER(HSIZE_T) :: strides(1:5)
    INTEGER(HSIZE_T) :: counts(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1/)
    counts = (/1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r5_parallel

  SUBROUTINE WriteArray_HDF5_int64_r6_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:6)
    TYPE(hfint64_r6),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:6)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:6)
    INTEGER(HSIZE_T) :: strides(1:6)
    INTEGER(HSIZE_T) :: counts(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1,1/)
    counts = (/1,1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r6_parallel

  SUBROUTINE WriteArray_HDF5_int64_r7_parallel(fileId,arrayName,hfArray,offset,globalDims)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:7)
    TYPE(hfint64_r7),INTENT(in) :: hfArray
    INTEGER(HID_T),INTENT(in) :: globalDims(1:7)
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HSIZE_T) :: dims(1:7)
    INTEGER(HSIZE_T) :: strides(1:7)
    INTEGER(HSIZE_T) :: counts(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,globalDims,filespace,error)
    CALL h5screate_simple_f(aRank,dims,memspace,error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F,plistId,error)
    CALL h5pset_chunk_f(plistId,aRank,dims,error)

    CALL h5dcreate_f(fileId,TRIM(arrayName),H5T_STD_I64LE,filespace, &
                     dsetId,error,plistId)

    strides = (/1,1,1,1,1,1,1/)
    counts = (/1,1,1,1,1,1,1/)
    CALL h5sselect_hyperslab_f(filespace,&
                               H5S_SELECT_SET_F,&
                               offset,&
                               counts,&
                               error,&
                               strides,&
                               dims)

    CALL h5dwrite_f(dsetId,H5T_STD_I64LE, &
                    hfArray % hostData,dims,error,memspace,filespace)
    IF( error /= 0 )THEN
      PRINT*, 'Failure to write dataset'
      STOP
    ENDIF

    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dSetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE WriteArray_HDF5_int64_r7_parallel

  SUBROUTINE ReadArray_HDF5_real_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r1),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r1_serial

  SUBROUTINE ReadArray_HDF5_real_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r2),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r2_serial

  SUBROUTINE ReadArray_HDF5_real_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r3),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:3)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r3_serial

  SUBROUTINE ReadArray_HDF5_real_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r4),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:4)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r4_serial

  SUBROUTINE ReadArray_HDF5_real_r5_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r5),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:5)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r5_serial

  SUBROUTINE ReadArray_HDF5_real_r6_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r6),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:6)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r6_serial

  SUBROUTINE ReadArray_HDF5_real_r7_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfReal_r7),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:7)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(prec,H5_REAL_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_real_r7_serial

  SUBROUTINE ReadArray_HDF5_int32_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt32_r1),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r1_serial

  SUBROUTINE ReadArray_HDF5_int32_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt32_r2),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r2_serial

  SUBROUTINE ReadArray_HDF5_int32_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt32_r3),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:3)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r3_serial

  SUBROUTINE ReadArray_HDF5_int32_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt32_r4),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:4)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r4_serial

  SUBROUTINE ReadArray_HDF5_int32_r5_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt32_r5),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:5)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r5_serial

  SUBROUTINE ReadArray_HDF5_int32_r6_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt32_r6),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:6)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r6_serial

  SUBROUTINE ReadArray_HDF5_int32_r7_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt32_r7),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:7)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT32,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int32_r7_serial

  SUBROUTINE ReadArray_HDF5_int64_r1_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt64_r1),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT64,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int64_r1_serial

  SUBROUTINE ReadArray_HDF5_int64_r2_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt64_r2),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT64,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int64_r2_serial

  SUBROUTINE ReadArray_HDF5_int64_r3_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt64_r3),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:3)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT64,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int64_r3_serial

  SUBROUTINE ReadArray_HDF5_int64_r4_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt64_r4),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: dims(1:4)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT64,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int64_r4_serial

  SUBROUTINE ReadArray_HDF5_int64_r5_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt64_r5),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:5)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT64,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int64_r5_serial

  SUBROUTINE ReadArray_HDF5_int64_r6_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt64_r6),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:6)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT64,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int64_r6_serial

  SUBROUTINE ReadArray_HDF5_int64_r7_serial(fileId,arrayName,hfArray)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    TYPE(hfInt64_r7),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dims(1:7)
    INTEGER :: error

    dims = SHAPE(hfArray % hostData)

    CALL h5dopen_f(fileId,arrayName,dsetId,error)

    CALL h5dread_f(dsetId,h5kind_to_type(INT64,H5_INTEGER_KIND),hfArray % hostData,dims,error)

    CALL h5dclose_f(dsetId,error)

  END SUBROUTINE ReadArray_HDF5_int64_r7_serial

  SUBROUTINE ReadArray_HDF5_real_r1_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    TYPE(hfReal_r1),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
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
    TYPE(hfReal_r2),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
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
    TYPE(hfReal_r3),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
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
    TYPE(hfReal_r4),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
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
    TYPE(hfReal_r5),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
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
    TYPE(hfReal_r6),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_real_r6_parallel

  SUBROUTINE ReadArray_HDF5_real_r7_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:7)
    TYPE(hfReal_r7),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_real_r7_parallel

  SUBROUTINE ReadArray_HDF5_int32_r1_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    TYPE(hfInt32_r1),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(1,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
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
    TYPE(hfInt32_r2),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int32_r2_parallel

  SUBROUTINE ReadArray_HDF5_int32_r3_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:3)
    TYPE(hfInt32_r3),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int32_r3_parallel

  SUBROUTINE ReadArray_HDF5_int32_r4_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:4)
    TYPE(hfInt32_r4),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int32_r4_parallel

  SUBROUTINE ReadArray_HDF5_int32_r5_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:5)
    TYPE(hfInt32_r5),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int32_r5_parallel

  SUBROUTINE ReadArray_HDF5_int32_r6_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:6)
    TYPE(hfInt32_r6),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int32_r6_parallel

  SUBROUTINE ReadArray_HDF5_int32_r7_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:7)
    TYPE(hfInt32_r7),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int32_r7_parallel

  SUBROUTINE ReadArray_HDF5_int64_r1_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1)
    TYPE(hfInt64_r1),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int64_r1_parallel

  SUBROUTINE ReadArray_HDF5_int64_r2_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:2)
    TYPE(hfInt64_r2),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:2)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int64_r2_parallel

  SUBROUTINE ReadArray_HDF5_int64_r3_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:3)
    TYPE(hfInt64_r3),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:3)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int64_r3_parallel

  SUBROUTINE ReadArray_HDF5_int64_r4_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:4)
    TYPE(hfInt64_r4),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:4)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int64_r4_parallel

  SUBROUTINE ReadArray_HDF5_int64_r5_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:5)
    TYPE(hfInt64_r5),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:5)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int64_r5_parallel

  SUBROUTINE ReadArray_HDF5_int64_r6_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:6)
    TYPE(hfInt64_r6),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:6)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int64_r6_parallel

  SUBROUTINE ReadArray_HDF5_int64_r7_parallel(fileId,arrayName,hfArray,offset)
    IMPLICIT NONE
    INTEGER(HID_T),INTENT(in) :: fileId
    CHARACTER(*),INTENT(in) :: arrayName
    INTEGER(HID_T),INTENT(in) :: offset(1:7)
    TYPE(hfInt64_r7),INTENT(inout) :: hfArray
    ! Local
    INTEGER(HID_T) :: plistId
    INTEGER(HID_T) :: dsetId
    INTEGER(HID_T) :: dtypeId
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: dims(1:7)
    INTEGER :: error
    INTEGER :: aRank

    aRank = RANK(hfArray % hostData)

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(aRank,dims,memspace,error)
    CALL h5dopen_f(fileId,arrayName,dsetId,error)
    CALL h5dget_space_f(dsetId,filespace,error)
    CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,dims,error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F,plistId,error)
    CALL h5pset_dxpl_mpio_f(plistId,H5FD_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId,dtypeId,error)

    CALL h5dread_f(dsetId,dtypeId,hfArray % hostData,dims, &
                   error,memspace,filespace,plistId)

    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

  END SUBROUTINE ReadArray_HDF5_int64_r7_parallel

END MODULE SELF_HDF5
