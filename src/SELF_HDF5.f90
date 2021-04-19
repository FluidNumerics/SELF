MODULE SELF_HDF5

USE SELF_Constants
USE SELF_Memory
USE HDF5

#ifdef DOUBLE_PRECISION
#define HDF5_IO_PREC H5T_IEEE_F64LE
#else
#define HDF5_IO_PREC H5T_IEEE_F32LE
#endif

  INTERFACE ReadArray_HDF5
    MODULE PROCEDURE :: ReadArray_HDF5_real_r1
    MODULE PROCEDURE :: ReadArray_HDF5_real_r2
    MODULE PROCEDURE :: ReadArray_HDF5_real_r3
    MODULE PROCEDURE :: ReadArray_HDF5_real_r4
    MODULE PROCEDURE :: ReadArray_HDF5_real_r5
    MODULE PROCEDURE :: ReadArray_HDF5_real_r6
    MODULE PROCEDURE :: ReadArray_HDF5_real_r7

    MODULE PROCEDURE :: ReadArray_HDF5_int32_r1
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r2
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r3
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r4
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r5
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r6
    MODULE PROCEDURE :: ReadArray_HDF5_int32_r7

    MODULE PROCEDURE :: ReadArray_HDF5_int64_r1
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r2
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r3
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r4
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r5
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r6
    MODULE PROCEDURE :: ReadArray_HDF5_int64_r7
  END INTERFACE

CONTAINS

SUBROUTINE Open_HDF5( fileName, accessFlag, fileId, mpiComm )
  IMPLICIT NONE
  CHARACTER(*), INTENT(in) :: fileName
  INTEGER, INTENT(in) :: accessFlag
  INTEGER(HID_T), INTENT(out) :: fileId
  INTEGER, OPTIONAL, INTENT(in) :: mpiComm
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER :: error

#ifdef MPI

    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistId, error)
    CALL h5pset_fapl_mpio_f(plistId, mpiComm, MPI_INFO_NULL, error)
    CALL h5fopen_f(TRIM(fileName), accessFlag, fileId, error, plistId)
    CALL h5pclose_f(plistId, error)

#else

    CALL h5fopen_f(TRIM(fileName), accessFlag, fileId, error)

#endif

END SUBROUTINE Open_HDF5

SUBROUTINE Close_HDF5( fileId )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  ! Local
  INTEGER :: error
  
    CALL h5fclose_f( fileId, error )

END SUBROUTINE Close_HDF5

SUBROUTINE ReadAttribute_HDF5( fileId, attributeName, attribute )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: attributeName
  INTEGER, INTENT(out) :: attribute
  ! Local
  INTEGER(HID_T) :: attrId
  INTEGER(HID_T) :: typeId
  INTEGER(HSIZE_T) :: dims(1:1)
  INTEGER :: error

    dims(1) = 1
    CALL h5aopen_f(fileId, TRIM(attributeName), attrId, error)
    CALL h5aget_type_f(attrId, typeId, error)

    CALL h5aread_f(attrId, typeId, attribute, dims, error) 

    CALL h5tclose_f(typeId, error)
    CALL h5aclose_f(attrId, error) 

END SUBROUTINE ReadAttribute_HDF5

SUBROUTINE ReadArray_HDF5_real_r1( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfReal_r1), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_real_r1

SUBROUTINE ReadArray_HDF5_real_r2( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfReal_r2), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_real_r2

SUBROUTINE ReadArray_HDF5_real_r3( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfReal_r3), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_real_r3

SUBROUTINE ReadArray_HDF5_real_r4( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfReal_r4), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_real_r4

SUBROUTINE ReadArray_HDF5_real_r5( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfReal_r5), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_real_r5

SUBROUTINE ReadArray_HDF5_real_r6( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfReal_r6), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_real_r6

SUBROUTINE ReadArray_HDF5_real_r7( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfReal_r7), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_real_r7

SUBROUTINE ReadArray_HDF5_int32_r1( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt32_r1), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int32_r1

SUBROUTINE ReadArray_HDF5_int32_r2( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt32_r2), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int32_r2

SUBROUTINE ReadArray_HDF5_int32_r3( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt32_r3), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int32_r3

SUBROUTINE ReadArray_HDF5_int32_r4( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt32_r4), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int32_r4

SUBROUTINE ReadArray_HDF5_int32_r5( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt32_r5), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int32_r5

SUBROUTINE ReadArray_HDF5_int32_r6( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt32_r6), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int32_r6

SUBROUTINE ReadArray_HDF5_int32_r7( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt32_r7), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int32_r7

SUBROUTINE ReadArray_HDF5_int64_r1( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt64_r1), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int64_r1

SUBROUTINE ReadArray_HDF5_int64_r2( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt64_r2), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int64_r2

SUBROUTINE ReadArray_HDF5_int64_r3( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt64_r3), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int64_r3

SUBROUTINE ReadArray_HDF5_int64_r4( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt64_r4), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int64_r4

SUBROUTINE ReadArray_HDF5_int64_r5( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt64_r5), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int64_r5

SUBROUTINE ReadArray_HDF5_int64_r6( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt64_r6), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int64_r6

SUBROUTINE ReadArray_HDF5_int64_r7( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfInt64_r7), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
  INTEGER(HID_T) :: dtypeId
  INTEGER(HID_T) :: filespace
  INTEGER(HID_T) :: memspace
  INTEGER(HID_T) :: dims(1:3)
  INTEGER :: error

#ifdef MPI

    dims = SHAPE(hfArray % hostData)
    CALL h5screate_simple_f(3, dims, memspace, error)
    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_space_f(dsetId, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, error)
    CALL h5pset_dxpl_mpio_f(plistId, HDF5D_MPIO_COLLECTIVE_F,error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, &
      error, memspace, filespace, plistId)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5pclose_f(plistId,error)
    CALL h5sclose_f(filespace,error)
    CALL h5dclose_f(dsetId,error)
    CALL h5sclose_f(memspace,error)

#else

    CALL h5dopen_f(fileId, arrayName, dsetId, error)
    CALL h5dget_type_f(dsetId, dtypeId, error)

    CALL h5dread_f(dsetId, dtypeId, hfArray % hostData, dims, error)
    
    CALL h5tclose_f(dtypeId,error)
    CALL h5dclose_f(dsetId,error)

#endif

END SUBROUTINE ReadArray_HDF5_int64_r7

END MODULE SELF_HDF5
