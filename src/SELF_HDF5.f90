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
    MODULE PROCEDURE :: ReadArray_HDF5_r3
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

SUBROUTINE ReadArray_HDF5_r3( fileId, arrayName, offset, hfArray )
  IMPLICIT NONE
  INTEGER(HID_T), INTENT(in) :: fileId
  CHARACTER(*), INTENT(in) :: arrayName
  INTEGER(HID_T), INTENT(in) :: offset
  TYPE(hfReal_r3), INTENT(inout) :: hfArray
  ! Local
  INTEGER(HID_T) :: plistId
  INTEGER(HID_T) :: dsetId
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
    CALL h5dget_type

#else

#endif

END SUBROUTINE ReadArray_HDF5_r3

END MODULE SELF_HDF5
