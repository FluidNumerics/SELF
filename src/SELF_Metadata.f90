! SELF_Metadata.F90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Metadata

  USE SELF_HDF5
  USE HDF5

  INTEGER,PARAMETER,PUBLIC :: SELF_MTD_NameLength = 250
  INTEGER,PARAMETER,PUBLIC :: SELF_MTD_DescriptionLength = 1000
  INTEGER,PARAMETER,PUBLIC :: SELF_MTD_UnitsLength = 20

  ! A class for storing metadata information, intended for file IO
  TYPE Metadata
    CHARACTER(SELF_MTD_NameLength) :: name
    CHARACTER(SELF_MTD_DescriptionLength) :: description
    CHARACTER(SELF_MTD_UnitsLength) :: units

  CONTAINS

    PROCEDURE,PUBLIC :: SetName => SetName_Metadata
    PROCEDURE,PUBLIC :: SetDescription => SetDescription_Metadata
    PROCEDURE,PUBLIC :: SetUnits => SetUnits_Metadata
    PROCEDURE,PUBLIC :: WriteHDF5 => WriteHDF5_Metadata

  END TYPE Metadata

CONTAINS

  SUBROUTINE SetName_Metadata(mtd,name)
    IMPLICIT NONE
    CLASS(Metadata),INTENT(inout) :: mtd
    CHARACTER(*),INTENT(in) :: name

    mtd % name = name

  END SUBROUTINE SetName_Metadata

  SUBROUTINE SetDescription_Metadata(mtd,description)
    IMPLICIT NONE
    CLASS(Metadata),INTENT(inout) :: mtd
    CHARACTER(*),INTENT(in) :: description

    mtd % description = description

  END SUBROUTINE SetDescription_Metadata

  SUBROUTINE SetUnits_Metadata(mtd,units)
    IMPLICIT NONE
    CLASS(Metadata),INTENT(inout) :: mtd
    CHARACTER(*),INTENT(in) :: units

    mtd % units = units

  END SUBROUTINE SetUnits_Metadata

  SUBROUTINE WriteHDF5_Metadata(mtd,group,fileId)
  !! Writes the metadata to a HDF5 file using the 
  !! fields :
  !!  * `/metadata/{group}/{name}/`
  !!  * `/metadata/{group}/{name}/description`
  !!  * `/metadata/{group}/{name}/units`
  !!
  !! This method assumes that an HDF5 file is already
  !! open for writing and is associated with the `fileId` 
  !! input.
    CLASS(Metadata), INTENT(in) :: mtd
    CHARACTER(*), INTENT(in) :: group
    INTEGER(HID_T), INTENT(in) :: fileId
    ! Local
    CHARACTER(4) :: varNumber

    ! Add variable names to the file
    CALL CreateGroup_HDF5(fileId,TRIM(group)//"/metadata")
    CALL CreateGroup_HDF5(fileId,TRIM(group)//"/metadata/"//TRIM(mtd % name))
    ! CALL CreateGroup_HDF5(fileId,TRIM(group)//"/metadata/"//TRIM(mtd % name)//"/description")
    ! CALL CreateGroup_HDF5(fileId,TRIM(group)//"/metadata/"//TRIM(mtd % name)//"/units")
  
    WRITE (varNumber,"(I0)") varid
    CALL WriteCharacter_HDF5(fileId, TRIM(group)//"/metadata/name/"//TRIM(varnumber), &
                                   TRIM(mtd % name))
    CALL WriteCharacter_HDF5(fileId,TRIM(group)//"/metadata/"//TRIM(mtd % name)//"/description", &
                                TRIM(mtd % description))
    CALL WriteCharacter_HDF5(fileId, TRIM(group)//"/metadata/"//TRIM(mtd % name)//"/units", &
                                TRIM(mtd % units))
  END SUBROUTINE WriteHDF5_Metadata

END MODULE SELF_Metadata
