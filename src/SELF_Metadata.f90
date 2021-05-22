! SELF_Metadata.F90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE SELF_Metadata

  USE SELF_HDF5

  INTEGER,PARAMETER,PUBLIC :: SELF_MTD_NameLength = 250
  INTEGER,PARAMETER,PUBLIC :: SELF_MTD_DescriptionLength = 1000
  INTEGER,PARAMETER,PUBLIC :: SELF_MTD_UnitsLength = 20

  ! A class for storing metadata information, intended for file IO
  TYPE Metadata
    CHARACTER(SELF_MTD_NameLength) :: name
    CHARACTER(SELF_MTD_DescriptionLength) :: description
    CHARACTER(SELF_MTD_UnitsLength) :: units

    CONTAINS

    PROCEDURE,PUBLIC :: Set_Name => Set_Name_Metadata
    PROCEDURE,PUBLIC :: Set_Description => Set_Description_Metadata
    PROCEDURE,PUBLIC :: Set_Units => Set_Unit_Metadata
    PROCEDURE,PUBLIC :: Write_HDF5 => Write_HDF5_Metadata

  END TYPE Metadata


CONTAINS

SUBROUTINE Set_Name_Metadata(mtd, name)
  IMPLICIT NONE
  CLASS(Metadata), INTENT(inout) :: mtd
  CHARACTER(*), INTENT(in) :: name

    mtd % name = name

SUBROUTINE Set_Name_Metadata
  
SUBROUTINE Set_Description_Metadata(mtd, description)
  IMPLICIT NONE
  CLASS(Metadata), INTENT(inout) :: mtd
  CHARACTER(*), INTENT(in) :: description

    mtd % description = description

SUBROUTINE Set_Description_Metadata

SUBROUTINE Set_Units_Metadata(mtd,  units)
  IMPLICIT NONE
  CLASS(Metadata), INTENT(inout) :: mtd
  CHARACTER(*), INTENT(in) :: units

    mtd % units = units

SUBROUTINE Set_Units_Metadata

SUBROUTINE Write_HDF5_Metadata(mtd, name)
  IMPLICIT NONE
  CLASS(Metadata), INTENT(inout) :: mtd
  CHARACTER(*), INTENT(in) :: name

    mtd % name = name

SUBROUTINE Write_HDF5_Metadata

END MODULE SELF_Metadata