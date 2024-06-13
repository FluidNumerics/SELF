! SELF_Metadata.F90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Metadata

  use SELF_HDF5
  use HDF5

  integer,parameter,public :: SELF_MTD_NameLength = 250
  integer,parameter,public :: SELF_MTD_DescriptionLength = 1000
  integer,parameter,public :: SELF_MTD_UnitsLength = 20

  ! A class for storing metadata information, intended for file IO
  type Metadata
    character(SELF_MTD_NameLength) :: name
    character(SELF_MTD_DescriptionLength) :: description
    character(SELF_MTD_UnitsLength) :: units

  contains

    procedure,public :: SetName => SetName_Metadata
    procedure,public :: SetDescription => SetDescription_Metadata
    procedure,public :: SetUnits => SetUnits_Metadata
    procedure,public :: WriteHDF5 => WriteHDF5_Metadata

  endtype Metadata

contains

  subroutine SetName_Metadata(mtd,name)
    implicit none
    class(Metadata),intent(inout) :: mtd
    character(*),intent(in) :: name

    mtd%name = name

  endsubroutine SetName_Metadata

  subroutine SetDescription_Metadata(mtd,description)
    implicit none
    class(Metadata),intent(inout) :: mtd
    character(*),intent(in) :: description

    mtd%description = description

  endsubroutine SetDescription_Metadata

  subroutine SetUnits_Metadata(mtd,units)
    implicit none
    class(Metadata),intent(inout) :: mtd
    character(*),intent(in) :: units

    mtd%units = units

  endsubroutine SetUnits_Metadata

  subroutine WriteHDF5_Metadata(mtd,group,varid,fileId)
  !! Writes the metadata to a HDF5 file using the
  !! fields :
  !!  * `/metadata/{group}/name/{varid}`
  !!  * `/metadata/{group}/description/{varid}`
  !!  * `/metadata/{group}/units/{varid}`
  !!
  !! This method assumes that an HDF5 file is already
  !! open for writing and is associated with the `fileId`
  !! input.
    class(Metadata),intent(in) :: mtd
    character(*),intent(in) :: group
    integer,intent(in) :: varid
    integer(HID_T),intent(in) :: fileId
    ! Local
    character(4) :: varNumber

    ! Add variable names to the file
    call CreateGroup_HDF5(fileId,trim(group)//"/metadata")
    call CreateGroup_HDF5(fileId,trim(group)//"/metadata/name")
    call CreateGroup_HDF5(fileId,trim(group)//"/metadata/description")
    call CreateGroup_HDF5(fileId,trim(group)//"/metadata/units")

    write(varNumber,"(I0)") varid
    call WriteCharacter_HDF5(fileId,trim(group)//"/metadata/name/"//trim(varnumber), &
                             trim(mtd%name))
    call WriteCharacter_HDF5(fileId,trim(group)//"/metadata/description/"//trim(varnumber), &
                             trim(mtd%description))
    call WriteCharacter_HDF5(fileId,trim(group)//"/metadata/units/"//trim(varnumber), &
                             trim(mtd%units))
  endsubroutine WriteHDF5_Metadata

endmodule SELF_Metadata
