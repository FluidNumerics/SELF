! SELF_Data.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Data

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5

  use HDF5
  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,public :: SELF_DataObj
  !! The SELF_DataObj class is a base class for all data objects in SELF.
  !! A data object in SELF is a multidimensional array of data, represented
  !! on both host and device, that is associated with an interpolant, metadata,
  !! and (optionally) an equation string.
  !! Type extensions of the SELF_DataObj include scalars, vectors, and tensors
  !! in 1-D, 2-D, and 3-D using the storage patterns that are expected for
  !! derivative and interpolation operations defined in SELF_Lagrange.f90
  !! Additionally, each extended type has the necessary attributes to store
  !! information on element interiors and element boundaries, both of which
  !! are commonly used for spectral element solvers.

    integer :: nVar
    integer :: nElem
    type(Lagrange),pointer :: interp
    type(Metadata),allocatable :: meta(:)
    type(EquationParser),allocatable :: eqn(:)

  contains

    ! Procedures for setting metadata for
    procedure,public :: SetName => SetName_DataObj
    procedure,public :: SetDescription => SetDescription_DataObj
    procedure,public :: SetUnits => SetUnits_DataObj
    generic,public :: SetEquation => SetEquation_DataObj
    procedure,private :: SetEquation_DataObj

  endtype SELF_DataObj

  integer,parameter :: selfStrongForm = 0
  integer,parameter :: selfWeakDGForm = 1
  integer,parameter :: selfWeakCGForm = 2
  integer,parameter :: selfWeakBRForm = 3

contains

! -- DataObj -- !

  subroutine SetName_DataObj(this,ivar,name)
    !! Set the name of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: name

    call this%meta(ivar)%SetName(name)

  endsubroutine SetName_DataObj

  subroutine SetDescription_DataObj(this,ivar,description)
    !! Set the description of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: description

    call this%meta(ivar)%SetDescription(description)

  endsubroutine SetDescription_DataObj

  subroutine SetUnits_DataObj(this,ivar,units)
    !! Set the units of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: units

    call this%meta(ivar)%SetUnits(units)

  endsubroutine SetUnits_DataObj

  subroutine SetEquation_DataObj(this,ivar,eqnChar)
    !! Sets the equation parser for the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: eqnChar

    this%eqn(ivar) = EquationParser(trim(eqnChar), &
                                    (/'x','y','z','t'/))

  endsubroutine SetEquation_DataObj

endmodule SELF_Data
