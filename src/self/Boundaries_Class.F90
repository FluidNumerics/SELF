! Boundaries_CLASS.f90
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Boundaries_CLASS

! src/COMMON/
  USE ModelPrecision

  IMPLICIT NONE



! Boundaries
! The Boundaries CLASS provides a convenient package of attributes for implementing
! boundary conditions.
!
! This structure was motivated by the need for a robust means of implementing boundary conditions
! on an unstructured mesh. This CLASS makes it trivial to implement message-passing for
! MPI parallelism.
!

  TYPE Boundaries
    INTEGER              :: nBoundaries
    INTEGER, ALLOCATABLE :: extProcIDs(:)
    INTEGER, ALLOCATABLE :: boundaryIDs(:)
    INTEGER, ALLOCATABLE :: boundaryGlobalIDs(:)
    INTEGER, ALLOCATABLE :: boundaryCondition(:)

  CONTAINS

    PROCEDURE :: Build => Build_Boundaries
    PROCEDURE :: Trash => Trash_Boundaries


  END TYPE Boundaries

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
  SUBROUTINE Build_Boundaries( myBoundaries, nBe )

    IMPLICIT NONE
    CLASS(Boundaries), INTENT(inout) :: myBoundaries
    INTEGER, INTENT(in)                        :: nBe


    myBoundaries % nBoundaries = nBe

    ALLOCATE( myBoundaries % extProcIDs(1:nBe) )
    ALLOCATE( myBoundaries % boundaryIDs(1:nBe) )
    ALLOCATE( myBoundaries % boundaryGlobalIDs(1:nBe) )
    ALLOCATE( myBoundaries % boundaryCondition(1:nBe) )

    myBoundaries % extProcIDs        = 0
    myBoundaries % boundaryIDs       = -1
    myBoundaries % boundaryGlobalIDs = -1
    myBoundaries % boundaryCondition = 0
   


  END SUBROUTINE Build_Boundaries
!
  SUBROUTINE Trash_Boundaries( myBoundaries )

    IMPLICIT NONE
    CLASS(Boundaries), INTENT(inout) :: myBoundaries

    DEALLOCATE( myBoundaries % extProcIDs, &
                myBoundaries % boundaryIDs, &
                myBoundaries % boundaryGlobalIDs, &
                myBoundaries % boundaryCondition  )

  END SUBROUTINE Trash_Boundaries


END MODULE Boundaries_CLASS
