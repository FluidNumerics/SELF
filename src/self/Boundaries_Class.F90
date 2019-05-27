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

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: extProcIDs_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: boundaryIDs_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: boundaryGlobalIDs_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: boundaryCondition_dev(:)
#endif

  CONTAINS

    PROCEDURE :: Build => Build_Boundaries
    PROCEDURE :: Trash => Trash_Boundaries

#ifdef HAVE_CUDA    
    PROCEDURE :: UpdateDevice => UpdateDevice_Boundaries
#endif

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
   
#ifdef HAVE_CUDA
    ALLOCATE( myBoundaries % extProcIDs_dev(1:nBe) )
    ALLOCATE( myBoundaries % boundaryIDs_dev(1:nBe) )
    ALLOCATE( myBoundaries % boundaryGlobalIDs_dev(1:nBe) )
    ALLOCATE( myBoundaries % boundaryCondition_dev(1:nBe) )
#endif


  END SUBROUTINE Build_Boundaries
!
  SUBROUTINE Trash_Boundaries( myBoundaries )

    IMPLICIT NONE
    CLASS(Boundaries), INTENT(inout) :: myBoundaries

    IF( ALLOCATED( myBoundaries % extProcIDs ) ) DEALLOCATE( myBoundaries % extProcIDs )
    IF( ALLOCATED( myBoundaries % boundaryIDs ) ) DEALLOCATE( myBoundaries % boundaryIDs )
    IF( ALLOCATED( myBoundaries % boundaryGlobalIDs ) ) DEALLOCATE( myBoundaries % boundaryGlobalIDs )
    IF( ALLOCATED( myBoundaries % boundaryCondition ) ) DEALLOCATE( myBoundaries % boundaryCondition )

#ifdef HAVE_CUDA
    IF( ALLOCATED( myBoundaries % extProcIDs_dev ) ) DEALLOCATE( myBoundaries % extProcIDs_dev )
    IF( ALLOCATED( myBoundaries % boundaryIDs_dev ) ) DEALLOCATE( myBoundaries % boundaryIDs_dev )
    IF( ALLOCATED( myBoundaries % boundaryGlobalIDs_dev ) ) DEALLOCATE( myBoundaries % boundaryGlobalIDs_dev )
    IF( ALLOCATED( myBoundaries % boundaryCondition_dev ) ) DEALLOCATE( myBoundaries % boundaryCondition_dev )
#endif

  END SUBROUTINE Trash_Boundaries

#ifdef HAVE_CUDA
  SUBROUTINE UpdateDevice_Boundaries( myBoundaries )
    CLASS( Boundaries ), INTENT(inout) :: myBoundaries

    myBoundaries % extProcIDs_dev  = myBoundaries % extProcIDs
    myBoundaries % boundaryIDs_dev = myBoundaries % boundaryIDs
    myBoundaries % boundaryGlobalIDs_dev = myBoundaries % boundaryGlobalIDs
    myBoundaries % boundaryCondition_dev = myBoundaries % boundaryCondition

  END SUBROUTINE UpdateDevice_Boundaries
#endif


END MODULE Boundaries_CLASS
