! BoundaryCommunicator_CLASS.f90
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE BoundaryCommunicator_CLASS

! src/COMMON/
  USE ModelPrecision
  USE ConstantsDictionary
  USE CommonRoutines
#ifdef HAVE_CUDA
  USE cudafor
#endif


  IMPLICIT NONE

#ifdef HAVE_MPI
  INCLUDE 'mpif.h'
#endif


! BoundaryCommunicator
! The BoundaryCommunicator CLASS provides a convenient package of attributes for implementing
! boundary conditions.
!
! This structure was motivated by the need for a robust means of implementing boundary conditions
! on an unstructured mesh. This CLASS makes it trivial to implement message-passing for
! MPI parallelism.
!

  TYPE BoundaryCommunicator
    LOGICAL              :: setup
    INTEGER              :: nBoundaries, myRank, nProc
    INTEGER, ALLOCATABLE :: extProcIDs(:)
    INTEGER, ALLOCATABLE :: boundaryIDs(:)
    INTEGER, ALLOCATABLE :: boundaryGlobalIDs(:)
    INTEGER, ALLOCATABLE :: unPackMap(:)

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: extProcIDs_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: boundaryIDs_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: boundaryGlobalIDs_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: unPackMap_dev(:)
#endif

#ifdef HAVE_MPI
    INTEGER              :: nNeighbors,  maxBufferSize
    INTEGER              :: MPI_COMM, MPI_PREC, mpiErr
    INTEGER, ALLOCATABLE :: bufferMap(:)
    INTEGER, ALLOCATABLE :: neighborRank(:)
    INTEGER, ALLOCATABLE :: bufferSize(:)
    INTEGER, ALLOCATABLE :: rankTable(:)

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: bufferMap_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: neighborRank_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: bufferSize_dev(:)
    INTEGER, DEVICE, ALLOCATABLE :: rankTable_dev(:)
#endif

#endif

  CONTAINS

    PROCEDURE :: Build => Build_BoundaryCommunicator
    PROCEDURE :: Trash => Trash_BoundaryCommunicator
    PROCEDURE :: Finalize => Finalize_BoundaryCommunicator
    PROCEDURE :: SetRanks

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateDevice => UpdateDevice_BoundaryCommunicator
    PROCEDURE :: UpdateHost   => UpdateHost_BoundaryCommunicator
#endif

    PROCEDURE :: ReadPickup  => ReadPickup_BoundaryCommunicator
    PROCEDURE :: WritePickup => WritePickup_BoundaryCommunicator

#ifdef HAVE_MPI
    PROCEDURE :: ConstructCommTables
#endif

  END TYPE BoundaryCommunicator

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup BoundaryCommunicator_CLASS
!! @{
! ================================================================================================ !
! S/R Build
!
!> \fn Build_BoundaryCommunicator
!! Allocates space for the BoundaryCommunicator structure and initializes all array values to zero.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(BoundaryCommunicator) :: this <BR>
!! <B>INTEGER</B>                    :: nBe <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( nBe ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> out <th> myComm <td> BoundaryCommunicator <td>
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree for the solution storage
!!   <tr> <td> in <th> nBe <td> INTEGER <td> The number of boundary edges in the mesh
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE Build_BoundaryCommunicator( myComm, nBe )

    IMPLICIT NONE
    CLASS(BoundaryCommunicator), INTENT(inout) :: myComm
    INTEGER, INTENT(in)                        :: nBe


    myComm % nBoundaries = nBe

    ALLOCATE( myComm % extProcIDs(1:nBe) )
    ALLOCATE( myComm % boundaryIDs(1:nBe) )
    ALLOCATE( myComm % boundaryGlobalIDs(1:nBe) )
    ALLOCATE( myComm % unPackMap(1:nBe) )

    myComm % extProcIDs  = 0
    myComm % boundaryIDs = -1
    myComm % unPackMap   = 0

#ifdef HAVE_CUDA

    ALLOCATE( myComm % extProcIDs_dev(1:nBe) )
    ALLOCATE( myComm % boundaryIDs_dev(1:nBe) )
    ALLOCATE( myComm % boundaryGlobalIDs_dev(1:nBe) )
    ALLOCATE( myComm % unPackMap_dev(1:nBe) )


#endif


#ifdef HAVE_MPI
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, myComm % nProc, myComm % mpiErr )

    ALLOCATE( myComm % bufferMap(1:myComm % nBoundaries), &
              myComm % rankTable(0:myComm % nProc-1) )

#ifdef HAVE_CUDA

    ALLOCATE( myComm % bufferMap_dev(1:myComm % nBoundaries), &
              myComm % rankTable_dev(0:myComm % nProc-1) )

#endif

#endif

myComm % setup = .TRUE.

  END SUBROUTINE Build_BoundaryCommunicator
!
!> \addtogroup BoundaryCommunicator_CLASS
!! @{
! ================================================================================================ !
! S/R Trash
!
!> \fn Trash_BoundaryCommunicator
!! Frees memory associated with the attributes of the BoundaryCommunicator DATA structure.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(BoundaryCommunicator) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myComm <td> BoundaryCommunicator <td>
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE Trash_BoundaryCommunicator( myComm )

    IMPLICIT NONE
    CLASS(BoundaryCommunicator), INTENT(inout) :: myComm

    DEALLOCATE( myComm % unPackMap, myComm % extProcIDs, myComm % boundaryIDs, myComm % boundaryGlobalIDs  )

#ifdef HAVE_CUDA
    DEALLOCATE( myComm % unPackMap_dev, myComm % extProcIDs_dev, myComm % boundaryIDs_dev, myComm % boundaryGlobalIDs_dev  )
#endif

#ifdef HAVE_MPI

    IF( ALLOCATED( myComm % neighborRank ) )DEALLOCATE( myComm % neighborRank )
    IF( ALLOCATED( myComm % bufferSize ) )DEALLOCATE( myComm % bufferSize )
    IF( ALLOCATED( myComm % bufferMap ) )DEALLOCATE( myComm % bufferMap )
    IF( ALLOCATED( myComm % rankTable ) )DEALLOCATE( myComm % rankTable )

#ifdef HAVE_CUDA

    IF( ALLOCATED( myComm % neighborRank_dev ) )DEALLOCATE( myComm % neighborRank_dev )
    IF( ALLOCATED( myComm % bufferSize_dev ) )DEALLOCATE( myComm % bufferSize_dev )
    IF( ALLOCATED( myComm % bufferMap_dev ) )DEALLOCATE( myComm % bufferMap_dev )
    IF( ALLOCATED( myComm % rankTable_dev ) )DEALLOCATE( myComm % rankTable_dev )

#endif


#endif

  END SUBROUTINE Trash_BoundaryCommunicator

  SUBROUTINE Finalize_BoundaryCommunicator( myComm )
    IMPLICIT NONE
    CLASS( BoundaryCommunicator ), INTENT(inout) :: myComm

#ifdef HAVE_MPI
      CALL MPI_FINALIZE( myComm % mpiErr )
#endif

  END SUBROUTINE Finalize_BoundaryCommunicator

#ifdef HAVE_CUDA

SUBROUTINE UpdateDevice_BoundaryCommunicator( myComm )

IMPLICIT NONE
CLASS( BoundaryCommunicator ), INTENT(inout) :: myComm


myComm % unPackMap_dev = myComm % unPackMap
myComm % extProcIDs_dev = myComm % extProcIDs
myComm % boundaryIDs_dev = myComm % boundaryIDs
myComm % boundaryGlobalIDs_dev = myComm % boundaryGlobalIDs

#ifdef HAVE_MPI

myComm % neighborRank_dev = myComm % neighborRank
mycomm % bufferSize_dev   = myComm % bufferSize
myComm % rankTable_dev    = myComm % rankTable

#endif

END SUBROUTINE UpdateDevice_BoundaryCommunicator

SUBROUTINE UpdateHost_BoundaryCommunicator( myComm )

IMPLICIT NONE
CLASS( BoundaryCommunicator ), INTENT(inout) :: myComm


myComm % unPackMap   = myComm % unPackMap_dev
myComm % extProcIDs  = myComm % extProcIDs_dev
myComm % boundaryIDs = myComm % boundaryIDs_dev
myComm % boundaryGlobalIDs = myComm % boundaryGlobalIDs_dev

#ifdef HAVE_MPI

myComm % neighborRank  = myComm % neighborRank_dev
mycomm % bufferSize    = myComm % bufferSize_dev
myComm % rankTable     = myComm % rankTable_dev

#endif

END SUBROUTINE UpdateHost_BoundaryCommunicator

#endif
!
!
!==================================================================================================!
!-------------------------------------- FILE I/O ROUTINES -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup BoundaryCommunicator_CLASS
!! @{
! ================================================================================================ !
! S/R WritePickup
!
!> \fn WritePickup_BoundaryCommunicator
!! Writes pickup files for the BoundaryCommunicator DATA structure.
!!
!! Given a file-name base (e.g. "foo"), this routine generates "foo.bcm" (2)
!! The .bcm file CONTAINS the boundary edge, external element, and external process information.
!! The file is an ASCII file.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(BoundaryCommunicator) :: this <BR>
!! <B>CHARACTER</B>                  :: filename
!!         .... <BR>
!!     <B>CALL</B> this % WritePickup( filename ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in <th> myComm <td> BoundaryCommunicator <td> Previously constructed boundary-
!!                               communicator DATA structure
!!   <tr> <td> in <th> filename <td> CHARACTER <td> File base-name for the pickup files
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE WritePickup_BoundaryCommunicator( myComm, filename )

    IMPLICIT NONE
    CLASS( BoundaryCommunicator ), INTENT(in) :: myComm
    CHARACTER(*), INTENT(in)                     :: filename
    ! LOCAL
    INTEGER       :: i, fUnit


    OPEN( UNIT   = NEWUNIT(fUnit), &
      FILE   = TRIM(filename)//'.bcm', &
      FORM   ='FORMATTED',&
      ACCESS ='SEQUENTIAL',&
      STATUS ='REPLACE',&
      ACTION ='WRITE' )

    WRITE( fUnit, * ) myComm % nBoundaries

    DO i = 1, myComm % nBoundaries

      WRITE( fUnit, * ) myComm % boundaryIDs(i), &
                        myComm % extProcIDs(i), &
                        myComm % unPackMap(i)

    ENDDO

    CLOSE(fUnit)


  END SUBROUTINE WritePickup_BoundaryCommunicator
!
!> \addtogroup BoundaryCommunicator_CLASS
!! @{
! ================================================================================================ !
! S/R ReadPickup
!
!> \fn ReadPickup_BoundaryCommunicator
!! Reads pickup files for the BoundaryCommunicator DATA structure.
!!
!! Given a file-name base (e.g. "foo"), this routine reads "foo.bcm"
!! The .bcm file CONTAINS the boundary edge, external element, and external process information.
!! The file is an ASCII file.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(BoundaryCommunicator) :: this <BR>
!! <B>CHARACTER</B>                     :: filename
!!         .... <BR>
!!     <B>CALL</B> this % ReadPickup( filename ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> out <th> myComm <td> BoundaryCommunicator <td> Previously constructed boundary-
!!                               communicator DATA structure
!!   <tr> <td> in <th> filename <td> CHARACTER <td> File base-name for the pickup files
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ReadPickup_BoundaryCommunicator( myComm )

    IMPLICIT NONE
    CLASS( BoundaryCommunicator ), INTENT(inout) :: myComm
    ! LOCAL
    CHARACTER(4)  :: rankChar
    INTEGER       :: i
    INTEGER       :: fUnit
    INTEGER       :: nBe


    WRITE( rankChar, '(I4.4)' ) myComm % myRank

    OPEN( UNIT   = NEWUNIT(fUnit), &
      FILE   = 'ExtComm.'//rankChar//'.bcm', &
      FORM   ='FORMATTED',&
      ACCESS ='SEQUENTIAL',&
      STATUS ='OLD',&
      ACTION ='READ' )

    READ( fUnit, * ) nBe

    CALL myComm % Build( nBe )


    DO i = 1, myComm % nBoundaries

      READ( fUnit, * ) myComm % boundaryIDs(i), &
                       myComm % extProcIDs(i), &
                       myComm % unPackMap(i)

    ENDDO

    CLOSE(fUnit)

#ifdef HAVE_MPI
    CALL myComm % ConstructCommTables(  )
#endif

#ifdef HAVE_CUDA
    CALL myComm % UpdateDevice( )
#endif

  END SUBROUTINE ReadPickup_BoundaryCommunicator
!
  SUBROUTINE SetRanks( myComm )
    IMPLICIT NONE
    CLASS( BoundaryCommunicator ), INTENT(inout) :: myComm
    ! Local
#ifdef HAVE_CUDA
    INTEGER :: iStat, cudaDeviceNumber, nDevices
#endif

    IF( .NOT. myComm % setup )THEN
#ifdef HAVE_MPI
      myComm % MPI_COMM = MPI_COMM_WORLD

      IF( prec == sp )THEN
        myComm % MPI_PREC=MPI_FLOAT
      ELSE
        myComm % MPI_PREC=MPI_DOUBLE
      ENDIF

      CALL MPI_INIT( myComm % mpiErr )
      CALL MPI_COMM_RANK( myComm % MPI_COMM, myComm % myRank, myComm % mpiErr )
      CALL MPI_COMM_SIZE( myComm % MPI_COMM, myComm % nProc, myComm % mpiErr )

#ifdef HAVE_CUDA

      ! Assuming the number of GPU's and the number of ranks per node is unIForm,
      ! each rank is assigned to it's own GPU.
      iStat = cudaGetDeviceCount( nDevices )
      cudaDeviceNumber = MOD( myComm % myRank, nDevices )
      iStat = cudaSetDevice( cudaDeviceNumber )

#endif

#else
      myComm % myRank = 0
      myComm % nProc  = 1
#endif

      myComm % setup = .TRUE.

    ENDIF

  END SUBROUTINE SetRanks

#ifdef HAVE_MPI
  SUBROUTINE ConstructCommTables( myComm )

    IMPLICIT NONE
    CLASS( BoundaryCommunicator ), INTENT(inout) :: myComm
    ! Local
    INTEGER, ALLOCATABLE :: bufferCounter(:)
    INTEGER :: sharedFaceCount(0:myComm % nProc-1)
    INTEGER :: IFace, bID, iNeighbor
    INTEGER :: tag, ierror
    INTEGER :: e1, e2, s1, p2, nmsg, maxFaceCount
    INTEGER :: fUnit


    ! Count up the number of neighboring ranks
    myComm % rankTable = 0
    sharedFaceCount    = 0
    DO bID = 1, myComm % nBoundaries

      p2 = myComm % extProcIDS(bID)

      IF( p2 /= myComm % myRank )THEN

        myComm % rankTable(p2) = 1
        sharedFaceCount(p2) = sharedFaceCount(p2)+1

      ENDIF

    ENDDO

    myComm % nNeighbors = SUM( myComm % rankTable )

    ALLOCATE( myComm % neighborRank(1:myComm % nNeighbors), &
              myComm % bufferSize(1:myComm % nNeighbors), &
              bufferCounter(1:myComm % nNeighbors) )


    ! For each neighbor, set the neighbor's rank
    iNeighbor = 0
    DO p2 = 0, myComm % nProc-1

      IF( myComm % rankTable(p2) == 1 )THEN

        iNeighbor = iNeighbor + 1
        myComm % neighborRank(iNeighbor) = p2
        myComm % rankTable(p2) = iNeighbor

      ENDIF

    ENDDO


    maxFaceCount = MAXVAL( sharedFaceCount )
    DO iNeighbor = 1, myComm % nNeighbors

      p2 = myComm % neighborRank(iNeighbor)
      myComm % bufferSize(iNeighbor) = sharedFaceCount(p2)

    ENDDO


    myComm % maxBufferSize = maxFaceCount
    bufferCounter = 0

    myComm % bufferMap = 0

    DO bID = 1, myComm % nBoundaries

      p2 = myComm % extProcIDs(bID)

      ! In the event that the external process ID (p2) is identical to the current rank (p1),
      ! THEN this boundary edge involves a physical boundary condition and DOes not require a
      ! message exchange

      IF( p2 /= myComm % myRank )THEN

        iNeighbor = myComm % rankTable(p2)

        bufferCounter(iNeighbor) = bufferCounter(iNeighbor) + 1
        myComm % bufferMap(bID)   = bufferCounter(iNeighbor)

      ENDIF

    ENDDO

    DEALLOCATE( bufferCounter )

#ifdef HAVE_CUDA
    ALLOCATE( myComm % neighborRank_dev(1:myComm % nNeighbors), &
              myComm % bufferSize_dev(1:myComm % nNeighbors) )

    myComm % rankTable_dev    = myComm % rankTable
    myComm % neighborRank_dev = myComm % neighborRank
    myComm % bufferSize_dev   = myComm % bufferSize
    myComm % bufferMap_dev    = myComm % bufferMap
#endif

  END SUBROUTINE ConstructCommTables
#endif

END MODULE BoundaryCommunicator_CLASS
