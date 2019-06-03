! MPILayer_Class.f90
!
! Copyright 2019 Fluid Numerics LLC
! Author : Joseph Schoonover <joe@fluidnumerics.com>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE MPILayer_Class

  USE ModelPrecision
  USE NodalDGSolution_3D_Class
  USE Faces_Class

  IMPLICIT NONE


#include "self_macros.h"

  TYPE MPILayer

    INTEGER              :: nVars, N, nMessages
    INTEGER, ALLOCATABLE :: boundaryToFaceID(:) 
    INTEGER, ALLOCATABLE :: requestHandle(:), requestStats(:,:) 

    CONTAINS

      PROCEDURE :: Build => Build_MPILayer
      PROCEDURE :: Trash => Trash_MPILayer

      PROCEDURE :: MPI_Exchange
      PROCEDURE :: Finalize_MPI_Exchange

  END TYPE MPILayer

  INTEGER :: mpiComm, mpiPrec, myRank, nProc
CONTAINS

  SUBROUTINE Build_MPILayer( myMPI, N, nVars, nMessages )
#undef __FUNC__
#define __FUNC__ "Build_MPILayer"
    CLASS( MPILayer ), INTENT(out) :: myMPI
    INTEGER, INTENT(in)            :: N, nVars, nMessages

    INFO('Start')
#ifdef HAVE_MPI

    ALLOCATE( myMPI % requestHandle(1:2*nMessages), &
              myMPI % requestStats(MPI_STATUS_SIZE,1:2*nMessages) )

    myMPI % requestHandle = 0
    myMPI % requestStats  = 0
#endif

    myMPI % N         = N
    myMPI % nVars     = nVars
    myMPI % nMessages = nMessages

    INFO('End')

  END SUBROUTINE Build_MPILayer
!
  SUBROUTINE Trash_MPILayer( myMPI )
#undef __FUNC__
#define __FUNC__ "Trash_MPILayer"
    CLASS( MPILayer ), INTENT(inout) :: myMPI

#ifdef HAVE_MPI
    DEALLOCATE( myMPI % requestHandle, &
                myMPI % requestStats )
#endif

  END SUBROUTINE Trash_MPILayer

  SUBROUTINE Initialize_MPILayer( )
#undef __FUNC__
#define __FUNC__ "Initialize_MPILayer"
    INTEGER       :: ierror
    CHARACTER(30) :: msg

    INFO('Start')
    myRank = 0
    nProc  = 1
#ifdef HAVE_MPI
    mpiComm = MPI_COMM_WORLD
    CALL MPI_INIT( ierror )
    CALL MPI_COMM_RANK( mpiComm, myRank, ierror )
    CALL MPI_COMM_SIZE( mpiComm, nProc,  ierror )

    IF( prec == sp )THEN
      mpiPrec=MPI_FLOAT
    ELSE
      mpiPrec=MPI_DOUBLE
    ENDIF
#endif
    WRITE(msg,'(I5)')myRank
    msg="Greetings from rank "//TRIM(msg)//"."
    INFO(TRIM(msg))
    INFO('End')

  END SUBROUTINE Initialize_MPILayer

  SUBROUTINE Finalize_MPILayer( )

    INTEGER    :: ierror
#ifdef HAVE_MPI
      CALL MPI_FINALIZE( ierror )
#endif

  END SUBROUTINE Finalize_MPILayer

  SUBROUTINE MPI_Exchange( myMPI, state, meshfaces, boundaryToFaceID, boundaryToProcID )
#undef __FUNC__
#define __FUNC__ "MPI_Exchange"
    CLASS( MPILayer ), INTENT(inout)          :: myMPI
    TYPE( NodalDGSolution_3D ), INTENT(inout) :: state
    TYPE( Faces ), INTENT(in)                 :: meshFaces
    INTEGER, INTENT(in)                       :: boundaryToFaceID(1:myMPI % nMessages)
    INTEGER, INTENT(in)                       :: boundaryToProcID(1:myMPI % nMessages)
#ifdef HAVE_MPI
    ! Local
    INTEGER    :: iError
    INTEGER    :: local_face_id, bID, message_id
    INTEGER    :: global_face_id
    INTEGER    :: local_element_id, e2, side_id, external_proc_id
    
INFO('Start')
    message_id = 0
    DO bID = 1, state % nBoundaryFaces

      local_face_id    = boundaryToFaceID(bID)
      external_proc_id = boundaryToProcID(bID)

      IF( external_proc_id /= myRank )THEN
 
        message_id = message_id + 1
        local_element_id = meshFaces % elementIDs(1,local_face_id)
        side_id          = meshFaces % elementSides(1,local_face_id)
        global_face_id   = meshFaces % faceID(local_face_id) 

        CALL MPI_IRECV( state % externalState(:,:,:,bID), &
                        (myMPI % N+1)*(myMPI % N+1)*myMPI % nVars, &
                        mpiPrec,   &
                        external_proc_id, global_face_id,  &
                        mpiComm,   &
                        myMPI % requestHandle(message_id*2-1), iError )
  
        CALL MPI_ISEND( state % boundarySolution(:,:,:,side_id,local_element_id), &
                        (myMPI % N+1)*(myMPI % N+1)*myMPI % nVars, &
                        mpiPrec, &
                        external_proc_id, global_face_id, &
                        mpiComm, &
                        myMPI % requestHandle(message_id*2), iError)
  
      ENDIF
    ENDDO

INFO('End')
#endif

  END SUBROUTINE MPI_Exchange
!
  SUBROUTINE Finalize_MPI_Exchange( myMPI )
#undef __FUNC__
#define __FUNC__ "Finalize_MPI_Exchange"
    CLASS( MPILayer ), INTENT(inout)          :: myMPI
#ifdef HAVE_MPI
    ! Local
    INTEGER    :: ierror

    CALL MPI_WaitAll( myMPI % nMessages*2, &
                      myMPI % requestHandle, &
                      myMPI % requestStats, &
                      iError)

#endif

  END SUBROUTINE Finalize_MPI_Exchange
!
END MODULE MPILayer_Class


