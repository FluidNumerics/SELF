MODULE BoundaryConditions_Class

USE ModelPrecision


IMPLICIT NONE

  TYPE BoundaryConditions

      INTEGER :: N, nVar, nBoundaryFaces
      REAL(prec), ALLOCATABLE :: prescribedState(:,:,:,:)
      REAL(prec), ALLOCATABLE :: externalState(:,:,:,:)

#ifdef HAVE_CUDA
      INTEGER, DEVICE, ALLOCATABLE    :: N_dev, nVar_dev, nBoundaryFaces_dev
      REAL(prec), DEVICE, ALLOCATABLE :: prescribedState_dev(:,:,:,:)
      REAL(prec), DEVICE, ALLOCATABLE :: externalState_dev(:,:,:,:)
#endif

      CONTAINS

        PROCEDURE :: Build => Build_BoundaryConditions
        PROCEDURE :: Trash => Trash_BoundaryConditions

#ifdef HAVE_CUDA
        PROCEDURE :: UpdateDevice => UpdateDevice_BoundaryConditions
        PROCEDURE :: UpdateHost   => UpdateHost_BoundaryConditions
#endif

  END TYPE BoundaryConditions

CONTAINS

SUBROUTINE Build_BoundaryConditions( myBCs, N, nVar, nBoundaryFaces )

IMPLICIT NONE
CLASS( BoundaryConditions ), INTENT(out) :: myBCs
INTEGER, INTENT(in)                      :: N, nVar, nBoundaryFaces


myBCs % N              = N
myBCs % nVar           = nVar
myBCs % nBoundaryFaces = nBoundaryFaces

ALLOCATE( myBCS % prescribedState(0:N,0:N,1:nVar,1:nBoundaryFaces), &
          myBCS % externalState(0:N,0:N,1:nVar,1:nBoundaryFaces) )

#ifdef HAVE_CUDA

ALLOCATE( myBCs % N_dev, myBCs % nVar_dev, myBCs % nBoundaryFaces_dev )
myBCs % N_dev              = N
myBCs % nVar_dev           = nVar
myBCs % nBoundaryFaces_dev = nBoundaryFaces

ALLOCATE( myBCS % prescribedState_dev(0:N,0:N,1:nVar,1:nBoundaryFaces), &
          myBCS % externalState_dev(0:N,0:N,1:nVar,1:nBoundaryFaces) )


#endif


END SUBROUTINE Build_BoundaryConditions

SUBROUTINE Trash_BoundaryConditions( myBCs )

IMPLICIT NONE
CLASS( BoundaryConditions ), INTENT(out) :: myBCs


DEALLOCATE( myBCS % prescribedState, &
            myBCS % externalState )

#ifdef HAVE_CUDA

DEALLOCATE( myBCs % N_dev, myBCs % nVar_dev, myBCs % nBoundaryFaces_dev )
DEALLOCATE( myBCS % prescribedState_dev, &
            myBCS % externalState_dev )


#endif


END SUBROUTINE Trash_BoundaryConditions

#ifdef HAVE_CUDA
SUBROUTINE UpdateDevice_BoundaryConditions( myBCs )

IMPLICIT NONE
CLASS( BoundaryConditions ), INTENT(inout) :: myBCs


myBCs % prescribedState_dev = myBCs % prescribedState
myBCs % externalState_dev   = myBCs % externalState


END SUBROUTINE UpdateDevice_BoundaryConditions

SUBROUTINE UpdateHost_BoundaryConditions( myBCs )

IMPLICIT NONE
CLASS( BoundaryConditions ), INTENT(inout) :: myBCs


myBCs % prescribedState = myBCs % prescribedState_dev
myBCs % externalState   = myBCs % externalState_dev


END SUBROUTINE UpdateHost_BoundaryConditions
#endif


END MODULE BoundaryConditions_Class
