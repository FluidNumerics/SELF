! Facess_Class.f90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Faces_Class

! src/COMMON/
  USE ModelPrecision
  USE ConstantsDictionary

  IMPLICIT NONE

  TYPE Faces
    INTEGER              :: nFaces, N
    INTEGER, ALLOCATABLE :: faceID(:)            ! The global face ID
    INTEGER, ALLOCATABLE :: boundaryID(:)        ! IF the face is part of the mesh boundary, the face gets assigned a boundary face ID
    INTEGER, ALLOCATABLE :: nodeIDs(:,:)         ! Node IDs which start and terminate this face
    INTEGER, ALLOCATABLE :: elementIDs(:,:)      ! Neighboring elements IDs across the face
    INTEGER, ALLOCATABLE :: elementSides(:,:)    ! Local side IDs for the neighboring elements
    INTEGER, ALLOCATABLE :: iStart(:), iInc(:)   ! Loop start and increment for the secondary element side (1st computational direction)
    INTEGER, ALLOCATABLE :: jStart(:), jInc(:)   ! Loop start and increment for the secondary element side (2nd computational direction)
    INTEGER, ALLOCATABLE :: swapDimensions(:)    ! A flag USEd to swap the computational directions between the primary and secondary elements
    INTEGER, ALLOCATABLE :: iMap(:,:,:)
    INTEGER, ALLOCATABLE :: jMap(:,:,:)

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: nFaces_dev, N_dev
    INTEGER, DEVICE, ALLOCATABLE :: faceID_dev(:)            ! The face ID
    INTEGER, DEVICE, ALLOCATABLE :: boundaryID_dev(:)        ! IF the face is part of the mesh boundary, the face gets assigned a boundary face ID
    INTEGER, DEVICE, ALLOCATABLE :: nodeIDs_dev(:,:)      ! Node IDs which start and terminate this face
    INTEGER, DEVICE, ALLOCATABLE :: elementIDs_dev(:,:)   ! Neighboring elements IDs across the face
    INTEGER, DEVICE, ALLOCATABLE :: elementSides_dev(:,:) ! Local side IDs for the neighboring elements
    INTEGER, DEVICE, ALLOCATABLE :: iMap_dev(:,:,:)
    INTEGER, DEVICE, ALLOCATABLE :: jMap_dev(:,:,:)

#endif

  CONTAINS

    PROCEDURE :: Build => Build_Faces
    PROCEDURE :: Trash => Trash_Faces

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateHost => UpdateHost_Faces
    PROCEDURE :: UpdateDevice => UpdateDevice_Faces
#endif

    PROCEDURE :: nBoundaryFaces

  END TYPE Faces


CONTAINS

  SUBROUTINE Build_Faces( myFaces, nFaces, N )
    IMPLICIT NONE
    Class( Faces ), INTENT(out) :: myFaces
    INTEGER, INTENT(in)         :: nFaces, N


    myFaces % nFaces = nFaces
    myFaces % N      = N

    ALLOCATE( myFaces % faceID(1:nFaces), &
      myFaces % boundaryID(1:nFaces), &
      myFaces % nodeIDs(1:4,1:nFaces), &
      myFaces % elementIDs(1:2,1:nFaces), &
      myFaces % elementSides(1:2,1:nFaces), &
      myFaces % iStart(1:nFaces), &
      myFaces % iInc(1:nFaces), &
      myFaces % jStart(1:nFaces), &
      myFaces % jInc(1:nFaces), &
      myFaces % swapDimensions(1:nFaces), &
      myFaces % iMap(0:N,0:N,1:nFaces), &
      myFaces % jMap(0:N,0:N,1:nFaces) )

    myFaces % faceID         = 0
    myFaces % boundaryID     = 0
    myFaces % nodeIDs        = 0
    myFaces % elementIDs     = 0
    myFaces % elementSides   = 0
    myFaces % iStart         = 0
    myFaces % iInc           = 1
    myFaces % jStart         = 0
    myFaces % jInc           = 1
    myFaces % swapDimensions = 0
    myFaces % iMap           = 0
    myFaces % jMap           = 0

#ifdef HAVE_CUDA

    ALLOCATE( myFaces % nFaces_dev, myFaces % N_dev )

    myFaces % nFaces_dev = nFaces
    myFaces % N_dev      = N

    ALLOCATE( myFaces % faceID_dev(1:nFaces), &
      myFaces % boundaryID_dev(1:nFaces), &
      myFaces % nodeIDs_dev(1:4,1:nFaces), &
      myFaces % elementIDs_dev(1:2,1:nFaces), &
      myFaces % elementSides_dev(1:2,1:nFaces), &
      myFaces % iMap_dev(0:N,0:N,1:nFaces), &
      myFaces % jMap_dev(0:N,0:N,1:nFaces) )

    myFaces % faceID_dev       = 0
    myFaces % boundaryID_dev   = 0
    myFaces % nodeIDs_dev      = 0
    myFaces % elementIDs_dev   = 0
    myFaces % elementSides_dev = 0
    myFaces % iMap_dev         = 0
    myFaces % jMap_dev         = 0

#endif

  END SUBROUTINE Build_Faces
!
  SUBROUTINE Trash_Faces( myFaces )
    IMPLICIT NONE
    Class( Faces ), INTENT(inout) :: myFaces



    DEALLOCATE( myFaces % faceID, &
      myFaces % boundaryID, &
      myFaces % nodeIDs, &
      myFaces % elementIDs, &
      myFaces % elementSides, &
      myFaces % iStart, &
      myFaces % iInc, &
      myFaces % jStart, &
      myFaces % jInc, &
      myFaces % swapDimensions, &
      myFaces % iMap, &
      myFaces % jMap )

#ifdef HAVE_CUDA

    DEALLOCATE( myFaces % nFaces_dev, myFaces % N_dev )
    DEALLOCATE( myFaces % faceID_dev, &
      myFaces % boundaryID_dev, &
      myFaces % nodeIDs_dev, &
      myFaces % elementIDs_dev, &
      myFaces % elementSides_dev, &
      myFaces % iMap_dev, &
      myFaces % jMap_dev )

#endif

  END SUBROUTINE Trash_Faces

#ifdef HAVE_CUDA

  SUBROUTINE UpdateHost_Faces( myFaces )
    IMPLICIT NONE
    Class( Faces ), INTENT(inout) :: myFaces

    myFaces % faceID       = myFaces % faceID_dev
    myFaces % boundaryID   = myFaces % boundaryID_dev
    myFaces % nodeIDs      = myFaces % nodeIDs_dev
    myFaces % elementIDs   = myFaces % elementIDs_dev
    myFaces % elementSides = myFaces % elementSides_dev
    myFaces % iMap         = myFaces % iMap_dev
    myFaces % jMap         = myFaces % jMap_dev

  END SUBROUTINE UpdateHost_Faces

  SUBROUTINE UpdateDevice_Faces( myFaces )
    IMPLICIT NONE
    Class( Faces ), INTENT(inout) :: myFaces

    myFaces % faceID_dev       = myFaces % faceID
    myFaces % boundaryID_dev   = myFaces % boundaryID
    myFaces % nodeIDs_dev      = myFaces % nodeIDs
    myFaces % elementIDs_dev   = myFaces % elementIDs
    myFaces % elementSides_dev = myFaces % elementSides
    myFaces % iMap_dev         = myFaces % iMap
    myFaces % jMap_dev         = myFaces % jMap

  END SUBROUTINE UpdateDevice_Faces

#endif

  FUNCTION nBoundaryFaces( myFaces ) RESULT( nbf )
    CLASS( Faces ) :: myFaces
    INTEGER :: nbf
    INTEGER :: i

    nbf = 0
    DO i = 1, myFaces % nFaces
      IF( myFaces % elementIDs(2,i) < 0 )THEN
        nbf = nbf + 1
      ENDIF
    ENDDO

  END FUNCTION nBoundaryFaces

END MODULE Faces_Class
