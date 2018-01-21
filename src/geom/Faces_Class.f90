! Facess_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Faces_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary

IMPLICIT NONE


  TYPE Faces
    INTEGER              :: nFaces
    INTEGER, ALLOCATABLE :: faceID(:)            ! The face ID
    INTEGER, ALLOCATABLE :: boundaryID(:)        ! If the face is part of the mesh boundary, the face gets assigned a boundary face ID
    INTEGER, ALLOCATABLE :: nodeIDs(:,:)      ! Node IDs which start and terminate this face
    INTEGER, ALLOCATABLE :: elementIDs(:,:)   ! Neighboring elements IDs across the face
    INTEGER, ALLOCATABLE :: elementSides(:,:) ! Local side IDs for the neighboring elements
    INTEGER, ALLOCATABLE :: iStart(:), iInc(:)      ! Loop start and increment for the secondary element side (1st computational direction)
    INTEGER, ALLOCATABLE :: jStart(:), jInc(:)      ! Loop start and increment for the secondary element side (2nd computational direction)
    INTEGER, ALLOCATABLE :: swapDimensions(:)    ! A flag used to swap the computational directions
!                                         between the primary and secondary elements

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: nFaces_dev
    INTEGER, DEVICE, ALLOCATABLE :: faceID_dev(:)            ! The face ID
    INTEGER, DEVICE, ALLOCATABLE :: boundaryID_dev(:)        ! If the face is part of the mesh boundary, the face gets assigned a boundary face ID
    INTEGER, DEVICE, ALLOCATABLE :: nodeIDs_dev(:,:)      ! Node IDs which start and terminate this face
    INTEGER, DEVICE, ALLOCATABLE :: elementIDs_dev(:,:)   ! Neighboring elements IDs across the face
    INTEGER, DEVICE, ALLOCATABLE :: elementSides_dev(:,:) ! Local side IDs for the neighboring elements
    INTEGER, DEVICE, ALLOCATABLE :: iStart_dev(:), iInc_dev(:)      ! Loop start and increment for the secondary element side (1st computational direction)
    INTEGER, DEVICE, ALLOCATABLE :: jStart_dev(:), jInc_dev(:)      ! Loop start and increment for the secondary element side (2nd computational direction)
    INTEGER, DEVICE, ALLOCATABLE :: swapDimensions_dev(:)    ! A flag used to swap the computational directions

#endif

     CONTAINS

  END TYPE Faces

END MODULE Faces_Class
