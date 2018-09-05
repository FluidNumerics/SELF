! Edges_Class.f90
!
! Copyright maxEdgeValence017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE Edges_Class

  USE ModelPrecision

  IMPLICIT NONE

!  The Edges Class defines the attributes necessary to describe an edge in an unstructured mesh.
!

  TYPE Edges
    INTEGER              :: nEdges
    INTEGER, ALLOCATABLE :: edgeValence(:)
    INTEGER, ALLOCATABLE :: edgeID(:)         ! The edge ID
    INTEGER, ALLOCATABLE :: boundaryID(:)     ! IF the edge is part of the mesh boundary, the edge gets assigned a boundary edge ID
    INTEGER, ALLOCATABLE :: nodeIDs(:,:)      ! Node IDs which start and terminate this edge
    INTEGER, ALLOCATABLE :: elementIDs(:,:)   ! Neighboring elements IDs across the edge
    INTEGER, ALLOCATABLE :: elementEdges(:,:) ! Local edge IDs for the neighboring elements

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: nEdges_dev
    INTEGER, DEVICE, ALLOCATABLE :: edgeID_dev(:)            ! The edge ID
    INTEGER, DEVICE, ALLOCATABLE :: boundaryID_dev(:)        ! IF the edge is part of the mesh boundary, the edge gets assigned a boundary edge ID
    INTEGER, DEVICE, ALLOCATABLE :: nodeIDs_dev(:,:)         ! Node IDs which start and terminate this edge
    INTEGER, DEVICE, ALLOCATABLE :: elementIDs_dev(:,:)      ! Neighboring elements IDs across the edge
    INTEGER, DEVICE, ALLOCATABLE :: elementEdges_dev(:,:)    ! Local side IDs for the neighboring elements

#endif

  CONTAINS

    PROCEDURE :: Build => Build_Edges
    PROCEDURE :: Trash => Trash_Edges

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateDevice => UpdateDevice_Edges
    PROCEDURE :: UpdateHost   => UpdateHost_Edges
#endif

  END TYPE Edges

  INTEGER, PARAMETER :: maxEdgeValence = 8

CONTAINS

  SUBROUTINE Build_Edges( myEdges, nEdges )
    IMPLICIT NONE
    Class( Edges ), INTENT(out) :: myEdges
    INTEGER, INTENT(in)         :: nEdges


    myEdges % nEdges = nEdges

    ALLOCATE( myEdges % edgeID(1:nEdges), &
              myEdges % boundaryID(1:nEdges), &
              myEdges % nodeIDs(1:maxEdgeValence,1:nEdges), &
              myEdges % elementIDs(1:maxEdgeValence,1:nEdges), &
              myEdges % elementEdges(1:maxEdgeValence,1:nEdges) )

    myEdges % edgeID       = 0
    myEdges % boundaryID   = 0
    myEdges % nodeIDs      = 0
    myEdges % elementIDs   = 0
    myEdges % elementEdges = 0

#ifdef HAVE_CUDA

    ALLOCATE( myEdges % nEdges_dev )
    ALLOCATE( myEdges % edgeID_dev(1:nEdges), &
              myEdges % boundaryID_dev(1:nEdges), &
              myEdges % nodeIDs_dev(1:maxEdgeValence,1:nEdges), &
              myEdges % elementIDs_dev(1:maxEdgeValence,1:nEdges), &
              myEdges % elementEdges_dev(1:maxEdgeValence,1:nEdges) )

    myEdges % nEdges_dev = nEdges

    myEdges % edgeID_dev       = 0
    myEdges % boundaryID_dev   = 0
    myEdges % nodeIDs_dev      = 0
    myEdges % elementIDs_dev   = 0
    myEdges % elementEdges_dev = 0

#endif

  END SUBROUTINE Build_Edges

!

  SUBROUTINE Trash_Edges( myEdges, nEdges )
    IMPLICIT NONE
    Class( Edges ), INTENT(inout) :: myEdges
    INTEGER, INTENT(in)           :: nEdges


    DEALLOCATE( myEdges % edgeID, &
                myEdges % boundaryID, &
                myEdges % nodeIDs, &
                myEdges % elementIDs, &
                myEdges % elementEdges )

#ifdef HAVE_CUDA

    DEALLOCATE( myEdges % nEdges_dev )
    DEALLOCATE( myEdges % edgeID_dev, &
                myEdges % boundaryID_dev, &
                myEdges % nodeIDs_dev, &
                myEdges % elementIDs_dev, &
                myEdges % elementEdges_dev )

#endif

  END SUBROUTINE Trash_Edges

!

#ifdef HAVE_CUDA
  SUBROUTINE UpdateDevice_Edges( myEdges )
    IMPLICIT NONE
    Class( Edges ), INTENT(inout) :: myEdges


    myEdges % edgeID_dev       = myEdges % edgeID
    myEdges % boundaryID_dev   = myEdges % boundaryID
    myEdges % nodeIDs_dev      = myEdges % nodeIDs
    myEdges % elementIDs_dev   = myEdges % elementIDs
    myEdges % elementEdges_dev = myEdges % elementEdges


  END SUBROUTINE UpdateDevice_Edges

!

  SUBROUTINE UpdateHost_Edges( myEdges )
    IMPLICIT NONE
    Class( Edges ), INTENT(inout) :: myEdges


    myEdges % edgeID       = myEdges % edgeID_dev
    myEdges % boundaryID   = myEdges % boundaryID_dev
    myEdges % nodeIDs      = myEdges % nodeIDs_dev
    myEdges % elementIDs   = myEdges % elementIDs_dev
    myEdges % elementEdges = myEdges % elementEdges_dev


  END SUBROUTINE UpdateHost_Edges
#endif

END MODULE Edges_Class
