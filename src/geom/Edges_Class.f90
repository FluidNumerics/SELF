! Edges_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE Edges_Class

USE ModelPrecision
USE ConstantsDictionary

IMPLICIT NONE

!  The Edges class defines the attributes necessary to describe an edge in an unstructured mesh.
!

  TYPE Edges
    INTEGER              :: nEdges
    INTEGER, ALLOCATABLE :: edgeID(:)         ! The edge ID
    INTEGER, ALLOCATABLE :: boundaryID(:)     ! If the edge is part of the mesh boundary, the edge gets assigned a boundary edge ID
    INTEGER, ALLOCATABLE :: nodeIDs(:,:)      ! Node IDs which start and terminate this edge
    INTEGER, ALLOCATABLE :: elementIDs(:,:)   ! Neighboring elements IDs across the edge
    INTEGER, ALLOCATABLE :: elementSides(:,:) ! Local side IDs for the neighboring elements
    INTEGER, ALLOCATABLE :: start(:), inc(:)  ! Loop start and increment for the secondary element side

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE :: nEdges_dev
    INTEGER, DEVICE, ALLOCATABLE :: edgeID_dev(:)            ! The edge ID
    INTEGER, DEVICE, ALLOCATABLE :: boundaryID_dev(:)        ! If the edge is part of the mesh boundary, the edge gets assigned a boundary edge ID
    INTEGER, DEVICE, ALLOCATABLE :: nodeIDs_dev(:,:)         ! Node IDs which start and terminate this edge
    INTEGER, DEVICE, ALLOCATABLE :: elementIDs_dev(:,:)      ! Neighboring elements IDs across the edge
    INTEGER, DEVICE, ALLOCATABLE :: elementSides_dev(:,:)    ! Local side IDs for the neighboring elements
    INTEGER, DEVICE, ALLOCATABLE :: start_dev(:), inc_dev(:) ! Loop start and increment for the secondary element side

#endif

    CONTAINS

      PROCEDURE :: Build => Build_Edges
      PROCEDURE :: Trash => Trash_Edges

#ifdef HAVE_CUDA
      PROCEDURE :: UpdateDevice => UpdateDevice_Edges
      PROCEDURE :: UpdateHost   => UpdateHost_Edges
#endif

  END TYPE Edges
   

CONTAINS

  SUBROUTINE Build_Edges( myEdges, nEdges )
    IMPLICIT NONE
    CLASS( Edges ), INTENT(out) :: myEdges
    INTEGER, INTENT(in)         :: nEdges
  
    
      myEdges % nEdges = nEdges
  
      ALLOCATE( myEdges % edgeID(1:nEdges), &
                myEdges % boundaryID(1:nEdges), &
                myEdges % nodeIDs(1:2,1:nEdges), &
                myEdges % elementIDs(1:2,1:nEdges), &
                myEdges % elementSides(1:2,1:nEdges), &
                myEdges % start(1:nEdges), &
                myEdges % inc(1:nEdges) )
  
      myEdges % edgeID       = 0
      myEdges % boundaryID   = 0
      myEdges % nodeIDs      = 0
      myEdges % elementIDs   = 0
      myEdges % elementSides = 0
      myEdges % start        = 0
      myEdges % inc          = 0
  
#ifdef HAVE_CUDA
  
      ALLOCATE( myEdges % nEdges_dev )
      ALLOCATE( myEdges % edgeID_dev(1:nEdges), &
                myEdges % boundaryID_dev(1:nEdges), &
                myEdges % nodeIDs_dev(1:2,1:nEdges), &
                myEdges % elementIDs_dev(1:2,1:nEdges), &
                myEdges % elementSides_dev(1:2,1:nEdges), &
                myEdges % start_dev(1:nEdges), &
                myEdges % inc_dev(1:nEdges) )
  
      myEdges % nEdges_dev = nEdges
  
      myEdges % edgeID_dev       = 0
      myEdges % boundaryID_dev   = 0
      myEdges % nodeIDs_dev      = 0
      myEdges % elementIDs_dev   = 0
      myEdges % elementSides_dev = 0
      myEdges % start_dev        = 0
      myEdges % inc_dev          = 0
  
#endif
  
  END SUBROUTINE Build_Edges

!

  SUBROUTINE Trash_Edges( myEdges, nEdges )
    IMPLICIT NONE
    CLASS( Edges ), INTENT(inout) :: myEdges
    INTEGER, INTENT(in)           :: nEdges
  
    
      DEALLOCATE( myEdges % edgeID, &
                  myEdges % boundaryID, &
                  myEdges % nodeIDs, &
                  myEdges % elementIDs, &
                  myEdges % elementSides, &
                  myEdges % start, &
                  myEdges % inc )
  
#ifdef HAVE_CUDA
  
      DEALLOCATE( myEdges % nEdges_dev )
      DEALLOCATE( myEdges % edgeID_dev, &
                  myEdges % boundaryID_dev, &
                  myEdges % nodeIDs_dev, &
                  myEdges % elementIDs_dev, &
                  myEdges % elementSides_dev, &
                  myEdges % start_dev, &
                  myEdges % inc_dev )
  
#endif
  
  END SUBROUTINE Trash_Edges

!

#ifdef HAVE_CUDA  
  SUBROUTINE UpdateDevice_Edges( myEdges )
    IMPLICIT NONE
    CLASS( Edges ), INTENT(inout) :: myEdges


      myEdges % edgeID_dev       = myEdges % edgeID
      myEdges % boundaryID_dev   = myEdges % boundaryID
      myEdges % nodeIDs_dev      = myEdges % nodeIDs
      myEdges % elementIDs_dev   = myEdges % elementIDs
      myEdges % elementSides_dev = myEdges % elementSides
      myEdges % start_dev        = myEdges % start
      myEdges % inc_dev          = myEdges % inc   
  

  END SUBROUTINE UpdateDevice_Edges

!

  SUBROUTINE UpdateHost_Edges( myEdges )
    IMPLICIT NONE
    CLASS( Edges ), INTENT(inout) :: myEdges


      myEdges % edgeID_dev       = myEdges % edgeID
      myEdges % boundaryID_dev   = myEdges % boundaryID
      myEdges % nodeIDs_dev      = myEdges % nodeIDs
      myEdges % elementIDs_dev   = myEdges % elementIDs
      myEdges % elementSides_dev = myEdges % elementSides
      myEdges % start_dev        = myEdges % start
      myEdges % inc_dev          = myEdges % inc   
  

  END SUBROUTINE UpdateHost_Edges
#endif

END MODULE Edges_Class
