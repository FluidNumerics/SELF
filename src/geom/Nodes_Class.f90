! Nodes_CLASS.f90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Nodes_CLASS

  USE ModelPrecision

  IMPLICIT NONE


!  The Nodes DATA structure defines attributes and TYPE-bound procedures for working with the "node"
!  mesh primitive in an unstructured mesh.
!
!  Nodess, elements, edges, and faces form the foundation of describing an unstructured mesh. The
!  relationship betweens nodes and elements and nodes and edges (or faces in 3-D) define the
!  connectivity in an unstructured mesh. In this DATA structure a node is defined through an
!  INTEGER ID, its TYPE (INTERIOR or BOUNDARY), and its position.

  TYPE Nodes
    INTEGER                 :: nNodes
    INTEGER, ALLOCATABLE    :: nodeID(:)
    INTEGER, ALLOCATABLE    :: nodeTYPE(:)
    REAL(prec), ALLOCATABLE :: x(:,:)

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE    :: nNodes_dev
    INTEGER, DEVICE, ALLOCATABLE    :: nodeID_dev(:)
    INTEGER, DEVICE, ALLOCATABLE    :: nodeTYPE_dev(:)
    REAL(prec), DEVICE, ALLOCATABLE :: x_dev(:,:)
#endif

  CONTAINS

    PROCEDURE :: Build => Build_Nodes
    PROCEDURE :: Trash => Trash_Nodes

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateDevice => UpdateDevice_Nodes
    PROCEDURE :: UpdateHost   => UpdateHost_Nodes
#endif

    PROCEDURE :: ScalePosition => ScalePosition_Nodes

  END TYPE Nodes

CONTAINS

  SUBROUTINE Build_Nodes( myNodes, nNodes )
    IMPLICIT NONE
    CLASS( Nodes ), INTENT(out) :: myNodes
    INTEGER, INTENT(in)         :: nNodes

    myNodes % nNodes = nNodes

    ALLOCATE( myNodes % nodeID(1:nNodes), &
      myNodes % nodeTYPE(1:nNodes), &
      myNodes % x(1:3,1:nNodes) )

    myNodes % nodeID   = 0
    myNodes % nodeTYPE = 0
    myNodes % x        = 0.0_prec

#ifdef HAVE_CUDA

    ALLOCATE( myNodes % nNodes_dev )
    ALLOCATE( myNodes % nodeID_dev(1:nNodes), &
      myNodes % nodeTYPE_dev(1:nNodes), &
      myNodes % x_dev(1:3,1:nNodes) )

    myNodes % nNodes_dev = nNodes

    myNodes % nodeID_dev   = 0
    myNodes % nodeTYPE_dev = 0
    myNodes % x_dev        = 0.0_prec

#endif

  END SUBROUTINE Build_Nodes

!

  SUBROUTINE Trash_Nodes( myNodes )
    IMPLICIT NONE
    CLASS( Nodes ), INTENT(inout) :: myNodes


    DEALLOCATE( myNodes % nodeID, &
      myNodes % nodeTYPE, &
      myNodes % x )

#ifdef HAVE_CUDA

    DEALLOCATE( myNodes % nNodes_dev )
    DEALLOCATE( myNodes % nodeID_dev, &
      myNodes % nodeTYPE_dev, &
      myNodes % x_dev )

#endif

  END SUBROUTINE Trash_Nodes

!
#ifdef HAVE_CUDA
  SUBROUTINE UpdateDevice_Nodes( myNodes )
    IMPLICIT NONE
    CLASS( Nodes ), INTENT(inout) :: myNodes


    myNodes % nodeID_dev   = myNodes % nodeID
    myNodes % nodeTYPE_dev = myNodes % nodeTYPE
    myNodes % x_dev        = myNodes % x


  END SUBROUTINE UpdateDevice_Nodes

!

  SUBROUTINE UpdateHost_Nodes( myNodes )
    IMPLICIT NONE
    CLASS( Nodes ), INTENT(inout) :: myNodes


    myNodes % nodeID   = myNodes % nodeID_dev
    myNodes % nodeTYPE = myNodes % nodeTYPE_dev
    myNodes % x        = myNodes % x_dev


  END SUBROUTINE UpdateHost_Nodes
#endif
!

  SUBROUTINE ScalePosition_Nodes( myNodes, xScale, yScale, zScale )
    IMPLICIT NONE
    CLASS( Nodes ), INTENT(inout) :: myNodes
    REAL(prec), INTENT(in)        :: xScale, yScale, zScale
    ! Local
    INTEGER :: i

    DO i = 1, myNodes % nNodes

      myNodes % x(1,i) = xScale*myNodes % x(1,i)
      myNodes % x(2,i) = yScale*myNodes % x(2,i)
      myNodes % x(3,i) = zScale*myNodes % x(3,i)

    ENDDO

#ifdef HAVE_CUDA

    CALL myNodes % UpdateDevice( )

#endif

  END SUBROUTINE ScalePosition_Nodes

END MODULE Nodes_CLASS
