MODULE SEMMesh

USE ModelPrecision
USE ConstantsDictionary
USE NodalSEMData
USE Lagrange_Class


IMPLICIT NONE

  TYPE, PUBLIC :: SEMGeometry2D
    INTEGER :: qType ! Quadrature Type
    TYPE(Lagrange) :: interp ! Interpolant for the Geometry
    TYPE(SEMVector2D) :: x ! Physical positions
    TYPE(SEMTensor2D) :: dxds ! Covariant basis vectors
    TYPE(SEMTensor2D) :: dsdx ! Contavariant basis vectors
    TYPE(SEMScalar2D) :: J ! Jacobian of the transformation

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMGeometry2D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMGeometry2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMGeometry2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMGeometry2D
#endif

  END TYPE SEMGeometry2D

  TYPE, PUBLIC :: SEMGeometry3D
    INTEGER :: qType ! Quadrature Type
    TYPE(Lagrange) :: interp ! Interpolant
    TYPE(SEMVector3D) :: x ! Physical positions
    TYPE(SEMTensor3D) :: dxds ! Covariant basis vectors
    TYPE(SEMTensor3D) :: dsdx ! Contavariant basis vectors
    TYPE(SEMScalar3D) :: J ! Jacobian of the transformation

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_SEMGeometry3D
      PROCEDURE, PUBLIC :: Trash => Trash_SEMGeometry3D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_SEMGeometry3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_SEMGeometry3D
#endif

  END TYPE SEMGeometry3D

  TYPE, PUBLIC :: Nodes
    INTEGER :: objCount
    INTEGER :: ndim
    INTEGER, POINTER :: id(:)
    INTEGER, POINTER :: nodeFlag(:)
    INTEGER, POINTER :: eValence(:)
    INTEGER, POINTER :: elements(:,:)
    INTEGER, POINTER :: elementNodes(:,:)
    REAL(prec), ALLOCATABLE :: x(:,:)

!    CONTAINS
!
!      PROCEDURE, PUBLIC :: Build => Build_Nodes
!      PROCEDURE, PUBLIC :: Trash => Trash_Nodes

  END TYPE Nodes

  TYPE, PUBLIC :: Edges
    INTEGER :: objCount
    INTEGER, POINTER :: id(:)
    INTEGER, POINTER :: eValence(:)
    INTEGER, POINTER :: boundaryID(:)
    INTEGER, POINTER :: nodes(:,:)
    INTEGER, POINTER :: elements(:,:)
    INTEGER, POINTER :: elementEdges(:,:)

!    CONTAINS

  END TYPE Edges

  TYPE, PUBLIC :: Faces
    INTEGER :: N
    INTEGER :: objCount
    INTEGER, POINTER :: id(:)
    INTEGER, POINTER :: boundaryID(:)
    INTEGER, POINTER :: nodes(:,:)
    INTEGER, POINTER :: elements(:,:)
    INTEGER, POINTER :: elementFaces(:,:)
    INTEGER, POINTER :: iMap(:,:,:)
    INTEGER, POINTER :: jMap(:,:,:)

!    CONTAINS

  END TYPE Faces

  TYPE, PUBLIC :: Elements
    INTEGER :: objCount
    INTEGER :: objType
    INTEGER, POINTER :: id(:)
    INTEGER, POINTER :: nodes(:,:)
    INTEGER, POINTER :: neighbors(:,:)

!    CONTAINS

  END TYPE Elements

!  TYPE, PUBLIC :: MeshObjHandles
!    INTEGER :: nNodes, nEdges, nFaces, nElements
!    INTEGER, POINTER :: nodeids(:)
!    INTEGER, POINTER :: edgeids(:)
!    INTEGER, POINTER :: faceids(:)
!    INTEGER, POINTER :: elementids(:)
!
!    CONTAINS
!
!  END TYPE MeshObjHandles
!
!  TYPE, PUBLIC :: DomainDecomposition
!    INTEGER :: nBlocks, nGlobalElements, nMPIMessages, nBoundaryFaces
!    INTEGER, POINTER :: blockID(:) 
!    INTEGER, POINTER :: localID(:) 
!    TYPE(MeshObjectList), ALLOCATABLE :: mesh_obj(:) 
!
!    CONTAINS
!
!  END TYPE DomainDecomposition

  TYPE, PUBLIC :: SEMMesh2D
    TYPE( SEMGeometry2D ) :: geometry
    TYPE( Elements ) :: elements
    TYPE( Nodes  ) :: nodes
    TYPE( Edges ) :: edges 
    !TYPE( DomainDecomposition ) :: decomp 
    INTEGER :: cornerMap(1:2,1:4)
    INTEGER :: sideMap(1:4)
    INTEGER :: edgeMap(1:2,1:4)

!    CONTAINS

  END TYPE SEMMesh2D

  TYPE, PUBLIC :: SEMMesh3D
    TYPE( SEMGeometry3D ) :: geometry
    TYPE( Elements ) :: elements
    TYPE( Nodes  ) :: nodes
    TYPE( Edges ) :: edges 
    TYPE( Faces ) :: faces
    !TYPE( DomainDecomposition ) :: decomp 
    INTEGER :: cornerMap(1:3,1:8)
    INTEGER :: sideMap(1:6)
    INTEGER :: faceMap(1:4,1:6)
    INTEGER :: edgeMap(1:2,1:12)
    INTEGER :: edgeFaceMap(1:2,1:4)
  
!    CONTAINS

!      PROCEDURE, PUBLIC :: Build => Build_SEMMesh3D
!      PROCEDURE, PUBLIC :: Trash => Trash_SEMMesh3D
!
!      PROCEDURE, PUBLIC :: Read_CGNSMesh
!      PROCEDURE, PUBLIC :: Read_UCDMesh
!      PROCEDURE, PUBLIC :: Read_TrellisUCDMesh
!
!      PROCEDURE, PUBLIC :: GenerateConnectivity => GenerateConnectivity_SEMMesh3D

  END TYPE SEMMesh3D

CONTAINS

SUBROUTINE Build_SEMGeometry2D( myGeom, quadrature, polyDegree, nPlotPoints, nElem ) 
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(out) :: myGeom
  INTEGER, INTENT(in) :: quadrature
  INTEGER, INTENT(in) :: polyDegree
  INTEGER, INTENT(in) :: nPlotPoints
  INTEGER, INTENT(in) :: nElem

    myGeom % qType = quadrature
    CALL myGeom % interp % Build( N = polyDegree, &
                                  controlNodeType = quadrature, &
                                  M = nPlotPoints, &
                                  targetNodeType = UNIFORM )

    CALL myGeom % x %  Build( N = polyDegree, &
                              nVar = 2, &
                              nElem = nElem ) 

    CALL myGeom % dxds %  Build( N = polyDegree, &
                                 nVar = 2, &
                                 nElem = nElem ) 

    CALL myGeom % dsdx %  Build( N = polyDegree, &
                                 nVar = 2, &
                                 nElem = nElem ) 

    CALL myGeom % J %  Build( N = polyDegree, &
                                 nVar = 2, &
                                 nElem = nElem ) 
    
END SUBROUTINE Build_SEMGeometry2D

SUBROUTINE Trash_SEMGeometry2D( myGeom ) 
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(inout) :: myGeom

    CALL myGeom % interp % Trash()
    CALL myGeom % x % Trash()
    CALL myGeom % dxds % Trash()
    CALL myGeom % dsdx % Trash()
    CALL myGeom % J % Trash()
  
END SUBROUTINE Trash_SEMGeometry2D 
#ifdef GPU
SUBROUTINE UpdateHost_SEMGeometry2D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(inout) :: myGeom
 
    CALL myGeom % interp % UpdateHost()
    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % J % UpdateHost()

END SUBROUTINE UpdateHost_SEMGeometry2D

SUBROUTINE UpdateDevice_SEMGeometry2D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(inout) :: myGeom
 
    CALL myGeom % interp % UpdateDevice()
    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMGeometry2D
#endif

SUBROUTINE Build_SEMGeometry3D( myGeom, quadrature, polyDegree, nPlotPoints, nElem ) 
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(out) :: myGeom
  INTEGER, INTENT(in) :: quadrature
  INTEGER, INTENT(in) :: polyDegree
  INTEGER, INTENT(in) :: nPlotPoints
  INTEGER, INTENT(in) :: nElem

    myGeom % qType = quadrature
    CALL myGeom % interp % Build( N = polyDegree, &
                                  controlNodeType = quadrature, &
                                  M = nPlotPoints, &
                                  targetNodeType = UNIFORM )

    CALL myGeom % x %  Build( N = polyDegree, &
                              nVar = 2, &
                              nElem = nElem ) 

    CALL myGeom % dxds %  Build( N = polyDegree, &
                                 nVar = 2, &
                                 nElem = nElem ) 

    CALL myGeom % dsdx %  Build( N = polyDegree, &
                                 nVar = 2, &
                                 nElem = nElem ) 

    CALL myGeom % J %  Build( N = polyDegree, &
                                 nVar = 2, &
                                 nElem = nElem ) 
    
END SUBROUTINE Build_SEMGeometry3D

SUBROUTINE Trash_SEMGeometry3D( myGeom ) 
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(inout) :: myGeom

    CALL myGeom % interp % Trash()
    CALL myGeom % x % Trash()
    CALL myGeom % dxds % Trash()
    CALL myGeom % dsdx % Trash()
    CALL myGeom % J % Trash()
  
END SUBROUTINE Trash_SEMGeometry3D 
#ifdef GPU
SUBROUTINE UpdateHost_SEMGeometry3D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(inout) :: myGeom
 
    CALL myGeom % interp % UpdateHost()
    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % J % UpdateHost()

END SUBROUTINE UpdateHost_SEMGeometry3D

SUBROUTINE UpdateDevice_SEMGeometry3D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(inout) :: myGeom
 
    CALL myGeom % interp % UpdateDevice()
    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMGeometry3D
#endif

END MODULE SEMMesh
