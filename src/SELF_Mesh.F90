MODULE SELF_Mesh

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Data

USE hipfort
USE ISO_C_BINDING



IMPLICIT NONE

#include "SELF_Macros.h"
  
  INTEGER, PARAMETER :: selfMinNodalValence2D = 4
  INTEGER, PARAMETER :: selfMinNodalValence3D = 8
  INTEGER, PARAMETER :: selfMaxNodalValence2D = 6
  INTEGER, PARAMETER :: selfMaxNodalValence3D = 10


  TYPE, PUBLIC :: Geometry2D
    INTEGER :: qType ! Quadrature Type
    INTEGER :: nElem
    TYPE(Vector2D) :: x ! Physical positions
    TYPE(Tensor2D) :: dxds ! Covariant basis vectors
    TYPE(Tensor2D) :: dsdx ! Contavariant basis vectors
    TYPE(Scalar2D) :: J ! Jacobian of the transformation

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Geometry2D
      PROCEDURE, PUBLIC :: Trash => Trash_Geometry2D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Geometry2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Geometry2D
#endif
      PROCEDURE, PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_Geometry2D
      PROCEDURE, PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_Geometry2D

  END TYPE Geometry2D

  TYPE, PUBLIC :: Geometry3D
    INTEGER :: qType ! Quadrature Type
    INTEGER :: nElem
    TYPE(Vector3D) :: x ! Physical positions
    TYPE(Tensor3D) :: dxds ! Covariant basis vectors
    TYPE(Tensor3D) :: dsdx ! Contavariant basis vectors
    TYPE(Scalar3D) :: J ! Jacobian of the transformation

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_Geometry3D
      PROCEDURE, PUBLIC :: Trash => Trash_Geometry3D
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_Geometry3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_Geometry3D
#endif

      PROCEDURE, PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_Geometry3D
      PROCEDURE, PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_Geometry3D

  END TYPE Geometry3D

  TYPE, PUBLIC :: MeshNodes
    INTEGER :: objCount
    INTEGER :: ndim
    TYPE(hfInt32_r1) :: id
    TYPE(hfInt32_r1) :: flag
    TYPE(hfInt32_r1) :: eValence
    TYPE(hfInt32_r2) :: elements
    TYPE(hfInt32_r2) :: elementNodes
    TYPE(hfReal_r2) :: x

    CONTAINS

      PROCEDURE, PUBLIC :: Build => Build_MeshNodes
      PROCEDURE, PUBLIC :: Trash => Trash_MeshNodes
#ifdef GPU
      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_MeshNodes
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_MeshNodes
#endif

  END TYPE MeshNodes

  TYPE, PUBLIC :: Edges
    INTEGER :: objCount
    TYPE(hfInt32_r1) :: id
    TYPE(hfInt32_r1) :: eValence
    TYPE(hfInt32_r1) :: boundaryID
    TYPE(hfInt32_r2) :: nodes
    TYPE(hfInt32_r2) :: elements
    TYPE(hfInt32_r2) :: elementEdges

!    CONTAINS

  END TYPE Edges

  TYPE, PUBLIC :: Faces
    INTEGER :: N
    INTEGER :: objCount
    TYPE(hfInt32_r1) :: id
    TYPE(hfInt32_r1) :: boundaryID
    TYPE(hfInt32_r2) :: nodes
    TYPE(hfInt32_r2) :: elements
    TYPE(hfInt32_r2) :: elementFaces
    TYPE(hfInt32_r3) :: iMap
    TYPE(hfInt32_r3) :: jMap

!    CONTAINS

  END TYPE Faces

  TYPE, PUBLIC :: Elements
    INTEGER :: objCount
    INTEGER :: objType
    TYPE(hfInt32_r1) :: id
    TYPE(hfInt32_r2) :: nodes
    TYPE(hfInt32_r2) :: neighbors

!    CONTAINS

  END TYPE Elements

!  TYPE, PUBLIC :: MeshObjHandles
!    INTEGER :: nMeshNodes, nEdges, nFaces, nElements
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

  TYPE, PUBLIC :: Mesh2D
    TYPE( Geometry2D ) :: geometry
    TYPE( Elements ) :: elements
    TYPE( MeshNodes  ) :: nodes
    TYPE( Edges ) :: edges 
    !TYPE( DomainDecomposition ) :: decomp 
    INTEGER :: cornerMap(1:2,1:4)
    INTEGER :: sideMap(1:4)
    INTEGER :: edgeMap(1:2,1:4)

!    CONTAINS

  END TYPE Mesh2D

  TYPE, PUBLIC :: Mesh3D
    TYPE( Geometry3D ) :: geometry
    TYPE( Elements ) :: elements
    TYPE( MeshNodes  ) :: nodes
    TYPE( Edges ) :: edges 
    TYPE( Faces ) :: faces
    !TYPE( DomainDecomposition ) :: decomp 
    INTEGER :: cornerMap(1:3,1:8)
    INTEGER :: sideMap(1:6)
    INTEGER :: faceMap(1:4,1:6)
    INTEGER :: edgeMap(1:2,1:12)
    INTEGER :: edgeFaceMap(1:2,1:4)
  
!    CONTAINS

!      PROCEDURE, PUBLIC :: Build => Build_Mesh3D
!      PROCEDURE, PUBLIC :: Trash => Trash_Mesh3D
!
!      PROCEDURE, PUBLIC :: Read_CGNSMesh
!      PROCEDURE, PUBLIC :: Read_UCDMesh
!      PROCEDURE, PUBLIC :: Read_TrellisUCDMesh
!
!      PROCEDURE, PUBLIC :: GenerateConnectivity => GenerateConnectivity_Mesh3D

  END TYPE Mesh3D

CONTAINS

SUBROUTINE Build_Geometry2D( myGeom, quadrature, polyDegree, nPlotPoints, nElem ) 
  IMPLICIT NONE
  CLASS(Geometry2D), INTENT(out) :: myGeom
  INTEGER, INTENT(in) :: quadrature
  INTEGER, INTENT(in) :: polyDegree
  INTEGER, INTENT(in) :: nPlotPoints
  INTEGER, INTENT(in) :: nElem

    myGeom % qType = quadrature
    myGeom % nElem = nElem

    CALL myGeom % x %  Build( N = polyDegree, &
                              quadratureType = GAUSS, &
                              M = 2*polyDegree, &
                              targetNodeType = UNIFORM, &
                              nVar = 2, &
                              nElem = nElem ) 

    CALL myGeom % dxds %  Build( N = polyDegree, &
                                 quadratureType = GAUSS, &
                                 M = 2*polyDegree, &
                                 targetNodeType = UNIFORM, &
                                 nVar = 2, &
                                 nElem = nElem ) 

    CALL myGeom % dsdx %  Build( N = polyDegree, &
                                 quadratureType = GAUSS, &
                                 M = 2*polyDegree, &
                                 targetNodeType = UNIFORM, &
                                 nVar = 2, &
                                 nElem = nElem ) 

    CALL myGeom % J %  Build( N = polyDegree, &
                              quadratureType = GAUSS, &
                              M = 2*polyDegree, &
                              targetNodeType = UNIFORM, &
                              nVar = 2, &
                              nElem = nElem ) 
    
END SUBROUTINE Build_Geometry2D

SUBROUTINE Build_MeshNodes( nodes, nNodes, nDim, maxValence )
#undef __FUNC__
#define __FUNC__ "Build_MeshNodes"
  IMPLICIT NONE
  CLASS( MeshNodes ), INTENT(out) :: nodes
  INTEGER, INTENT(in) :: nNodes
  INTEGER, INTENT(in) :: nDim
  INTEGER, INTENT(in), OPTIONAL :: maxValence
  ! LOGICAL
  INTEGER :: eValence

    nodes % objCount = nNodes
    nodes % nDim = nDim

    ! set the element valence size
    IF( PRESENT(maxValence) )THEN

      INFO("User provided optional maxValence")
      IF( nDim == 2 )THEN

        IF(maxValence >= selfMinNodalValence2D)THEN
          eValence = maxValence
        ELSE
          ERROR("Max nodal valence is too small") 
          STOP
        ENDIF

      ELSE IF( nDim == 3 )THEN

        IF(maxValence >= selfMinNodalValence3D)THEN
          eValence = maxValence
        ELSE
          ERROR("Max nodal valence is too small") 
          STOP
        ENDIF

      ENDIF

    ELSE

      INFO("Setting nodal valence to SELF default")
      IF( nDim == 2 )THEN
        eValence = selfMaxNodalValence2D
      ELSE IF( nDim == 3 )THEN
        eValence = selfMaxNodalValence2D
      ENDIF

    ENDIF

    ! Allocate space
    CALL nodes % id % Alloc( loBound = 1, &
                             upBound = nNodes )

    CALL nodes % flag % Alloc( loBound = 1, &
                               upBound = nNodes )

    CALL nodes % eValence % Alloc( loBound = 1, &
                                   upBound = nNodes )

    CALL nodes % elements % Alloc( loBound = (/1, 1/), &
                                   upBound = (/eValence, nNodes/) )

    CALL nodes % elementNodes % Alloc( loBound = (/1, 1/), &
                                       upBound = (/eValence, nNodes/) )

    CALL nodes % x % Alloc( loBound = (/1, 1/), &
                            upBound = (/nDim, nNodes/) )

END SUBROUTINE Build_MeshNodes

SUBROUTINE Trash_Geometry2D( myGeom ) 
  IMPLICIT NONE
  CLASS(Geometry2D), INTENT(inout) :: myGeom

    CALL myGeom % x % Trash()
    CALL myGeom % dxds % Trash()
    CALL myGeom % dsdx % Trash()
    CALL myGeom % J % Trash()
  
END SUBROUTINE Trash_Geometry2D 

SUBROUTINE Trash_MeshNodes( nodes )
#undef __FUNC__
#define __FUNC__ "Trash_MeshNodes"
  IMPLICIT NONE
  CLASS( MeshNodes ), INTENT(inout) :: nodes

    CALL nodes % id % Free()
    CALL nodes % flag % Free()
    CALL nodes % eValence % Free()
    CALL nodes % elements % Free()
    CALL nodes % elementNodes % Free()
    CALL nodes % x % Free()

END SUBROUTINE Trash_MeshNodes

#ifdef GPU
SUBROUTINE UpdateHost_Geometry2D( myGeom )
  IMPLICIT NONE
  CLASS(Geometry2D), INTENT(inout) :: myGeom
 
    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % J % UpdateHost()

END SUBROUTINE UpdateHost_Geometry2D

SUBROUTINE UpdateHost_MeshNodes( nodes )
#undef __FUNC__
#define __FUNC__ "UpdateHost_MeshNodes"
  IMPLICIT NONE
  CLASS( MeshNodes ), INTENT(inout) :: nodes

    CALL nodes % id % UpdateHost()
    CALL nodes % flag % UpdateHost()
    CALL nodes % eValence % UpdateHost()
    CALL nodes % elements % UpdateHost()
    CALL nodes % elementNodes % UpdateHost()
    CALL nodes % x % UpdateHost()

END SUBROUTINE UpdateHost_MeshNodes

SUBROUTINE UpdateDevice_Geometry2D( myGeom )
  IMPLICIT NONE
  CLASS(Geometry2D), INTENT(inout) :: myGeom
 
    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

END SUBROUTINE UpdateDevice_Geometry2D

SUBROUTINE UpdateDevice_MeshNodes( nodes )
#undef __FUNC__
#define __FUNC__ "UpdateDevice_MeshNodes"
  IMPLICIT NONE
  CLASS( MeshNodes ), INTENT(inout) :: nodes

    CALL nodes % id % UpdateDevice()
    CALL nodes % flag % UpdateDevice()
    CALL nodes % eValence % UpdateDevice()
    CALL nodes % elements % UpdateDevice()
    CALL nodes % elementNodes % UpdateDevice()
    CALL nodes % x % UpdateDevice()

END SUBROUTINE UpdateDevice_MeshNodes
#endif

SUBROUTINE CalculateContravariantBasis_Geometry2D( myGeom )
  IMPLICIT NONE
  CLASS(Geometry2D), INTENT(inout) :: myGeom
  ! Local
  INTEGER :: iEl, i, j, k
  REAL(prec) :: fac


    ! Now calculate the contravariant basis vectors
    DO iEl = 1, myGeom % nElem
      DO j = 0, myGeom % dxds % N
        DO i = 0, myGeom % dxds % N

          myGeom % dsdx % interior % hostData(1,1,i,j,1,iEl) = myGeom % dxds % interior % hostData(2,2,i,j,1,iEl)
          myGeom % dsdx % interior % hostData(2,1,i,j,1,iEl) = -myGeom % dxds % interior % hostData(1,2,i,j,1,iEl)
          myGeom % dsdx % interior % hostData(1,2,i,j,1,iEl) = -myGeom % dxds % interior % hostData(2,1,i,j,1,iEl)
          myGeom % dsdx % interior % hostData(2,2,i,j,1,iEl) = myGeom % dxds % interior % hostData(1,1,i,j,1,iEl)

        ENDDO
      ENDDO
    ENDDO

    ! Interpolate the contravariant tensor to the boundaries
    myGeom % dsdx = myGeom % dsdx % BoundaryInterp( )

    ! Now, modify the sign of dsdx so that
    ! myGeom % dsdx % boundary is equal to the outward pointing normal vector
    DO iEl = 1, myGeom % nElem
      DO k = 1, 4
        DO i = 0, myGeom % J % N
          IF( k == TOP .OR. k == EAST .OR. k == NORTH )THEN
            fac = SIGN(1.0_prec, myGeom % J % boundary % hostData(i,1,k,iEl))
          ELSE
            fac = -SIGN(1.0_prec, myGeom % J % boundary % hostData(i,1,k,iEl))
          ENDIF

          myGeom % dsdx % boundary % hostData(1:2,1:2,i,1,k,iEl) = fac*myGeom % dsdx % boundary % hostData(1:2,1:2,i,1,k,iEl)
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE CalculateContravariantBasis_Geometry2D

SUBROUTINE CalculateMetricTerms_Geometry2D( myGeom )
  IMPLICIT NONE
  CLASS(Geometry2D), INTENT(inout) :: myGeom

    myGeom % dxds = myGeom % x % Gradient( )
    myGeom % dxds = myGeom % dxds % BoundaryInterp( )

    ! Calculate the Jacobian = determinant of the covariant matrix at each point
    myGeom % J = myGeom % dxds % Determinant()
    myGeom % J = myGeom % J % BoundaryInterp( )

    CALL myGeom % CalculateContravariantBasis()

    CALL myGeom % UpdateDevice()

END SUBROUTINE CalculateMetricTerms_Geometry2D

SUBROUTINE Build_Geometry3D( myGeom, quadrature, polyDegree, nPlotPoints, nElem ) 
  IMPLICIT NONE
  CLASS(Geometry3D), INTENT(out) :: myGeom
  INTEGER, INTENT(in) :: quadrature
  INTEGER, INTENT(in) :: polyDegree
  INTEGER, INTENT(in) :: nPlotPoints
  INTEGER, INTENT(in) :: nElem

    myGeom % qType = quadrature
    myGeom % nElem = nElem

    CALL myGeom % x %  Build( N = polyDegree, &
                              quadratureType = GAUSS, &
                              M = 2*polyDegree, &
                              targetNodeType = UNIFORM, &
                              nVar = 2, &
                              nElem = nElem ) 

    CALL myGeom % dxds %  Build( N = polyDegree, &
                                 quadratureType = GAUSS, &
                                 M = 2*polyDegree, &
                                 targetNodeType = UNIFORM, &
                                 nVar = 2, &
                                 nElem = nElem ) 

    CALL myGeom % dsdx %  Build( N = polyDegree, &
                                 quadratureType = GAUSS, &
                                 M = 2*polyDegree, &
                                 targetNodeType = UNIFORM, &
                                 nVar = 2, &
                                 nElem = nElem ) 

    CALL myGeom % J %  Build( N = polyDegree, &
                              quadratureType = GAUSS, &
                              M = 2*polyDegree, &
                              targetNodeType = UNIFORM, &
                              nVar = 2, &
                              nElem = nElem ) 
    
END SUBROUTINE Build_Geometry3D

SUBROUTINE Trash_Geometry3D( myGeom ) 
  IMPLICIT NONE
  CLASS(Geometry3D), INTENT(inout) :: myGeom

    CALL myGeom % x % Trash()
    CALL myGeom % dxds % Trash()
    CALL myGeom % dsdx % Trash()
    CALL myGeom % J % Trash()
  
END SUBROUTINE Trash_Geometry3D 
#ifdef GPU
SUBROUTINE UpdateHost_Geometry3D( myGeom )
  IMPLICIT NONE
  CLASS(Geometry3D), INTENT(inout) :: myGeom
 
    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % J % UpdateHost()

END SUBROUTINE UpdateHost_Geometry3D

SUBROUTINE UpdateDevice_Geometry3D( myGeom )
  IMPLICIT NONE
  CLASS(Geometry3D), INTENT(inout) :: myGeom
 
    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

END SUBROUTINE UpdateDevice_Geometry3D
#endif

SUBROUTINE CalculateContravariantBasis_Geometry3D( myGeom )
  IMPLICIT NONE
  CLASS(Geometry3D), INTENT(inout) :: myGeom
  ! Local
  INTEGER :: iEl, i, j, k
  REAL(prec) :: fac


    ! Now calculate the contravariant basis vectors
    DO iEl = 1, myGeom % nElem
      DO k = 0, myGeom % dxds % N
        DO j = 0, myGeom % dxds % N
          DO i = 0, myGeom % dxds % N

            myGeom % dsdx % interior % hostData(1,1,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(2,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(3,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(2,1,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(3,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(1,3,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(1,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(3,3,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(3,1,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(1,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(2,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(1,3,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(1,2,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(3,1,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(3,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(2,1,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(2,2,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(3,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(1,1,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(1,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(3,1,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(3,2,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(1,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(2,1,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(2,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(1,1,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(1,3,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(2,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(3,2,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(3,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(2,2,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(2,3,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(3,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(1,2,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(1,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(3,2,i,j,k,1,iEl)

            myGeom % dsdx % interior % hostData(3,3,i,j,k,1,iEl) = myGeom % dxds % interior % hostData(1,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(2,2,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior % hostData(2,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior % hostData(1,2,i,j,k,1,iEl)


          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Interpolate the contravariant tensor to the boundaries
    myGeom % dsdx = myGeom % dsdx % BoundaryInterp( )

    ! Now, modify the sign of dsdx so that
    ! myGeom % dsdx % boundary is equal to the outward pointing normal vector
    DO iEl = 1, myGeom % nElem
      DO k = 1, 6
        DO j = 0, myGeom % J % N
          DO i = 0, myGeom % J % N
            IF( k == TOP .OR. k == EAST .OR. k == NORTH )THEN
              fac = SIGN(1.0_prec, myGeom % J % boundary % hostData(i,j,1,k,iEl))
            ELSE
              fac = -SIGN(1.0_prec, myGeom % J % boundary % hostData(i,j,1,k,iEl))
            ENDIF

            myGeom % dsdx % boundary % hostData(1:3,1:3,i,j,1,k,iEl) = fac*myGeom % dsdx % boundary % hostData(1:3,1:3,i,j,1,k,iEl)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE CalculateContravariantBasis_Geometry3D

SUBROUTINE CalculateMetricTerms_Geometry3D( myGeom )
  IMPLICIT NONE
  CLASS(Geometry3D), INTENT(inout) :: myGeom

    myGeom % dxds = myGeom % x % Gradient()
    myGeom % dxds = myGeom % dxds % BoundaryInterp()

    ! Calculate the Jacobian = determinant of the covariant matrix at each point
    myGeom % J = myGeom % dxds % Determinant()
    myGeom % J = myGeom % J % BoundaryInterp()

    CALL myGeom % CalculateContravariantBasis()

    CALL myGeom % UpdateDevice()

END SUBROUTINE CalculateMetricTerms_Geometry3D

END MODULE SELF_Mesh
