MODULE SEMMesh

USE ModelPrecision
USE ConstantsDictionary
USE NodalSEMData
USE Lagrange_Class

USE hipfort
USE ISO_C_BINDING

INCLUDE "SELF_Macros.h"


IMPLICIT NONE
  
  INTEGER, PARAMETER :: selfMinNodalValence2D = 4
  INTEGER, PARAMETER :: selfMinNodalValence3D = 8
  INTEGER, PARAMETER :: selfMaxNodalValence2D = 6
  INTEGER, PARAMETER :: selfMaxNodalValence3D = 10


  TYPE, PUBLIC :: SEMGeometry2D
    INTEGER :: qType ! Quadrature Type
    INTEGER :: nElem
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
      PROCEDURE, PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_SEMGeometry2D
      PROCEDURE, PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_SEMGeometry2D

  END TYPE SEMGeometry2D

  TYPE, PUBLIC :: SEMGeometry3D
    INTEGER :: qType ! Quadrature Type
    INTEGER :: nElem
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

      PROCEDURE, PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_SEMGeometry3D
      PROCEDURE, PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_SEMGeometry3D

  END TYPE SEMGeometry3D

  TYPE, PUBLIC :: MeshNodes
    INTEGER :: objCount
    INTEGER :: ndim
    INTEGER, POINTER :: id(:)
    INTEGER, POINTER :: flag(:)
    INTEGER, POINTER :: eValence(:)
    INTEGER, POINTER :: elements(:,:)
    INTEGER, POINTER :: elementNodes(:,:)
    REAL(prec), ALLOCATABLE :: x(:,:)

    TYPE(c_ptr) :: id_dev
    TYPE(c_ptr) :: flag_dev
    TYPE(c_ptr) :: eValence_dev
    TYPE(c_ptr) :: elements_dev
    TYPE(c_ptr) :: elementNodes_dev
    TYPE(c_ptr) :: x_dev

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

  TYPE, PUBLIC :: SEMMesh2D
    TYPE( SEMGeometry2D ) :: geometry
    TYPE( Elements ) :: elements
    TYPE( MeshNodes  ) :: nodes
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
    
END SUBROUTINE Build_SEMGeometry2D

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
    ALLOCATE( nodes % id(1:nNodes), &
              nodes % flag(1:nNodes), &
              nodes % eValence(1:nNodes), &
              nodes % elements(1:eValence,1:nNodes), &
              nodes % elementNodes(1:eValence,1:nNodes), &
              nodes % x(1:nDim,1:nNodes) )

    nodes % id = 0
    nodes % flag = 0
    nodes % eValence = 0
    nodes % elements = 0
    nodes % elementNodes = 0
    nodes % x = 0.0_prec

#ifdef GPU
    CALL hipCheck(hipMalloc(nodes % id_dev, SIZEOF(nodes % id)))
    CALL hipCheck(hipMalloc(nodes % flag_dev, SIZEOF(nodes % flag)))
    CALL hipCheck(hipMalloc(nodes % eValence_dev, SIZEOF(nodes % eValence)))
    CALL hipCheck(hipMalloc(nodes % elements_dev, SIZEOF(nodes % elements)))
    CALL hipCheck(hipMalloc(nodes % elementNodes_dev, SIZEOF(nodes % elementNodes)))
    CALL hipCheck(hipMalloc(nodes % x_dev, SIZEOF(nodes % x)))
#endif
  

END SUBROUTINE Build_MeshNodes

SUBROUTINE Trash_SEMGeometry2D( myGeom ) 
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(inout) :: myGeom

    CALL myGeom % x % Trash()
    CALL myGeom % dxds % Trash()
    CALL myGeom % dsdx % Trash()
    CALL myGeom % J % Trash()
  
END SUBROUTINE Trash_SEMGeometry2D 

SUBROUTINE Trash_MeshNodes( nodes )
#undef __FUNC__
#define __FUNC__ "Trash_MeshNodes"
  IMPLICIT NONE
  CLASS( MeshNodes ), INTENT(inout) :: nodes

    DEALLOCATE( nodes % id, &
                nodes % flag, &
                nodes % eValence, &
                nodes % elements, &
                nodes % elementNodes, &
                nodes % x )

#ifdef GPU
    CALL hipCheck(hipFree(nodes % id_dev))
    CALL hipCheck(hipFree(nodes % flag_dev))
    CALL hipCheck(hipFree(nodes % eValence_dev))
    CALL hipCheck(hipFree(nodes % elements_dev))
    CALL hipCheck(hipFree(nodes % elementNodes_dev))
    CALL hipCheck(hipFree(nodes % x_dev))
#endif

END SUBROUTINE Trash_MeshNodes

#ifdef GPU
SUBROUTINE UpdateHost_SEMGeometry2D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(inout) :: myGeom
 
    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % J % UpdateHost()

END SUBROUTINE UpdateHost_SEMGeometry2D

SUBROUTINE UpdateHost_MeshNodes( nodes )
#undef __FUNC__
#define __FUNC__ "UpdateHost_MeshNodes"
  IMPLICIT NONE
  CLASS( MeshNodes ), INTENT(inout) :: nodes

    CALL hipCheck(hipMemcpy(c_loc(nodes % id), &
                            nodes % id_dev, &
                            SIZEOF(nodes % id), &
                            hipMemcpyDeviceToHost))

    CALL hipCheck(hipMemcpy(c_loc(nodes % flag ), &
                            nodes % flag_dev, &
                            SIZEOF(nodes % flag), &
                            hipMemcpyDeviceToHost))

    CALL hipCheck(hipMemcpy(c_loc(nodes % eValence), &
                            nodes % eValence_dev, &
                            SIZEOF(nodes % eValence), &
                            hipMemcpyDeviceToHost))

    CALL hipCheck(hipMemcpy(c_loc(nodes % elements), &
                            nodes % elements_dev, &
                            SIZEOF(nodes % elements), &
                            hipMemcpyDeviceToHost))

    CALL hipCheck(hipMemcpy(c_loc(nodes % elementNodes), &
                            nodes % elementNodes_dev, &
                            SIZEOF(nodes % elementNodes), &
                            hipMemcpyDeviceToHost))

    CALL hipCheck(hipMemcpy(c_loc(nodes % x), &
                            nodes % x_dev, &
                            SIZEOF(nodes % x), &
                            hipMemcpyDeviceToHost))

END SUBROUTINE UpdateHost_MeshNodes

SUBROUTINE UpdateDevice_SEMGeometry2D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(inout) :: myGeom
 
    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMGeometry2D

SUBROUTINE UpdateDevice_MeshNodes( nodes )
#undef __FUNC__
#define __FUNC__ "UpdateDevice_MeshNodes"
  IMPLICIT NONE
  CLASS( MeshNodes ), INTENT(inout) :: nodes

    CALL hipCheck(hipMemcpy(nodes % id_dev, &
                            c_loc(nodes % id), &
                            SIZEOF(nodes % id), &
                            hipMemcpyHostToDevice))

    CALL hipCheck(hipMemcpy(nodes % flag_dev, &
                            c_loc(nodes % flag), &
                            SIZEOF(nodes % flag), &
                            hipMemcpyHostToDevice))

    CALL hipCheck(hipMemcpy(nodes % eValence_dev, &
                            c_loc(nodes % eValence), &
                            SIZEOF(nodes % eValence), &
                            hipMemcpyHostToDevice))

    CALL hipCheck(hipMemcpy(nodes % elements_dev, &
                            c_loc(nodes % elements), &
                            SIZEOF(nodes % elements), &
                            hipMemcpyHostToDevice))

    CALL hipCheck(hipMemcpy(nodes % elementNodes_dev, &
                            c_loc(nodes % elementNodes), &
                            SIZEOF(nodes % elementNodes), &
                            hipMemcpyHostToDevice))

    CALL hipCheck(hipMemcpy(nodes % x_dev, &
                            c_loc(nodes % x), &
                            SIZEOF(nodes % x), &
                            hipMemcpyHostToDevice))

END SUBROUTINE UpdateDevice_MeshNodes
#endif

SUBROUTINE CalculateContravariantBasis_SEMGeometry2D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(inout) :: myGeom
  ! Local
  INTEGER :: iEl, i, j, k
  REAL(prec) :: fac


    ! Now calculate the contravariant basis vectors
    DO iEl = 1, myGeom % nElem
      DO j = 0, myGeom % dxds % N
        DO i = 0, myGeom % dxds % N

          myGeom % dsdx % interior(1,1,i,j,1,iEl) = myGeom % dxds % interior(2,2,i,j,1,iEl)
          myGeom % dsdx % interior(2,1,i,j,1,iEl) = -myGeom % dxds % interior(1,2,i,j,1,iEl)
          myGeom % dsdx % interior(1,2,i,j,1,iEl) = -myGeom % dxds % interior(2,1,i,j,1,iEl)
          myGeom % dsdx % interior(2,2,i,j,1,iEl) = myGeom % dxds % interior(1,1,i,j,1,iEl)

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
            fac = SIGN(1.0_prec, myGeom % J % boundary(i,1,k,iEl))
          ELSE
            fac = -SIGN(1.0_prec, myGeom % J % boundary(i,1,k,iEl))
          ENDIF

          myGeom % dsdx % boundary(1:2,1:2,i,1,k,iEl) = fac*myGeom % dsdx % boundary(1:2,1:2,i,1,k,iEl)
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE CalculateContravariantBasis_SEMGeometry2D

SUBROUTINE CalculateMetricTerms_SEMGeometry2D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry2D), INTENT(inout) :: myGeom

    myGeom % dxds = myGeom % x % Gradient( )
    myGeom % dxds = myGeom % dxds % BoundaryInterp( )

    ! Calculate the Jacobian = determinant of the covariant matrix at each point
    myGeom % J = myGeom % dxds % Determinant()
    myGeom % J = myGeom % J % BoundaryInterp( )

    CALL myGeom % CalculateContravariantBasis()

    CALL myGeom % UpdateDevice()

END SUBROUTINE CalculateMetricTerms_SEMGeometry2D

SUBROUTINE Build_SEMGeometry3D( myGeom, quadrature, polyDegree, nPlotPoints, nElem ) 
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(out) :: myGeom
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
    
END SUBROUTINE Build_SEMGeometry3D

SUBROUTINE Trash_SEMGeometry3D( myGeom ) 
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(inout) :: myGeom

    CALL myGeom % x % Trash()
    CALL myGeom % dxds % Trash()
    CALL myGeom % dsdx % Trash()
    CALL myGeom % J % Trash()
  
END SUBROUTINE Trash_SEMGeometry3D 
#ifdef GPU
SUBROUTINE UpdateHost_SEMGeometry3D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(inout) :: myGeom
 
    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % J % UpdateHost()

END SUBROUTINE UpdateHost_SEMGeometry3D

SUBROUTINE UpdateDevice_SEMGeometry3D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(inout) :: myGeom
 
    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

END SUBROUTINE UpdateDevice_SEMGeometry3D
#endif

SUBROUTINE CalculateContravariantBasis_SEMGeometry3D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(inout) :: myGeom
  ! Local
  INTEGER :: iEl, i, j, k
  REAL(prec) :: fac


    ! Now calculate the contravariant basis vectors
    DO iEl = 1, myGeom % nElem
      DO k = 0, myGeom % dxds % N
        DO j = 0, myGeom % dxds % N
          DO i = 0, myGeom % dxds % N

            myGeom % dsdx % interior(1,1,i,j,k,1,iEl) = myGeom % dxds % interior(2,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(2,3,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(3,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(2,3,i,j,k,1,iEl)

            myGeom % dsdx % interior(2,1,i,j,k,1,iEl) = myGeom % dxds % interior(3,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(1,3,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(1,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(3,3,i,j,k,1,iEl)

            myGeom % dsdx % interior(3,1,i,j,k,1,iEl) = myGeom % dxds % interior(1,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(2,3,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(2,2,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(1,3,i,j,k,1,iEl)

            myGeom % dsdx % interior(1,2,i,j,k,1,iEl) = myGeom % dxds % interior(2,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(3,1,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(3,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(2,1,i,j,k,1,iEl)

            myGeom % dsdx % interior(2,2,i,j,k,1,iEl) = myGeom % dxds % interior(3,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(1,1,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(1,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(3,1,i,j,k,1,iEl)

            myGeom % dsdx % interior(3,2,i,j,k,1,iEl) = myGeom % dxds % interior(1,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(2,1,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(2,3,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(1,1,i,j,k,1,iEl)

            myGeom % dsdx % interior(1,3,i,j,k,1,iEl) = myGeom % dxds % interior(2,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(3,2,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(3,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(2,2,i,j,k,1,iEl)

            myGeom % dsdx % interior(2,3,i,j,k,1,iEl) = myGeom % dxds % interior(3,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(1,2,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(1,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(3,2,i,j,k,1,iEl)

            myGeom % dsdx % interior(3,3,i,j,k,1,iEl) = myGeom % dxds % interior(1,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(2,2,i,j,k,1,iEl)-&
                                             myGeom % dxds % interior(2,1,i,j,k,1,iEl)*&
                                             myGeom % dxds % interior(1,2,i,j,k,1,iEl)


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
              fac = SIGN(1.0_prec, myGeom % J % boundary(i,j,1,k,iEl))
            ELSE
              fac = -SIGN(1.0_prec, myGeom % J % boundary(i,j,1,k,iEl))
            ENDIF

            myGeom % dsdx % boundary(1:3,1:3,i,j,1,k,iEl) = fac*myGeom % dsdx % boundary(1:3,1:3,i,j,1,k,iEl)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE CalculateContravariantBasis_SEMGeometry3D

SUBROUTINE CalculateMetricTerms_SEMGeometry3D( myGeom )
  IMPLICIT NONE
  CLASS(SEMGeometry3D), INTENT(inout) :: myGeom

    myGeom % dxds = myGeom % x % Gradient()
    myGeom % dxds = myGeom % dxds % BoundaryInterp()

    ! Calculate the Jacobian = determinant of the covariant matrix at each point
    myGeom % J = myGeom % dxds % Determinant()
    myGeom % J = myGeom % J % BoundaryInterp()

    CALL myGeom % CalculateContravariantBasis()

    CALL myGeom % UpdateDevice()

END SUBROUTINE CalculateMetricTerms_SEMGeometry3D

END MODULE SEMMesh
