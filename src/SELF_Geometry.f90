!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Geometry

  USE SELF_Constants
  USE SELF_Lagrange
  USE SELF_Data
  USE SELF_SupportRoutines
  USE SELF_Mesh

  IMPLICIT NONE

#include "SELF_Macros.h"
  TYPE,PUBLIC :: SEMGeometry
  INTEGER :: nElem
  END TYPE SEMGeometry

  TYPE,EXTENDS(SEMGeometry),PUBLIC :: Geometry1D
    TYPE(Scalar1D) :: x ! Physical Positions
    TYPE(Scalar1D) :: dxds ! Conversion from computational to physical space

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_Geometry1D
    PROCEDURE,PUBLIC :: Free => Free_Geometry1D
    PROCEDURE,PUBLIC :: GenerateFromMesh => GenerateFromMesh_Geometry1D
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_Geometry1D
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_Geometry1D
    PROCEDURE,PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_Geometry1D

    PROCEDURE :: Write => Write_Geometry1D

  END TYPE Geometry1D

  TYPE,EXTENDS(SEMGeometry),PUBLIC :: SEMQuad
    TYPE(Vector2D) :: x ! Physical positions
    TYPE(Tensor2D) :: dxds ! Covariant basis vectors
    TYPE(Tensor2D) :: dsdx ! Contavariant basis vectors
    TYPE(Vector2D) :: nHat ! Normal Vectors pointing across coordinate lines
    TYPE(Scalar2D) :: nScale ! Boundary scale
    TYPE(Scalar2D) :: J ! Jacobian of the transformation

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_SEMQuad
    PROCEDURE,PUBLIC :: Free => Free_SEMQuad
    PROCEDURE,PUBLIC :: GenerateFromMesh => GenerateFromMesh_SEMQuad
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_SEMQuad
    PROCEDURE,PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_SEMQuad
    PROCEDURE,PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_SEMQuad

    PROCEDURE,PUBLIC :: CovariantArcMin => CovariantArcMin_SEMQuad
    !PROCEDURE :: Write => Write_SEMQuad

  END TYPE SEMQuad

  TYPE,EXTENDS(SEMGeometry),PUBLIC :: SEMHex
    TYPE(Vector3D) :: x ! Physical positions
    TYPE(Tensor3D) :: dxds ! Covariant basis vectors
    TYPE(Tensor3D) :: dsdx ! Contavariant basis vectors
    TYPE(Vector3D) :: nHat ! Normal Vectors pointing across coordinate lines
    TYPE(Scalar3D) :: nScale ! Boundary scale
    TYPE(Scalar3D) :: J ! Jacobian of the transformation

  CONTAINS

    PROCEDURE,PUBLIC :: Init => Init_SEMHex
    PROCEDURE,PUBLIC :: Free => Free_SEMHex
    PROCEDURE,PUBLIC :: GenerateFromMesh => GenerateFromMesh_SEMHex
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_SEMHex
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_SEMHex
    PROCEDURE,PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_SEMHex
    PROCEDURE,PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_SEMHex
    PROCEDURE,PRIVATE :: CheckSides => CheckSides_SEMHex

    !PROCEDURE :: Write => Write_SEMHex

  END TYPE SEMHex

CONTAINS

  SUBROUTINE Init_Geometry1D(myGeom,interp,nElem)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(out) :: myGeom
    TYPE(Lagrange),POINTER,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nElem

    myGeom % nElem = nElem

    CALL myGeom % x % Init(interp=interp,&
                           nVar=1, &
                           nElem=nElem)

    CALL myGeom % dxds % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

  END SUBROUTINE Init_Geometry1D

  SUBROUTINE Free_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % Free()
    CALL myGeom % dxds % Free()

  END SUBROUTINE Free_Geometry1D

  SUBROUTINE UpdateHost_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()

  END SUBROUTINE UpdateHost_Geometry1D

  SUBROUTINE GenerateFromMesh_Geometry1D(myGeom,mesh)
    ! Generates the geometry for a 1-D mesh ( set of line segments )
    ! Assumes that mesh is using Gauss-Lobatto quadrature and the degree is given by mesh % nGeo
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom
    TYPE(Mesh1D),INTENT(in) :: mesh
    ! Local
    INTEGER :: iel,i,nid
    TYPE(Lagrange),TARGET :: meshToModel
    TYPE(Scalar1D) :: xMesh

    CALL meshToModel % Init(mesh % nGeo, mesh % quadrature,&
            myGeom % x % interp % N, &
            myGeom % x % interp % controlNodeType )

    CALL xMesh % Init(meshToModel,&
                      1,mesh % nElem)

    ! Set the element internal mesh locations
    nid = 1
    DO iel = 1,mesh % nElem
      DO i = 0,mesh % nGeo
        xMesh % interior(i,iel,1) = mesh % nodeCoords(nid)
        nid = nid + 1
      END DO
    END DO

    ! Interpolate from the mesh hopr_nodeCoords to the geometry (Possibly not gauss_lobatto quadrature)
    CALL xMesh % GridInterp(myGeom % x)

    CALL myGeom % x % BoundaryInterp()

    CALL myGeom % CalculateMetricTerms()

    CALL myGeom % UpdateDevice()

    CALL xMesh % Free()

    CALL meshToModel % Free() 

  END SUBROUTINE GenerateFromMesh_Geometry1D

  SUBROUTINE CalculateMetricTerms_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % Derivative(myGeom % dxds)
    CALL myGeom % dxds % BoundaryInterp()
    CALL myGeom % UpdateDevice()

  END SUBROUTINE CalculateMetricTerms_Geometry1D

  SUBROUTINE Write_Geometry1D(myGeom,fileName)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(in) :: myGeom
    CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    ! Local
    CHARACTER(LEN=self_FileNameLength) :: pickupFile

    IF( PRESENT(filename) )THEN
      pickupFile = filename
    ELSE
      pickupFile = 'mesh.h5'
    ENDIF

    CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

    CALL CreateGroup_HDF5(fileId,'/quadrature')

    CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
                         myGeom % x % interp % controlPoints)

    CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
                         myGeom % x % interp % qWeights)

    CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
                         myGeom % x % interp % dgMatrix)

    CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
                         myGeom % x % interp % dMatrix)

    CALL CreateGroup_HDF5(fileId,'/mesh')

    CALL CreateGroup_HDF5(fileId,'/mesh/interior')

    CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

    CALL WriteArray_HDF5(fileId,'/mesh/interior/x',myGeom % x % interior)

    CALL WriteArray_HDF5(fileId,'/mesh/interior/dxds',myGeom % dxds % interior)

    CALL WriteArray_HDF5(fileId,'/mesh/boundary/x',myGeom % x % boundary)

    CALL WriteArray_HDF5(fileId,'/mesh/boundary/dxds',myGeom % dxds % boundary)

    CALL Close_HDF5(fileId)

  END SUBROUTINE Write_Geometry1D

  SUBROUTINE Init_SEMQuad(myGeom,interp,nElem)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(out) :: myGeom
    TYPE(Lagrange),POINTER,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nElem

    myGeom % nElem = nElem

    CALL myGeom % x % Init(interp=interp,&
                           nVar=1, &
                           nElem=nElem)

    CALL myGeom % dxds % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % dsdx % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % nHat % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % nScale % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % J % Init(interp=interp,&
                           nVar=1, &
                           nElem=nElem)

  END SUBROUTINE Init_SEMQuad

  SUBROUTINE Free_SEMQuad(myGeom)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(inout) :: myGeom

    CALL myGeom % x % Free()
    CALL myGeom % dxds % Free()
    CALL myGeom % dsdx % Free()
    CALL myGeom % nHat % Free()
    CALL myGeom % nScale % Free()
    CALL myGeom % J % Free()

  END SUBROUTINE Free_SEMQuad

  SUBROUTINE UpdateDevice_SEMQuad(myGeom)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % nHat % UpdateDevice()
    CALL myGeom % nScale % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

  END SUBROUTINE UpdateDevice_SEMQuad

  SUBROUTINE GenerateFromMesh_SEMQuad(myGeom,mesh)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(inout) :: myGeom
    TYPE(Mesh2D),INTENT(in) :: mesh
    ! Local
    INTEGER :: iel
    INTEGER :: i,j,nid
    TYPE(Lagrange),TARGET :: meshToModel
    TYPE(Vector2D) :: xMesh

    CALL meshToModel % Init(mesh % nGeo, &
            mesh % quadrature, &
            myGeom % x % interp % N, &
            myGeom % x % interp % controlNodeType )

    CALL xMesh % Init(meshToModel,1,mesh % nElem)

    ! Set the element internal mesh locations
    DO iel = 1, mesh % nElem
      DO j = 0, mesh % nGeo
        DO i = 0, mesh % nGeo
          xMesh % interior(i,j,iel,1,1:2) = mesh % nodeCoords(1:2,i,j,iel)
        END DO
      END DO
    END DO

    CALL xMesh % GridInterp(myGeom % x)
    CALL myGeom % x % BoundaryInterp()
    CALL myGeom % CalculateMetricTerms()

    CALL myGeom % UpdateDevice()
    CALL xMesh % Free()
    CALL meshToModel % Free()

  END SUBROUTINE GenerateFromMesh_SEMQuad

  SUBROUTINE CalculateContravariantBasis_SEMQuad(myGeom)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(inout) :: myGeom
    ! Local
    INTEGER :: iEl,i,j,k
    REAL(prec) :: fac
    REAL(prec) :: mag

    ! Now calculate the contravariant basis vectors
    ! In this convention, dsdx(j,i) is contravariant vector i, component j
    ! To project onto contravariant vector i, dot vector along the first dimension
    DO iEl = 1,myGeom % nElem
      DO j = 0,1,myGeom % dxds % interp % N+1
        DO i = 0,1,myGeom % dxds % interp % N+1

          myGeom % dsdx % interior(i,j,iel,1,1,1) =  myGeom % dxds % interior(i,j,iel,1,2,2)
          myGeom % dsdx % interior(i,j,iel,1,2,1) = -myGeom % dxds % interior(i,j,iel,1,1,2)
          myGeom % dsdx % interior(i,j,iel,1,1,2) = -myGeom % dxds % interior(i,j,iel,1,2,1)
          myGeom % dsdx % interior(i,j,iel,1,2,2) =  myGeom % dxds % interior(i,j,iel,1,1,1)

        END DO
      END DO
    END DO

    ! Interpolate the contravariant tensor to the boundaries
    CALL myGeom % dsdx % BoundaryInterp()

    ! Now, modify the sign of dsdx so that
    ! myGeom % dsdx % boundary is equal to the outward pointing normal vector
    DO iEl = 1,myGeom % nElem
      DO k = 1,4
        DO i = 1,myGeom % J % interp % N+1
          IF (k == selfSide2D_East .OR. k == selfSide2D_North) THEN
            fac = SIGN(1.0_prec,myGeom % J % boundary(i,k,iEl,1))
          ELSE
            fac = -SIGN(1.0_prec,myGeom % J % boundary(i,k,iEl,1))
          END IF

          IF( k == 1 )THEN ! South

            mag = SQRT( myGeom % dsdx % boundary(i,k,iEl,1,1,1,2)**2 +&
                        myGeom % dsdx % boundary(i,k,iEl,1,1,2,2)**2 )
 
            myGeom % nScale % boundary(i,k,iEl,1) = mag

            myGeom % nHat % boundary(i,k,iEl,1,1:2) = &
              fac*myGeom % dsdx % boundary(i,k,iEl,1,1:2,2)/mag


          ELSEIF( k == 2 )THEN ! East

            mag = SQRT( myGeom % dsdx % boundary(i,k,iEl,1,1,1)**2 +&
                        myGeom % dsdx % boundary(i,k,iEl,1,2,1)**2 )
 
            myGeom % nScale % boundary(i,k,iEl,1) = mag

            myGeom % nHat % boundary(i,k,iEl,1,1:2) = &
              fac*myGeom % dsdx % boundary(i,k,iEl,1,1:2)/mag

          ELSEIF( k == 3 )THEN ! North

            mag = SQRT( myGeom % dsdx % boundary(i,k,iEl,1,1,2)**2 +&
                        myGeom % dsdx % boundary(i,k,iEl,1,2,2)**2 )
 
            myGeom % nScale % boundary(i,k,iEl,1) = mag

            myGeom % nHat % boundary(i,k,iEl,1,1:2) = &
              fac*myGeom % dsdx % boundary(i,k,iEl,1,1:2,2)/mag

          ELSEIF( k == 4 )THEN ! West

            mag = SQRT( myGeom % dsdx % boundary(i,k,iEl,1,1,1)**2 +&
                        myGeom % dsdx % boundary(i,k,iEl,1,2,1)**2 )
 
            myGeom % nScale % boundary(i,k,iEl,1) = mag

            myGeom % nHat % boundary(i,k,iEl,1:2) = &
              fac*myGeom % dsdx % boundary(i,k,iEl,1,1:2,1)/mag

          ENDIF

          ! Set the directionality for dsdx on the boundaries
          myGeom % dsdx % boundary(i,k,iEl,1,1:2,1:2) = & 
              myGeom % dsdx % boundary(i,k,iEl,1,1:2,1:2)*fac

        END DO
      END DO
    END DO

  END SUBROUTINE CalculateContravariantBasis_SEMQuad

  SUBROUTINE CalculateMetricTerms_SEMQuad(myGeom)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(inout) :: myGeom

    CALL myGeom % x % Gradient(myGeom % dxds)
    CALL myGeom % dxds % BoundaryInterp()
    CALL myGeom % dxds % Determinant(myGeom % J)
    CALL myGeom % J % BoundaryInterp()

    CALL myGeom % CalculateContravariantBasis()

    IF (GPUAvailable()) THEN
      CALL myGeom % dsdx % UpdateDevice()
      CALL myGeom % nHat % UpdateDevice()
      CALL myGeom % nScale % UpdateDevice()
    ENDIF

  END SUBROUTINE CalculateMetricTerms_SEMQuad

  ! FUNCTION CovariantArcMin_SEMQuad(myGeom) RESULT(dxMin)
  !   IMPLICIT NONE
  !   CLASS(SEMQuad) :: myGeom
  !   REAL(prec) :: dxMin
  !   ! Local
  !   INTEGER :: i, j, iEl, N
  !   REAL(prec) :: dx, dy
  !   REAL(prec) :: dxds(1:2,1:2)
  !   REAL(prec) :: ds(0:1,myGeom % dxds % interp % N+1,&
  !                    0:1,myGeom % dxds % interp % N+1,&
  !                    1:myGeom % nElem)

  !   N = 1,myGeom % dxds % interp % N+1
  !   DO iEl = 1,myGeom % nElem
  !     DO j = 0, N
  !       DO i = 0, N

  !         dxds =  myGeom % dxds % interior(1:2,1:2,i,j,iel,1)
  !         dx = SQRT(dxds(1,1)**2 + dxds(1,2)**2)
  !         dy = SQRT(dxds(2,1)**2 + dxds(2,2)**2)
  !         ds(i,j,iEl) = 2.0_prec*MIN(dx,dy)/(REAL(N,prec)**2)

  !       ENDDO
  !     ENDDO
  !   ENDDO

  !   dxMin = MINVAL(ds) 

  ! END FUNCTION CovariantArcMin_SEMQuad

  ! SUBROUTINE Write_SEMQuad(myGeom,fileName)
  !   IMPLICIT NONE
  !   CLASS(SEMQuad),INTENT(in) :: myGeom
  !   CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
  !   ! Local
  !   INTEGER(HID_T) :: fileId
  !   ! Local
  !   CHARACTER(LEN=self_FileNameLength) :: pickupFile

  !   IF( PRESENT(filename) )THEN
  !     pickupFile = filename
  !   ELSE
  !     pickupFile = 'mesh.h5'
  !   ENDIF

  !   CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

  !   CALL CreateGroup_HDF5(fileId,'/quadrature')

  !   CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
  !                        myGeom % x % interp % controlPoints)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
  !                        myGeom % x % interp % qWeights)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
  !                        myGeom % x % interp % dgMatrix)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
  !                        myGeom % x % interp % dMatrix)

  !   CALL CreateGroup_HDF5(fileId,'/mesh')

  !   CALL CreateGroup_HDF5(fileId,'/mesh/interior')

  !   CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/x',myGeom % x % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/dxds',myGeom % dxds % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/dsdx',myGeom % dsdx % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/J',myGeom % J % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/x',myGeom % x % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/dxds',myGeom % dxds % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/dsdx',myGeom % dsdx % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/nHat',myGeom % nHat % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/nScale',myGeom % nScale % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/J',myGeom % J % boundary)

  !   CALL Close_HDF5(fileId)

  ! END SUBROUTINE Write_SEMQuad

  SUBROUTINE Init_SEMHex(myGeom,interp,nElem)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(out) :: myGeom
    TYPE(Lagrange),POINTER,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nElem

    myGeom % nElem = nElem

    CALL myGeom % x % Init(interp=interp,&
                           nVar=1, &
                           nElem=nElem)

    CALL myGeom % dxds % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % dsdx % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % nHat % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % nScale % Init(interp=interp,&
                              nVar=1, &
                              nElem=nElem)

    CALL myGeom % J % Init(interp=interp,&
                           nVar=1, &
                           nElem=nElem)

  END SUBROUTINE Init_SEMHex

  SUBROUTINE Free_SEMHex(myGeom)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(inout) :: myGeom

    CALL myGeom % x % Free()
    CALL myGeom % dxds % Free()
    CALL myGeom % dsdx % Free()
    CALL myGeom % nHat % Free()
    CALL myGeom % nScale % Free()
    CALL myGeom % J % Free()

  END SUBROUTINE Free_SEMHex

  SUBROUTINE UpdateDevice_SEMHex(myGeom)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()
    CALL myGeom % dsdx % UpdateDevice()
    CALL myGeom % nHat % UpdateDevice()
    CALL myGeom % nScale % UpdateDevice()
    CALL myGeom % J % UpdateDevice()

  END SUBROUTINE UpdateDevice_SEMHex

  SUBROUTINE GenerateFromMesh_SEMHex(myGeom,mesh)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(inout) :: myGeom
    TYPE(Mesh3D),INTENT(in) :: mesh
    ! Local
    INTEGER :: iel
    INTEGER :: i,j,k,nid
    TYPE(Lagrange),TARGET :: meshToModel
    TYPE(Vector3D) :: xMesh

    CALL meshToModel % Init(mesh % nGeo, mesh % quadrature,&
            myGeom % x % interp % N, &
            myGeom % x % interp % controlNodeType )

    CALL xMesh % Init(meshToModel,&
                      1,mesh % nElem)

    ! Set the element internal mesh locations
    DO iel = 1,mesh % nElem
      DO k = 0,mesh % nGeo
        DO j = 0,mesh % nGeo
          DO i = 0,mesh % nGeo
            xMesh % interior(i,j,k,iel,1,1:3) = mesh % nodeCoords(1:3,i,j,k,iel)
          END DO
        END DO
      END DO
    END DO

    CALL xMesh % GridInterp(myGeom % x)
    CALL myGeom % x % BoundaryInterp()
    CALL myGeom % CalculateMetricTerms()

    CALL myGeom % UpdateDevice()

    CALL xMesh % Free()
    CALL meshToModel % Free()

  END SUBROUTINE GenerateFromMesh_SEMHex

  SUBROUTINE CheckSides_SEMHex(myGeom,mesh)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(in) :: myGeom
    TYPE(Mesh3D),INTENT(in) :: mesh
    ! 
    INTEGER :: e1, s1
    INTEGER :: e2, s2
    INTEGER :: i1, j1
    INTEGER :: i2, j2
    INTEGER :: flip, bcid
    REAL(prec) :: rms

      DO e1 = 1,mesh % nElem
        DO s1 = 1,6

          e2 = mesh % sideInfo(3,s1,e1)
          s2 = mesh % sideInfo(4,s1,e1)/10
          flip = mesh % sideInfo(4,s1,e1) - s2*10
          bcid = mesh % sideInfo(5,s1,e1)

          IF (bcid == 0) THEN ! Interior

            rms = 0.0_prec

            IF (flip == 0) THEN

                DO j1 = 1,myGeom % x % interp % N+1
                  DO i1 = 1,myGeom % x % interp % N+1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary(i1,j1,s1,e1,1,1)-&
                                 myGeom % x % boundary(i1,j1,s2,e2,1,1))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,2)-&
                                 myGeom % x % boundary(i1,j1,s2,e2,1,2))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,3)-&
                                 myGeom % x % boundary(i1,j1,s2,e2,1,3))**2 )
                  END DO
                END DO

            ELSEIF (flip == 1) THEN

                DO j1 = 1,myGeom % x % interp % N+1
                  DO i1 = 1,myGeom % x % interp % N+1
                    i2 = j1
                    j2 = myGeom % x % interp % N - i1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary(i1,j1,s1,e1,1,1)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,1))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,2)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,2))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,3)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,3))**2 )
                  END DO
                END DO

            ELSEIF (flip == 2) THEN

                DO j1 = 1,myGeom % x % interp % N+1
                  DO i1 = 1,myGeom % x % interp % N+1
                    i2 = myGeom % x % interp % N - i1
                    j2 = myGeom % x % interp % N - j1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary(i1,j1,s1,e1,1,1)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,1))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,2)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,2))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,3)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,3))**2 )
                  END DO
                END DO

            ELSEIF (flip == 3) THEN

                DO j1 = 1,myGeom % x % interp % N+1
                  DO i1 = 1,myGeom % x % interp % N+1
                    i2 = myGeom % x % interp % N - j1
                    j2 = i1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary(i1,j1,s1,e1,1,1)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,1))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,2)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,2))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,3)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,3))**2 )
                  END DO
                END DO

            ELSEIF (flip == 4) THEN

                DO j1 = 1,myGeom % x % interp % N+1
                  DO i1 = 1,myGeom % x % interp % N+1
                    i2 = j1
                    j2 = i1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary(i1,j1,s1,e1,1,1)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,1))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,2)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,2))**2+&
                                (myGeom % x % boundary(i1,j1,s1,e1,1,3)-&
                                 myGeom % x % boundary(i2,j2,s2,e2,1,3))**2 )
                  END DO
                END DO
            END IF

          END IF

        END DO
      END DO

  END SUBROUTINE CheckSides_SEMHex

  SUBROUTINE CalculateContravariantBasis_SEMHex(myGeom)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(inout) :: myGeom
    ! Local
    INTEGER :: iEl,i,j,k
    REAL(prec) :: fac
    REAL(prec) :: mag

    ! Now calculate the contravariant basis vectors
    ! In this convention, dsdx(j,i) is contravariant vector i, component j
    ! To project onto contravariant vector i, dot vector along the first dimension
    ! TO DO : Curl Invariant Form
    DO iEl = 1,myGeom % nElem
      DO k = 1,myGeom % dxds % interp % N+1
        DO j = 1,myGeom % dxds % interp % N+1
          DO i = 1,myGeom % dxds % interp % N+1

            ! Ja1
              myGeom % dsdx % interior(i,j,k,iel,1,1,1) = &
              myGeom % dxds % interior(i,j,k,iel,1,2,2)* &
              myGeom % dxds % interior(i,j,k,iel,1,3,3) - &
              myGeom % dxds % interior(i,j,k,iel,1,3,2)* &
              myGeom % dxds % interior(i,j,k,iel,1,2,3)

              myGeom % dsdx % interior(i,j,k,iel,1,2,1) = &
              myGeom % dxds % interior(i,j,k,iel,1,1,3)* &
              myGeom % dxds % interior(i,j,k,iel,1,3,2) - &
              myGeom % dxds % interior(i,j,k,iel,1,3,3)* &
              myGeom % dxds % interior(i,j,k,iel,1,1,2)

              myGeom % dsdx % interior(i,j,k,iel,1,3,1) = &
              myGeom % dxds % interior(i,j,k,iel,1,1,2)* &
              myGeom % dxds % interior(i,j,k,iel,1,2,3) - &
              myGeom % dxds % interior(i,j,k,iel,1,2,2)* &
              myGeom % dxds % interior(i,j,k,iel,1,1,3)

            ! Ja2
              myGeom % dsdx % interior(i,j,k,iel,1,1,2) = &
              myGeom % dxds % interior(i,j,k,iel,1,2,3)* &
              myGeom % dxds % interior(i,j,k,iel,1,3,1) - &
              myGeom % dxds % interior(i,j,k,iel,1,3,3)* &
              myGeom % dxds % interior(i,j,k,iel,1,2,1)

              myGeom % dsdx % interior(i,j,k,iel,1,2,2) = &
              myGeom % dxds % interior(i,j,k,iel,1,1,1)* &
              myGeom % dxds % interior(i,j,k,iel,1,3,3) - &
              myGeom % dxds % interior(i,j,k,iel,1,3,1)* &
              myGeom % dxds % interior(i,j,k,iel,1,1,3)

              myGeom % dsdx % interior(i,j,k,iel,1,3,2) = &
              myGeom % dxds % interior(i,j,k,iel,1,1,3)* &
              myGeom % dxds % interior(i,j,k,iel,1,2,1) - &
              myGeom % dxds % interior(i,j,k,iel,1,2,3)* &
              myGeom % dxds % interior(i,j,k,iel,1,1,1)

            ! Ja3
              myGeom % dsdx % interior(i,j,k,iel,1,1,3) = &
              myGeom % dxds % interior(i,j,k,iel,1,2,1)* &
              myGeom % dxds % interior(i,j,k,iel,1,3,2) - &
              myGeom % dxds % interior(i,j,k,iel,1,3,1)* &
              myGeom % dxds % interior(i,j,k,iel,1,2,2)

              myGeom % dsdx % interior(i,j,k,iel,1,2,3) = &
              myGeom % dxds % interior(i,j,k,iel,1,1,2)* &
              myGeom % dxds % interior(i,j,k,iel,1,3,1) - &
              myGeom % dxds % interior(i,j,k,iel,1,3,2)* &
              myGeom % dxds % interior(i,j,k,iel,1,1,1)

              myGeom % dsdx % interior(i,j,k,iel,1,3,3) = &
              myGeom % dxds % interior(i,j,k,iel,1,1,1)* &
              myGeom % dxds % interior(i,j,k,iel,1,2,2) - &
              myGeom % dxds % interior(i,j,k,iel,1,2,1)* &
              myGeom % dxds % interior(i,j,k,iel,1,1,2)

          END DO
        END DO
      END DO
    END DO

    ! Interpolate the contravariant tensor to the boundaries
    CALL myGeom % dsdx % BoundaryInterp(gpuAccel=.FALSE.)

    ! Now, calculate nHat (outward pointing normal)
    DO iEl = 1,myGeom % nElem
      DO k = 1,6
        DO j = 1,myGeom % J % interp % N+1
          DO i = 1,myGeom % J % interp % N+1
            IF (k == selfSide3D_Top .OR. k == selfSide3D_East .OR. k == selfSide3D_North) THEN
              fac = SIGN(1.0_prec,myGeom % J % boundary(i,j,k,iEl,1))
            ELSE
              fac = -SIGN(1.0_prec,myGeom % J % boundary(i,j,k,iEl,1))
            END IF

            IF( k == 1 )THEN ! Bottom

              mag = SQRT( myGeom % dsdx % boundary(i,j,k,iEl,1,1,3)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,2,3)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,3,3)**2 )
 
              myGeom % nScale % boundary(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom % dsdx % boundary(i,j,k,iEl,1,1:3,3)/mag

            ELSEIF( k == 2 )THEN  ! South

              mag = SQRT( myGeom % dsdx % boundary(i,j,k,iEl,1,1,2)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,2,2)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,3,2)**2 )

              myGeom % nScale % boundary(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom % dsdx % boundary(i,j,k,iEl,1,1:3,2)/mag

            ELSEIF( k == 3 )THEN  ! East

              mag = SQRT( myGeom % dsdx % boundary(i,j,k,iEl,1,1,1)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,2,1)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,3,1)**2 )

              myGeom % nScale % boundary(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom % dsdx % boundary(i,j,k,iEl,1,1:3,1)/mag

            ELSEIF( k == 4 )THEN  ! North

              mag = SQRT( myGeom % dsdx % boundary(i,j,k,iEl,1,1,2)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,2,2)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,3,2)**2 )

              myGeom % nScale % boundary(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom % dsdx % boundary(i,j,k,iEl,1,1:3,2)/mag

            ELSEIF( k == 5 )THEN  ! West

              mag = SQRT( myGeom % dsdx % boundary(i,j,k,iEl,1,1,1)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,2,1)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,3,1)**2 )

              myGeom % nScale % boundary(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom % dsdx % boundary(i,j,k,iEl,1,1:3,1)/mag

            ELSEIF( k == 6 )THEN  ! Top

              mag = SQRT( myGeom % dsdx % boundary(i,j,k,iEl,1,1,3)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,2,3)**2 +&
                          myGeom % dsdx % boundary(i,j,k,iEl,1,3,3)**2 )
 
              myGeom % nScale % boundary(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom % dsdx % boundary(i,j,k,iEl,1,1:3,3)/mag

            ENDIF

            ! Set the directionality for dsdx on the boundaries
            ! This is primarily used for DG gradient calculations,
            ! which do not use nHat for the boundary terms.
            myGeom % dsdx % boundary(i,j,k,iEl,1,1:3,1:3) = & 
                myGeom % dsdx % boundary(i,j,k,iEl,1,1:3,1:3)*fac

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE CalculateContravariantBasis_SEMHex

  SUBROUTINE CalculateMetricTerms_SEMHex(myGeom)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(inout) :: myGeom

    CALL myGeom % x % Gradient(myGeom % dxds)
    CALL myGeom % dxds % BoundaryInterp()
    CALL myGeom % dxds % Determinant(myGeom % J)
    CALL myGeom % J % BoundaryInterp()

    CALL myGeom % CalculateContravariantBasis()
    IF (GPUAvailable()) THEN
      CALL myGeom % UpdateDevice()
    ENDIF

  END SUBROUTINE CalculateMetricTerms_SEMHex

  ! SUBROUTINE Write_SEMHex(myGeom,fileName)
  !   IMPLICIT NONE
  !   CLASS(SEMHex),INTENT(in) :: myGeom
  !   CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
  !   ! Local
  !   INTEGER(HID_T) :: fileId
  !   ! Local
  !   CHARACTER(LEN=self_FileNameLength) :: pickupFile

  !   IF( PRESENT(filename) )THEN
  !     pickupFile = filename
  !   ELSE
  !     pickupFile = 'mesh.h5'
  !   ENDIF

  !   CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

  !   CALL CreateGroup_HDF5(fileId,'/quadrature')

  !   CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
  !                        myGeom % x % interp % controlPoints)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
  !                        myGeom % x % interp % qWeights)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
  !                        myGeom % x % interp % dgMatrix)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
  !                        myGeom % x % interp % dMatrix)

  !   CALL CreateGroup_HDF5(fileId,'/mesh')

  !   CALL CreateGroup_HDF5(fileId,'/mesh/interior')

  !   CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/x',myGeom % x % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/dxds',myGeom % dxds % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/dsdx',myGeom % dsdx % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/J',myGeom % J % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/x',myGeom % x % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/dxds',myGeom % dxds % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/dsdx',myGeom % dsdx % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/nHat',myGeom % nHat % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/nScale',myGeom % nScale % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/J',myGeom % J % boundary)

  !   CALL Close_HDF5(fileId)

  ! END SUBROUTINE Write_SEMHex

END MODULE SELF_Geometry
