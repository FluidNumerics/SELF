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
    PROCEDURE,PUBLIC :: UpdateHost => UpdateHost_SEMQuad
    PROCEDURE,PUBLIC :: UpdateDevice => UpdateDevice_SEMQuad
    PROCEDURE,PUBLIC :: CalculateMetricTerms => CalculateMetricTerms_SEMQuad
    PROCEDURE,PRIVATE :: CalculateContravariantBasis => CalculateContravariantBasis_SEMQuad

    PROCEDURE,PUBLIC :: CovariantArcMin => CovariantArcMin_SEMQuad

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

  END TYPE SEMHex

CONTAINS

  SUBROUTINE Init_Geometry1D(myGeom,interp,nElem)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(out) :: myGeom
    TYPE(Lagrange),POINTER,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nElem

    myGeom % nElem = nElem

    CALL myGeom % x % Init(interp=interp,&
                           nElem=nElem,&
                           nVar=1)

    CALL myGeom % dxds % Init(interp=interp,&
                              nElem=nElem,&
                              nVar=1)

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

  SUBROUTINE UpdateDevice_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateDevice()
    CALL myGeom % dxds % UpdateDevice()

  END SUBROUTINE UpdateDevice_Geometry1D

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
        xMesh % interior % hostData(i,iel,1) = mesh % nodeCoords % hostData(nid)
        nid = nid + 1
      END DO
    END DO

    ! Interpolate from the mesh hopr_nodeCoords to the geometry (Possibly not gauss_lobatto quadrature)
    CALL xMesh % GridInterp(myGeom % x,.FALSE.)

    CALL myGeom % x % BoundaryInterp(gpuAccel=.FALSE.)

    CALL myGeom % CalculateMetricTerms()

    CALL myGeom % UpdateDevice()

    CALL xMesh % Free()

    CALL meshToModel % Free() 

  END SUBROUTINE GenerateFromMesh_Geometry1D

  SUBROUTINE CalculateMetricTerms_Geometry1D(myGeom)
    IMPLICIT NONE
    CLASS(Geometry1D),INTENT(inout) :: myGeom

    CALL myGeom % x % Derivative(myGeom % dxds,gpuAccel=.FALSE.)
    CALL myGeom % dxds % BoundaryInterp(gpuAccel=.FALSE.)
    CALL myGeom % UpdateDevice()

  END SUBROUTINE CalculateMetricTerms_Geometry1D

  SUBROUTINE Init_SEMQuad(myGeom,interp,nElem)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(out) :: myGeom
    TYPE(Lagrange),POINTER,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nElem

    myGeom % nElem = nElem

    CALL myGeom % x % Init(interp=interp,&
                           nElem=nElem,&
                           nVar=1)

    CALL myGeom % dxds % Init(interp=interp,&
                              nElem=nElem,&
                              nVar=1)

    CALL myGeom % dsdx % Init(interp=interp,&
                              nElem=nElem,&
                              nVar=1)

    CALL myGeom % nHat % Init(interp=interp,&
                              nElem=nElem,&
                              nVar=1)

    CALL myGeom % nScale % Init(interp=interp,&
                              nElem=nElem,&
                              nVar=1)

    CALL myGeom % J % Init(interp=interp,&
                           nElem=nElem,&
                              nVar=1)

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

  SUBROUTINE UpdateHost_SEMQuad(myGeom)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % nHat % UpdateHost()
    CALL myGeom % nScale % UpdateHost()
    CALL myGeom % J % UpdateHost()

  END SUBROUTINE UpdateHost_SEMQuad

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
          xMesh % interior % hostData(1:2,i,j,iel,1) = mesh % nodeCoords % hostData(1:2,i,j,iel)
        END DO
      END DO
    END DO

    CALL xMesh % GridInterp(myGeom % x,.FALSE.)
    CALL myGeom % x % BoundaryInterp(gpuAccel=.FALSE.)
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
      DO j = 0,myGeom % dxds % interp % N
        DO i = 0,myGeom % dxds % interp % N

          myGeom % dsdx % interior % hostData(1,1,i,j,iEl,1) = myGeom % dxds % interior % hostData(2,2,i,j,iEl,1)
          myGeom % dsdx % interior % hostData(2,1,i,j,iEl,1) = -myGeom % dxds % interior % hostData(1,2,i,j,iEl,1)
          myGeom % dsdx % interior % hostData(1,2,i,j,iEl,1) = -myGeom % dxds % interior % hostData(2,1,i,j,iEl,1)
          myGeom % dsdx % interior % hostData(2,2,i,j,iEl,1) = myGeom % dxds % interior % hostData(1,1,i,j,iEl,1)

        END DO
      END DO
    END DO

    ! Interpolate the contravariant tensor to the boundaries
    CALL myGeom % dsdx % BoundaryInterp(gpuAccel=.FALSE.)

    ! Now, modify the sign of dsdx so that
    ! myGeom % dsdx % boundary is equal to the outward pointing normal vector
    DO iEl = 1,myGeom % nElem
      DO k = 1,4
        DO i = 0,myGeom % J % interp % N
          IF (k == selfSide2D_East .OR. k == selfSide2D_North) THEN
            fac = SIGN(1.0_prec,myGeom % J % boundary % hostData(i,k,iEl,1))
          ELSE
            fac = -SIGN(1.0_prec,myGeom % J % boundary % hostData(i,k,iEl,1))
          END IF

          IF( k == 1 )THEN ! South

            mag = SQRT( myGeom % dsdx % boundary % hostData(1,2,i,k,iEl,1)**2 +&
                        myGeom % dsdx % boundary % hostData(2,2,i,k,iEl,1)**2 )
 
            myGeom % nScale % boundary % hostData(i,k,iEl,1) = mag

            myGeom % nHat % boundary % hostData(1:2,i,k,iEl,1) = &
              fac*myGeom % dsdx % boundary % hostData(1:2,2,i,k,iEl,1)/mag


          ELSEIF( k == 2 )THEN ! East

            mag = SQRT( myGeom % dsdx % boundary % hostData(1,1,i,k,iEl,1)**2 +&
                        myGeom % dsdx % boundary % hostData(2,1,i,k,iEl,1)**2 )
 
            myGeom % nScale % boundary % hostData(i,k,iEl,1) = mag

            myGeom % nHat % boundary % hostData(1:2,i,k,iEl,1) = &
              fac*myGeom % dsdx % boundary % hostData(1:2,1,i,k,iEl,1)/mag

          ELSEIF( k == 3 )THEN ! North

            mag = SQRT( myGeom % dsdx % boundary % hostData(1,2,i,k,iEl,1)**2 +&
                        myGeom % dsdx % boundary % hostData(2,2,i,k,iEl,1)**2 )
 
            myGeom % nScale % boundary % hostData(i,k,iEl,1) = mag

            myGeom % nHat % boundary % hostData(1:2,i,k,iEl,1) = &
              fac*myGeom % dsdx % boundary % hostData(1:2,2,i,k,iEl,1)/mag

          ELSEIF( k == 4 )THEN ! West

            mag = SQRT( myGeom % dsdx % boundary % hostData(1,1,i,k,iEl,1)**2 +&
                        myGeom % dsdx % boundary % hostData(2,1,i,k,iEl,1)**2 )
 
            myGeom % nScale % boundary % hostData(i,k,iEl,1) = mag

            myGeom % nHat % boundary % hostData(1:2,i,k,iEl,1) = &
              fac*myGeom % dsdx % boundary % hostData(1:2,1,i,k,iEl,1)/mag

          ENDIF

          ! Set the directionality for dsdx on the boundaries
          myGeom % dsdx % boundary % hostData(1:2,1:2,i,k,iEl,1) = & 
              myGeom % dsdx % boundary % hostData(1:2,1:2,i,k,iEl,1)*fac

        END DO
      END DO
    END DO

  END SUBROUTINE CalculateContravariantBasis_SEMQuad

  SUBROUTINE CalculateMetricTerms_SEMQuad(myGeom)
    IMPLICIT NONE
    CLASS(SEMQuad),INTENT(inout) :: myGeom

    CALL myGeom % x % Gradient(myGeom % dxds,gpuAccel=.FALSE.)
    CALL myGeom % dxds % BoundaryInterp(gpuAccel=.FALSE.)
    CALL myGeom % dxds % Determinant(myGeom % J,gpuAccel=.FALSE.)
    CALL myGeom % J % BoundaryInterp(gpuAccel=.FALSE.)

    CALL myGeom % CalculateContravariantBasis()
    IF (GPUAvailable()) THEN
      CALL myGeom % dsdx % UpdateDevice()
      CALL myGeom % nHat % UpdateDevice()
      CALL myGeom % nScale % UpdateDevice()
    ENDIF

  END SUBROUTINE CalculateMetricTerms_SEMQuad

  FUNCTION CovariantArcMin_SEMQuad(myGeom) RESULT(dxMin)
    IMPLICIT NONE
    CLASS(SEMQuad) :: myGeom
    REAL(prec) :: dxMin
    ! Local
    INTEGER :: i, j, iEl, N
    REAL(prec) :: dx, dy
    REAL(prec) :: dxds(1:2,1:2)
    REAL(prec) :: ds(0:myGeom % dxds % interp % N,&
                     0:myGeom % dxds % interp % N,&
                     1:myGeom % nElem)

    N = myGeom % dxds % interp % N
    DO iEl = 1,myGeom % nElem
      DO j = 0, N
        DO i = 0, N

          dxds =  myGeom % dxds % interior % hostData(1:2,1:2,i,j,iEl,1)
          dx = SQRT(dxds(1,1)**2 + dxds(1,2)**2)
          dy = SQRT(dxds(2,1)**2 + dxds(2,2)**2)
          ds(i,j,iEl) = 2.0_prec*MIN(dx,dy)/(REAL(N,prec)**2)

        ENDDO
      ENDDO
    ENDDO

    dxMin = MINVAL(ds) 

  END FUNCTION CovariantArcMin_SEMQuad

  SUBROUTINE Init_SEMHex(myGeom,interp,nElem)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(out) :: myGeom
    TYPE(Lagrange),POINTER,INTENT(in) :: interp
    INTEGER,INTENT(in) :: nElem

    myGeom % nElem = nElem

    CALL myGeom % x % Init(interp=interp,&
                           nElem=nElem, &
                           nVar=1)

    CALL myGeom % dxds % Init(interp=interp,&
                              nElem=nElem, &
                              nVar=1)

    CALL myGeom % dsdx % Init(interp=interp,&
                              nElem=nElem, &
                              nVar=1)

    CALL myGeom % nHat % Init(interp=interp,&
                              nElem=nElem, &
                              nVar=1)

    CALL myGeom % nScale % Init(interp=interp,&
                              nElem=nElem, &
                              nVar=1)

    CALL myGeom % J % Init(interp=interp,&
                           nElem=nElem, &
                              nVar=1)

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

  SUBROUTINE UpdateHost_SEMHex(myGeom)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(inout) :: myGeom

    CALL myGeom % x % UpdateHost()
    CALL myGeom % dxds % UpdateHost()
    CALL myGeom % dsdx % UpdateHost()
    CALL myGeom % nHat % UpdateHost()
    CALL myGeom % nScale % UpdateHost()
    CALL myGeom % J % UpdateHost()

  END SUBROUTINE UpdateHost_SEMHex

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
            xMesh % interior % hostData(1:3,i,j,k,iEl,1) = mesh % nodeCoords % hostData(1:3,i,j,k,iel)
          END DO
        END DO
      END DO
    END DO

    ! Interpolate from the mesh hopr_nodeCoords to the geometry (Possibly not gauss_lobatto quadrature)
    CALL xMesh % GridInterp(myGeom % x,.FALSE.)
    CALL myGeom % x % BoundaryInterp(gpuAccel=.FALSE.)
    CALL myGeom % CalculateMetricTerms()
!    CALL myGeom % CheckSides(mesh)

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

          e2 = mesh % sideInfo % hostData(3,s1,e1)
          s2 = mesh % sideInfo % hostData(4,s1,e1)/10
          flip = mesh % sideInfo % hostData(4,s1,e1) - s2*10
          bcid = mesh % sideInfo % hostData(5,s1,e1)

          IF (bcid == 0) THEN ! Interior

            rms = 0.0_prec

            IF (flip == 0) THEN

                DO j1 = 0,myGeom % x % interp % N
                  DO i1 = 0,myGeom % x % interp % N
                    rms = rms + &
                          sqrt( (myGeom % x % boundary % hostdata(1,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(1,i1,j1,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(2,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(2,i1,j1,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(3,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(3,i1,j1,1,s2,e2))**2 )
                  END DO
                END DO

            ELSEIF (flip == 1) THEN

                DO j1 = 0,myGeom % x % interp % N
                  DO i1 = 0,myGeom % x % interp % N
                    i2 = j1
                    j2 = myGeom % x % interp % N - i1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary % hostdata(1,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(1,i2,j2,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(2,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(2,i2,j2,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(3,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(3,i2,j2,1,s2,e2))**2 )
                  END DO
                END DO

            ELSEIF (flip == 2) THEN

                DO j1 = 0,myGeom % x % interp % N
                  DO i1 = 0,myGeom % x % interp % N
                    i2 = myGeom % x % interp % N - i1
                    j2 = myGeom % x % interp % N - j1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary % hostdata(1,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(1,i2,j2,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(2,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(2,i2,j2,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(3,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(3,i2,j2,1,s2,e2))**2 )
                  END DO
                END DO

            ELSEIF (flip == 3) THEN

                DO j1 = 0,myGeom % x % interp % N
                  DO i1 = 0,myGeom % x % interp % N
                    i2 = myGeom % x % interp % N - j1
                    j2 = i1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary % hostdata(1,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(1,i2,j2,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(2,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(2,i2,j2,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(3,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(3,i2,j2,1,s2,e2))**2 )
                  END DO
                END DO

            ELSEIF (flip == 4) THEN

                DO j1 = 0,myGeom % x % interp % N
                  DO i1 = 0,myGeom % x % interp % N
                    i2 = j1
                    j2 = i1
                    rms = rms + &
                          sqrt( (myGeom % x % boundary % hostdata(1,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(1,i2,j2,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(2,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(2,i2,j2,1,s2,e2))**2+&
                                (myGeom % x % boundary % hostdata(3,i1,j1,1,s1,e1)-&
                                 myGeom % x % boundary % hostdata(3,i2,j2,1,s2,e2))**2 )
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
      DO k = 0,myGeom % dxds % interp % N
        DO j = 0,myGeom % dxds % interp % N
          DO i = 0,myGeom % dxds % interp % N

            ! Ja1
            myGeom % dsdx % interior % hostData(1,1,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(2,2,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(3,3,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(3,2,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(2,3,i,j,k,iEl,1)

            myGeom % dsdx % interior % hostData(2,1,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(1,3,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(3,2,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(3,3,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(1,2,i,j,k,iEl,1)

            myGeom % dsdx % interior % hostData(3,1,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(1,2,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(2,3,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(2,2,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(1,3,i,j,k,iEl,1)

            ! Ja2
            myGeom % dsdx % interior % hostData(1,2,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(2,3,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(3,1,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(3,3,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(2,1,i,j,k,iEl,1)

            myGeom % dsdx % interior % hostData(2,2,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(1,1,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(3,3,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(3,1,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(1,3,i,j,k,iEl,1)

            myGeom % dsdx % interior % hostData(3,2,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(1,3,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(2,1,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(2,3,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(1,1,i,j,k,iEl,1)

            ! Ja3
            myGeom % dsdx % interior % hostData(1,3,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(2,1,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(3,2,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(3,1,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(2,2,i,j,k,iEl,1)

            myGeom % dsdx % interior % hostData(2,3,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(1,2,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(3,1,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(3,2,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(1,1,i,j,k,iEl,1)

            myGeom % dsdx % interior % hostData(3,3,i,j,k,iEl,1) = &
              myGeom % dxds % interior % hostData(1,1,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(2,2,i,j,k,iEl,1) - &
              myGeom % dxds % interior % hostData(2,1,i,j,k,iEl,1)* &
              myGeom % dxds % interior % hostData(1,2,i,j,k,iEl,1)

          END DO
        END DO
      END DO
    END DO

    ! Interpolate the contravariant tensor to the boundaries
    CALL myGeom % dsdx % BoundaryInterp(gpuAccel=.FALSE.)

    ! Now, calculate nHat (outward pointing normal)
    DO iEl = 1,myGeom % nElem
      DO k = 1,6
        DO j = 0,myGeom % J % interp % N
          DO i = 0,myGeom % J % interp % N
            IF (k == selfSide3D_Top .OR. k == selfSide3D_East .OR. k == selfSide3D_North) THEN
              fac = SIGN(1.0_prec,myGeom % J % boundary % hostData(i,j,k,iEl,1))
            ELSE
              fac = -SIGN(1.0_prec,myGeom % J % boundary % hostData(i,j,k,iEl,1))
            END IF

            IF( k == 1 )THEN ! Bottom

              mag = SQRT( myGeom % dsdx % boundary % hostData(1,3,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(2,3,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(3,3,i,j,k,iEl,1)**2 )
 
              myGeom % nScale % boundary % hostData(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary % hostData(1:3,i,j,k,iEl,1) = &
                fac*myGeom % dsdx % boundary % hostData(1:3,3,i,j,k,iEl,1)/mag

            ELSEIF( k == 2 )THEN  ! South

              mag = SQRT( myGeom % dsdx % boundary % hostData(1,2,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(2,2,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(3,2,i,j,k,iEl,1)**2 )

              myGeom % nScale % boundary % hostData(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary % hostData(1:3,i,j,k,iEl,1) = &
                fac*myGeom % dsdx % boundary % hostData(1:3,2,i,j,k,iEl,1)/mag

            ELSEIF( k == 3 )THEN  ! East

              mag = SQRT( myGeom % dsdx % boundary % hostData(1,1,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(2,1,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(3,1,i,j,k,iEl,1)**2 )

              myGeom % nScale % boundary % hostData(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary % hostData(1:3,i,j,k,iEl,1) = &
                fac*myGeom % dsdx % boundary % hostData(1:3,1,i,j,k,iEl,1)/mag

            ELSEIF( k == 4 )THEN  ! North

              mag = SQRT( myGeom % dsdx % boundary % hostData(1,2,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(2,2,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(3,2,i,j,k,iEl,1)**2 )

              myGeom % nScale % boundary % hostData(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary % hostData(1:3,i,j,k,iEl,1) = &
                fac*myGeom % dsdx % boundary % hostData(1:3,2,i,j,k,iEl,1)/mag

            ELSEIF( k == 5 )THEN  ! West

              mag = SQRT( myGeom % dsdx % boundary % hostData(1,1,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(2,1,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(3,1,i,j,k,iEl,1)**2 )

              myGeom % nScale % boundary % hostData(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary % hostData(1:3,i,j,k,iEl,1) = &
                fac*myGeom % dsdx % boundary % hostData(1:3,1,i,j,k,iEl,1)/mag

            ELSEIF( k == 6 )THEN  ! Top

              mag = SQRT( myGeom % dsdx % boundary % hostData(1,3,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(2,3,i,j,k,iEl,1)**2 +&
                          myGeom % dsdx % boundary % hostData(3,3,i,j,k,iEl,1)**2 )
 
              myGeom % nScale % boundary % hostData(i,j,k,iEl,1) = mag

              myGeom % nHat % boundary % hostData(1:3,i,j,k,iEl,1) = &
                fac*myGeom % dsdx % boundary % hostData(1:3,3,i,j,k,iEl,1)/mag

            ENDIF

            ! Set the directionality for dsdx on the boundaries
            ! This is primarily used for DG gradient calculations,
            ! which do not use nHat for the boundary terms.
            myGeom % dsdx % boundary % hostData(1:3,1:3,i,j,k,iEl,1) = & 
                myGeom % dsdx % boundary % hostData(1:3,1:3,i,j,k,iEl,1)*fac

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE CalculateContravariantBasis_SEMHex

  SUBROUTINE CalculateMetricTerms_SEMHex(myGeom)
    IMPLICIT NONE
    CLASS(SEMHex),INTENT(inout) :: myGeom


    CALL myGeom % x % Gradient(myGeom % dxds,gpuAccel=.FALSE.)
    CALL myGeom % dxds % BoundaryInterp(gpuAccel=.FALSE.)
    CALL myGeom % dxds % Determinant(myGeom % J,gpuAccel=.FALSE.)
    CALL myGeom % J % BoundaryInterp(gpuAccel=.FALSE.)

    CALL myGeom % CalculateContravariantBasis()
    IF (GPUAvailable()) THEN
      CALL myGeom % UpdateDevice()
    ENDIF

  END SUBROUTINE CalculateMetricTerms_SEMHex

END MODULE SELF_Geometry
