MODULE SELF_Tests

  USE SELF_Constants
  USE SELF_Memory
  USE SELF_SupportRoutines
  USE SELF_Quadrature
  USE SELF_Lagrange
  USE SELF_Data
  USE SELF_Mesh
  USE SELF_MappedData
  USE FEQParse

  IMPLICIT NONE

#include "SELF_Macros.h"

CONTAINS

  SUBROUTINE BlockMesh1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,tolerance,error)
#undef __FUNC__
#define __FUNC__ "BlockMesh1D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh1D) :: mesh
    TYPE(Geometry1D) :: geometry
    REAL(prec) :: expect_dxds,dxds_error,expect_boundx,boundx_error
    INTEGER :: iel,i

    error = 0
    INFO('Number of elements : '//Int2Str(nElem))
    INFO('Control point degree : '//Int2Str(cqDegree))
    INFO('Target point degree : '//Int2Str(tqDegree))
    INFO('Error tolerance : '//Float2Str(tolerance))

    CALL mesh % UniformBlockMesh(cqDegree,nElem, (/0.0_prec,1.0_prec/))

    ! Create the geometry
    CALL geometry % GenerateFromMesh(mesh,cqType,tqType,cqDegree,tqDegree)

    ! Verify the mesh
    expect_dxds = (1.0_prec/REAL(nElem,prec))/2.0_prec

    ! Calculate error in metric terms
    dxds_error = 0.0_prec
    DO iel = 1,nElem
      DO i = 0,cqDegree
        dxds_error = MAX(dxds_error,ABS(geometry % dxds % interior % hostData(i,1,iel) - expect_dxds))
      END DO
    END DO

    ! Calculate error in boundary interpolation
    boundx_error = 0.0_prec
    DO iel = 1,nElem
      i = mesh % eleminfo % hostData(3,iel)
      expect_boundx = mesh % nodeCoords % hostData(i)
      boundx_error = MAX(boundx_error,ABS(geometry % x % boundary % hostData(1,1,iel) - expect_boundx))

      i = mesh % eleminfo % hostData(4,iel)
      expect_boundx = mesh % nodeCoords % hostData(i)
      boundx_error = MAX(boundx_error,ABS(geometry % x % boundary % hostData(1,2,iel) - expect_boundx))
    END DO

    CALL mesh % Free()
    CALL geometry % Free()

    msg = "Max dx/ds error : "//Float2Str(dxds_error)
    IF (dxds_error > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] Metric Terms Test")
    ELSE
      INFO(TRIM(msg))
      INFO("Max dx/ds error : "//Float2Str(dxds_error))
      INFO("[PASS] Metric Terms Test Pass")
    END IF

    msg = "Max boundx error : "//Float2Str(boundx_error)
    IF (boundx_error > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] Boundary Interpolation Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] Boundary Interpolation Test")
    END IF

  END SUBROUTINE BlockMesh1D_Test

  SUBROUTINE BlockMesh2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,tolerance,error)
#undef __FUNC__
#define __FUNC__ "BlockMesh2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: mesh
    TYPE(SEMQuad) :: geometry
    REAL(prec) :: expect_dxds(1:2,1:2),dxds_error(1:2,1:2)
    REAL(prec) :: expect_J,J_error
    INTEGER :: iel,jel,i,j
    INTEGER :: row,col

    error = 0
    INFO('Number of elements : '//Int2Str(nElem*nElem))
    INFO('Control point degree : '//Int2Str(cqDegree))
    INFO('Target point degree : '//Int2Str(tqDegree))
    INFO('Error tolerance : '//Float2Str(tolerance))

    CALL mesh % UniformBlockMesh(cqDegree, &
                                 (/nElem,nElem/), &
                                 (/0.0_prec,1.0_prec, &
                                   0.0_prec,1.0_prec/))

    ! Create the geometry
    CALL geometry % GenerateFromMesh(mesh,cqType,tqType,cqDegree,tqDegree)

    ! Verify the mesh
    expect_dxds(1,1) = (1.0_prec/REAL(nElem,prec))/2.0_prec
    expect_dxds(1,2) = 0.0_prec
    expect_dxds(2,1) = 0.0_prec
    expect_dxds(2,2) = (1.0_prec/REAL(nElem,prec))/2.0_prec

    expect_J = expect_dxds(1,1)*expect_dxds(2,2)

    ! Calculate error in metric terms
    dxds_error = 0.0_prec
    DO iel = 1,nElem
      DO j = 0,cqDegree
        DO i = 0,cqDegree
          DO col = 1,2
            DO row = 1,2
              dxds_error(row,col) = MAX(dxds_error(row,col), &
                                        ABS(geometry % dxds % &
                                            interior % hostData(row,col,i,j,1,iel) - &
                                            expect_dxds(row,col)))
            END DO
          END DO
          J_error = MAX(J_error,ABS(geometry % J % &
                                    interior % hostData(i,j,1,iel) - &
                                    expect_J))
        END DO
      END DO
    END DO

    CALL mesh % Free()
    CALL geometry % Free()

    DO col = 1,2
      DO row = 1,2
        msg = "Max dx/ds error ("// &
              TRIM(Int2Str(row))//","// &
              TRIM(Int2Str(col))//") : "// &
              Float2Str(dxds_error(row,col))
        IF (dxds_error(row,col) > tolerance) THEN
          error = error + 1
          ERROR(TRIM(msg))
          ERROR("[FAIL] Covariant Tensor Test")
        ELSE
          INFO(TRIM(msg))
          INFO("[PASS] Covariant Tensor Test")
        END IF
      END DO
    END DO

    IF (J_error > tolerance) THEN
      error = error + 1
      ERROR("Max Jacobian error : "//Float2Str(J_error))
      ERROR("[FAIL] Jacobian Test")
    ELSE
      INFO("Max Jacobian error : "//Float2Str(J_error))
      INFO("[PASS] Jacobian Test")
    END IF

  END SUBROUTINE BlockMesh2D_Test

  SUBROUTINE BlockMesh3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,tolerance,error)
#undef __FUNC__
#define __FUNC__ "BlockMesh3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: mesh
    TYPE(SEMHex) :: geometry
    REAL(prec) :: expect_dxds(1:3,1:3),dxds_error(1:3,1:3)
    REAL(prec) :: expect_J,J_error
    INTEGER :: iel,jel,kel,i,j,k
    INTEGER :: row,col

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem*nElem*nElem)
    INFO(TRIM(msg))
    INFO('Control point degree : '//Int2Str(cqDegree))
    INFO('Target point degree : '//Int2Str(tqDegree))
    INFO('Error tolerance : '//Float2Str(tolerance))

    CALL mesh % UniformBlockMesh(cqDegree, &
                                 (/nElem,nElem,nElem/), &
                                 (/0.0_prec,1.0_prec, &
                                   0.0_prec,1.0_prec, &
                                   0.0_prec,1.0_prec/))

    ! Create the geometry
    CALL geometry % GenerateFromMesh(mesh,cqType,tqType,cqDegree,tqDegree)

    ! Verify the mesh
    expect_dxds(1,1) = (1.0_prec/REAL(nElem,prec))/2.0_prec
    expect_dxds(1,2) = 0.0_prec
    expect_dxds(1,3) = 0.0_prec
    expect_dxds(2,1) = 0.0_prec
    expect_dxds(2,2) = (1.0_prec/REAL(nElem,prec))/2.0_prec
    expect_dxds(2,3) = 0.0_prec
    expect_dxds(3,1) = 0.0_prec
    expect_dxds(3,2) = 0.0_prec
    expect_dxds(3,3) = (1.0_prec/REAL(nElem,prec))/2.0_prec

    expect_J = expect_dxds(1,1)*expect_dxds(2,2)*expect_dxds(3,3)

    ! Calculate error in metric terms
    dxds_error = 0.0_prec
    DO iel = 1,nElem
      DO k = 0,cqDegree
        DO j = 0,cqDegree
          DO i = 0,cqDegree
            DO col = 1,3
              DO row = 1,3
                dxds_error(row,col) = MAX(dxds_error(row,col), &
                                          ABS(geometry % dxds % &
                                              interior % hostData(row,col,i,j,k,1,iel) - &
                                              expect_dxds(row,col)))
              END DO
            END DO
            J_error = MAX(J_error,ABS(geometry % J % &
                                      interior % hostData(i,j,k,1,iel) - &
                                      expect_J))
          END DO
        END DO
      END DO
    END DO

    CALL mesh % Free()
    CALL geometry % Free()

    DO col = 1,3
      DO row = 1,3
        msg = "Max dx/ds error ("// &
              TRIM(Int2Str(row))//","// &
              TRIM(Int2Str(col))//") : "// &
              Float2Str(dxds_error(row,col))
        IF (dxds_error(row,col) > tolerance) THEN
          error = error + 1
          ERROR(TRIM(msg))
          ERROR("[FAIL] Covariant Tensor Test")
        ELSE
          INFO(TRIM(msg))
          INFO("[PASS] Covariant Tensor Test")
        END IF
      END DO
    END DO

    IF (J_error > tolerance) THEN
      error = error + 1
      ERROR("Max Jacobian error : "//Float2Str(J_error))
      ERROR("[FAIL] Jacobian Test")
    ELSE
      INFO("Max Jacobian error : "//Float2Str(J_error))
      INFO("[PASS] Jacobian Test")
    END IF

  END SUBROUTINE BlockMesh3D_Test

  SUBROUTINE ScalarInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,functionChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarInterp1D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: functionChar
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh1D) :: controlMesh,targetMesh
    TYPE(Geometry1D) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar1D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: iel,i,ivar

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    INFO('Error tolerance : '//Float2Str(tolerance))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree,nElem,(/0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the target mesh and geometry
    CALL targetMesh % UniformBlockMesh(tqDegree,nElem,(/0.0_prec,1.0_prec/))
    CALL targetGeometry % GenerateFromMesh(targetMesh,tqType,tqType,tqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)
    CALL fInterp % Init(tqDegree,tqType,tqDegree,tqType,nvar,nElem)
    CALL fActual % Init(tqDegree,tqType,tqDegree,tqType,nvar,nElem)
    CALL fError % Init(tqDegree,tqType,tqDegree,tqType,nvar,nElem)

    ! Create the equation parser object
    feq = EquationParser(functionChar, (/'x'/))

    ! Load the control function
     DO iel = 1, controlGeometry % nElem
       DO ivar = 1, nvar
         DO i = 0, cqDegree
           f % interior % hostData(i,ivar,iel) = &
             feq % Evaluate( (/controlGeometry % x % interior % hostData(i,1,iel)/) )
         ENDDO
       ENDDO
     ENDDO

    ! Load the target function
     DO iel = 1, targetGeometry % nElem
       DO ivar = 1, nvar
         DO i = 0, tqDegree
           fActual % interior % hostData(i,ivar,iel) = &
             feq % Evaluate( (/targetGeometry % x % interior % hostData(i,1,iel)/) )
         ENDDO
       ENDDO
     ENDDO

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,.FALSE.)
    fError = fActual - fInterp

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxInterior( )

    msg = "Max ScalarGridInterp_1D Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] ScalarGridInterp_1D Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] ScalarGridInterp_1D Test")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE ScalarInterp1D_Test

  SUBROUTINE ScalarBoundaryInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,functionChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarBoundaryInterp1D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: functionChar
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh1D) :: controlMesh
    TYPE(Geometry1D) :: controlGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar1D) :: f,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:2)
    INTEGER :: iel,i,ivar,iSide

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    msg = 'Error tolerance : '//Float2Str(tolerance)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree,nElem,(/0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)
    CALL fActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)
    CALL fError % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)

    ! Create the equation parser object
    feq = EquationParser(functionChar, (/'x'/))

    ! Load the control function
     DO iel = 1, controlGeometry % nElem
       DO ivar = 1, nvar
         DO i = 0, cqDegree
           f % interior % hostData(i,ivar,iel) = &
             feq % Evaluate( (/controlGeometry % x % interior % hostData(i,1,iel)/) )
         ENDDO
         ! Right Boundary
         fActual % boundary % hostData(ivar,1,iel) = &
             feq % Evaluate( (/controlGeometry % x % boundary % hostData(1,1,iel)/) )
         ! Right boundary
         fActual % boundary % hostData(ivar,2,iel) = &
             feq % Evaluate( (/controlGeometry % x % boundary % hostData(1,1,iel)/) )
       ENDDO
     ENDDO

    ! Run the grid interpolation
    CALL f % BoundaryInterp(.FALSE.)
    fError = fActual - f

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxBoundary( )

    DO iSide = 1,2
      msg = "Max ScalarBoundaryInterp_1D Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "[FAIL] ScalarBoundaryInterp_1D Test"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "[PASS] ScalarBoundaryInterp_1D Test"
        INFO(TRIM(msg))
      END IF
    ENDDO

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE ScalarBoundaryInterp1D_Test

  SUBROUTINE ScalarDerivative1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,fChar,dfChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarDerivative1D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: fChar
    CHARACTER(*),INTENT(in) :: dfChar
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh1D) :: controlMesh
    TYPE(Geometry1D) :: controlGeometry
    TYPE(EquationParser)  :: feq,dfeq
    TYPE(MappedScalar1D) :: f,dfInterp,dfActual,dfError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: iel,i,ivar

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    msg = 'Error tolerance : '//Float2Str(tolerance)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree,nElem,(/0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)
    CALL dfInterp % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)
    CALL dfActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)
    CALL dfError % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)

    ! Create the equation parser object
    feq = EquationParser(fChar, (/'x'/))
    dfeq = EquationParser(dfChar, (/'x'/))

    ! Load the control function
     DO iel = 1, controlGeometry % nElem
       DO ivar = 1, nvar
         DO i = 0, cqDegree
           f % interior % hostData(i,ivar,iel) = &
             feq % Evaluate( (/controlGeometry % x % interior % hostData(i,1,iel)/) )
           dfActual % interior % hostData(i,ivar,iel) = &
             dfeq % Evaluate( (/controlGeometry % x % interior % hostData(i,1,iel)/) )
         ENDDO
       ENDDO
     ENDDO

    ! Run the grid interpolation
    CALL f % Derivative(controlGeometry,dfInterp,.FALSE.)
    dfError = dfActual - dfInterp

    ! Calculate Absolute Maximum Error
    maxErrors = dfError % AbsMaxInterior( )

    msg = "Max ScalarDerivative_1D Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] ScalarDerivative_1D Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] ScalarDerivative_1D Test")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL dfInterp % Free()
    CALL dfActual % Free()
    CALL dfError % Free()

  END SUBROUTINE ScalarDerivative1D_Test

  SUBROUTINE ScalarInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,functionChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarInterp2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: functionChar
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh,targetMesh
    TYPE(SEMQuad) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar2D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar

    nel = nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    INFO('Error tolerance : '//Float2Str(tolerance))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the target mesh and geometry
    CALL targetMesh % UniformBlockMesh(tqDegree, &
                                       (/nElem,nElem/), &
                                       (/0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec/))
    CALL targetGeometry % GenerateFromMesh(targetMesh,tqType,tqType,tqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fInterp % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    feq = EquationParser(functionChar, (/'x','y'/))

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
             f % interior % hostData(i,j,ivar,iel) = &
               feq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Load the target function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, tqDegree
          DO i = 0, tqDegree
           fActual % interior % hostData(i,j,ivar,iel) = &
             feq % Evaluate( targetGeometry % x % interior % hostData(1:2,i,j,1,iel) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,.FALSE.)
    fError = fActual - fInterp

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxInterior( )

    msg = "Max ScalarGridInterp_2D Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] ScalarGridInterp_2D Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] ScalarGridInterp_2D Test")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE ScalarInterp2D_Test

  SUBROUTINE ScalarBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,functionChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarBoundaryInterp2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: functionChar
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar2D) :: f,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:4)
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,iside

    nel = nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    msg = 'Error tolerance : '//Float2Str(tolerance)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    feq = EquationParser(functionChar, (/'x','y'/))

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
             f % interior % hostData(i,j,ivar,iel) = &
               feq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
          ENDDO
          DO iside = 1,4
            fActual % boundary % hostData(j,ivar,iside,iel) = &
               feq % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % BoundaryInterp(.FALSE.)
    fError = fActual - f

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxBoundary( )

    DO iside = 1,4
      msg = "Max ScalarBoundaryInterp_2D Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "[FAIL] ScalarBoundaryInterp_2D Test"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "[PASS] ScalarBoundaryInterp_2D Test"
        INFO(TRIM(msg))
      END IF
    ENDDO

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE ScalarBoundaryInterp2D_Test

  SUBROUTINE VectorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vx,vy,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorInterp2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: vx,vy
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh,targetMesh
    TYPE(SEMQuad) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: vxeq,vyeq
    TYPE(Vector2D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar

    nel = nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    INFO('Error tolerance : '//Float2Str(tolerance))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the target mesh and geometry
    CALL targetMesh % UniformBlockMesh(tqDegree, &
                                       (/nElem,nElem/), &
                                       (/0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec/))
    CALL targetGeometry % GenerateFromMesh(targetMesh,tqType,tqType,tqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fInterp % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    vxeq = EquationParser(vx, (/'x','y'/))
    vyeq = EquationParser(vy, (/'x','y'/))

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
             f % interior % hostData(1,i,j,ivar,iel) = &
               vxeq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
             f % interior % hostData(2,i,j,ivar,iel) = &
               vyeq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Load the target function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, tqDegree
          DO i = 0, tqDegree
           fActual % interior % hostData(1,i,j,ivar,iel) = &
             vxeq % Evaluate( targetGeometry % x % interior % hostData(1:2,i,j,1,iel) )
           fActual % interior % hostData(2,i,j,ivar,iel) = &
             vyeq % Evaluate( targetGeometry % x % interior % hostData(1:2,i,j,1,iel) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,.FALSE.)
    fError = fActual - fInterp

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxInterior( )

    msg = "Max VectorGridInterp_2D Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] VectorGridInterp_2D Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] VectorGridInterp_2D Test")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE VectorInterp2D_Test

  SUBROUTINE VectorBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vx,vy,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorBoundaryInterp2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: vx,vy
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: vxeq,vyeq
    TYPE(Vector2D) :: f,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:4)
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,iside

    nel = nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    msg = 'Error tolerance : '//Float2Str(tolerance)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    vxeq = EquationParser(vx, (/'x','y'/))
    vyeq = EquationParser(vy, (/'x','y'/))

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
             f % interior % hostData(1,i,j,ivar,iel) = &
               vxeq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
             f % interior % hostData(2,i,j,ivar,iel) = &
               vyeq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
          ENDDO
          DO iside = 1,4
           fActual % boundary % hostData(1,j,ivar,iside,iel) = &
             vxeq % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
           fActual % boundary % hostData(2,j,ivar,iside,iel) = &
             vyeq % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % BoundaryInterp(.FALSE.)
    fError = fActual - f

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxBoundary( )

    DO iside = 1,4
      msg = "Max VectorBoundaryInterp_2D Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "[FAIL] VectorBoundaryInterp_2D Test"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "[PASS] VectorBoundaryInterp_2D Test"
        INFO(TRIM(msg))
      END IF
    ENDDO

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE VectorBoundaryInterp2D_Test

  SUBROUTINE TensorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,tensorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "TensorInterp2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: tensorChar(1:2,1:2)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh,targetMesh
    TYPE(SEMQuad) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: tensorEq(1:2,1:2)
    TYPE(Tensor2D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar
    INTEGER :: row,col

    nel = nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    INFO('Error tolerance : '//Float2Str(tolerance))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the target mesh and geometry
    CALL targetMesh % UniformBlockMesh(tqDegree, &
                                       (/nElem,nElem/), &
                                       (/0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec/))
    CALL targetGeometry % GenerateFromMesh(targetMesh,tqType,tqType,tqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fInterp % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    DO col = 1,2
      DO row = 1,2
        tensorEq(row,col) = EquationParser(TRIM(tensorChar(row,col)), (/'x','y'/))
      ENDDO
    ENDDO

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
            DO col = 1,2
              DO row = 1,2
                f % interior % hostData(row,col,i,j,ivar,iel) = &
                  tensorEq(row,col) % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Load the target function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, tqDegree
          DO i = 0, tqDegree
            DO col = 1,2
              DO row = 1,2
               fActual % interior % hostData(row,col,i,j,ivar,iel) = &
                 tensorEq(row,col) % Evaluate( targetGeometry % x % interior % hostData(1:2,i,j,1,iel) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,.FALSE.)
    fError = fActual - fInterp

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxInterior( )

    msg = "Max TensorGridInterp_2D Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] TensorGridInterp_2D Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] TensorGridInterp_2D Test")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE TensorInterp2D_Test

  SUBROUTINE TensorBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,tensorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "TensorBoundaryInterp2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: tensorChar(1:2,1:2)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: tensorEq(1:2,1:2)
    TYPE(Tensor2D) :: f,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:4)
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,iside
    INTEGER :: row,col

    nel = nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    msg = 'Error tolerance : '//Float2Str(tolerance)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    DO col = 1,2
      DO row = 1,2
        tensorEq(row,col) = EquationParser(TRIM(tensorChar(row,col)), (/'x','y'/))
      ENDDO
    ENDDO

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
            DO col = 1,2
              DO row = 1,2
                f % interior % hostData(row,col,i,j,ivar,iel) = &
                  tensorEq(row,col) % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
              ENDDO
            ENDDO
          ENDDO
          DO iside = 1,4
            DO col = 1,2
              DO row = 1,2
               fActual % boundary % hostData(row,col,j,ivar,iside,iel) = &
                 tensorEq(row,col) % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % BoundaryInterp(.FALSE.)
    fError = fActual - f

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxBoundary( )

    DO iside = 1,4
      msg = "Max TensorBoundaryInterp_2D Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "[FAIL] TensorBoundaryInterp_2D Test"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "[PASS] TensorBoundaryInterp_2D Test"
        INFO(TRIM(msg))
      END IF
    ENDDO

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE TensorBoundaryInterp2D_Test

  SUBROUTINE ScalarInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,functionChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarInterp3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: functionChar
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar3D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar

    nel = nElem*nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    INFO('Error tolerance : '//Float2Str(tolerance))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the target mesh and geometry
    CALL targetMesh % UniformBlockMesh(tqDegree, &
                                       (/nElem,nElem,nElem/), &
                                       (/0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec/))
    CALL targetGeometry % GenerateFromMesh(targetMesh,tqType,tqType,tqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fInterp % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    feq = EquationParser(functionChar, (/'x','y','z'/))

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree
               f % interior % hostData(i,j,k,ivar,iel) = &
                 feq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Load the target function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, tqDegree
          DO j = 0, tqDegree
            DO i = 0, tqDegree
             fActual % interior % hostData(i,j,k,ivar,iel) = &
               feq % Evaluate( targetGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,.FALSE.)
    fError = fActual - fInterp

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxInterior( )

    msg = "Max GridInterp_3D Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] GridInterp_3D Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] GridInterp_3D Test")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE ScalarInterp3D_Test

  SUBROUTINE ScalarBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,functionChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarBoundaryInterp3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: functionChar
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar3D) :: f,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:6)
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,iside

    nel = nElem*nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    msg = 'Error tolerance : '//Float2Str(tolerance)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    feq = EquationParser(functionChar, (/'x','y','z'/))

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree
               f % interior % hostData(i,j,k,ivar,iel) = &
                 feq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
            ENDDO
            DO iside = 1,6
             fActual % boundary % hostData(j,k,ivar,iside,iel) = &
               feq % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % BoundaryInterp(.FALSE.)
    fError = fActual - f

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxBoundary( )

    DO iside = 1,6
      msg = "Max ScalarBoundaryInterp_3D Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "[FAIL] ScalarBoundaryInterp_3D Test"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "[PASS] ScalarBoundaryInterp_3D Test"
        INFO(TRIM(msg))
      END IF
    ENDDO

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE ScalarBoundaryInterp3D_Test

  SUBROUTINE VectorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vx,vy,vz,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorInterp3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: vx,vy,vz
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: vxeq,vyeq,vzeq
    TYPE(Vector3D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar

    nel = nElem*nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    INFO('Error tolerance : '//Float2Str(tolerance))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the target mesh and geometry
    CALL targetMesh % UniformBlockMesh(tqDegree, &
                                       (/nElem,nElem,nElem/), &
                                       (/0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec/))
    CALL targetGeometry % GenerateFromMesh(targetMesh,tqType,tqType,tqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fInterp % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    vxeq = EquationParser(vx, (/'x','y','z'/))
    vyeq = EquationParser(vy, (/'x','y','z'/))
    vzeq = EquationParser(vz, (/'x','y','z'/))

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree
               f % interior % hostData(1,i,j,k,ivar,iel) = &
                 vxeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
               f % interior % hostData(2,i,j,k,ivar,iel) = &
                 vyeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
               f % interior % hostData(3,i,j,k,ivar,iel) = &
                 vzeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Load the target function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, tqDegree
          DO j = 0, tqDegree
            DO i = 0, tqDegree
             fActual % interior % hostData(1,i,j,k,ivar,iel) = &
               vxeq % Evaluate( targetGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
             fActual % interior % hostData(2,i,j,k,ivar,iel) = &
               vyeq % Evaluate( targetGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
             fActual % interior % hostData(3,i,j,k,ivar,iel) = &
               vzeq % Evaluate( targetGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,.FALSE.)
    fError = fActual - fInterp

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxInterior( )

    msg = "Max GridInterp_3D Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] VectorGridInterp_3D Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] VectorGridInterp_3D Test")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE VectorInterp3D_Test

  SUBROUTINE VectorBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vx,vy,vz,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorBoundaryInterp3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: vx,vy,vz
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser)  :: vxeq,vyeq,vzeq
    TYPE(Vector3D) :: f,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:6)
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,iside

    nel = nElem*nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    msg = 'Error tolerance : '//Float2Str(tolerance)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    vxeq = EquationParser(vx, (/'x','y','z'/))
    vyeq = EquationParser(vy, (/'x','y','z'/))
    vzeq = EquationParser(vz, (/'x','y','z'/))

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree
              f % interior % hostData(1,i,j,k,ivar,iel) = &
                vxeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
              f % interior % hostData(2,i,j,k,ivar,iel) = &
                vyeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
              f % interior % hostData(3,i,j,k,ivar,iel) = &
                vzeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
            ENDDO
            DO iside = 1,6
              fActual % boundary % hostData(1,j,k,ivar,iside,iel) = &
                vxeq % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
              fActual % boundary % hostData(2,j,k,ivar,iside,iel) = &
                vyeq % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
              fActual % boundary % hostData(3,j,k,ivar,iside,iel) = &
                vzeq % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % BoundaryInterp(.FALSE.)
    fError = fActual - f

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxBoundary( )

    DO iside = 1,6
      msg = "Max VectorBoundaryInterp_3D Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "[FAIL] VectorBoundaryInterp_3D Test"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "[PASS] VectorBoundaryInterp_3D Test"
        INFO(TRIM(msg))
      END IF
    ENDDO

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE VectorBoundaryInterp3D_Test

  SUBROUTINE TensorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,tensorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "TensorInterp3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: tensorChar(1:3,1:3)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: tensorEq(1:3,1:3)
    TYPE(Tensor3D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar
    INTEGER :: row,col

    nel = nElem*nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    INFO('Error tolerance : '//Float2Str(tolerance))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the target mesh and geometry
    CALL targetMesh % UniformBlockMesh(tqDegree, &
                                       (/nElem,nElem,nElem/), &
                                       (/0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec, &
                                         0.0_prec,1.0_prec/))
    CALL targetGeometry % GenerateFromMesh(targetMesh,tqType,tqType,tqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fInterp % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(tqDegree,tqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    DO col = 1,3
      DO row = 1,3
        tensorEq(row,col) = EquationParser(TRIM(tensorChar(row,col)), (/'x','y','z'/))
      ENDDO
    ENDDO

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree
              DO col = 1,3
                DO row = 1,3
                  f % interior % hostData(row,col,i,j,k,ivar,iel) = &
                    tensorEq(row,col) % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Load the target function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, tqDegree
          DO j = 0, tqDegree
            DO i = 0, tqDegree
              DO col = 1,3
                DO row = 1,3
                 fActual % interior % hostData(row,col,i,j,k,ivar,iel) = &
                   tensorEq(row,col) % Evaluate( targetGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,.FALSE.)
    fError = fActual - fInterp

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxInterior( )

    msg = "Max TensorGridInterp_3D Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("[FAIL] TensorGridInterp_3D Test")
    ELSE
      INFO(TRIM(msg))
      INFO("[PASS] TensorGridInterp_3D Test")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE TensorInterp3D_Test

  SUBROUTINE TensorBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,tensorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "TensorBoundaryInterp3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: tensorChar(1:3,1:3)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser)  :: tensorEq(1:3,1:3)
    TYPE(Tensor3D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:6)
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,iside
    INTEGER :: row,col

    nel = nElem*nElem*nElem
    error = 0
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))
    msg = 'Error tolerance : '//Float2Str(tolerance)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)
    CALL fError % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

    ! Create the equation parser object
    DO col = 1,3
      DO row = 1,3
        tensorEq(row,col) = EquationParser(TRIM(tensorChar(row,col)), (/'x','y','z'/))
      ENDDO
    ENDDO

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree
              DO col = 1,3
                DO row = 1,3
                  f % interior % hostData(row,col,i,j,k,ivar,iel) = &
                    tensorEq(row,col) % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
                ENDDO
              ENDDO
            ENDDO
            DO iside = 1,6
              DO col = 1,3
                DO row = 1,3
                  fActual % boundary % hostData(row,col,j,k,ivar,iside,iel) = &
                    tensorEq(row,col) % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % BoundaryInterp(.FALSE.)
    fError = fActual - f

    ! Calculate Absolute Maximum Error
    maxErrors = fError % AbsMaxBoundary( )

    DO iside = 1,6
      msg = "Max TensorBoundaryInterp_3D Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "[FAIL] TensorBoundaryInterp_3D Test"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "[PASS] TensorBoundaryInterp_3D Test"
        INFO(TRIM(msg))
      END IF
    ENDDO

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL fActual % Free()
    CALL fError % Free()

  END SUBROUTINE TensorBoundaryInterp3D_Test

END MODULE SELF_Tests
