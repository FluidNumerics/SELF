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

    msg = "Numerical Error (dx/ds) : "//Float2Str(dxds_error)
    IF (dxds_error > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
    END IF

    msg = "Numerical Error (xBound) : "//Float2Str(boundx_error)
    IF (boundx_error > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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
        msg = "Numerical Error (dx/ds) ("// &
              TRIM(Int2Str(row))//","// &
              TRIM(Int2Str(col))//") : "// &
              Float2Str(dxds_error(row,col))
        IF (dxds_error(row,col) > tolerance) THEN
          error = error + 1
          ERROR(TRIM(msg))
          ERROR("Status : [FAIL]")
        ELSE
          INFO(TRIM(msg))
          INFO("Status : [PASS]")
        END IF
      END DO
    END DO

    msg = "Numerical Error (Jacobian) : "//Float2Str(J_error)
    IF (J_error > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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
    J_error = 0.0_prec
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
        msg = "Numerical Error (dx/ds) ("// &
              TRIM(Int2Str(row))//","// &
              TRIM(Int2Str(col))//") : "// &
              Float2Str(dxds_error(row,col))
        IF (dxds_error(row,col) > tolerance) THEN
          error = error + 1
          ERROR(TRIM(msg))
          ERROR("Status : [FAIL]")
        ELSE
          INFO(TRIM(msg))
          INFO("Status : [PASS]")
        END IF
      END DO
    END DO

    msg = "Numerical Error (Jacobian) : "//Float2Str(J_error)
    IF (J_error > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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
      msg = "Numerical Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "Status : [FAIL]"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "Status : [PASS]"
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

  SUBROUTINE ScalarDerivative1D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,fChar,dfChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarDerivative1D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: dForm
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
         ! Right Boundary
         f % boundary % hostData(ivar,1,iel) = &
             feq % Evaluate( (/controlGeometry % x % boundary % hostData(1,2,iel)/) )
         ! Right boundary
         f % boundary % hostData(ivar,2,iel) = &
             feq % Evaluate( (/controlGeometry % x % boundary % hostData(1,1,iel)/) )
       ENDDO
     ENDDO

    ! Run the grid interpolation
    CALL f % Derivative(controlGeometry,dfInterp,dForm,.FALSE.)
    dfError = dfActual - dfInterp

    ! Calculate Absolute Maximum Error
    maxErrors = dfError % AbsMaxInterior( )

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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
      msg = "Numerical Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "Status : [FAIL]"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "Status : [PASS]"
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

  SUBROUTINE ScalarGradient2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,fChar,gradientChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarGradient2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: dForm
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: fChar
    CHARACTER(240),INTENT(in) :: gradientChar(1:2)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: feq,gxeq,gyeq
    TYPE(MappedScalar2D) :: f
    TYPE(MappedTensor2D) :: workTensor
    TYPE(MappedVector2D) :: dfInterp,dfActual,dfError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: iel,i,j,ivar,iside

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem*nElem)
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
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL workTensor % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfInterp % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfError % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)

    ! Create the equation parser object
    feq = EquationParser(fChar, (/'x','y'/))
    gxeq = EquationParser(gradientChar(1), (/'x','y'/))
    gyeq = EquationParser(gradientChar(2), (/'x','y'/))

    ! Load the control function
    DO iel = 1, controlGeometry % nElem
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
            f % interior % hostData(i,j,ivar,iel) = &
              feq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )

            dfActual % interior % hostData(1,i,j,ivar,iel) = &
              gxeq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )

            dfActual % interior % hostData(2,i,j,ivar,iel) = &
              gyeq % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
          ENDDO
          DO iside = 1,4
            f % boundary % hostData(j,ivar,iside,iel) = &
               feq % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % Gradient(workTensor,controlGeometry,dfInterp,dForm,.FALSE.)
    dfError = dfActual - dfInterp

    ! Calculate Absolute Maximum Error
    maxErrors = dfError % AbsMaxInterior( )

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workTensor % Free()
    CALL dfInterp % Free()
    CALL dfActual % Free()
    CALL dfError % Free()

  END SUBROUTINE ScalarGradient2D_Test

  SUBROUTINE VectorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vectorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorInterp2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:2)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh,targetMesh
    TYPE(SEMQuad) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: vEq(1:2)
    TYPE(Vector2D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,idir

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
    DO idir = 1,2
      vEq(idir) = EquationParser(vectorChar(idir), (/'x','y'/))
    ENDDO

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
            DO idir=1,2
              f % interior % hostData(idir,i,j,ivar,iel) = &
                vEq(idir) % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
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
            DO idir = 1,2
              fActual % interior % hostData(idir,i,j,ivar,iel) = &
                vEq(idir) % Evaluate( targetGeometry % x % interior % hostData(1:2,i,j,1,iel) )
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

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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

  SUBROUTINE VectorBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vectorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorBoundaryInterp2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:2)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: vEq(1:2)
    TYPE(Vector2D) :: f,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:4)
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,iside,idir

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
    DO idir = 1,2
      vEq(idir) = EquationParser(vectorChar(idir), (/'x','y'/))
    ENDDO

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree
            DO idir = 1,2
              f % interior % hostData(idir,i,j,ivar,iel) = &
                vEq(idir) % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
            ENDDO
          ENDDO
          DO iside = 1,4
            DO idir = 1,2
              fActual % boundary % hostData(idir,j,ivar,iside,iel) = &
                vEq(idir) % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
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
      msg = "Numerical Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "Status : [FAIL]"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "Status : [PASS]"
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

  SUBROUTINE VectorGradient2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,vectorChar,tensorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorGradient2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: dForm
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:2)
    CHARACTER(240),INTENT(in) :: tensorChar(1:2,1:2)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: feq(1:2),dfeqChar(1:2,1:2)
    TYPE(MappedVector2D) :: f
    TYPE(MappedScalar2D) :: workScalar
    TYPE(MappedVector2D) :: workVector
    TYPE(MappedTensor2D) :: workTensor
    TYPE(MappedTensor2D) :: dfInterp,dfActual,dfError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: iel,i,j,ivar,row,col,iside

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem*nElem)
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
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfInterp % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfError % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)

    ! Create work objects
    CALL workScalar % Init(cqDegree,cqType,tqDegree,tqType,2*nvar,controlGeometry % nElem)
    CALL workVector % Init(cqDegree,cqType,tqDegree,tqType,2*nvar,controlGeometry % nElem)
    CALL workTensor % Init(cqDegree,cqType,tqDegree,tqType,2*nvar,controlGeometry % nElem)

    ! Create the equation parser object
    DO row = 1, 2
      feq(row) = EquationParser(vectorChar(row), (/'x','y'/))
    ENDDO

    DO col = 1,2
      DO row = 1,2
        dfeqChar(row,col) = EquationParser(tensorChar(row,col), (/'x','y'/))
      ENDDO
    ENDDO

    ! Load the control function
    DO iel = 1, controlGeometry % nElem
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree

            DO row = 1,2
              f % interior % hostData(row,i,j,ivar,iel) = &
                feq(row) % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
            ENDDO

            DO col = 1,2
              DO row = 1,2
                dfActual % interior % hostData(row,col,i,j,ivar,iel) = &
                  dfeqChar(row,col) % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
              ENDDO
            ENDDO

          ENDDO

          DO iside = 1,4
            DO row = 1,2
              f % boundary % hostData(row,j,ivar,iside,iel) = &
                fEq(row) % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
            ENDDO
          ENDDO

        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % Gradient(workScalar,workVector,workTensor,controlGeometry,dfInterp,dForm,.FALSE.)
    dfError = dfActual - dfInterp

    ! Calculate Absolute Maximum Error
    maxErrors = dfError % AbsMaxInterior( )

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workScalar % Free()
    CALL workVector % Free()
    CALL workTensor % Free()
    CALL dfInterp % Free()
    CALL dfActual % Free()
    CALL dfError % Free()

  END SUBROUTINE VectorGradient2D_Test

  SUBROUTINE VectorDivergence2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,vectorChar,scalarChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorDivergence2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: dForm
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:2)
    CHARACTER(240),INTENT(in) :: scalarChar
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: feq(1:2),dfeqChar
    TYPE(MappedVector2D) :: f
    TYPE(MappedVector2D) :: workVector
    TYPE(MappedScalar2D) :: dfInterp,dfActual,dfError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: iel,i,j,ivar,row,col,iside

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem*nElem)
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
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfInterp % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfError % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)

    ! Create work objects
    CALL workVector % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)

    ! Create the equation parser object
    DO row = 1, 2
      feq(row) = EquationParser(vectorChar(row), (/'x','y'/))
    ENDDO

    dfeqChar = EquationParser(scalarChar, (/'x','y'/))

    ! Load the control function
    DO iel = 1, controlGeometry % nElem
      DO ivar = 1, nvar
        DO j = 0, cqDegree
          DO i = 0, cqDegree

            DO row = 1,2
              f % interior % hostData(row,i,j,ivar,iel) = &
                feq(row) % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )
            ENDDO

            dfActual % interior % hostData(i,j,ivar,iel) = &
              dfeqChar % Evaluate( controlGeometry % x % interior % hostData(1:2,i,j,1,iel) )

          ENDDO

          DO iside = 1,4
            DO row = 1,2
              f % boundary % hostData(row,j,ivar,iside,iel) = &
                feq(row) % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
            ENDDO
          ENDDO

        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % Divergence(workVector,controlGeometry,dfInterp,dForm,.FALSE.)
!      physVector,       compVector,geometry,       divVector,dForm,gpuAccel)
    dfError = dfActual - dfInterp

    ! Calculate Absolute Maximum Error
    maxErrors = dfError % AbsMaxInterior( )

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workVector % Free()
    CALL dfInterp % Free()
    CALL dfActual % Free()
    CALL dfError % Free()

  END SUBROUTINE VectorDivergence2D_Test

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

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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
      msg = "Numerical Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "Status : [FAIL]"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "Status : [PASS]"
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

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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
      msg = "Numerical Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "Status : [FAIL]"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "Status : [PASS]"
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

  SUBROUTINE ScalarGradient3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,fChar,gradientChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "ScalarGradient3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: fChar
    CHARACTER(240),INTENT(in) :: gradientChar(1:3)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser) :: feq,gxeq,gyeq,gzeq
    TYPE(MappedScalar3D) :: f
    TYPE(MappedTensor3D) :: workTensor
    TYPE(MappedVector3D) :: dfInterp,dfActual,dfError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: iel,i,j,k,ivar

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem*nElem*nElem)
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
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL workTensor % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfInterp % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfError % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)

    ! Create the equation parser object
    feq = EquationParser(fChar, (/'x','y','z'/))
    gxeq = EquationParser(gradientChar(1), (/'x','y','z'/))
    gyeq = EquationParser(gradientChar(2), (/'x','y','z'/))
    gzeq = EquationParser(gradientChar(3), (/'x','y','z'/))

    ! Load the control function
     DO iel = 1, controlGeometry % nElem
       DO ivar = 1, nvar
         DO k = 0, cqDegree
           DO j = 0, cqDegree
             DO i = 0, cqDegree
               f % interior % hostData(i,j,k,ivar,iel) = &
                 feq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )

               dfActual % interior % hostData(1,i,j,k,ivar,iel) = &
                 gxeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )

               dfActual % interior % hostData(2,i,j,k,ivar,iel) = &
                 gyeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )

               dfActual % interior % hostData(3,i,j,k,ivar,iel) = &
                 gzeq % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
             ENDDO
           ENDDO
         ENDDO
       ENDDO
     ENDDO

    ! Run the grid interpolation
    CALL f % Gradient(workTensor,controlGeometry,dfInterp,.FALSE.)
    dfError = dfActual - dfInterp

    ! Calculate Absolute Maximum Error
    maxErrors = dfError % AbsMaxInterior( )

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workTensor % Free()
    CALL dfInterp % Free()
    CALL dfActual % Free()
    CALL dfError % Free()

  END SUBROUTINE ScalarGradient3D_Test

  SUBROUTINE VectorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vectorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorInterp3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:3)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: vEq(1:3)
    TYPE(Vector3D) :: f,fInterp,fActual,fError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,idir

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
    DO idir = 1,3
      vEq(idir) = EquationParser(vectorChar(idir), (/'x','y','z'/))
    ENDDO

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree
              DO idir = 1,3
                f % interior % hostData(idir,i,j,k,ivar,iel) = &
                  vEq(idir) % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
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
              DO idir = 1,3
                fActual % interior % hostData(idir,i,j,k,ivar,iel) = &
                  vEq(idir) % Evaluate( targetGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
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

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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

  SUBROUTINE VectorBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vectorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorBoundaryInterp3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:3)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser)  :: vEq(1:3)
    TYPE(Vector3D) :: f,fActual,fError
    REAL(prec) :: maxErrors(1:nvar,1:6)
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,iside,idir

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
    DO idir = 1,3
      vEq(idir) = EquationParser(vectorChar(idir), (/'x','y','z'/))
    ENDDO

    ! Load the control function
    DO iel = 1, nel
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree
              DO idir = 1,3
                f % interior % hostData(idir,i,j,k,ivar,iel) = &
                  vEq(idir) % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
              ENDDO
            ENDDO
            DO iside = 1,6
              DO idir = 1,3
                fActual % boundary % hostData(idir,j,k,ivar,iside,iel) = &
                  vEq(idir) % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
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
      msg = "Numerical Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "Status : [FAIL]"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "Status : [PASS]"
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

  SUBROUTINE VectorGradient3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,vectorChar,tensorChar,tolerance,error)
#undef __FUNC__
#define __FUNC__ "VectorGradient3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:3)
    CHARACTER(240),INTENT(in) :: tensorChar(1:3,1:3)
    REAL(prec),INTENT(in) :: tolerance
    INTEGER,INTENT(out) :: error
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser)  :: feq(1:3),dfeqChar(1:3,1:3)
    TYPE(MappedVector3D) :: f
    TYPE(MappedScalar3D) :: workScalar
    TYPE(MappedVector3D) :: workVector
    TYPE(MappedTensor3D) :: workTensor
    TYPE(MappedTensor3D) :: dfInterp,dfActual,dfError
    REAL(prec) :: maxErrors(1:nvar)
    INTEGER :: iel,i,j,k,ivar,row,col

    error = 0
    msg = 'Number of elements : '//Int2Str(nElem*nElem*nElem)
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
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfInterp % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfActual % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfError % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)

    ! Create work objects
    CALL workScalar % Init(cqDegree,cqType,tqDegree,tqType,3*nvar,controlGeometry % nElem)
    CALL workVector % Init(cqDegree,cqType,tqDegree,tqType,3*nvar,controlGeometry % nElem)
    CALL workTensor % Init(cqDegree,cqType,tqDegree,tqType,3*nvar,controlGeometry % nElem)

    ! Create the equation parser object
    DO row = 1,3
      feq(row) = EquationParser(vectorChar(row), (/'x','y','z'/))
    ENDDO

    DO col = 1,3
      DO row = 1,3
        dfeqChar(row,col) = EquationParser(tensorChar(row,col), (/'x','y','z'/))
      ENDDO
    ENDDO

    ! Load the control function
    DO iel = 1, controlGeometry % nElem
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree

              DO row = 1,3
                f % interior % hostData(row,i,j,k,ivar,iel) = &
                  feq(row) % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
              ENDDO

              DO col = 1,3
                DO row = 1,3
                  dfActual % interior % hostData(row,col,i,j,k,ivar,iel) = &
                    dfeqChar(row,col) % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
                ENDDO
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
    CALL f % Gradient(workScalar,workVector,workTensor,controlGeometry,dfInterp,.FALSE.)
    dfError = dfActual - dfInterp

    ! Calculate Absolute Maximum Error
    maxErrors = dfError % AbsMaxInterior( )

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
    END IF

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workScalar % Free()
    CALL workVector % Free()
    CALL workTensor % Free()
    CALL dfInterp % Free()
    CALL dfActual % Free()
    CALL dfError % Free()

  END SUBROUTINE VectorGradient3D_Test

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

    msg = "Numerical Error : "//Float2Str(maxErrors(1))
    IF (maxErrors(1) > tolerance) THEN
      error = error + 1
      ERROR(TRIM(msg))
      ERROR("Status : [FAIL]")
    ELSE
      INFO(TRIM(msg))
      INFO("Status : [PASS]")
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
      msg = "Numerical Error : "//TRIM(Int2Str(iSide))//Float2Str(maxErrors(1,iSide))
      IF (maxErrors(1,iSide) > tolerance) THEN
        error = error + 1
        ERROR(TRIM(msg))
        msg = "Status : [FAIL]"
        ERROR(TRIM(msg))
      ELSE
        INFO(TRIM(msg))
        msg = "Status : [PASS]"
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
