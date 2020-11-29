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

  REAL(prec),PRIVATE,PARAMETER :: exactTolerance = 10.0_prec**(-5)

CONTAINS

  SUBROUTINE BlockMesh1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,error)
#undef __FUNC__
#define __FUNC__ "BlockMesh1D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(out) :: error
    ! Local
    TYPE(Mesh1D) :: mesh
    TYPE(Geometry1D) :: geometry
    REAL(prec) :: expect_dxds,dxds_error,expect_boundx,boundx_error
    INTEGER :: iel,i

    error = 0
    INFO('Number of elements : '//Int2Str(nElem))
    INFO('Control point degree : '//Int2Str(cqDegree))
    INFO('Target point degree : '//Int2Str(tqDegree))

    CALL mesh % UniformBlockMesh(cqDegree,nElem,(/0.0_prec,1.0_prec/))

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

    IF (dxds_error > exactTolerance) THEN
      error = error + 1
      ERROR("Max dx/ds error : "//Float2Str(dxds_error))
      ERROR("[FAIL] Metric Terms Test")
    ELSE
      INFO("Max dx/ds error : "//Float2Str(dxds_error))
      INFO("[PASS] Metric Terms Test Pass")
    END IF

    IF (boundx_error > exactTolerance) THEN
      error = error + 1
      ERROR("Max boundx error : "//Float2Str(boundx_error))
      ERROR("[FAIL] Boundary Interpolation Test")
    ELSE
      INFO("Max boundx error : "//Float2Str(boundx_error))
      INFO("[PASS] Boundary Interpolation Test")
    END IF

  END SUBROUTINE BlockMesh1D_Test

  SUBROUTINE BlockMesh2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,error)
#undef __FUNC__
#define __FUNC__ "BlockMesh2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(out) :: error
    ! Local
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

    CALL mesh % UniformBlockMesh(cqDegree,&
                                 (/nElem,nElem/),&
                                 (/0.0_prec,1.0_prec,&
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
              dxds_error(row,col) = MAX(dxds_error(row,col),&
                      ABS(geometry % dxds % &
                          interior % hostData(row,col,i,j,1,iel) -&
                          expect_dxds(row,col)))
            ENDDO
          ENDDO
          J_error = MAX(J_error,ABS(geometry % J % &
                                    interior % hostData(i,j,1,iel) -&
                                    expect_J)) 
        ENDDO
      ENDDO
    ENDDO

    CALL mesh % Free()

    DO col = 1,2
      DO row = 1,2
        IF (dxds_error(row,col) > exactTolerance) THEN
          error = error + 1
          ERROR("Max dx/ds error ("//TRIM(Int2Str(row))//","//TRIM(Int2Str(col))//") : "//Float2Str(dxds_error(row,col)))
          ERROR("[FAIL] Covariant Tensor Test")
        ELSE
          INFO("Max dx/ds error ("//TRIM(Int2Str(row))//","//TRIM(Int2Str(col))//") : "//Float2Str(dxds_error(row,col)))
          INFO("[PASS] Covariant Tensor Test")
        END IF
      ENDDO
    ENDDO

    IF (J_error > exactTolerance) THEN
      error = error + 1
      ERROR("Max Jacobian error : "//Float2Str(J_error))
      ERROR("[FAIL] Jacobian Test")
    ELSE
      INFO("Max Jacobian error : "//Float2Str(J_error))
      INFO("[PASS] Jacobian Test")
    END IF

  END SUBROUTINE BlockMesh2D_Test

!  SUBROUTINE ScalarInterp1D_Test(cqType,tqType,nControlPoints,nTargetPoints,nElem,nvar,functionChar,error)
!#undef __FUNC__
!#define __FUNC__ "ScalarInterp1D_Test"
!    IMPLICIT NONE
!    INTEGER,INTENT(in) :: cqType
!    INTEGER,INTENT(in) :: tqType
!    INTEGER,INTENT(in) :: nControlPoints
!    INTEGER,INTENT(in) :: nTargetPoints
!    INTEGER,INTENT(in) :: nElem
!    INTEGER,INTENT(in) :: nVar
!    CHARACTER(*),INTENT(in) :: functionChar
!    INTEGER,INTENT(out) :: error
!    ! Local
!    TYPE(Mesh1D) :: controlMesh,targetMesh
!    TYPE(EquationParser)  :: feq
!    TYPE(Scalar1D) :: f,fInterp,fActual,fError
!    REAL(prec) :: expect_dxds,dxds_error,expect_boundx,boundx_error
!    INTEGER :: iel,i
!
!    error = 0
!    INFO('Number of elements : '//Int2Str(nElem))
!    INFO('Number of control points : '//Int2Str(nControlPoints))
!    INFO('Number of target points : '//Int2Str(nTargetPoints))
!    INFO('Number of variables : '//Int2Str(nvar))
!
!    ! Create the mesh
!    CALL controlMesh % UniformBlockMesh(cqType,tqType,nControlPoints,nTargetPoints,nElem,(/0.0_prec,1.0_prec/))
!    CALL targetMesh % UniformBlockMesh(tqType,tqType,nTargetPoints,nTargetPoints,nElem,(/0.0_prec,1.0_prec/))
!
!    ! Create the scalar1d objects
!    CALL f % Init(nControlPoints,cqType,nTargetPoints,tqType,nvar,nElem)
!    CALL fInterp % Init(nTargetPoints,tqType,nTargetPoints,tqType,nvar,nElem)
!    CALL fActual % Init(nTargetPoints,tqType,nTargetPoints,tqType,nvar,nElem)
!    CALL fError % Init(nTargetPoints,tqType,nTargetPoints,tqType,nvar,nElem)
!
!    ! Create the equation parser object
!    feq = EquationParser(functionChar, (/'x'/))
!
!    ! Load the control function
!    ! Load the target function
!  END SUBROUTINE ScalarInterp1D_Test

END MODULE SELF_Tests
