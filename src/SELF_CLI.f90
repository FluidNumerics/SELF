MODULE SELF_CLI

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

  SUBROUTINE BlockMesh1D_Test(cqType,tqType,cqDegree,tqDegree,nElem)
#undef __FUNC__
#define __FUNC__ "BlockMesh1D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh1D) :: mesh
    TYPE(Geometry1D) :: geometry
    INTEGER :: iel,i

    
    INFO('Number of elements : '//Int2Str(nElem))
    INFO('Control point degree : '//Int2Str(cqDegree))
    INFO('Target point degree : '//Int2Str(tqDegree))

    CALL mesh % UniformBlockMesh(cqDegree,nElem, (/0.0_prec,1.0_prec/))

    ! Create the geometry
    CALL geometry % GenerateFromMesh(mesh,cqType,tqType,cqDegree,tqDegree)

    ! To do : File IO for mesh and geometry

    CALL mesh % Free()
    CALL geometry % Free()

  END SUBROUTINE BlockMesh1D_Test

  SUBROUTINE BlockMesh2D_Test(cqType,tqType,cqDegree,tqDegree,nElem)
#undef __FUNC__
#define __FUNC__ "BlockMesh2D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: mesh
    TYPE(SEMQuad) :: geometry
    INTEGER :: iel,jel,i,j
    INTEGER :: row,col

    
    INFO('Number of elements : '//Int2Str(nElem*nElem))
    INFO('Control point degree : '//Int2Str(cqDegree))
    INFO('Target point degree : '//Int2Str(tqDegree))

    CALL mesh % UniformBlockMesh(cqDegree, &
                                 (/nElem,nElem/), &
                                 (/0.0_prec,1.0_prec, &
                                   0.0_prec,1.0_prec/))

    ! Create the geometry
    CALL geometry % GenerateFromMesh(mesh,cqType,tqType,cqDegree,tqDegree)

    ! To do : file io for mesh and geometry

    CALL mesh % Free()
    CALL geometry % Free()

  END SUBROUTINE BlockMesh2D_Test

  SUBROUTINE BlockMesh3D_Test(cqType,tqType,cqDegree,tqDegree,nElem)
#undef __FUNC__
#define __FUNC__ "BlockMesh3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nElem
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: mesh
    TYPE(SEMHex) :: geometry
    INTEGER :: iel,jel,kel,i,j,k
    INTEGER :: row,col

    
    msg = 'Number of elements : '//Int2Str(nElem*nElem*nElem)
    INFO(TRIM(msg))
    INFO('Control point degree : '//Int2Str(cqDegree))
    INFO('Target point degree : '//Int2Str(tqDegree))

    CALL mesh % UniformBlockMesh(cqDegree, &
                                 (/nElem,nElem,nElem/), &
                                 (/0.0_prec,1.0_prec, &
                                   0.0_prec,1.0_prec, &
                                   0.0_prec,1.0_prec/))

    ! Create the geometry
    CALL geometry % GenerateFromMesh(mesh,cqType,tqType,cqDegree,tqDegree)

    ! To do : file IO for mesh and geometry

    CALL mesh % Free()
    CALL geometry % Free()

  END SUBROUTINE BlockMesh3D_Test

  SUBROUTINE ScalarInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,functionChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh1D) :: controlMesh,targetMesh
    TYPE(Geometry1D) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar1D) :: f,fInterp
    INTEGER :: iel,i,ivar

    
    msg = 'Number of elements : '//Int2Str(nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree,nElem,(/0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the target mesh and geometry
    CALL targetMesh % UniformBlockMesh(tqDegree,nElem,(/0.0_prec,1.0_prec/))
    CALL targetGeometry % GenerateFromMesh(targetMesh,tqType,tqType,tqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)
    CALL fInterp % Init(tqDegree,tqType,tqDegree,tqType,nvar,nElem)

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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL fInterp % UpdateHost()
    END IF
#endif
    
    ! To do : file IO for fInterp, targetMesh, targetGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()

  END SUBROUTINE ScalarInterp1D_Test

  SUBROUTINE ScalarBoundaryInterp1D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                         functionChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh1D) :: controlMesh
    TYPE(Geometry1D) :: controlGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar1D) :: f
    INTEGER :: iel,i,ivar,iSide

    
    msg = 'Number of elements : '//Int2Str(nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree,nElem,(/0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)

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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % BoundaryInterp(gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateHost()
    END IF
#endif

    ! To do : file IO for f

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()

  END SUBROUTINE ScalarBoundaryInterp1D_Test

  SUBROUTINE ScalarDerivative1D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,&
                                     fChar,dfChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh1D) :: controlMesh
    TYPE(Geometry1D) :: controlGeometry
    TYPE(EquationParser)  :: feq,dfeq
    TYPE(MappedScalar1D) :: f,dfInterp
    INTEGER :: iel,i,ivar

    
    IF (dForm == selfStrongForm ) THEN
      msg = 'Formulation Type : Strong Form'
    ELSEIF (dForm == selfWeakDGForm ) THEN
      msg = 'Formulation Type : Weak DG Form'
    ELSEIF (dForm == selfWeakCGForm ) THEN
      msg = 'Formulation Type : Weak CG Form'
    ENDIF
    INFO(TRIM(msg))
    msg = 'Number of elements : '//Int2Str(nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree,nElem,(/0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)
    CALL dfInterp % Init(cqDegree,cqType,tqDegree,tqType,nvar,nElem)

    ! Create the equation parser object
    feq = EquationParser(fChar, (/'x'/))
    dfeq = EquationParser(dfChar, (/'x'/))

    ! Load the control function
     DO iel = 1, controlGeometry % nElem
       DO ivar = 1, nvar
         DO i = 0, cqDegree
           f % interior % hostData(i,ivar,iel) = &
             feq % Evaluate( (/controlGeometry % x % interior % hostData(i,1,iel)/) )
         ENDDO
         ! Left Boundary
         f % boundary % hostData(ivar,1,iel) = &
             feq % Evaluate( (/controlGeometry % x % boundary % hostData(1,1,iel)/) )
         ! Right boundary
         f % boundary % hostData(ivar,2,iel) = &
             feq % Evaluate( (/controlGeometry % x % boundary % hostData(1,2,iel)/) )
       ENDDO
     ENDDO

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % Derivative(controlGeometry,dfInterp,dForm,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL dfInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for dfInterp

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL dfInterp % Free()

  END SUBROUTINE ScalarDerivative1D_Test

  SUBROUTINE ScalarInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                 functionChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh,targetMesh
    TYPE(SEMQuad) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar2D) :: f,fInterp
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar

    nel = nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

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


#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL fInterp % UpdateHost()
    END IF
#endif

    ! To do : file IO for fInterp, targetMesh, targetGeometry
    
    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()

  END SUBROUTINE ScalarInterp2D_Test

  SUBROUTINE ScalarBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                         functionChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar2D) :: f
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,iside

    nel = nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % BoundaryInterp(gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateHost()
    END IF
#endif

    ! To do : file io for f

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()

  END SUBROUTINE ScalarBoundaryInterp2D_Test

  SUBROUTINE ScalarGradient2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,&
                                   fChar,gradientChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: feq,gxeq,gyeq
    TYPE(MappedScalar2D) :: f
    TYPE(MappedTensor2D) :: workTensor
    TYPE(MappedVector2D) :: dfInterp
    INTEGER :: iel,i,j,ivar,iside

    
    IF (dForm == selfStrongForm ) THEN
      msg = 'Formulation Type : Strong Form'
    ELSEIF (dForm == selfWeakDGForm ) THEN
      msg = 'Formulation Type : Weak DG Form'
    ELSEIF (dForm == selfWeakCGForm ) THEN
      msg = 'Formulation Type : Weak CG Form'
    ENDIF
    INFO(TRIM(msg))
    msg = 'Number of elements : '//Int2Str(nElem)
    msg = 'Number of elements : '//Int2Str(nElem*nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
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
          ENDDO
          DO iside = 1,4
            f % boundary % hostData(j,ivar,iside,iel) = &
               feq % Evaluate( controlGeometry % x % boundary % hostData(1:2,j,1,iside,iel) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % Gradient(workTensor,controlGeometry,dfInterp,dForm,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL dfInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for dfInterp

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workTensor % Free()
    CALL dfInterp % Free()

  END SUBROUTINE ScalarGradient2D_Test

  SUBROUTINE VectorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                 vectorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh,targetMesh
    TYPE(SEMQuad) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: vEq(1:2)
    TYPE(Vector2D) :: f,fInterp
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,idir

    nel = nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

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


#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL fInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for fInterp

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()

  END SUBROUTINE VectorInterp2D_Test

  SUBROUTINE VectorBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                         vectorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: vEq(1:2)
    TYPE(Vector2D) :: f
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,iside,idir

    nel = nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

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
        ENDDO
      ENDDO
    ENDDO

    ! Run the grid interpolation
#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    CALL f % BoundaryInterp(gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateHost()
    END IF
#endif

    ! To do : file io for f

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()

  END SUBROUTINE VectorBoundaryInterp2D_Test

  SUBROUTINE VectorGradient2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,&
                                   vectorChar,tensorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: feq(1:2),dfeqChar(1:2,1:2)
    TYPE(MappedVector2D) :: f
    TYPE(MappedScalar2D) :: workScalar
    TYPE(MappedVector2D) :: workVector
    TYPE(MappedTensor2D) :: workTensor
    TYPE(MappedTensor2D) :: dfInterp
    INTEGER :: iel,i,j,ivar,row,col,iside

    
    IF (dForm == selfStrongForm ) THEN
      msg = 'Formulation Type : Strong Form'
    ELSEIF (dForm == selfWeakDGForm ) THEN
      msg = 'Formulation Type : Weak DG Form'
    ELSEIF (dForm == selfWeakCGForm ) THEN
      msg = 'Formulation Type : Weak CG Form'
    ENDIF
    INFO(TRIM(msg))
    msg = 'Number of elements : '//Int2Str(nElem)
    msg = 'Number of elements : '//Int2Str(nElem*nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % Gradient(workScalar,workVector,workTensor,controlGeometry,dfInterp,dForm,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL dfInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for dfInterp

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workScalar % Free()
    CALL workVector % Free()
    CALL workTensor % Free()
    CALL dfInterp % Free()

  END SUBROUTINE VectorGradient2D_Test

  SUBROUTINE VectorDivergence2D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,&
                                     vectorChar,scalarChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: feq(1:2),dfeqChar
    TYPE(MappedVector2D) :: f
    TYPE(MappedVector2D) :: workVector
    TYPE(MappedScalar2D) :: dfInterp
    INTEGER :: iel,i,j,ivar,row,col,iside

    
    IF (dForm == selfStrongForm ) THEN
      msg = 'Formulation Type : Strong Form'
    ELSEIF (dForm == selfWeakDGForm ) THEN
      msg = 'Formulation Type : Weak DG Form'
    ELSEIF (dForm == selfWeakCGForm ) THEN
      msg = 'Formulation Type : Weak CG Form'
    ENDIF
    INFO(TRIM(msg))
    msg = 'Number of elements : '//Int2Str(nElem)
    msg = 'Number of elements : '//Int2Str(nElem*nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % Divergence(workVector,controlGeometry,dfInterp,dForm,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL dfInterp % UpdateHost( )
    END IF
#endif

    ! To do : file io for dfInterp

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workVector % Free()
    CALL dfInterp % Free()

  END SUBROUTINE VectorDivergence2D_Test

  SUBROUTINE TensorInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                 tensorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh,targetMesh
    TYPE(SEMQuad) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: tensorEq(1:2,1:2)
    TYPE(Tensor2D) :: f,fInterp
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar
    INTEGER :: row,col

    nel = nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

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


#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL fInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for fInterp

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()

  END SUBROUTINE TensorInterp2D_Test

  SUBROUTINE TensorBoundaryInterp2D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                         tensorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh2D) :: controlMesh
    TYPE(SEMQuad) :: controlGeometry
    TYPE(EquationParser)  :: tensorEq(1:2,1:2)
    TYPE(Tensor2D) :: f
    INTEGER :: nel,iel,jel
    INTEGER :: i,j,ivar,iside
    INTEGER :: row,col

    nel = nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the scalar1d objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,nEl)

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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % BoundaryInterp(gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateHost()
    END IF
#endif

    ! To do : file io for f

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()

  END SUBROUTINE TensorBoundaryInterp2D_Test

  SUBROUTINE ScalarInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                 functionChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar3D) :: f,fInterp
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar

    nel = nElem*nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL fInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for fInterp, targetMesh, targetGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()

  END SUBROUTINE ScalarInterp3D_Test

  SUBROUTINE ScalarBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                         functionChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: feq
    TYPE(Scalar3D) :: f
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,iside

    nel = nElem*nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % BoundaryInterp(gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateHost()
    END IF
#endif

    ! To do : file io for f

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()

  END SUBROUTINE ScalarBoundaryInterp3D_Test

  SUBROUTINE ScalarGradient3D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,&
                                   fChar,gradientChar,gpuAccel)
#undef __FUNC__
#define __FUNC__ "ScalarGradient3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: dForm
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(*),INTENT(in) :: fChar
    CHARACTER(240),INTENT(in) :: gradientChar(1:3)
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser) :: feq,gxeq,gyeq,gzeq
    TYPE(MappedScalar3D) :: f
    TYPE(MappedTensor3D) :: workTensor
    TYPE(MappedVector3D) :: dfInterp
    INTEGER :: iel,i,j,k,ivar,iside

    
    IF (dForm == selfStrongForm ) THEN
      msg = 'Formulation Type : Strong Form'
    ELSEIF (dForm == selfWeakDGForm ) THEN
      msg = 'Formulation Type : Weak DG Form'
    ELSEIF (dForm == selfWeakCGForm ) THEN
      msg = 'Formulation Type : Weak CG Form'
    ENDIF
    INFO(TRIM(msg))
    msg = 'Number of elements : '//Int2Str(nElem)
    msg = 'Number of elements : '//Int2Str(nElem*nElem*nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
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
            ENDDO
            DO iside = 1,6
             f % boundary % hostData(j,k,ivar,iside,iel) = &
               feq % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % Gradient(workTensor,controlGeometry,dfInterp,dForm,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL dfInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for dfInterp and controlGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workTensor % Free()
    CALL dfInterp % Free()

  END SUBROUTINE ScalarGradient3D_Test

  SUBROUTINE VectorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                 vectorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: vEq(1:3)
    TYPE(Vector3D) :: f,fInterp
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,idir

    nel = nElem*nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL fInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for fInterp, targetMesh, targetGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()

  END SUBROUTINE VectorInterp3D_Test

  SUBROUTINE VectorBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                         vectorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser)  :: vEq(1:3)
    TYPE(Vector3D) :: f
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,iside,idir

    nel = nElem*nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % BoundaryInterp(gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateHost()
    END IF
#endif

    ! To do : file io for f, controlGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()

  END SUBROUTINE VectorBoundaryInterp3D_Test

  SUBROUTINE VectorGradient3D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,&
                                   vectorChar,tensorChar,gpuAccel)
#undef __FUNC__
#define __FUNC__ "VectorGradient3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: dForm
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:3)
    CHARACTER(240),INTENT(in) :: tensorChar(1:3,1:3)
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser)  :: feq(1:3),dfeqChar(1:3,1:3)
    TYPE(MappedVector3D) :: f
    TYPE(MappedScalar3D) :: workScalar
    TYPE(MappedVector3D) :: workVector
    TYPE(MappedTensor3D) :: workTensor
    TYPE(MappedTensor3D) :: dfInterp
    INTEGER :: iel,i,j,k,ivar,row,col,iside

    
    IF (dForm == selfStrongForm ) THEN
      msg = 'Formulation Type : Strong Form'
    ELSEIF (dForm == selfWeakDGForm ) THEN
      msg = 'Formulation Type : Weak DG Form'
    ELSEIF (dForm == selfWeakCGForm ) THEN
      msg = 'Formulation Type : Weak CG Form'
    ENDIF
    INFO(TRIM(msg))
    msg = 'Number of elements : '//Int2Str(nElem)
    msg = 'Number of elements : '//Int2Str(nElem*nElem*nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
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

            ENDDO
            DO iside = 1,6
              DO row = 1,3
                f % boundary % hostData(row,j,k,ivar,iside,iel) = &
                  feq(row) % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % Gradient(workScalar,workVector,workTensor,controlGeometry,dfInterp,dForm,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL dfInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for dfInterp, controlGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workScalar % Free()
    CALL workVector % Free()
    CALL workTensor % Free()
    CALL dfInterp % Free()

  END SUBROUTINE VectorGradient3D_Test

  SUBROUTINE VectorDivergence3D_Test(cqType,tqType,cqDegree,tqDegree,dForm,nElem,nvar,&
                                     vectorChar,scalarChar,gpuAccel)
#undef __FUNC__
#define __FUNC__ "VectorDivergence3D_Test"
    IMPLICIT NONE
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: dForm
    INTEGER,INTENT(in) :: nElem
    INTEGER,INTENT(in) :: nVar
    CHARACTER(240),INTENT(in) :: vectorChar(1:3)
    CHARACTER(240),INTENT(in) :: scalarChar
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser)  :: feq(1:3),dfeqChar
    TYPE(MappedVector3D) :: f
    TYPE(MappedVector3D) :: workVector
    TYPE(MappedScalar3D) :: dfInterp
    INTEGER :: iel,i,j,k,ivar,idir,iside

    
    IF (dForm == selfStrongForm ) THEN
      msg = 'Formulation Type : Strong Form'
    ELSEIF (dForm == selfWeakDGForm ) THEN
      msg = 'Formulation Type : Weak DG Form'
    ELSEIF (dForm == selfWeakCGForm ) THEN
      msg = 'Formulation Type : Weak CG Form'
    ENDIF
    INFO(TRIM(msg))
    msg = 'Number of elements : '//Int2Str(nElem*nElem*nElem)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

    ! Create the control mesh and geometry
    CALL controlMesh % UniformBlockMesh(cqDegree, &
                                        (/nElem,nElem,nElem/), &
                                        (/0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec, &
                                          0.0_prec,1.0_prec/))
    CALL controlGeometry % GenerateFromMesh(controlMesh,cqType,tqType,cqDegree,tqDegree)

    ! Create the objects
    CALL f % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)
    CALL dfInterp % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)

    ! Create work objects
    CALL workVector % Init(cqDegree,cqType,tqDegree,tqType,nvar,controlGeometry % nElem)

    ! Create the equation parser object
    DO idir = 1,3
      feq(idir) = EquationParser(vectorChar(idir), (/'x','y','z'/))
    ENDDO

    dfeqChar = EquationParser(scalarChar, (/'x','y','z'/))

    ! Load the control function
    DO iel = 1, controlGeometry % nElem
      DO ivar = 1, nvar
        DO k = 0, cqDegree
          DO j = 0, cqDegree
            DO i = 0, cqDegree

              DO idir = 1,3
                f % interior % hostData(idir,i,j,k,ivar,iel) = &
                  feq(idir) % Evaluate( controlGeometry % x % interior % hostData(1:3,i,j,k,1,iel) )
              ENDDO

            ENDDO

            DO iside = 1,6
              DO idir = 1,3
                f % boundary % hostData(idir,j,k,ivar,iside,iel) = &
                  feq(idir) % Evaluate( controlGeometry % x % boundary % hostData(1:3,j,k,1,iside,iel) )
              ENDDO
            ENDDO

          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % Divergence(workVector,controlGeometry,dfInterp,dForm,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL dfInterp % UpdateHost()
    END IF
#endif
 
    ! To do : file io for dfInterp, controlGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()
    CALL workVector % Free()
    CALL dfInterp % Free()

  END SUBROUTINE VectorDivergence3D_Test

  SUBROUTINE TensorInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                 tensorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh,targetMesh
    TYPE(SEMHex) :: controlGeometry,targetGeometry
    TYPE(EquationParser)  :: tensorEq(1:3,1:3)
    TYPE(Tensor3D) :: f,fInterp
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar
    INTEGER :: row,col

    nel = nElem*nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
    INFO(TRIM(msg))

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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
    CALL f % GridInterp(fInterp,gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL fInterp % UpdateHost()
    END IF
#endif

    ! To do : file io for finterp, targetMesh, targetGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL targetMesh % Free()
    CALL targetGeometry % Free()
    CALL f % Free()
    CALL fInterp % Free()

  END SUBROUTINE TensorInterp3D_Test

  SUBROUTINE TensorBoundaryInterp3D_Test(cqType,tqType,cqDegree,tqDegree,nElem,nvar,&
                                         tensorChar,gpuAccel)
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
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    CHARACTER(240) :: msg
    TYPE(Mesh3D) :: controlMesh
    TYPE(SEMHex) :: controlGeometry
    TYPE(EquationParser)  :: tensorEq(1:3,1:3)
    TYPE(Tensor3D) :: f,fInterp
    INTEGER :: nel,iel
    INTEGER :: i,j,k,ivar,iside
    INTEGER :: row,col

    nel = nElem*nElem*nElem
    
    msg = 'Number of elements : '//Int2Str(nEl)
    INFO(TRIM(msg))
    msg = 'Number of control points : '//Int2Str(cqDegree)
    INFO(TRIM(msg))
    msg = 'Number of target points : '//Int2Str(tqDegree)
    INFO(TRIM(msg))
    msg = 'Number of variables : '//Int2Str(nvar)
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

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateDevice()
    END IF
#endif

    ! Run the grid interpolation
     CALL f % BoundaryInterp(gpuAccel)

#ifdef GPU     
    IF (gpuAccel) THEN
      CALL f % UpdateHost()
    END IF
#endif

    ! To do : file io for f and controlGeometry

    ! Clean up
    CALL controlMesh % Free()
    CALL controlGeometry % Free()
    CALL f % Free()

  END SUBROUTINE TensorBoundaryInterp3D_Test

END MODULE SELF_CLI
