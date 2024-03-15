!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_DGModel3D

    USE SELF_SupportRoutines
    USE SELF_Metadata
    USE SELF_Mesh
    USE SELF_MappedData
    USE SELF_HDF5
    USE HDF5
    USE FEQParse
    USE SELF_Model
  
    IMPLICIT NONE
  
#include "SELF_Macros.h"
  
  
    TYPE,EXTENDS(Model) :: DGModel3D
      TYPE(MappedScalar3D)   :: solution
      TYPE(MappedVector3D)   :: solutionGradient
      TYPE(MappedVector3D)   :: flux
      TYPE(MappedScalar3D)   :: source
      TYPE(MappedScalar3D)   :: fluxDivergence
      TYPE(MappedScalar3D)   :: dSdt
      TYPE(MappedScalar3D)   :: workSol
      TYPE(MappedScalar3D)   :: prevSol
      TYPE(Mesh3D),POINTER   :: mesh
      TYPE(SEMHex),POINTER  :: geometry
  
    CONTAINS
  
      PROCEDURE :: Init => Init_DGModel3D
      PROCEDURE :: Free => Free_DGModel3D
  
      PROCEDURE :: UpdateDevice => UpdateDevice_DGModel3D
  
      PROCEDURE :: UpdateSolution => UpdateSolution_DGModel3D
  
      PROCEDURE :: ResizePrevSol => ResizePrevSol_DGModel3D
  
      PROCEDURE :: UpdateGAB2 => UpdateGAB2_DGModel3D
      PROCEDURE :: UpdateGAB3 => UpdateGAB3_DGModel3D
      PROCEDURE :: UpdateGAB4 => UpdateGAB4_DGModel3D
  
      PROCEDURE :: UpdateGRK2 => UpdateGRK2_DGModel3D
      PROCEDURE :: UpdateGRK3 => UpdateGRK3_DGModel3D
      PROCEDURE :: UpdateGRK4 => UpdateGRK4_DGModel3D
  
      PROCEDURE :: CalculateTendency => CalculateTendency_DGModel3D
      PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_DGModel3D
  
      GENERIC :: SetSolution => SetSolutionFromChar_DGModel3D, &
        SetSolutionFromEqn_DGModel3D
      PROCEDURE,PRIVATE :: SetSolutionFromChar_DGModel3D
      PROCEDURE,PRIVATE :: SetSolutionFromEqn_DGModel3D
    
      PROCEDURE :: ReadModel => Read_DGModel3D
      PROCEDURE :: WriteModel => Write_DGModel3D
      PROCEDURE :: WriteTecplot => WriteTecplot_DGModel3D
  
    END TYPE DGModel3D
  
    INTERFACE
      SUBROUTINE UpdateSolution_DGModel3D_gpu(solution,dSdt,dt,N,nVar,nEl) &
        BIND(c,name="UpdateSolution_DGModel3D_gpu")
        USE ISO_C_BINDING
        USE SELF_Constants
        IMPLICIT NONE
        TYPE(C_PTR), value :: solution,dSdt
        INTEGER(C_INT),VALUE :: N,nVar,nEl
        REAL(c_prec),VALUE :: dt
      END SUBROUTINE UpdateSolution_DGModel3D_gpu
    END INTERFACE
  
    INTERFACE
      SUBROUTINE UpdateGAB2_DGModel3D_gpu(prevsol,solution,m,nPrev,N,nVar,nEl) &
        BIND(c,name="UpdateGAB2_DGModel3D_gpu")
        USE ISO_C_BINDING
        USE SELF_Constants
        IMPLICIT NONE
        TYPE(C_PTR), value :: prevsol,solution
        INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
      END SUBROUTINE UpdateGAB2_DGModel3D_gpu
    END INTERFACE
  
    INTERFACE
      SUBROUTINE UpdateGAB3_DGModel3D_gpu(prevsol,solution,m,nPrev,N,nVar,nEl) &
        BIND(c,name="UpdateGAB3_DGModel3D_gpu")
        USE ISO_C_BINDING
        USE SELF_Constants
        IMPLICIT NONE
        TYPE(C_PTR), value :: prevsol,solution
        INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
      END SUBROUTINE UpdateGAB3_DGModel3D_gpu
    END INTERFACE
  
    INTERFACE
      SUBROUTINE UpdateGAB4_DGModel3D_gpu(prevsol,solution,m,nPrev,N,nVar,nEl) &
        BIND(c,name="UpdateGAB4_DGModel3D_gpu")
        USE ISO_C_BINDING
        USE SELF_Constants
        IMPLICIT NONE
        TYPE(C_PTR), value :: prevsol,solution
        INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
      END SUBROUTINE UpdateGAB4_DGModel3D_gpu
    END INTERFACE
  
    INTERFACE
      SUBROUTINE UpdateGRK_DGModel3D_gpu(grk,solution,dSdt,rk_a,rk_g,dt,nWork,N,nVar,nEl) &
        BIND(c,name="UpdateGRK_DGModel3D_gpu")
        USE ISO_C_BINDING
        USE SELF_Constants
        IMPLICIT NONE
        TYPE(C_PTR), value :: grk,solution,dSdt
        INTEGER(C_INT),VALUE :: nWork,N,nVar,nEl
        REAL(c_prec),VALUE :: rk_a,rk_g,dt
      END SUBROUTINE UpdateGRK_DGModel3D_gpu
    END INTERFACE
  
    INTERFACE
      SUBROUTINE CalculateDSDt_DGModel3D_gpu(fluxDivergence,source,dSdt,N,nVar,nEl) &
        BIND(c,name="CalculateDSDt_DGModel3D_gpu")
        USE ISO_C_BINDING
        IMPLICIT NONE
        TYPE(C_PTR), value :: fluxDivergence,source,dSdt
        INTEGER(C_INT),VALUE :: N,nVar,nEl
      END SUBROUTINE CalculateDSDt_DGModel3D_gpu
    END INTERFACE
  
  CONTAINS
  
    SUBROUTINE Init_DGModel3D(this,nvar,mesh,geometry,decomp)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(out) :: this
      INTEGER,INTENT(in) :: nvar
      TYPE(Mesh3D),INTENT(in),TARGET :: mesh
      TYPE(SEMHex),INTENT(in),TARGET :: geometry
      TYPE(MPILayer),INTENT(in),TARGET :: decomp
      ! Local
      INTEGER :: ivar
      CHARACTER(LEN=3) :: ivarChar
      CHARACTER(LEN=25) :: varname
  
      this % decomp => decomp
      this % mesh => mesh
      this % geometry => geometry
      this % GPUBackend = .FALSE.
  
      CALL this % solution % Init(geometry % x % interp,nVar,this % mesh % nElem)
      CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
      CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
      CALL this % dSdt % Init(geometry % x % interp,nVar,this % mesh % nElem)
      CALL this % solutionGradient % Init(geometry % x % interp,nVar,this % mesh % nElem)
      CALL this % flux % Init(geometry % x % interp,nVar,this % mesh % nElem)
      CALL this % source % Init(geometry % x % interp,nVar,this % mesh % nElem)
      CALL this % fluxDivergence % Init(geometry % x % interp,nVar,this % mesh % nElem)
  
      ! set default metadata
      DO ivar = 1,nvar
        WRITE (ivarChar,'(I3.3)') ivar
        varname = "solution"//TRIM(ivarChar)
        CALL this % solution % SetName(ivar,varname)
        CALL this % solution % SetUnits(ivar,"[null]")
      END DO
  
    END SUBROUTINE Init_DGModel3D
  
    SUBROUTINE Free_DGModel3D(this)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
  
      CALL this % solution % Free()
      CALL this % workSol % Free()
      CALL this % prevSol % Free()
      CALL this % dSdt % Free()
      CALL this % solutionGradient % Free()
      CALL this % flux % Free()
      CALL this % source % Free()
      CALL this % fluxDivergence % Free()
  
    END SUBROUTINE Free_DGModel3D
  
    SUBROUTINE ResizePrevSol_DGModel3D(this,m)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
      ! Local
      INTEGER :: nVar
  
      ! Free space, if necessary
      CALL this % prevSol % Free()
  
      ! Reallocate with increased variable dimension for
      ! storing "m" copies of solution data
      nVar = this % solution % nVar
      CALL this % prevSol % Init(this % geometry % x % interp,m*nVar,this % mesh % nElem)
  
    END SUBROUTINE ResizePrevSol_DGModel3D
  
    SUBROUTINE UpdateDevice_DGModel3D(this)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
  
      CALL this % mesh % UpdateDevice()
      CALL this % geometry % UpdateDevice()
      CALL this % dSdt % UpdateDevice()
      CALL this % solution % UpdateDevice()
      CALL this % solutionGradient % UpdateDevice()
      CALL this % flux % UpdateDevice()
      CALL this % source % UpdateDevice()
      CALL this % fluxDivergence % UpdateDevice()
  
    END SUBROUTINE UpdateDevice_DGModel3D
  
    SUBROUTINE SetSolutionFromEqn_DGModel3D(this,eqn)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      TYPE(EquationParser),INTENT(in) :: eqn(1:this % solution % nVar)
      ! Local
      INTEGER :: iVar
  
      ! Copy the equation parser
      DO iVar = 1,this % solution % nVar
        CALL this % solution % SetEquation(ivar,eqn(iVar) % equation)
      END DO
  
      CALL this % solution % SetInteriorFromEquation(this % geometry,this % t)
      if ( this % GPUBackend ) then
        CALL this % solution % BoundaryInterp(this % hipblas_handle)
      else 
        CALL this % solution % BoundaryInterp()
      endif
  
    END SUBROUTINE SetSolutionFromEqn_DGModel3D
  
    SUBROUTINE SetSolutionFromChar_DGModel3D(this,eqnChar)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      CHARACTER(*),INTENT(in) :: eqnChar(1:this % solution % nVar)
      ! Local
      INTEGER :: iVar
  
      DO iVar = 1,this % solution % nVar
        CALL this % solution % SetEquation(ivar,TRIM(eqnChar(iVar)))
      END DO
  
      CALL this % solution % SetInteriorFromEquation(this % geometry,this % t)
      if ( this % GPUBackend ) then
        CALL this % solution % BoundaryInterp(this % hipblas_handle)
      else 
        CALL this % solution % BoundaryInterp()
      endif
  
    END SUBROUTINE SetSolutionFromChar_DGModel3D
  
    SUBROUTINE UpdateSolution_DGModel3D(this,dt)
      !! Computes a solution update as , where dt is either provided through the interface
      !! or taken as the Model's stored time step size (model % dt)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      REAL(prec),OPTIONAL,INTENT(in) :: dt
      ! Local
      REAL(prec) :: dtLoc
      INTEGER :: i,j,k,iVar,iEl
  
      IF (PRESENT(dt)) THEN
        dtLoc = dt
      ELSE
        dtLoc = this % dt
      END IF
  
      IF (this % GPUBackend) THEN
  
        CALL UpdateSolution_DGModel3D_gpu(c_loc(this % solution % interior), &
                                          c_loc(this % dSdt % interior), &
                                                dtLoc, &
                                                this % solution % interp % N, &
                                                this % solution % nVar, &
                                                this % solution % nElem)
  
      ELSE
  
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO k = 1,this % solution % interp % N+1

            DO j = 1,this % solution % interp % N+1
              DO i = 1,this % solution % interp % N+1
  
                this % solution % interior(i,j,k,iEl,iVar) = &
                  this % solution % interior(i,j,k,iEl,iVar) + &
                  dtLoc*this % dSdt % interior(i,j,k,iEl,iVar)
  
              END DO
            END DO
        ENDDO
          END DO
        END DO
  
      END IF
  
    END SUBROUTINE UpdateSolution_DGModel3D
  
    SUBROUTINE UpdateGAB2_DGModel3D(this,m)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
      ! Local
      INTEGER :: i,j,k,nVar,iVar,iEl
  
      IF (this % GPUBackend) THEN
  
        CALL UpdateGAB2_DGModel3D_gpu(c_loc(this % prevSol % interior), &
                                      c_loc(this % solution % interior), &
                                            m, &
                                            this % prevsol % nVar, &
                                            this % solution % interp % N, &
                                            this % solution % nVar, &
                                            this % solution % nElem)
  
      ELSE
  
        ! ab2_weight
        IF (m == 0) THEN ! Initialization step - store the solution in the prevSol
  
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  this % prevSol % interior(i,j,k,iEl,iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        ELSEIF (m == 1) THEN ! Reset solution
  
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  this % solution % interior(i,j,k,iEl,iVar) = this % prevSol % interior(i,j,k,iEl,iVar)
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        ELSE ! Main looping section - nVar the previous solution, store the new solution, and
          ! create an interpolated solution to use for tendency calculation
  
          nVar = this % solution % nVar
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  ! Bump the last solution
            this % prevSol % interior(i,j,k,iEl,nVar + iVar) = this % prevSol % interior(i,j,k,iEl,iVar)
  
                  ! Store the new solution
                  this % prevSol % interior(i,j,k,iEl,iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                  this % solution % interior(i,j,k,iEl,iVar) = &
                    1.5_PREC*this % prevSol % interior(i,j,k,iEl,iVar) - &
                    0.5_PREC*this % prevSol % interior(i,j,k,iEl,nVar + iVar)
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        END IF
  
      END IF
  
    END SUBROUTINE UpdateGAB2_DGModel3D
  
    SUBROUTINE UpdateGAB3_DGModel3D(this,m)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
      ! Local
      INTEGER :: i,j,k,nVar,iVar,iEl
  
      IF (this % GPUBackend) THEN
  
        CALL UpdateGAB3_DGModel3D_gpu(c_loc(this % prevSol % interior), &
                                            c_loc(this % solution % interior), &
                                            m, &
                                            this % prevsol % nVar, &
                                            this % solution % interp % N, &
                                            this % solution % nVar, &
                                            this % solution % nElem)
  
      ELSE
  
        IF (m == 0) THEN ! Initialization step - store the solution in the prevSol at nvar+ivar
  
          nVar = this % solution % nVar
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
           this % prevSol % interior(i,j,k,iEl,nVar + iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        ELSEIF (m == 1) THEN ! Initialization step - store the solution in the prevSol at ivar
  
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  this % prevSol % interior(i,j,k,iEl,iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        ELSEIF (m == 2) THEN ! Copy the solution back from the most recent prevsol
  
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  this % solution % interior(i,j,k,iEl,iVar) = this % prevSol % interior(i,j,k,iEl,iVar)
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        ELSE ! Main looping section - nVar the previous solution, store the new solution, and
          ! create an interpolated solution to use for tendency calculation
  
          nVar = this % solution % nVar
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  ! Bump the last two stored solutions
   this % prevSol % interior(i,j,k,iEl,2*nVar + iVar) = this % prevSol % interior(i,j,k,iEl,nVar + iVar)
            this % prevSol % interior(i,j,k,iEl,nVar + iVar) = this % prevSol % interior(i,j,k,iEl,iVar)
  
                  ! Store the new solution
                  this % prevSol % interior(i,j,k,iEl,iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                  this % solution % interior(i,j,k,iEl,iVar) = &
                    (23.0_PREC*this % prevSol % interior(i,j,k,iEl,iVar) - &
                     16.0_PREC*this % prevSol % interior(i,j,k,iEl,nVar + iVar) + &
                     5.0_PREC*this % prevSol % interior(i,j,k,iEl,2*nVar + iVar))/12.0_PREC
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
        END IF
  
      END IF
  
    END SUBROUTINE UpdateGAB3_DGModel3D
  
    SUBROUTINE UpdateGAB4_DGModel3D(this,m)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
      ! Local
      INTEGER :: i,j,k,nVar,iVar,iEl
  
      IF (this % GPUBackend) THEN
  
        CALL UpdateGAB4_DGModel3D_gpu(c_loc(this % prevSol % interior), &
                                            c_loc(this % solution % interior), &
                                            m, &
                                            this % prevsol % nVar, &
                                            this % solution % interp % N, &
                                            this % solution % nVar, &
                                            this % solution % nElem)
  
      ELSE
  
        IF (m == 0) THEN ! Initialization step - store the solution in the prevSol at nvar+ivar
  
          nVar = this % solution % nVar
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
         this % prevSol % interior(i,j,k,iEl,2*nVar + iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        ELSEIF (m == 1) THEN ! Initialization step - store the solution in the prevSol at ivar
  
          nVar = this % solution % nVar
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
           this % prevSol % interior(i,j,k,iEl,nVar + iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        ELSEIF (m == 2) THEN ! Initialization step - store the solution in the prevSol at ivar
  
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  this % prevSol % interior(i,j,k,iEl,iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        ELSEIF (m == 3) THEN ! Copy the solution back from the most recent prevsol
  
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  this % solution % interior(i,j,k,iEl,iVar) = this % prevSol % interior(i,j,k,iEl,iVar)
  
                END DO
            ENDDO
              END DO
            END DO
          END DO
  
        ELSE ! Main looping section - nVar the previous solution, store the new solution, and
          ! create an interpolated solution to use for tendency calculation
  
          nVar = this % solution % nVar
          DO iEl = 1,this % solution % nElem
            DO iVar = 1,this % solution % nVar
                DO k = 1,this % solution % interp % N+1

              DO j = 1,this % solution % interp % N+1
                DO i = 1,this % solution % interp % N+1
  
                  ! Bump the last two stored solutions
     this % prevSol % interior(i,j,k,iEl,3*nVar + iVar) = this % prevSol % interior(i,j,k,2*nVar+iVar,iEl)
   this % prevSol % interior(i,j,k,iEl,2*nVar + iVar) = this % prevSol % interior(i,j,k,iEl,nVar + iVar)
            this % prevSol % interior(i,j,k,iEl,nVar + iVar) = this % prevSol % interior(i,j,k,iEl,iVar)
  
                  ! Store the new solution
                  this % prevSol % interior(i,j,k,iEl,iVar) = this % solution % interior(i,j,k,iEl,iVar)
  
                  this % solution % interior(i,j,k,iEl,iVar) = &
                    (55.0_PREC*this % prevSol % interior(i,j,k,iEl,iVar) - &
                     59.0_PREC*this % prevSol % interior(i,j,k,iEl,nVar + iVar) + &
                     37.0_PREC*this % prevSol % interior(i,j,k,iEl,2*nVar + iVar) - &
                     9.0_PREC*this % prevSol % interior(i,j,k,iEl,3*nVar + iVar))/24.0_PREC
  
                END DO
              END DO
            ENDDO
            END DO
          END DO
  
        END IF
  
      END IF
  
    END SUBROUTINE UpdateGAB4_DGModel3D
  
    SUBROUTINE UpdateGRK2_DGModel3D(this,m)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
      ! Local
      INTEGER :: i,j,k,iVar,iEl
  
      IF (this % GPUBackend) THEN
  
        CALL UpdateGRK_DGModel3D_gpu(c_loc(this % workSol % interior), &
                                           c_loc(this % solution % interior), &
                                           c_loc(this % dSdt % interior), &
                                           rk2_a(m),rk2_g(m),this % dt, &
                                           this % worksol % nVar, &
                                           this % solution % interp % N, &
                                           this % solution % nVar, &
                                           this % solution % nElem)
  
      ELSE
  
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO k = 1,this % solution % interp % N+1

            DO j = 1,this % solution % interp % N+1
              DO i = 1,this % solution % interp % N+1
  
                this % workSol % interior(i,j,k,iEl,iVar) = rk2_a(m)* &
                                                                    this % workSol % interior(i,j,k,iEl,iVar) + &
                                                                     this % dSdt % interior(i,j,k,iEl,iVar)
  
                this % solution % interior(i,j,k,iEl,iVar) = &
                  this % solution % interior(i,j,k,iEl,iVar) + &
                  rk2_g(m)*this % dt*this % workSol % interior(i,j,k,iEl,iVar)
  
              END DO
            ENDDO
            END DO
          END DO
        END DO
  
      END IF
  
    END SUBROUTINE UpdateGRK2_DGModel3D
  
    SUBROUTINE UpdateGRK3_DGModel3D(this,m)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
      ! Local
      INTEGER :: i,j,k,iVar,iEl
  
      IF (this % GPUBackend) THEN
  
        CALL UpdateGRK_DGModel3D_gpu(c_loc(this % workSol % interior), &
                                           c_loc(this % solution % interior), &
                                           c_loc(this % dSdt % interior), &
                                           rk3_a(m),rk3_g(m),this % dt, &
                                           this % worksol % nVar, &
                                           this % solution % interp % N, &
                                           this % solution % nVar, &
                                           this % solution % nElem)
  
      ELSE
  
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO k = 1,this % solution % interp % N+1

            DO j = 1,this % solution % interp % N+1
              DO i = 1,this % solution % interp % N+1
  
                this % workSol % interior(i,j,k,iEl,iVar) = rk3_a(m)* &
                                                                    this % workSol % interior(i,j,k,iEl,iVar) + &
                                                                     this % dSdt % interior(i,j,k,iEl,iVar)
  
                this % solution % interior(i,j,k,iEl,iVar) = &
                  this % solution % interior(i,j,k,iEl,iVar) + &
                  rk3_g(m)*this % dt*this % workSol % interior(i,j,k,iEl,iVar)
  
              END DO
            END DO
        ENDDO
          END DO
        END DO
  
      END IF
  
    END SUBROUTINE UpdateGRK3_DGModel3D
  
    SUBROUTINE UpdateGRK4_DGModel3D(this,m)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
      ! Local
      INTEGER :: i,j,k,iVar,iEl
  
      IF (this % GPUBackend) THEN
  
        CALL UpdateGRK_DGModel3D_gpu(c_loc(this % workSol % interior), &
                                           c_loc(this % solution % interior), &
                                           c_loc(this % dSdt % interior), &
                                           rk4_a(m),rk4_g(m),this % dt, &
                                           this % workSol % nVar, &
                                           this % solution % interp % N, &
                                           this % solution % nVar, &
                                           this % solution % nElem)
  
      ELSE
  
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO k = 1,this % solution % interp % N+1

            DO j = 1,this % solution % interp % N+1
              DO i = 1,this % solution % interp % N+1
  
                this % workSol % interior(i,j,k,iEl,iVar) = rk4_a(m)* &
                                                                    this % workSol % interior(i,j,k,iEl,iVar) + &
                                                                     this % dSdt % interior(i,j,k,iEl,iVar)
  
                this % solution % interior(i,j,k,iEl,iVar) = &
                  this % solution % interior(i,j,k,iEl,iVar) + &
                  rk4_g(m)*this % dt*this % workSol % interior(i,j,k,iEl,iVar)
  
              END DO
            END DO
        ENDDO
          END DO
        END DO
  
      END IF
  
    END SUBROUTINE UpdateGRK4_DGModel3D
  
    SUBROUTINE CalculateFluxDivergence_DGModel3D(this)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
  
      if( this % GPUBackend )then
        CALL this % flux % DGDivergence(this % geometry, this % fluxDivergence, this % hipblas_handle)
      else
        CALL this % flux % DGDivergence(this % geometry, this % fluxDivergence)
      end if
  
  
    END SUBROUTINE CalculateFluxDivergence_DGModel3D
  
    SUBROUTINE CalculateTendency_DGModel3D(this)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      ! Local
      INTEGER :: i,j,k,iVar,iEl
  
      CALL this % PreTendency()
      if (this % GPUBackend)then
        CALL this % solution % BoundaryInterp(this % hipblas_handle)
        CALL this % solution % SideExchange(this % mesh,this % decomp,this % hipblas_handle)
      else
        CALL this % solution % BoundaryInterp()
        CALL this % solution % SideExchange(this % mesh,this % decomp)
      endif
      CALL this % SetBoundaryCondition()
      CALL this % SourceMethod()
      CALL this % RiemannSolver()
      CALL this % FluxMethod()
      CALL this % CalculateFluxDivergence()
  
      IF (this % GPUBackend) THEN
  
        CALL CalculateDSDt_DGModel3D_gpu(c_loc(this % fluxDivergence % interior), &
                                               c_loc(this % source % interior), &
                                               c_loc(this % dSdt % interior), &
                                               this % solution % interp % N, &
                                               this % solution % nVar, &
                                               this % solution % nElem)
  
      ELSE
  
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO k = 1,this % solution % interp % N+1

            DO j = 1,this % solution % interp % N+1
              DO i = 1,this % solution % interp % N+1
  
                this % dSdt % interior(i,j,k,iEl,iVar) = &
                  this % source % interior(i,j,k,iEl,iVar) - &
                  this % fluxDivergence % interior(i,j,k,iEl,iVar)
  
              END DO
            END DO
        ENDDO
          END DO
        END DO
  
      END IF
  
    END SUBROUTINE CalculateTendency_DGModel3D
  
    SUBROUTINE Write_DGModel3D(this,fileName)
#undef __FUNC__
#define __FUNC__ "Write_DGModel3D"
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
      ! Local
      INTEGER(HID_T) :: fileId
      TYPE(Scalar3D) :: solution
      TYPE(Vector3D) :: x
      TYPE(Lagrange),TARGET :: interp
      CHARACTER(LEN=self_FileNameLength) :: pickupFile
      CHARACTER(13) :: timeStampString
  
      IF (PRESENT(filename)) THEN
        pickupFile = filename
      ELSE
        WRITE (timeStampString,'(I13.13)') this % ioIterate
        pickupFile = 'solution.'//timeStampString//'.h5'
      END IF
  
      INFO("Writing pickup file : "//TRIM(pickupFile))
  
      IF (this % decomp % mpiEnabled) THEN
  
        CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId,this % decomp % mpiComm)
  
        ! Write the interpolant to the file
        INFO("Writing interpolant data to file")
        CALL this % solution % interp % WriteHDF5( fileId )
  
        ! In this section, we write the solution and geometry on the control (quadrature) grid
        ! which can be used for model pickup runs or post-processing
        ! Write the model state to file
        INFO("Writing control grid solution to file")
        CALL CreateGroup_HDF5(fileId,'/controlgrid')
        CALL this % solution % WriteHDF5( fileId, '/controlgrid/solution', &
        this % decomp % offsetElem(this % decomp % rankId), this % decomp % nElem )
  
        ! Write the geometry to file
        INFO("Writing control grid geometry to file")
        CALL CreateGroup_HDF5(fileId,'/controlgrid/geometry')
        CALL this % geometry % x % WriteHDF5( fileId, '/controlgrid/geometry/x', &
        this % decomp % offsetElem(this % decomp % rankId), this % decomp % nElem )
  
        ! -- END : writing solution on control grid -- !
  
        ! Interpolate the solution to a grid for plotting results
        ! Create an interpolant for the uniform grid
        CALL interp % Init(this % solution % interp % M, &
                            this % solution % interp % targetNodeType, &
                            this % solution % interp % N, &
                            this % solution % interp % controlNodeType)
  
        CALL solution % Init(interp, &
                              this % solution % nVar,this % solution % nElem)
  
        CALL x % Init(interp,1,this % solution % nElem)
  
        ! Map the mesh positions to the target grid
        CALL this % geometry % x % GridInterp(x)
  
        ! Map the solution to the target grid
        CALL this % solution % GridInterp(solution)
  
        ! Write the model state to file
        CALL CreateGroup_HDF5(fileId,'/targetgrid')
        CALL solution % WriteHDF5( fileId, '/targetgrid/solution', &
        this % decomp % offsetElem(this % decomp % rankId), this % decomp % nElem )
  
        ! Write the geometry to file
        CALL CreateGroup_HDF5(fileId,'/targetgrid/mesh')
        CALL x % WriteHDF5( fileId, '/targetgrid/mesh/coords', &
        this % decomp % offsetElem(this % decomp % rankId), this % decomp % nElem )
  
        CALL Close_HDF5(fileId)
  
      ELSE
  
        CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)
  
        ! Write the interpolant to the file
        INFO("Writing interpolant data to file")
        CALL this % solution % interp % WriteHDF5( fileId )
  
        ! In this section, we write the solution and geometry on the control (quadrature) grid
        ! which can be used for model pickup runs or post-processing
  
        ! Write the model state to file
        INFO("Writing control grid solution to file")
        CALL CreateGroup_HDF5(fileId,'/controlgrid')
        CALL this % solution % WriteHDF5( fileId, '/controlgrid/solution' )
  
        ! Write the geometry to file
        INFO("Writing control grid  geometry to file")
        CALL CreateGroup_HDF5(fileId,'/controlgrid/geometry')
        CALL this % geometry % x % WriteHDF5(fileId,'/controlgrid/geometry/x')
        ! -- END : writing solution on control grid -- !
  
        ! Interpolate the solution to a grid for plotting results
        ! Create an interpolant for the uniform grid
        CALL interp % Init(this % solution % interp % M, &
                            this % solution % interp % targetNodeType, &
                            this % solution % interp % N, &
                            this % solution % interp % controlNodeType)
  
        CALL solution % Init(interp, &
                              this % solution % nVar,this % solution % nElem)
  
        CALL x % Init(interp,1,this % solution % nElem)
  
        ! Map the mesh positions to the target grid
        CALL this % geometry % x % GridInterp(x)
  
        ! Map the solution to the target grid
        CALL this % solution % GridInterp(solution)
  
        ! Write the model state to file
        INFO("Writing target grid solution to file")
        CALL CreateGroup_HDF5(fileId,'/targetgrid')
        CALL solution % WriteHDF5(fileId, '/targetgrid/solution')
  
        ! Write the geometry to file
        INFO("Writing target grid geometry to file")
        CALL CreateGroup_HDF5(fileId,'/targetgrid/geometry')
        CALL x % WriteHDF5(fileId,'/targetgrid/geometry/x')
  
        CALL Close_HDF5(fileId)
  
      END IF
  
      CALL x % Free()
      CALL solution % Free()
      CALL interp % Free()
  
    END SUBROUTINE Write_DGModel3D
  
    SUBROUTINE Read_DGModel3D(this,fileName)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      CHARACTER(*),INTENT(in) :: fileName
      ! Local
      INTEGER(HID_T) :: fileId
      INTEGER(HID_T) :: solOffset(1:5)
      INTEGER :: firstElem
      INTEGER :: N
  
      IF (this % decomp % mpiEnabled) THEN
        CALL Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId, &
                       this % decomp % mpiComm)
      ELSE
        CALL Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId)
      END IF
  
      ! CALL ReadAttribute_HDF5(fileId,'N',N)
  
      ! IF (this % solution % interp % N /= N) THEN
      !   STOP 'Error : Solution polynomial degree does not match input file'
      ! END IF
  
      IF (this % decomp % mpiEnabled) THEN
        firstElem = this % decomp % offsetElem(this % decomp % rankId) + 1
        solOffset(1:5) = (/0,0,0,1,firstElem/)
        CALL ReadArray_HDF5(fileId,'/controlgrid/solution/interior', &
                            this % solution % interior,solOffset)
      ELSE
        CALL ReadArray_HDF5(fileId,'/controlgrid/solution/interior',this % solution % interior)
      END IF
  
      CALL Close_HDF5(fileId)
  
    END SUBROUTINE Read_DGModel3D
  
    SUBROUTINE WriteTecplot_DGModel3D(this,filename)
      IMPLICIT NONE
      CLASS(DGModel3D),INTENT(inout) :: this
      CHARACTER(*),INTENT(in),OPTIONAL :: filename
      ! Local
      CHARACTER(8) :: zoneID
      INTEGER :: fUnit
      INTEGER :: iEl,i,j,k,iVar
      CHARACTER(LEN=self_FileNameLength) :: tecFile
      CHARACTER(LEN=self_TecplotHeaderLength) :: tecHeader
      CHARACTER(LEN=self_FormatLength) :: fmat
      CHARACTER(13) :: timeStampString
      CHARACTER(5) :: rankString
      TYPE(Scalar3D) :: solution
      TYPE(Vector3D) :: solutionGradient
      TYPE(Vector3D) :: x
      TYPE(Lagrange),TARGET :: interp
  
      IF (PRESENT(filename)) THEN
        tecFile = filename
      ELSE
        WRITE (timeStampString,'(I13.13)') this % ioIterate
  
        IF (this % decomp % mpiEnabled) THEN
          WRITE (rankString,'(I5.5)') this % decomp % rankId
          tecFile = 'solution.'//rankString//'.'//timeStampString//'.tec'
        ELSE
          tecFile = 'solution.'//timeStampString//'.tec'
        END IF
  
      END IF
  
      ! Create an interpolant for the uniform grid
      CALL interp % Init(this % solution % interp % M, &
                         this % solution % interp % targetNodeType, &
                         this % solution % interp % N, &
                         this % solution % interp % controlNodeType)
  
      CALL solution % Init(interp, &
                           this % solution % nVar,this % solution % nElem)
  
      CALL solutionGradient % Init(interp, &
                                   this % solution % nVar,this % solution % nElem)
  
      CALL x % Init(interp,1,this % solution % nElem)
  
      ! Map the mesh positions to the target grid
      CALL this % geometry % x % GridInterp(x)
  
      ! Map the solution to the target grid
      CALL this % solution % GridInterp(solution)
  
      ! Map the solution to the target grid
      CALL this % solutionGradient % GridInterp(solutionGradient)
  
      OPEN (UNIT=NEWUNIT(fUnit), &
            FILE=TRIM(tecFile), &
            FORM='formatted', &
            STATUS='replace')
  
      tecHeader = 'VARIABLES = "X", "Y", "Z"'
      DO iVar = 1,this % solution % nVar
        tecHeader = TRIM(tecHeader)//', "'//TRIM(this % solution % meta(iVar) % name)//'"'
      END DO
  
      DO iVar = 1,this % solution % nVar
        tecHeader = TRIM(tecHeader)//', "d/dx('//TRIM(this % solution % meta(iVar) % name)//')"'
      END DO
  
      DO iVar = 1,this % solution % nVar
        tecHeader = TRIM(tecHeader)//', "d/dy('//TRIM(this % solution % meta(iVar) % name)//')"'
      END DO
  
      WRITE (fUnit,*) TRIM(tecHeader)
  
      ! Create format statement
      WRITE (fmat,*) 3*this % solution % nvar + 3
      fmat = '('//TRIM(fmat)//'(ES16.7E3,1x))'
  
      DO iEl = 1,this % solution % nElem
  
        ! TO DO :: Get the global element ID
        WRITE (zoneID,'(I8.8)') iEl
        WRITE (fUnit,*) 'ZONE T="el'//TRIM(zoneID)//'", I=',this % solution % interp % M + 1, &
          ', J=',this % solution % interp % M + 1
  
        DO k = 1,this % solution % interp % M+1
        DO j = 1,this % solution % interp % M+1
          DO i = 1,this % solution % interp % M+1
  
            WRITE (fUnit,fmat) x % interior(1,i,j,k,iEl,1), &
              x % interior(2,i,j,k,iEl,1), &
              x % interior(3,i,j,k,iEl,1), &
              solution % interior(i,j,k,iEl,1:this % solution % nvar), &
              solutionGradient % interior(1,i,j,k,iEl,1:this % solution % nvar), &
              solutionGradient % interior(2,i,j,k,iEl,1:this % solution % nvar), &
              solutionGradient % interior(3,i,j,k,iEl,1:this % solution % nvar)
  
          END DO
        END DO
        ENDDO
  
      END DO
  
      CLOSE (UNIT=fUnit)
  
      CALL x % Free()
      CALL solution % Free()
      CALL interp % Free()
  
    END SUBROUTINE WriteTecplot_DGModel3D
  
  END MODULE SELF_DGModel3D
  