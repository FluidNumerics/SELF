!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Model1D

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

  TYPE,EXTENDS(Model) :: Model1D
    TYPE(MappedScalar1D) :: solution
    TYPE(MappedScalar1D) :: solutionGradient
    TYPE(MappedScalar1D) :: velocity
    TYPE(MappedScalar1D) :: flux
    TYPE(MappedScalar1D) :: source
    TYPE(MappedScalar1D) :: fluxDivergence
    TYPE(MappedScalar1D) :: dSdt
    TYPE(MappedScalar1D) :: workSol
    TYPE(MappedScalar1D) :: prevSol
    TYPE(Mesh1D),POINTER :: mesh
    TYPE(Geometry1D),POINTER :: geometry

  CONTAINS

    PROCEDURE :: Init => Init_Model1D
    PROCEDURE :: Free => Free_Model1D

    PROCEDURE :: UpdateHost => UpdateHost_Model1D
    PROCEDURE :: UpdateDevice => UpdateDevice_Model1D

    PROCEDURE :: UpdateSolution => UpdateSolution_Model1D

    PROCEDURE :: ResizePrevSol => ResizePrevSol_Model1D

    PROCEDURE :: UpdateGAB2 => UpdateGAB2_Model1D
    PROCEDURE :: UpdateGAB3 => UpdateGAB3_Model1D
    PROCEDURE :: UpdateGAB4 => UpdateGAB4_Model1D

    PROCEDURE :: UpdateGRK2 => UpdateGRK2_Model1D
    PROCEDURE :: UpdateGRK3 => UpdateGRK3_Model1D
    PROCEDURE :: UpdateGRK4 => UpdateGRK4_Model1D
    PROCEDURE :: CalculateTendency => CalculateTendency_Model1D
    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_Model1D

    GENERIC :: SetSolution => SetSolutionFromChar_Model1D, &
      SetSolutionFromEqn_Model1D
    PROCEDURE,PRIVATE :: SetSolutionFromChar_Model1D
    PROCEDURE,PRIVATE :: SetSolutionFromEqn_Model1D

    GENERIC :: SetVelocityField => SetVelocityFieldFromChar_Model1D, &
      SetVelocityFieldFromEqn_Model1D
    PROCEDURE,PRIVATE :: SetVelocityFieldFromChar_Model1D
    PROCEDURE,PRIVATE :: SetVelocityFieldFromEqn_Model1D

    PROCEDURE :: ReadModel => Read_Model1D
    PROCEDURE :: WriteModel => Write_Model1D
    PROCEDURE :: WriteTecplot => WriteTecplot_Model1D

  END TYPE Model1D

  INTERFACE
    SUBROUTINE UpdateSolution_Model1D_gpu_wrapper(solution,dSdt,dt,N,nVar,nEl) &
      BIND(c,name="UpdateSolution_Model1D_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: solution,dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: dt
    END SUBROUTINE UpdateSolution_Model1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE UpdateGAB2_Model1D_gpu_wrapper(prevsol,solution,m,nPrev,N,nVar,nEl) &
      BIND(c,name="UpdateGAB2_Model1D_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: prevsol,solution
      INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
    END SUBROUTINE UpdateGAB2_Model1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE UpdateGAB3_Model1D_gpu_wrapper(prevsol,solution,m,nPrev,N,nVar,nEl) &
      BIND(c,name="UpdateGAB3_Model1D_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: prevsol,solution
      INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
    END SUBROUTINE UpdateGAB3_Model1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE UpdateGAB4_Model1D_gpu_wrapper(prevsol,solution,m,nPrev,N,nVar,nEl) &
      BIND(c,name="UpdateGAB4_Model1D_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: prevsol,solution
      INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
    END SUBROUTINE UpdateGAB4_Model1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE UpdateGRK_Model1D_gpu_wrapper(grk,solution,dSdt,rk_a,rk_g,dt,nWork,N,nVar,nEl) &
      BIND(c,name="UpdateGRK_Model1D_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: grk,solution,dSdt
      INTEGER(C_INT),VALUE :: nWork,N,nVar,nEl
      REAL(c_prec),VALUE :: rk_a,rk_g,dt
    END SUBROUTINE UpdateGRK_Model1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE CalculateDSDt_Model1D_gpu_wrapper(fluxDivergence,source,dSdt,N,nVar,nEl) &
      BIND(c,name="CalculateDSDt_Model1D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(C_PTR) :: fluxDivergence,source,dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE CalculateDSDt_Model1D_gpu_wrapper
  END INTERFACE

CONTAINS

  SUBROUTINE Init_Model1D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh1D),INTENT(in),TARGET :: mesh
    TYPE(Geometry1D),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.

    CALL this % solution % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % velocity % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nVar,this % mesh % nElem)

  END SUBROUTINE Init_Model1D

  SUBROUTINE Free_Model1D(this)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this

    CALL this % solution % Free()
    CALL this % workSol % Free()
    CALL this % prevSol % Free()
    CALL this % velocity % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_Model1D

  SUBROUTINE ResizePrevSol_Model1D(this,m)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: nVar

    ! Free space, if necessary
    CALL this % prevSol % Free()

    ! Reallocate with increased variable dimension for
    ! storing "m" copies of solution data
    nVar = this % solution % nVar
    CALL this % prevSol % Init(this % geometry % x % interp,m*nVar,this % mesh % nElem)

  END SUBROUTINE ResizePrevSol_Model1D

  SUBROUTINE UpdateHost_Model1D(this)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this

    CALL this % mesh % UpdateHost()
    CALL this % geometry % UpdateHost()
    CALL this % velocity % UpdateHost()
    CALL this % solution % UpdateHost()
    CALL this % dSdt % UpdateHost()
    CALL this % solutionGradient % UpdateHost()
    CALL this % flux % UpdateHost()
    CALL this % source % UpdateHost()
    CALL this % fluxDivergence % UpdateHost()

  END SUBROUTINE UpdateHost_Model1D

  SUBROUTINE UpdateDevice_Model1D(this)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this

    CALL this % mesh % UpdateDevice()
    CALL this % geometry % UpdateDevice()
    CALL this % dSdt % UpdateDevice()
    CALL this % solution % UpdateDevice()
    CALL this % velocity % UpdateDevice()
    CALL this % solutionGradient % UpdateDevice()
    CALL this % flux % UpdateDevice()
    CALL this % source % UpdateDevice()
    CALL this % fluxDivergence % UpdateDevice()

  END SUBROUTINE UpdateDevice_Model1D

  SUBROUTINE SetSolutionFromEqn_Model1D(this,eqn)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

    ! Copy the equation parser
    DO iVar = 1,this % solution % nVar
      CALL this % solution % SetEquation(ivar,eqn(iVar) % equation)
    END DO

    CALL this % solution % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % solution % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateDevice()
    END IF

  END SUBROUTINE SetSolutionFromEqn_Model1D

  SUBROUTINE SetSolutionFromChar_Model1D(this,eqnChar)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

    DO iVar = 1,this % solution % nVar
      PRINT *, iVar,eqnChar(iVar)
      CALL this % solution % SetEquation(ivar,eqnChar(iVar))
    END DO

    CALL this % solution % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % solution % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateDevice()
    END IF

  END SUBROUTINE SetSolutionFromChar_Model1D

  SUBROUTINE SetVelocityFieldFromEqn_Model1D(this,eqn)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn

    ! Copy the equation parser
    ! Set the x-component of the velocity
    CALL this % velocity % SetEquation(1,eqn % equation)

    ! Set the velocity values using the equation parser
    CALL this % velocity % SetInteriorFromEquation(this % geometry,this % t)

    CALL this % velocity % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % velocity % UpdateDevice()
    END IF

  END SUBROUTINE SetVelocityFieldFromEqn_Model1D

  SUBROUTINE SetVelocityFieldFromChar_Model1D(this,eqnChar)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar

    ! Set the x-component of the velocity
    CALL this % velocity % SetEquation(1,eqnChar)

    ! Set the velocity values using the equation parser
    CALL this % velocity % SetInteriorFromEquation(this % geometry,this % t)

    CALL this % velocity % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % velocity % UpdateDevice()
    END IF

  END SUBROUTINE SetVelocityFieldFromChar_Model1D

  SUBROUTINE UpdateSolution_Model1D(this,dt)
    !! Computes a solution update as , where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    REAL(prec),OPTIONAL,INTENT(in) :: dt
    ! Local
    REAL(prec) :: dtLoc
    INTEGER :: i,iVar,iEl

    IF (PRESENT(dt)) THEN
      dtLoc = dt
    ELSE
      dtLoc = this % dt
    END IF

    IF (this % gpuAccel) THEN

      CALL UpdateSolution_Model1D_gpu_wrapper(this % solution % interior % deviceData, &
                                              this % dSdt % interior % deviceData, &
                                              dtLoc, &
                                              this % solution % interp % N, &
                                              this % solution % nVar, &
                                              this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 0,this % solution % interp % N

            this % solution % interior % hostData(i,iVar,iEl) = &
              this % solution % interior % hostData(i,iVar,iEl) + &
              dtLoc*this % dSdt % interior % hostData(i,iVar,iEl)

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE UpdateSolution_Model1D

  SUBROUTINE UpdateGAB2_Model1D(this,m)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,nVar,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGAB2_Model1D_gpu_wrapper(this % prevSol % interior % deviceData, &
                                          this % solution % interior % deviceData, &
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
            DO i = 0,this % solution % interp % N

              this % prevSol % interior % hostData(i,iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSEIF (m == 1) THEN ! Copy the solution back from prevsol

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

              this % solution % interior % hostData(i,iVar,iEl) = this % prevSol % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSE ! Main looping section - nVar the previous solution, store the new solution, and
        ! create an interpolated solution to use for tendency calculation

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

              ! Bump the last solution
              nVar = this % solution % nVar
              this % prevSol % interior % hostData(i,nVar + iVar,iEl) = this % prevSol % interior % hostData(i,iVar,iEl)

              ! Store the new solution
              this % prevSol % interior % hostData(i,iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

              this % solution % interior % hostData(i,iVar,iEl) = &
                1.5_PREC*this % prevSol % interior % hostData(i,iVar,iEl) - &
                0.5_PREC*this % prevSol % interior % hostData(i,nVar + iVar,iEl)
            END DO
          END DO
        END DO
      END IF

    END IF

  END SUBROUTINE UpdateGAB2_Model1D

  SUBROUTINE UpdateGAB3_Model1D(this,m)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,nVar,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGAB3_Model1D_gpu_wrapper(this % prevSol % interior % deviceData, &
                                          this % solution % interior % deviceData, &
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
            DO i = 0,this % solution % interp % N

             this % prevSol % interior % hostData(i,nVar + iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSEIF (m == 1) THEN ! Initialization step - store the solution in the prevSol at ivar

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

              this % prevSol % interior % hostData(i,iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSEIF (m == 2) THEN ! Copy the solution back from the most recent prevsol

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

              this % solution % interior % hostData(i,iVar,iEl) = this % prevSol % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSE ! Main looping section - nVar the previous solution, store the new solution, and
        ! create an interpolated solution to use for tendency calculation

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

              ! Bump the last two stored solutions
              nVar = this % solution % nVar
     this % prevSol % interior % hostData(i,2*nVar + iVar,iEl) = this % prevSol % interior % hostData(i,nVar + iVar,iEl)
              this % prevSol % interior % hostData(i,nVar + iVar,iEl) = this % prevSol % interior % hostData(i,iVar,iEl)

              ! Store the new solution
              this % prevSol % interior % hostData(i,iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

              this % solution % interior % hostData(i,iVar,iEl) = &
                (23.0_PREC*this % prevSol % interior % hostData(i,iVar,iEl) - &
                 16.0_PREC*this % prevSol % interior % hostData(i,nVar + iVar,iEl) + &
                 5.0_PREC*this % prevSol % interior % hostData(i,2*nVar + iVar,iEl))/12.0_PREC

            END DO
          END DO
        END DO

      END IF

    END IF

  END SUBROUTINE UpdateGAB3_Model1D

  SUBROUTINE UpdateGAB4_Model1D(this,m)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,nVar,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGAB4_Model1D_gpu_wrapper(this % prevSol % interior % deviceData, &
                                          this % solution % interior % deviceData, &
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
            DO i = 0,this % solution % interp % N

           this % prevSol % interior % hostData(i,2*nVar + iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSEIF (m == 1) THEN ! Initialization step - store the solution in the prevSol at ivar

        nVar = this % solution % nVar
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

             this % prevSol % interior % hostData(i,nVar + iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSEIF (m == 2) THEN ! Initialization step - store the solution in the prevSol at ivar

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

              this % prevSol % interior % hostData(i,iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSEIF (m == 3) THEN ! Copy the solution back from the most recent prevsol

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

              this % solution % interior % hostData(i,iVar,iEl) = this % prevSol % interior % hostData(i,iVar,iEl)

            END DO
          END DO
        END DO

      ELSE ! Main looping section - nVar the previous solution, store the new solution, and
        ! create an interpolated solution to use for tendency calculation

        nVar = this % solution % nVar
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 0,this % solution % interp % N

              ! Bump the last two stored solutions
   this % prevSol % interior % hostData(i,3*nVar + iVar,iEl) = this % prevSol % interior % hostData(i,2*nVar + iVar,iEl)
     this % prevSol % interior % hostData(i,2*nVar + iVar,iEl) = this % prevSol % interior % hostData(i,nVar + iVar,iEl)
              this % prevSol % interior % hostData(i,nVar + iVar,iEl) = this % prevSol % interior % hostData(i,iVar,iEl)

              ! Store the new solution
              this % prevSol % interior % hostData(i,iVar,iEl) = this % solution % interior % hostData(i,iVar,iEl)

              this % solution % interior % hostData(i,iVar,iEl) = &
                (55.0_PREC*this % prevSol % interior % hostData(i,iVar,iEl) - &
                 59.0_PREC*this % prevSol % interior % hostData(i,nVar + iVar,iEl) + &
                 37.0_PREC*this % prevSol % interior % hostData(i,2*nVar + iVar,iEl) - &
                 9.0_PREC*this % prevSol % interior % hostData(i,3*nVar + iVar,iEl))/24.0_PREC

            END DO
          END DO
        END DO

      END IF

    END IF

  END SUBROUTINE UpdateGAB4_Model1D

  SUBROUTINE UpdateGRK2_Model1D(this,m)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGRK_Model1D_gpu_wrapper(this % workSol % interior % deviceData, &
                                         this % solution % interior % deviceData, &
                                         this % dSdt % interior % deviceData, &
                                         rk2_a(m),rk2_g(m),this % dt, &
                                         this % worksol % nVar, &
                                         this % solution % interp % N, &
                                         this % solution % nVar, &
                                         this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 0,this % solution % interp % N

            this % workSol % interior % hostData(i,iVar,iEl) = rk2_a(m)* &
                                                               this % workSol % interior % hostData(i,iVar,iEl) + &
                                                               this % dSdt % interior % hostData(i,iVar,iEl)

            this % solution % interior % hostData(i,iVar,iEl) = &
              this % solution % interior % hostData(i,iVar,iEl) + &
              rk2_g(m)*this % dt*this % workSol % interior % hostData(i,iVar,iEl)

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE UpdateGRK2_Model1D

  SUBROUTINE UpdateGRK3_Model1D(this,m)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGRK_Model1D_gpu_wrapper(this % workSol % interior % deviceData, &
                                         this % solution % interior % deviceData, &
                                         this % dSdt % interior % deviceData, &
                                         rk3_a(m),rk3_g(m),this % dt, &
                                         this % worksol % nVar, &
                                         this % solution % interp % N, &
                                         this % solution % nVar, &
                                         this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 0,this % solution % interp % N

            this % workSol % interior % hostData(i,iVar,iEl) = rk3_a(m)* &
                                                               this % workSol % interior % hostData(i,iVar,iEl) + &
                                                               this % dSdt % interior % hostData(i,iVar,iEl)

            this % solution % interior % hostData(i,iVar,iEl) = &
              this % solution % interior % hostData(i,iVar,iEl) + &
              rk3_g(m)*this % dt*this % workSol % interior % hostData(i,iVar,iEl)

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE UpdateGRK3_Model1D

  SUBROUTINE UpdateGRK4_Model1D(this,m)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGRK_Model1D_gpu_wrapper(this % workSol % interior % deviceData, &
                                         this % solution % interior % deviceData, &
                                         this % dSdt % interior % deviceData, &
                                         rk4_a(m),rk4_g(m),this % dt, &
                                         this % worksol % nVar, &
                                         this % solution % interp % N, &
                                         this % solution % nVar, &
                                         this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 0,this % solution % interp % N

            this % workSol % interior % hostData(i,iVar,iEl) = rk4_a(m)* &
                                                               this % workSol % interior % hostData(i,iVar,iEl) + &
                                                               this % dSdt % interior % hostData(i,iVar,iEl)

            this % solution % interior % hostData(i,iVar,iEl) = &
              this % solution % interior % hostData(i,iVar,iEl) + &
              rk4_g(m)*this % dt*this % workSol % interior % hostData(i,iVar,iEl)

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE UpdateGRK4_Model1D

  SUBROUTINE CalculateFluxDivergence_Model1D(this)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this

    CALL this % flux % Derivative(this % geometry, &
                                  this % fluxDivergence, &
                                  selfWeakDGForm, &
                                  this % gpuAccel)

  END SUBROUTINE CalculateFluxDivergence_Model1D

  SUBROUTINE CalculateTendency_Model1D(this)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iVar,iEl

    CALL this % PreTendency()
    CALL this % solution % BoundaryInterp(this % gpuAccel)
    CALL this % solution % SideExchange(this % mesh,this % decomp,this % gpuAccel)
    CALL this % SetBoundaryCondition()
    CALL this % SourceMethod()
    CALL this % RiemannSolver()
    CALL this % FluxMethod()
    CALL this % CalculateFluxDivergence()

    IF (this % gpuAccel) THEN

      CALL CalculateDSDt_Model1D_gpu_wrapper(this % fluxDivergence % interior % deviceData, &
                                             this % source % interior % deviceData, &
                                             this % dSdt % interior % deviceData, &
                                             this % solution % interp % N, &
                                             this % solution % nVar, &
                                             this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 0,this % solution % interp % N

            this % dSdt % interior % hostData(i,iVar,iEl) = &
              this % source % interior % hostData(i,iVar,iEl) - &
              this % fluxDivergence % interior % hostData(i,iVar,iEl)

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE CalculateTendency_Model1D

  SUBROUTINE Write_Model1D(this,fileName)
#undef __FUNC__
#define __FUNC__ "Write_Model1D"
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    TYPE(Scalar1D) :: solution
    TYPE(Scalar1D) :: x
    TYPE(Lagrange),TARGET :: interp
    CHARACTER(LEN=self_FileNameLength) :: pickupFile
    CHARACTER(13) :: timeStampString

    IF (PRESENT(filename)) THEN
      pickupFile = filename
    ELSE
      WRITE (timeStampString,'(I13.13)') this % ioIterate
      pickupFile = 'solution.'//timeStampString//'.h5'
    END IF

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateHost()
      CALL this % solutionGradient % UpdateHost()
    END IF

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
      this % decomp % offsetElem % hostData(this % decomp % rankId), this % decomp % nElem )

      ! Write the geometry to file
      INFO("Writing control grid geometry to file")
      CALL CreateGroup_HDF5(fileId,'/controlgrid/geometry')
      CALL this % geometry % x % WriteHDF5( fileId, '/controlgrid/geometry/x', &
      this % decomp % offsetElem % hostData(this % decomp % rankId), this % decomp % nElem )

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
      CALL this % geometry % x % GridInterp(x,gpuAccel=.FALSE.)

      ! Map the solution to the target grid
      CALL this % solution % GridInterp(solution,gpuAccel=.FALSE.)

      ! Write the model state to file
      CALL CreateGroup_HDF5(fileId,'/targetgrid')
      CALL solution % WriteHDF5( fileId, '/targetgrid/solution', &
      this % decomp % offsetElem % hostData(this % decomp % rankId), this % decomp % nElem )

      ! Write the geometry to file
      CALL CreateGroup_HDF5(fileId,'/targetgrid/mesh')
      CALL x % WriteHDF5( fileId, '/targetgrid/mesh/coords', &
      this % decomp % offsetElem % hostData(this % decomp % rankId), this % decomp % nElem )

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
      CALL this % geometry % x % GridInterp(x,gpuAccel=.FALSE.)

      ! Map the solution to the target grid
      CALL this % solution % GridInterp(solution,gpuAccel=.FALSE.)

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

  END SUBROUTINE Write_Model1D

  SUBROUTINE Read_Model1D(this,fileName)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:3)
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
      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId) + 1
      solOffset(1:3) = (/0,1,firstElem/)
      CALL ReadArray_HDF5(fileId,'/controlgrid/solution/interior', &
                          this % solution % interior,solOffset)
    ELSE
      CALL ReadArray_HDF5(fileId,'/controlgrid/solution/interior',this % solution % interior)
    END IF

    CALL Close_HDF5(fileId)
    IF (this % gpuAccel) THEN
      CALL this % solution % interior % UpdateDevice()
    END IF

  END SUBROUTINE Read_Model1D

  SUBROUTINE WriteTecplot_Model1D(this,filename)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in),OPTIONAL :: filename
    ! Local
    CHARACTER(8) :: zoneID
    INTEGER :: fUnit
    INTEGER :: iEl,i,iVar
    CHARACTER(LEN=self_FileNameLength) :: tecFile
    CHARACTER(LEN=self_TecplotHeaderLength) :: tecHeader
    CHARACTER(LEN=self_FormatLength) :: fmat
    CHARACTER(13) :: timeStampString
    CHARACTER(5) :: rankString
    TYPE(Scalar1D) :: solution
    TYPE(Scalar1D) :: x
    TYPE(Lagrange),TARGET :: interp

    IF (PRESENT(filename)) THEN
      tecFile = filename
    ELSE
      ! Create a 0-padded integer for the output iterate
      WRITE (timeStampString,'(I13.13)') this % ioIterate
      ! Increment the ioIterate
      this % ioIterate = this % ioIterate + 1
      !timeStampString = TimeStamp(this % t,'s')
      ! IF (this % decomp % mpiEnabled) THEN
      !   WRITE (rankString,'(I5.5)') this % decomp % rankId
      !   tecFile = 'solution.'//rankString//'.'//timeStampString//'.curve'
      ! ELSE
      tecFile = 'solution.'//timeStampString//'.curve'
      ! END IF

    END IF

    ! Create an interpolant for the uniform grid
    CALL interp % Init(this % solution % interp % M, &
                       this % solution % interp % targetNodeType, &
                       this % solution % interp % N, &
                       this % solution % interp % controlNodeType)

    CALL solution % Init(interp, &
                         this % solution % nVar,this % solution % nElem)

    CALL x % Init(interp,1,this % solution % nElem)

    ! Map the mesh positions to the target grid
    CALL this % geometry % x % GridInterp(x,gpuAccel=.FALSE.)

    ! Map the solution to the target grid
    CALL this % solution % GridInterp(solution,gpuAccel=.FALSE.)

    fmat = '(2(ES16.7E3,1x))'
    ! Let's write some tecplot!!
    OPEN (UNIT=NEWUNIT(fUnit), &
          FILE=TRIM(tecFile), &
          FORM='formatted', &
          STATUS='replace')

    DO iVar = 1,this % solution % nVar
      WRITE (tecHeader,'(E15.6)') this % t
      tecHeader = "#TIME "//TRIM(tecHeader)
      WRITE (fUnit,*) TRIM(tecHeader)

      tecHeader = "#"//TRIM(this % solution % meta(iVar) % name)//" vs position"
      WRITE (fUnit,*) TRIM(tecHeader)
      DO iEl = 1,this % solution % nElem

        !WRITE (zoneID,'(I8.8)') iEl
        !WRITE (fUnit,*) 'ZONE T="el'//TRIM(zoneID)//'", I=',this % solution % interp % M + 1

        DO i = 0,this % solution % interp % M

          WRITE (fUnit,fmat) x % interior % hostData(i,1,iEl), &
            solution % interior % hostData(i,ivar,iEl)

        END DO

      END DO
    END DO

    CLOSE (UNIT=fUnit)

    CALL x % Free()
    CALL solution % Free()
    CALL interp % Free()

  END SUBROUTINE WriteTecplot_Model1D

END MODULE SELF_Model1D
