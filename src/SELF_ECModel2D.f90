!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_ECModel2D

  USE SELF_SupportRoutines
  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_HDF5
  USE HDF5
  USE FEQParse
  USE SELF_Model
  USE SELF_Model2D

  TYPE,EXTENDS(Model) :: ECModel2D
    TYPE(MappedScalar2D)   :: solution
    TYPE(MappedVector2D)   :: solutionGradient
    TYPE(MappedP2Vector2D) :: flux
    TYPE(MappedScalar2D)   :: source
    TYPE(MappedScalar2D)   :: fluxDivergence
    TYPE(MappedScalar2D)   :: dSdt
    TYPE(MappedScalar2D)   :: workSol
    TYPE(MappedScalar2D)   :: prevSol
    TYPE(Mesh2D),POINTER   :: mesh
    TYPE(SEMQuad),POINTER  :: geometry

  CONTAINS

    PROCEDURE :: Init => Init_ECModel2D
    PROCEDURE :: Free => Free_ECModel2D

    PROCEDURE :: UpdateHost => UpdateHost_ECModel2D
    PROCEDURE :: UpdateDevice => UpdateDevice_ECModel2D

    PROCEDURE :: UpdateSolution => UpdateSolution_ECModel2D

    PROCEDURE :: ResizePrevSol => ResizePrevSol_ECModel2D

    PROCEDURE :: UpdateGAB2 => UpdateGAB2_ECModel2D
    PROCEDURE :: UpdateGAB3 => UpdateGAB3_ECModel2D
    PROCEDURE :: UpdateGAB4 => UpdateGAB4_ECModel2D

    PROCEDURE :: UpdateGRK2 => UpdateGRK2_ECModel2D
    PROCEDURE :: UpdateGRK3 => UpdateGRK3_ECModel2D
    PROCEDURE :: UpdateGRK4 => UpdateGRK4_ECModel2D

    PROCEDURE :: CalculateTendency => CalculateTendency_ECModel2D
    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_ECModel2D

    GENERIC :: SetSolution => SetSolutionFromChar_ECModel2D, &
      SetSolutionFromEqn_ECModel2D
    PROCEDURE,PRIVATE :: SetSolutionFromChar_ECModel2D
    PROCEDURE,PRIVATE :: SetSolutionFromEqn_ECModel2D

    PROCEDURE :: UpdateBoundary => UpdateBoundary_ECModel2D
    PROCEDURE :: ReprojectFlux => ReprojectFlux_ECModel2D

    PROCEDURE :: ReadModel => Read_ECModel2D
    PROCEDURE :: WriteModel => Write_ECModel2D
    PROCEDURE :: WriteTecplot => WriteTecplot_ECModel2D

  END TYPE ECModel2D

  ! INTERFACE
  !   SUBROUTINE UpdateSolution_Model2D_gpu_wrapper(solution, dSdt, dt, N, nVar, nEl) &
  !     bind(c,name="UpdateSolution_Model2D_gpu_wrapper")
  !     USE iso_c_binding
  !     USE SELF_Constants
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: solution, dSdt
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !     REAL(c_prec),VALUE :: dt
  !   END SUBROUTINE UpdateSolution_Model2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE UpdateGAB2_Model2D_gpu_wrapper(prevsol, solution, m, nPrev, N, nVar, nEl) &
  !     bind(c,name="UpdateGAB2_Model2D_gpu_wrapper")
  !     USE iso_c_binding
  !     USE SELF_Constants
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: prevsol, solution
  !     INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
  !   END SUBROUTINE UpdateGAB2_Model2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE UpdateGAB3_Model2D_gpu_wrapper(prevsol, solution, m, nPrev, N, nVar, nEl) &
  !     bind(c,name="UpdateGAB3_Model2D_gpu_wrapper")
  !     USE iso_c_binding
  !     USE SELF_Constants
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: prevsol, solution
  !     INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
  !   END SUBROUTINE UpdateGAB3_Model2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE UpdateGAB4_Model2D_gpu_wrapper(prevsol, solution, m, nPrev, N, nVar, nEl) &
  !     bind(c,name="UpdateGAB4_Model2D_gpu_wrapper")
  !     USE iso_c_binding
  !     USE SELF_Constants
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: prevsol, solution
  !     INTEGER(C_INT),VALUE :: m,nPrev,N,nVar,nEl
  !   END SUBROUTINE UpdateGAB4_Model2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE UpdateGRK_Model2D_gpu_wrapper(grk, solution, dSdt, rk_a, rk_g, dt, nWork, N, nVar, nEl) &
  !     bind(c,name="UpdateGRK_Model2D_gpu_wrapper")
  !     USE iso_c_binding
  !     USE SELF_Constants
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: grk, solution, dSdt
  !     INTEGER(C_INT),VALUE :: nWork,N,nVar,nEl
  !     REAL(c_prec),VALUE :: rk_a, rk_g, dt
  !   END SUBROUTINE UpdateGRK_Model2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE CalculateDSDt_Model2D_gpu_wrapper(fluxDivergence, source, dSdt, N, nVar, nEl) &
  !     bind(c,name="CalculateDSDt_Model2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: fluxDivergence, source, dSdt
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE CalculateDSDt_Model2D_gpu_wrapper
  ! END INTERFACE

CONTAINS

  SUBROUTINE Init_ECModel2D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh2D),INTENT(in),TARGET :: mesh
    TYPE(SEMQuad),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp
    ! Local
    INTEGER :: ivar
    CHARACTER(LEN=3) :: ivarChar
    CHARACTER(LEN=25) :: varname

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.

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

  END SUBROUTINE Init_ECModel2D

  SUBROUTINE Free_ECModel2D(this)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this

    CALL this % solution % Free()
    CALL this % workSol % Free()
    CALL this % prevSol % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_ECModel2D

  SUBROUTINE ResizePrevSol_ECModel2D(this,m)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: nVar

    ! Free space, if necessary
    CALL this % prevSol % Free()

    ! Reallocate with increased variable dimension for
    ! storing "m" copies of solution data
    nVar = this % solution % nVar
    CALL this % prevSol % Init(this % geometry % x % interp,m*nVar,this % mesh % nElem)

  END SUBROUTINE ResizePrevSol_ECModel2D

  SUBROUTINE UpdateHost_ECModel2D(this)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this

    CALL this % mesh % UpdateHost()
    CALL this % geometry % UpdateHost()
    CALL this % solution % UpdateHost()
    CALL this % dSdt % UpdateHost()
    CALL this % solution % UpdateHost()
    CALL this % solutionGradient % UpdateHost()
    CALL this % flux % UpdateHost()
    CALL this % source % UpdateHost()
    CALL this % fluxDivergence % UpdateHost()

  END SUBROUTINE UpdateHost_ECModel2D

  SUBROUTINE UpdateDevice_ECModel2D(this)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this

    CALL this % mesh % UpdateDevice()
    CALL this % geometry % UpdateDevice()
    CALL this % dSdt % UpdateDevice()
    CALL this % solution % UpdateDevice()
    CALL this % solutionGradient % UpdateDevice()
    CALL this % flux % UpdateDevice()
    CALL this % source % UpdateDevice()
    CALL this % fluxDivergence % UpdateDevice()

  END SUBROUTINE UpdateDevice_ECModel2D

  SUBROUTINE SetSolutionFromEqn_ECModel2D(this,eqn)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
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

  END SUBROUTINE SetSolutionFromEqn_ECModel2D

  ! SUBROUTINE SetVelocityFieldFromEqn_ECModel2D(this, eqn)
  !   IMPLICIT NONE
  !   CLASS(ECModel2D),INTENT(inout) :: this
  !   TYPE(EquationParser),INTENT(in) :: eqn(1:2)

  !     ! Copy the equation parser
  !     ! Set the x-component of the velocity
  !     CALL this % velocity % SetEquation(1,1,eqn(1) % equation)

  !     ! Set the y-component of the velocity
  !     CALL this % velocity % SetEquation(2,1,eqn(2) % equation)

  !     ! Set the velocity values using the equation parser
  !     CALL this % velocity % SetInteriorFromEquation( this % geometry, this % t )

  !     CALL this % velocity % BoundaryInterp( gpuAccel = .FALSE. )

  !     IF( this % gpuAccel )THEN
  !       CALL this % velocity % UpdateDevice()
  !     ENDIF

  ! END SUBROUTINE SetVelocityFieldFromEqn_ECModel2D

  ! SUBROUTINE SetVelocityFieldFromChar_ECModel2D(this, eqnChar)
  !   IMPLICIT NONE
  !   CLASS(ECModel2D),INTENT(inout) :: this
  !   CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:2)

  !     ! Set the x-component of the velocity
  !     CALL this % velocity % SetEquation(1,1,eqnChar(1))

  !     ! Set the y-component of the velocity
  !     CALL this % velocity % SetEquation(2,1,eqnChar(2))

  !     ! Set the velocity values using the equation parser
  !     CALL this % velocity % SetInteriorFromEquation( this % geometry, this % t )

  !     CALL this % velocity % BoundaryInterp( gpuAccel = .FALSE. )

  !     IF( this % gpuAccel )THEN
  !       CALL this % velocity % UpdateDevice()
  !     ENDIF

  ! END SUBROUTINE SetVelocityFieldFromChar_ECModel2D

  SUBROUTINE SetSolutionFromChar_ECModel2D(this,eqnChar)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

    DO iVar = 1,this % solution % nVar
      CALL this % solution % SetEquation(ivar,eqnChar(iVar))
    END DO

    CALL this % solution % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % solution % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateDevice()
    END IF

  END SUBROUTINE SetSolutionFromChar_ECModel2D

  SUBROUTINE UpdateSolution_ECModel2D(this,dt)
    !! Computes a solution update as , where dt is either provided through the interface
    !! or taken as the ECModel's stored time step size (model % dt)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    REAL(prec),OPTIONAL,INTENT(in) :: dt
    ! Local
    REAL(prec) :: dtLoc
    INTEGER :: i,j,iVar,iEl

    IF (PRESENT(dt)) THEN
      dtLoc = dt
    ELSE
      dtLoc = this % dt
    END IF

    IF (this % gpuAccel) THEN

      CALL UpdateSolution_Model2D_gpu_wrapper(this % solution % interior % deviceData, &
                                              this % dSdt % interior % deviceData, &
                                              dtLoc, &
                                              this % solution % interp % N, &
                                              this % solution % nVar, &
                                              this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO j = 0,this % solution % interp % N
            DO i = 0,this % solution % interp % N

              this % solution % interior % hostData(i,j,iVar,iEl) = &
                this % solution % interior % hostData(i,j,iVar,iEl) + &
                dtLoc*this % dSdt % interior % hostData(i,j,iVar,iEl)

            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE UpdateSolution_ECModel2D

  SUBROUTINE UpdateGAB2_ECModel2D(this,m)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,j,nVar,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGAB2_Model2D_gpu_wrapper(this % prevSol % interior % deviceData, &
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
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                this % prevSol % interior % hostData(i,j,iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSEIF (m == 1) THEN ! Reset solution

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                this % solution % interior % hostData(i,j,iVar,iEl) = this % prevSol % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSE ! Main looping section - nVar the previous solution, store the new solution, and
        ! create an interpolated solution to use for tendency calculation

        nVar = this % solution % nVar
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                ! Bump the last solution
          this % prevSol % interior % hostData(i,j,nVar + iVar,iEl) = this % prevSol % interior % hostData(i,j,iVar,iEl)

                ! Store the new solution
                this % prevSol % interior % hostData(i,j,iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

                this % solution % interior % hostData(i,j,iVar,iEl) = &
                  1.5_PREC*this % prevSol % interior % hostData(i,j,iVar,iEl) - &
                  0.5_PREC*this % prevSol % interior % hostData(i,j,nVar + iVar,iEl)
              END DO
            END DO
          END DO
        END DO

      END IF

    END IF

  END SUBROUTINE UpdateGAB2_ECModel2D

  SUBROUTINE UpdateGAB3_ECModel2D(this,m)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,j,nVar,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGAB3_Model2D_gpu_wrapper(this % prevSol % interior % deviceData, &
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
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

         this % prevSol % interior % hostData(i,j,nVar + iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSEIF (m == 1) THEN ! Initialization step - store the solution in the prevSol at ivar

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                this % prevSol % interior % hostData(i,j,iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSEIF (m == 2) THEN ! Copy the solution back from the most recent prevsol

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                this % solution % interior % hostData(i,j,iVar,iEl) = this % prevSol % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSE ! Main looping section - nVar the previous solution, store the new solution, and
        ! create an interpolated solution to use for tendency calculation

        nVar = this % solution % nVar
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                ! Bump the last two stored solutions
 this % prevSol % interior % hostData(i,j,2*nVar + iVar,iEl) = this % prevSol % interior % hostData(i,j,nVar + iVar,iEl)
          this % prevSol % interior % hostData(i,j,nVar + iVar,iEl) = this % prevSol % interior % hostData(i,j,iVar,iEl)

                ! Store the new solution
                this % prevSol % interior % hostData(i,j,iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

                this % solution % interior % hostData(i,j,iVar,iEl) = &
                  (23.0_PREC*this % prevSol % interior % hostData(i,j,iVar,iEl) - &
                   16.0_PREC*this % prevSol % interior % hostData(i,j,nVar + iVar,iEl) + &
                   5.0_PREC*this % prevSol % interior % hostData(i,j,2*nVar + iVar,iEl))/12.0_PREC

              END DO
            END DO
          END DO
        END DO
      END IF

    END IF

  END SUBROUTINE UpdateGAB3_ECModel2D

  SUBROUTINE UpdateGAB4_ECModel2D(this,m)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,j,nVar,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGAB4_Model2D_gpu_wrapper(this % prevSol % interior % deviceData, &
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
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

       this % prevSol % interior % hostData(i,j,2*nVar + iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSEIF (m == 1) THEN ! Initialization step - store the solution in the prevSol at ivar

        nVar = this % solution % nVar
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

         this % prevSol % interior % hostData(i,j,nVar + iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSEIF (m == 2) THEN ! Initialization step - store the solution in the prevSol at ivar

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                this % prevSol % interior % hostData(i,j,iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSEIF (m == 3) THEN ! Copy the solution back from the most recent prevsol

        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                this % solution % interior % hostData(i,j,iVar,iEl) = this % prevSol % interior % hostData(i,j,iVar,iEl)

              END DO
            END DO
          END DO
        END DO

      ELSE ! Main looping section - nVar the previous solution, store the new solution, and
        ! create an interpolated solution to use for tendency calculation

        nVar = this % solution % nVar
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO j = 0,this % solution % interp % N
              DO i = 0,this % solution % interp % N

                ! Bump the last two stored solutions
   this % prevSol % interior % hostData(i,j,3*nVar+iVar,iEl) = this % prevSol % interior % hostData(i,j,2*nVar+iVar,iEl)
 this % prevSol % interior % hostData(i,j,2*nVar + iVar,iEl) = this % prevSol % interior % hostData(i,j,nVar + iVar,iEl)
          this % prevSol % interior % hostData(i,j,nVar + iVar,iEl) = this % prevSol % interior % hostData(i,j,iVar,iEl)

                ! Store the new solution
                this % prevSol % interior % hostData(i,j,iVar,iEl) = this % solution % interior % hostData(i,j,iVar,iEl)

                this % solution % interior % hostData(i,j,iVar,iEl) = &
                  (55.0_PREC*this % prevSol % interior % hostData(i,j,iVar,iEl) - &
                   59.0_PREC*this % prevSol % interior % hostData(i,j,nVar + iVar,iEl) + &
                   37.0_PREC*this % prevSol % interior % hostData(i,j,2*nVar + iVar,iEl) - &
                   9.0_PREC*this % prevSol % interior % hostData(i,j,3*nVar + iVar,iEl))/24.0_PREC

              END DO
            END DO
          END DO
        END DO

      END IF

    END IF

  END SUBROUTINE UpdateGAB4_ECModel2D

  SUBROUTINE UpdateGRK2_ECModel2D(this,m)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,j,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGRK_Model2D_gpu_wrapper(this % workSol % interior % deviceData, &
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
          DO j = 0,this % solution % interp % N
            DO i = 0,this % solution % interp % N

              this % workSol % interior % hostData(i,j,iVar,iEl) = rk2_a(m)* &
                                                                  this % workSol % interior % hostData(i,j,iVar,iEl) + &
                                                                   this % dSdt % interior % hostData(i,j,iVar,iEl)

              this % solution % interior % hostData(i,j,iVar,iEl) = &
                this % solution % interior % hostData(i,j,iVar,iEl) + &
                rk2_g(m)*this % dt*this % workSol % interior % hostData(i,j,iVar,iEl)

            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE UpdateGRK2_ECModel2D

  SUBROUTINE UpdateGRK3_ECModel2D(this,m)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,j,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGRK_Model2D_gpu_wrapper(this % workSol % interior % deviceData, &
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
          DO j = 0,this % solution % interp % N
            DO i = 0,this % solution % interp % N

              this % workSol % interior % hostData(i,j,iVar,iEl) = rk3_a(m)* &
                                                                  this % workSol % interior % hostData(i,j,iVar,iEl) + &
                                                                   this % dSdt % interior % hostData(i,j,iVar,iEl)

              this % solution % interior % hostData(i,j,iVar,iEl) = &
                this % solution % interior % hostData(i,j,iVar,iEl) + &
                rk3_g(m)*this % dt*this % workSol % interior % hostData(i,j,iVar,iEl)

            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE UpdateGRK3_ECModel2D

  SUBROUTINE UpdateGRK4_ECModel2D(this,m)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,j,iVar,iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGRK_Model2D_gpu_wrapper(this % workSol % interior % deviceData, &
                                         this % solution % interior % deviceData, &
                                         this % dSdt % interior % deviceData, &
                                         rk4_a(m),rk4_g(m),this % dt, &
                                         this % workSol % nVar, &
                                         this % solution % interp % N, &
                                         this % solution % nVar, &
                                         this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO j = 0,this % solution % interp % N
            DO i = 0,this % solution % interp % N

              this % workSol % interior % hostData(i,j,iVar,iEl) = rk4_a(m)* &
                                                                  this % workSol % interior % hostData(i,j,iVar,iEl) + &
                                                                   this % dSdt % interior % hostData(i,j,iVar,iEl)

              this % solution % interior % hostData(i,j,iVar,iEl) = &
                this % solution % interior % hostData(i,j,iVar,iEl) + &
                rk4_g(m)*this % dt*this % workSol % interior % hostData(i,j,iVar,iEl)

            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE UpdateGRK4_ECModel2D

  SUBROUTINE ReprojectFlux_ECModel2D(this)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this

    CALL this % flux % ContravariantProjection(this % geometry,this % gpuAccel)

  END SUBROUTINE ReprojectFlux_ECModel2D

  SUBROUTINE CalculateFluxDivergence_ECModel2D(this)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this

    CALL this % flux % Divergence(this % geometry, &
                                  this % fluxDivergence, &
                                  selfWeakDGForm, &
                                  this % gpuAccel)

  END SUBROUTINE CalculateFluxDivergence_ECModel2D

  SUBROUTINE UpdateBoundary_ECModel2D(this)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this

    CALL this % solution % BoundaryInterp(this % gpuAccel)
    CALL this % solution % SideExchange(this % mesh,this % decomp,this % gpuAccel)

  END SUBROUTINE UpdateBoundary_ECModel2D

  SUBROUTINE CalculateTendency_ECModel2D(this)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    INTEGER :: i,j,iVar,iEl

    CALL this % PreTendency()
    CALL this % UpdateBoundary()
    CALL this % SetBoundaryCondition()
    CALL this % SourceMethod()
    CALL this % RiemannSolver()
    CALL this % FluxMethod()
    CALL this % ReprojectFlux()
    CALL this % CalculateFluxDivergence()

    IF (this % gpuAccel) THEN

      CALL CalculateDSDt_Model2D_gpu_wrapper(this % fluxDivergence % interior % deviceData, &
                                             this % source % interior % deviceData, &
                                             this % dSdt % interior % deviceData, &
                                             this % solution % interp % N, &
                                             this % solution % nVar, &
                                             this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO j = 0,this % solution % interp % N
            DO i = 0,this % solution % interp % N

              this % dSdt % interior % hostData(i,j,iVar,iEl) = &
                this % source % interior % hostData(i,j,iVar,iEl) - &
                this % fluxDivergence % interior % hostData(i,j,iVar,iEl)

            END DO
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE CalculateTendency_ECModel2D

  SUBROUTINE Write_ECModel2D(this,fileName)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:4)
    INTEGER(HID_T) :: xOffset(1:5)
    INTEGER(HID_T) :: vOffset(1:5)
    INTEGER(HID_T) :: bOffset(1:4)
    INTEGER(HID_T) :: bxOffset(1:5)
    INTEGER(HID_T) :: bvOffset(1:5)
    INTEGER(HID_T) :: solGlobalDims(1:4)
    INTEGER(HID_T) :: xGlobalDims(1:5)
    INTEGER(HID_T) :: vGlobalDims(1:5)
    INTEGER(HID_T) :: bGlobalDims(1:4)
    INTEGER(HID_T) :: bxGlobalDims(1:5)
    INTEGER(HID_T) :: bvGlobalDims(1:5)
    TYPE(Scalar2D) :: solution
    TYPE(Vector2D) :: x
    TYPE(Lagrange),TARGET :: interp
    INTEGER :: firstElem
    INTEGER :: ivar
    CHARACTER(LEN=self_FileNameLength) :: pickupFile
    CHARACTER(13) :: timeStampString
    CHARACTER(4) :: varNumber


    IF (PRESENT(filename)) THEN
      pickupFile = filename
    ELSE
      timeStampString = TimeStamp(this % t,'s')
      pickupFile = 'solution.'//timeStampString//'.h5'
    END IF

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateHost()
      CALL this % solutionGradient % UpdateHost()
    END IF

    IF (this % decomp % mpiEnabled) THEN

      CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId,this % decomp % mpiComm)

      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)
      solOffset(1:4) = (/0,0,0,firstElem/)
      solGlobalDims(1:4) = (/this % solution % interp % N + 1, &
                             this % solution % interp % N + 1, &
                             this % solution % nVar, &
                             this % decomp % nElem/)

      xOffset(1:5) = (/0,0,0,0,firstElem/)
      xGlobalDims(1:5) = (/2, &
                           this % solution % interp % N + 1, &
                           this % solution % interp % N + 1, &
                           1, &
                           this % decomp % nElem/)

      vOffset(1:5) = (/0,0,0,0,firstElem/)
      vGlobalDims(1:5) = (/2, &
                           this % solution % interp % N + 1, &
                           this % solution % interp % N + 1, &
                           this % solution % nVar, &
                           this % decomp % nElem/)

      ! Offsets and dimensions for element boundary data
      bOffset(1:4) = (/0,0,0,firstElem/)
      bGlobalDims(1:4) = (/this % solution % interp % N + 1, &
                           this % solution % nVar, &
                           4, &
                           this % decomp % nElem/)

      bxOffset(1:5) = (/0,0,0,0,firstElem/)
      bxGlobalDims(1:5) = (/2, &
                            this % solution % interp % N + 1, &
                            1, &
                            4, &
                            this % decomp % nElem/)

      bvOffset(1:5) = (/0,0,0,0,firstElem/)
      bvGlobalDims(1:5) = (/2, &
                            this % solution % interp % N + 1, &
                            this % solution % nVar, &
                            4, &
                            this % decomp % nElem/)

      CALL CreateGroup_HDF5(fileId,'/quadrature')

      CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
                           this % solution % interp % controlPoints)

      CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
                           this % solution % interp % qWeights)

      CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
                           this % solution % interp % dgMatrix)

      CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
                           this % solution % interp % dMatrix)

      CALL WriteArray_HDF5(fileId,'/quadrature/bmatrix', &
                           this % solution % interp % bMatrix)

      CALL WriteArray_HDF5(fileId,'/quadrature/imatrix', &
                           this % solution % interp % iMatrix)

! Add variable names to the file
                           CALL CreateGroup_HDF5(fileId,'/metadata')
                           CALL CreateGroup_HDF5(fileId,'/metadata/vars')
                           DO ivar = 1, this % solution % nvar
                             WRITE (varNumber,'(I0)') ivar
                             CALL WriteCharacter_HDF5(fileId, "/metadata/vars/"//TRIM(varnumber), &
                                     TRIM(this % solution % meta(ivar) % name))
                           ENDDO
                     
      CALL CreateGroup_HDF5(fileId,'/state')
      CALL CreateGroup_HDF5(fileId,'/state/solution')

      CALL CreateGroup_HDF5(fileId,'/mesh')
      CALL CreateGroup_HDF5(fileId,'/mesh/coords')

      CALL WriteArray_HDF5(fileId,'/state/solution/interior', &
                           this % solution % interior,solOffset,solGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/solution/boundary', &
                           this % solution % boundary,bOffset,bGlobalDims)

      CALL WriteArray_HDF5(fileId,'/mesh/coords/interior', &
                           this % geometry % x % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/mesh/coords/boundary', &
                           this % geometry % x % boundary,bxOffset,bxGlobalDims)

      ! Write data on plotting mesh
      solOffset(1:4) = (/0,0,0,firstElem/)
      solGlobalDims(1:4) = (/this % solution % interp % M + 1, &
                             this % solution % interp % M + 1, &
                             this % solution % nVar, &
                             this % decomp % nElem/)

      xOffset(1:5) = (/0,0,0,0,firstElem/)
      xGlobalDims(1:5) = (/2, &
                           this % solution % interp % M + 1, &
                           this % solution % interp % M + 1, &
                           1, &
                           this % decomp % nElem/)

      vOffset(1:5) = (/0,0,0,0,firstElem/)
      vGlobalDims(1:5) = (/2, &
                           this % solution % interp % M + 1, &
                           this % solution % interp % M + 1, &
                           this % solution % nVar, &
                           this % decomp % nElem/)

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

      CALL CreateGroup_HDF5(fileId,'/plot')
      CALL CreateGroup_HDF5(fileId,'/plot/mesh')
      CALL CreateGroup_HDF5(fileId,'/plot/mesh/coords')
      CALL CreateGroup_HDF5(fileId,'/plot/state')
      CALL CreateGroup_HDF5(fileId,'/plot/state/solution')

      CALL WriteArray_HDF5(fileId,'/plot/mesh/coords/interior', &
                           x % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/plot/state/solution/interior', &
                           solution % interior,solOffset,solGlobalDims)

      CALL Close_HDF5(fileId)

    ELSE

      CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

      CALL CreateGroup_HDF5(fileId,'/quadrature')

      CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
                           this % solution % interp % controlPoints)

      CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
                           this % solution % interp % qWeights)

      CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
                           this % solution % interp % dgMatrix)

      CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
                           this % solution % interp % dMatrix)

      CALL WriteArray_HDF5(fileId,'/quadrature/bmatrix', &
                           this % solution % interp % bMatrix)

      CALL WriteArray_HDF5(fileId,'/quadrature/imatrix', &
                           this % solution % interp % iMatrix)

! Add variable names to the file
                           CALL CreateGroup_HDF5(fileId,'/metadata')
                           CALL CreateGroup_HDF5(fileId,'/metadata/vars')
                           DO ivar = 1, this % solution % nvar
                             WRITE (varNumber,'(I0)') ivar
                             CALL WriteCharacter_HDF5(fileId, "/metadata/vars/"//TRIM(varnumber), &
                                     TRIM(this % solution % meta(ivar) % name))
                           ENDDO

      CALL CreateGroup_HDF5(fileId,'/state')
      CALL CreateGroup_HDF5(fileId,'/state/solution')

      CALL CreateGroup_HDF5(fileId,'/mesh')
      CALL CreateGroup_HDF5(fileId,'/mesh/coords')

      CALL WriteArray_HDF5(fileId,'/state/solution/interior', &
                           this % solution % interior)

      CALL WriteArray_HDF5(fileId,'/state/solution/boundary', &
                           this % solution % boundary)

      CALL WriteArray_HDF5(fileId,'/mesh/coords/interior', &
                           this % geometry % x % interior)

      CALL WriteArray_HDF5(fileId,'/mesh/coords/boundary', &
                           this % geometry % x % boundary)
                     


      ! Write data on plotting mesh

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

      CALL CreateGroup_HDF5(fileId,'/plot')
      CALL CreateGroup_HDF5(fileId,'/plot/mesh')
      CALL CreateGroup_HDF5(fileId,'/plot/mesh/coords')
      CALL CreateGroup_HDF5(fileId,'/plot/state')
      CALL CreateGroup_HDF5(fileId,'/plot/state/solution')
      CALL WriteArray_HDF5(fileId,'/plot/mesh/coords/interior', x % interior)
      CALL WriteArray_HDF5(fileId,'/plot/state/solution/interior',solution % interior)

      CALL Close_HDF5(fileId)

    END IF

    CALL x % Free()
    CALL solution % Free()
    CALL interp % Free()

  END SUBROUTINE Write_ECModel2D

  SUBROUTINE Read_ECModel2D(this,fileName)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:4)
    INTEGER :: firstElem
    INTEGER :: N

    IF (this % decomp % mpiEnabled) THEN
      CALL Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId, &
                     this % decomp % mpiComm)
    ELSE
      CALL Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId)
    END IF

    CALL ReadAttribute_HDF5(fileId,'N',N)

    IF (this % solution % interp % N /= N) THEN
      STOP 'Error : Solution polynomial degree does not match input file'
    END IF

    IF (this % decomp % mpiEnabled) THEN
      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId) + 1
      solOffset(1:4) = (/0,0,1,firstElem/)
      CALL ReadArray_HDF5(fileId,'/state/solution/interior', &
                          this % solution % interior,solOffset)
    ELSE
      CALL ReadArray_HDF5(fileId,'/state/solution/interior',this % solution % interior)
    END IF

    CALL Close_HDF5(fileId)

    IF (this % gpuAccel) THEN
      CALL this % solution % interior % UpdateDevice()
    END IF

  END SUBROUTINE Read_ECModel2D

  SUBROUTINE WriteTecplot_ECModel2D(this,filename)
    IMPLICIT NONE
    CLASS(ECModel2D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in),OPTIONAL :: filename
    ! Local
    CHARACTER(8) :: zoneID
    INTEGER :: fUnit
    INTEGER :: iEl,i,j,iVar
    CHARACTER(LEN=self_FileNameLength) :: tecFile
    CHARACTER(LEN=self_TecplotHeaderLength) :: tecHeader
    CHARACTER(LEN=self_FormatLength) :: fmat
    CHARACTER(13) :: timeStampString
    CHARACTER(5) :: rankString
    TYPE(Scalar2D) :: solution
    TYPE(Vector2D) :: solutionGradient
    TYPE(Vector2D) :: x
    TYPE(Lagrange),TARGET :: interp

    IF (PRESENT(filename)) THEN
      tecFile = filename
    ELSE
      timeStampString = TimeStamp(this % t,'s')

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
    CALL this % geometry % x % GridInterp(x,gpuAccel=.FALSE.)

    ! Map the solution to the target grid
    CALL this % solution % GridInterp(solution,gpuAccel=.FALSE.)

    ! Map the solution to the target grid
    CALL this % solutionGradient % GridInterp(solutionGradient,gpuAccel=.FALSE.)

    OPEN (UNIT=NEWUNIT(fUnit), &
          FILE=TRIM(tecFile), &
          FORM='formatted', &
          STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y"'
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
    WRITE (fmat,*) 3*this % solution % nvar + 2
    fmat = '('//TRIM(fmat)//'(ES16.7E3,1x))'

    DO iEl = 1,this % solution % nElem

      ! TO DO :: Get the global element ID
      WRITE (zoneID,'(I8.8)') iEl
      WRITE (fUnit,*) 'ZONE T="el'//TRIM(zoneID)//'", I=',this % solution % interp % M + 1, &
        ', J=',this % solution % interp % M + 1

      DO j = 0,this % solution % interp % M
        DO i = 0,this % solution % interp % M

          WRITE (fUnit,fmat) x % interior % hostData(1,i,j,1,iEl), &
            x % interior % hostData(2,i,j,1,iEl), &
            solution % interior % hostData(i,j,1:this % solution % nvar,iEl), &
            solutionGradient % interior % hostData(1,i,j,1:this % solution % nvar,iEl), &
            solutionGradient % interior % hostData(2,i,j,1:this % solution % nvar,iEl)

        END DO
      END DO

    END DO

    CLOSE (UNIT=fUnit)

    CALL x % Free()
    CALL solution % Free()
    CALL solutionGradient % Free()
    CALL interp % Free()

  END SUBROUTINE WriteTecplot_ECModel2D

END MODULE SELF_ECModel2D
