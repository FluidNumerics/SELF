!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_DGModel1D

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

  TYPE,EXTENDS(Model) :: DGModel1D
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

    PROCEDURE :: Init => Init_DGModel1D
    PROCEDURE :: Free => Free_DGModel1D

    PROCEDURE :: UpdateSolution => UpdateSolution_DGModel1D

    PROCEDURE :: ResizePrevSol => ResizePrevSol_DGModel1D

    PROCEDURE :: UpdateGAB2 => UpdateGAB2_DGModel1D
    PROCEDURE :: UpdateGAB3 => UpdateGAB3_DGModel1D
    PROCEDURE :: UpdateGAB4 => UpdateGAB4_DGModel1D

    PROCEDURE :: UpdateGRK2 => UpdateGRK2_DGModel1D
    PROCEDURE :: UpdateGRK3 => UpdateGRK3_DGModel1D
    PROCEDURE :: UpdateGRK4 => UpdateGRK4_DGModel1D
    PROCEDURE :: CalculateTendency => CalculateTendency_DGModel1D
    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_DGModel1D

    GENERIC :: SetSolution => SetSolutionFromChar_DGModel1D, &
      SetSolutionFromEqn_DGModel1D
    PROCEDURE,PRIVATE :: SetSolutionFromChar_DGModel1D
    PROCEDURE,PRIVATE :: SetSolutionFromEqn_DGModel1D

    PROCEDURE :: ReadModel => Read_DGModel1D
    PROCEDURE :: WriteModel => Write_DGModel1D
    PROCEDURE :: WriteTecplot => WriteTecplot_DGModel1D

  END TYPE DGModel1D

CONTAINS

  SUBROUTINE Init_DGModel1D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh1D),INTENT(in),TARGET :: mesh
    TYPE(Geometry1D),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry

    CALL this % solution % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % velocity % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nVar,this % mesh % nElem)

  END SUBROUTINE Init_DGModel1D

  SUBROUTINE Free_DGModel1D(this)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this

    CALL this % solution % Free()
    CALL this % workSol % Free()
    CALL this % prevSol % Free()
    CALL this % velocity % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_DGModel1D

  SUBROUTINE ResizePrevSol_DGModel1D(this,m)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: nVar

    ! Free space, if necessary
    CALL this % prevSol % Free()

    ! Reallocate with increased variable dimension for
    ! storing "m" copies of solution data
    nVar = this % solution % nVar
    CALL this % prevSol % Init(this % geometry % x % interp,m*nVar,this % mesh % nElem)

  END SUBROUTINE ResizePrevSol_DGModel1D

  SUBROUTINE SetSolutionFromEqn_DGModel1D(this,eqn)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

    ! Copy the equation parser
    DO iVar = 1,this % solution % nVar
      CALL this % solution % SetEquation(ivar,eqn(iVar) % equation)
    END DO

    CALL this % solution % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % solution % BoundaryInterp()

  END SUBROUTINE SetSolutionFromEqn_DGModel1D

  SUBROUTINE SetSolutionFromChar_DGModel1D(this,eqnChar)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

    DO iVar = 1,this % solution % nVar
      PRINT *, iVar,eqnChar(iVar)
      CALL this % solution % SetEquation(ivar,eqnChar(iVar))
    END DO

    CALL this % solution % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % solution % BoundaryInterp()

  END SUBROUTINE SetSolutionFromChar_DGModel1D

  SUBROUTINE UpdateSolution_DGModel1D(this,dt)
    !! Computes a solution update as , where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    REAL(prec),OPTIONAL,INTENT(in) :: dt
    ! Local
    REAL(prec) :: dtLoc
    INTEGER :: i,iEl,iVar

    IF (PRESENT(dt)) THEN
      dtLoc = dt
    ELSE
      dtLoc = this % dt
    END IF

    !$omp target map(to:this % dsdt % interior) map(tofrom:this % solution)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    DO iEl = 1,this % solution % nElem
      DO iVar = 1,this % solution % nVar
        DO i = 1,this % solution % interp % N+1

          this % solution % interior(i,iEl,iVar) = &
            this % solution % interior(i,iEl,iVar) + &
            dtLoc*this % dSdt % interior(i,iEl,iVar)

        END DO
      END DO
    END DO
    !$omp end target

  END SUBROUTINE UpdateSolution_DGModel1D

  SUBROUTINE UpdateGAB2_DGModel1D(this,m)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,nVar,iEl,iVar


      ! ab2_weight
      IF (m == 0) THEN ! Initialization step - store the solution in the prevSol

        !$omp target map(tofrom: this % solution % interior) map(from:this % prevSol % interior)
        !$omp teams distribute parallel do collapse(3) num_threads(256)
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 1,this % solution % interp % N+1

              this % prevSol % interior(i,iEl,iVar) = this % solution % interior(i,iEl,iVar)

            END DO
          END DO
        END DO
        !$omp end target

      ELSEIF (m == 1) THEN ! Copy the solution back from prevsol

        !$omp target map(from: this % solution % interior) map(to:this % prevSol % interior)
        !$omp teams distribute parallel do collapse(3) num_threads(256)
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 1,this % solution % interp % N+1

              this % solution % interior(i,iEl,iVar) = this % prevSol % interior(i,iEl,iVar)

            END DO
          END DO
        END DO
        !$omp end target

      ELSE ! Main looping section - nVar the previous solution, store the new solution, and
        ! create an interpolated solution to use for tendency calculation

        nVar = this % solution % nVar
        !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
        !$omp teams distribute parallel do collapse(3) num_threads(256)
        DO iEl = 1,this % solution % nElem
          DO iVar = 1,this % solution % nVar
            DO i = 1,this % solution % interp % N+1

              ! Bump the last solution
              this % prevSol % interior(i,iEl, nVar + iVar) = this % prevSol % interior(i,iEl,iVar)

              ! Store the new solution
              this % prevSol % interior(i,iEl,iVar) = this % solution % interior(i,iEl,iVar)

              this % solution % interior(i,iEl,iVar) = &
                1.5_PREC*this % prevSol % interior(i,iEl,iVar) - &
                0.5_PREC*this % prevSol % interior(i,iEl, nVar + iVar)
            END DO
          END DO
        END DO
        !$omp end target
      END IF

  END SUBROUTINE UpdateGAB2_DGModel1D

  SUBROUTINE UpdateGAB3_DGModel1D(this,m)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,nVar,iEl,iVar

    IF (m == 0) THEN ! Initialization step - store the solution in the prevSol at nvar+ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      nVar = this % solution % nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

            this % prevSol % interior(i,iEl, nVar + iVar) = this % solution % interior(i,iEl,iVar)

          END DO
        END DO
      END DO

    ELSEIF (m == 1) THEN ! Initialization step - store the solution in the prevSol at ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

            this % prevSol % interior(i,iEl,iVar) = this % solution % interior(i,iEl,iVar)

          END DO
        END DO
      END DO
      !$omp end target

    ELSEIF (m == 2) THEN ! Copy the solution back from the most recent prevsol

      !$omp target map(from: this % solution % interior) map(to: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

            this % solution % interior(i,iEl,iVar) = this % prevSol % interior(i,iEl,iVar)

          END DO
        END DO
      END DO
      !$omp end target

    ELSE ! Main looping section - nVar the previous solution, store the new solution, and
      ! create an interpolated solution to use for tendency calculation

      !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
      nVar = this % solution % nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

            ! Bump the last two stored solutions
            nVar = this % solution % nVar
    this % prevSol % interior(i,iEl, 2*nVar + iVar) = this % prevSol % interior(i,iEl, nVar + iVar)
            this % prevSol % interior(i,iEl, nVar + iVar) = this % prevSol % interior(i,iEl,iVar)

            ! Store the new solution
            this % prevSol % interior(i,iEl,iVar) = this % solution % interior(i,iEl,iVar)

            this % solution % interior(i,iEl,iVar) = &
              (23.0_PREC*this % prevSol % interior(i,iEl,iVar) - &
                16.0_PREC*this % prevSol % interior(i,iEl, nVar + iVar) + &
                5.0_PREC*this % prevSol % interior(i,iEl, 2*nVar + iVar))/12.0_PREC

          END DO
        END DO
      END DO
      !$omp end target

    END IF

  END SUBROUTINE UpdateGAB3_DGModel1D

  SUBROUTINE UpdateGAB4_DGModel1D(this,m)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,nVar,iEl,iVar


    IF (m == 0) THEN ! Initialization step - store the solution in the prevSol at nvar+ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      nVar = this % solution % nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

          this % prevSol % interior(i,iEl, 2*nVar + iVar) = this % solution % interior(i,iEl,iVar)

          END DO
        END DO
      END DO
      !$omp end target

    ELSEIF (m == 1) THEN ! Initialization step - store the solution in the prevSol at ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      nVar = this % solution % nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

            this % prevSol % interior(i,iEl, nVar + iVar) = this % solution % interior(i,iEl,iVar)

          END DO
        END DO
      END DO
      !$omp end target

    ELSEIF (m == 2) THEN ! Initialization step - store the solution in the prevSol at ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

            this % prevSol % interior(i,iEl,iVar) = this % solution % interior(i,iEl,iVar)

          END DO
        END DO
      END DO
      !$omp end target

    ELSEIF (m == 3) THEN ! Copy the solution back from the most recent prevsol

      !$omp target map(from: this % solution % interior) map(to: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

            this % solution % interior(i,iEl,iVar) = this % prevSol % interior(i,iEl,iVar)

          END DO
        END DO
      END DO

    ELSE ! Main looping section - nVar the previous solution, store the new solution, and
      ! create an interpolated solution to use for tendency calculation

      !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
      nVar = this % solution % nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO i = 1,this % solution % interp % N+1

            ! Bump the last two stored solutions
  this % prevSol % interior(i,iEl, 3*nVar + iVar) = this % prevSol % interior(i,iEl, 2*nVar + iVar)
    this % prevSol % interior(i,iEl, 2*nVar + iVar) = this % prevSol % interior(i,iEl, nVar + iVar)
            this % prevSol % interior(i,iEl, nVar + iVar) = this % prevSol % interior(i,iEl,iVar)

            ! Store the new solution
            this % prevSol % interior(i,iEl,iVar) = this % solution % interior(i,iEl,iVar)

            this % solution % interior(i,iEl,iVar) = &
              (55.0_PREC*this % prevSol % interior(i,iEl,iVar) - &
                59.0_PREC*this % prevSol % interior(i,iEl, nVar + iVar) + &
                37.0_PREC*this % prevSol % interior(i,iEl, 2*nVar + iVar) - &
                9.0_PREC*this % prevSol % interior(i,iEl, 3*nVar + iVar))/24.0_PREC

          END DO
        END DO
      END DO
      !$omp end target

    END IF

  END SUBROUTINE UpdateGAB4_DGModel1D

  SUBROUTINE UpdateGRK2_DGModel1D(this,m)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,iEl,iVar

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    DO iEl = 1,this % solution % nElem
      DO iVar = 1,this % solution % nVar
        DO i = 1,this % solution % interp % N+1

          this % workSol % interior(i,iEl,iVar) = rk2_a(m)* &
                                                              this % workSol % interior(i,iEl,iVar) + &
                                                              this % dSdt % interior(i,iEl,iVar)

          this % solution % interior(i,iEl,iVar) = &
            this % solution % interior(i,iEl,iVar) + &
            rk2_g(m)*this % dt*this % workSol % interior(i,iEl,iVar)

        END DO
      END DO
    END DO
    !$omp end target

  END SUBROUTINE UpdateGRK2_DGModel1D

  SUBROUTINE UpdateGRK3_DGModel1D(this,m)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,iEl,iVar

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    DO iEl = 1,this % solution % nElem
      DO iVar = 1,this % solution % nVar
        DO i = 1,this % solution % interp % N+1

          this % workSol % interior(i,iEl,iVar) = rk3_a(m)* &
                                                              this % workSol % interior(i,iEl,iVar) + &
                                                              this % dSdt % interior(i,iEl,iVar)

          this % solution % interior(i,iEl,iVar) = &
            this % solution % interior(i,iEl,iVar) + &
            rk3_g(m)*this % dt*this % workSol % interior(i,iEl,iVar)

        END DO
      END DO
    END DO
    !$omp end target

  END SUBROUTINE UpdateGRK3_DGModel1D

  SUBROUTINE UpdateGRK4_DGModel1D(this,m)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    INTEGER,INTENT(in) :: m
    ! Local
    INTEGER :: i,iEl,iVar

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    DO iEl = 1,this % solution % nElem
      DO iVar = 1,this % solution % nVar
        DO i = 1,this % solution % interp % N+1

          this % workSol % interior(i,iEl,iVar) = rk4_a(m)* &
                                                              this % workSol % interior(i,iEl,iVar) + &
                                                              this % dSdt % interior(i,iEl,iVar)

          this % solution % interior(i,iEl,iVar) = &
            this % solution % interior(i,iEl,iVar) + &
            rk4_g(m)*this % dt*this % workSol % interior(i,iEl,iVar)

        END DO
      END DO
    END DO
    !$omp end target

  END SUBROUTINE UpdateGRK4_DGModel1D

  SUBROUTINE CalculateFluxDivergence_DGModel1D(this)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this

    CALL this % flux % DGDerivative(this % geometry,this % fluxDivergence)

  END SUBROUTINE CalculateFluxDivergence_DGModel1D

  SUBROUTINE CalculateTendency_DGModel1D(this)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iEl,iVar

    CALL this % solution % BoundaryInterp()
    CALL this % solution % SideExchange(this % mesh,this % decomp)
    CALL this % PreTendency()
    CALL this % SetBoundaryCondition()
    CALL this % SourceMethod()
    CALL this % RiemannSolver()
    CALL this % FluxMethod()
    CALL this % CalculateFluxDivergence()

    !$omp target map(to: this % source, this % fluxDivergence) map(from:this % dSdt)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    DO iEl = 1,this % solution % nElem
      DO iVar = 1,this % solution % nVar
        DO i = 1,this % solution % interp % N+1

          this % dSdt % interior(i,iEl,iVar) = &
            this % source % interior(i,iEl,iVar) - &
            this % fluxDivergence % interior(i,iEl,iVar)

        END DO
      END DO
    END DO
    !$omp end target

  END SUBROUTINE CalculateTendency_DGModel1D

  SUBROUTINE Write_DGModel1D(this,fileName)
#undef __FUNC__
#define __FUNC__ "Write_DGModel1D"
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
    CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    TYPE(Scalar1D) :: solution
    TYPE(Scalar1D) :: x
    TYPE(Lagrange),TARGET :: interp
    CHARACTER(LEN=self_FileNameLength) :: pickupFile
    CHARACTER(13) :: timeStampString

    WRITE (timeStampString,'(I13.13)') this % ioIterate
    IF (PRESENT(filename)) THEN
      pickupFile = TRIM(filename)//timeStampString//'.h5'
    ELSE
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

  END SUBROUTINE Write_DGModel1D

  SUBROUTINE Read_DGModel1D(this,fileName)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
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
      firstElem = this % decomp % offsetElem(this % decomp % rankId) + 1
      solOffset(1:3) = (/0,1,firstElem/)
      CALL ReadArray_HDF5(fileId,'/controlgrid/solution/interior', &
                          this % solution % interior,solOffset)
    ELSE
      CALL ReadArray_HDF5(fileId,'/controlgrid/solution/interior',this % solution % interior)
    END IF

    CALL Close_HDF5(fileId)

  END SUBROUTINE Read_DGModel1D

  SUBROUTINE WriteTecplot_DGModel1D(this,filename)
    IMPLICIT NONE
    CLASS(DGModel1D),INTENT(inout) :: this
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
      tecFile = 'solution.'//timeStampString//'.curve'

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
    CALL this % geometry % x % GridInterp(x)

    ! Map the solution to the target grid
    CALL this % solution % GridInterp(solution)

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

          WRITE (fUnit,fmat) x % interior(i,1,iEl), &
            solution % interior(i,iEl,iVar)

        END DO

      END DO
    END DO

    CLOSE (UNIT=fUnit)

    CALL x % Free()
    CALL solution % Free()
    CALL interp % Free()

  END SUBROUTINE WriteTecplot_DGModel1D

END MODULE SELF_DGModel1D
