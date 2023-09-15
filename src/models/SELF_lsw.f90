!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_lsw

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model2D

  IMPLICIT NONE
#include "../SELF_Macros.h"

  TYPE,EXTENDS(Model2D) :: lsw
    !! iVar = 1 ~> u velocity component
    !! iVar = 2 ~> v velocity component
    !! iVar = 3 ~> free surface height
    TYPE(MappedScalar2D) :: fCori ! Coriolis parameter ( 1/s )
    TYPE(MappedScalar2D) :: H ! bottom topography ( m )
    REAL(prec) :: g     ! gravity ( m/s^2)

  CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_lsw
    PROCEDURE :: Free => Free_lsw
    PROCEDURE :: PrintType => PrintType_lsw
    PROCEDURE :: CalculateEntropy => CalculateEntropy_lsw

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_lsw
    PROCEDURE :: FluxMethod => Flux_lsw
    PROCEDURE :: RiemannSolver => RiemannSolver_lsw
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_lsw

    ! New Methods
    PROCEDURE :: SetInitialConditions => SetInitialConditions_lsw

    GENERIC :: SetCoriolis => SetCoriolisFromChar_lsw, &
      SetCoriolisFromEqn_lsw
    PROCEDURE,PRIVATE :: SetCoriolisFromChar_lsw
    PROCEDURE,PRIVATE :: SetCoriolisFromEqn_lsw

    GENERIC :: SetBathymetry => SetBathymetryFromChar_lsw, &
      SetBathymetryFromEqn_lsw, &
      SetBathymetryFromConstant_lsw
    PROCEDURE,PRIVATE :: SetBathymetryFromChar_lsw
    PROCEDURE,PRIVATE :: SetBathymetryFromEqn_lsw
    PROCEDURE,PRIVATE :: SetBathymetryFromConstant_lsw

    PROCEDURE :: DiagnoseGeostrophicVelocity => DiagnoseGeostrophicVelocity_lsw

  END TYPE lsw

  INTERFACE
    SUBROUTINE SetBoundaryCondition_lsw_gpu_wrapper(solution,extBoundary,nHat,sideInfo,N,nVar,nEl) &
      BIND(c,name="SetBoundaryCondition_lsw_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: solution,extBoundary,nHat,sideInfo
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SetBoundaryCondition_lsw_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Source_lsw_gpu_wrapper(source,solution,f,N,nVar,nEl) &
      BIND(c,name="Source_lsw_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: source,solution,f
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE Source_lsw_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Flux_lsw_gpu_wrapper(flux,solution,H,g,N,nVar,nEl) &
      BIND(c,name="Flux_lsw_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: flux,solution,H
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: g
    END SUBROUTINE Flux_lsw_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE RiemannSolver_lsw_gpu_wrapper(flux,solution,extBoundary,H,nHat,nScale,g,N,nVar,nEl) &
      BIND(c,name="RiemannSolver_lsw_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: flux,solution,extBoundary,H,nHat,nScale
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: g
    END SUBROUTINE RiemannSolver_lsw_gpu_wrapper
  END INTERFACE

CONTAINS

  SUBROUTINE Init_lsw(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(lsw),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh2D),INTENT(in),TARGET :: mesh
    TYPE(SEMQuad),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp
    ! Local
    INTEGER :: ivar
    CHARACTER(LEN=3) :: ivarChar
    CHARACTER(LEN=25) :: varname
    INTEGER :: nvarloc

    ! Ensure that the number of variables is 3
    ! nvar is unused in this class extension
    nvarloc = 3

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.
    this % g = 1.0_PREC

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % fCori % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % H % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    ! First three variables are treated as u, v, eta
    ! Any additional are treated as passive tracers
    CALL this % solution % SetName(1,"u")
    CALL this % solution % SetUnits(1,"m/s")
    CALL this % solution % SetDescription(1,"x-component of the barotropic velocity field")

    CALL this % solution % SetName(2,"v")
    CALL this % solution % SetUnits(2,"m/s")
    CALL this % solution % SetDescription(2,"y-component of the barotropic velocity field")

    CALL this % solution % SetName(3,"eta")
    CALL this % solution % SetUnits(3,"m")
    CALL this % solution % SetDescription(3,"Free surface height anomaly")

  END SUBROUTINE Init_lsw

  SUBROUTINE PrintType_lsw(this)
#undef __FUNC__
#define __FUNC__ "Init"
    IMPLICIT NONE
    CLASS(lsw),INTENT(in) :: this

    INFO("Model : Linear Shallow Water")

  END SUBROUTINE PrintType_lsw

  SUBROUTINE Free_lsw(this)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this

    CALL this % solution % Free()
    CALL this % fCori % Free()
    CALL this % H % Free()
    CALL this % workSol % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_lsw

  SUBROUTINE SetInitialConditions_lsw(this, config)
#undef __FUNC__
#define __FUNC__ "SetInitialConditions"
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    TYPE(SELFConfig), INTENT(inout) :: config
    ! Local
    LOGICAL :: setStaticState
    LOGICAL :: hydrostaticAdjust
    LOGICAL :: found
    CHARACTER(LEN=self_EquationLength) :: u, v, eta
    CHARACTER(LEN=self_EquationLength) :: fcori ! coriolis parameter
    CHARACTER(LEN=self_EquationLength) :: bathymetry ! bathymetry (bottom topography)
    CHARACTER(LEN=self_EquationLength) :: initialCondition(1:3)
    CHARACTER(LEN=self_QuadratureTypeCharLength) :: integrator
    CHARACTER(LEN=SELF_EQUATION_LENGTH) :: velocity(1:2)
    INTEGER,PARAMETER :: ucs2 = SELECTED_CHAR_KIND('ISO_10646')
    CHARACTER(KIND=ucs2,len=20) :: tUSC
    CHARACTER(KIND=ucs2,len=20) :: rhoUSC
    CHARACTER(KIND=ucs2,len=20) :: CpUSC
    CHARACTER(KIND=ucs2,len=20) :: CvUSC
    CHARACTER(KIND=ucs2,len=:),ALLOCATABLE :: str
    CHARACTER(4) :: arrayCount
    REAL(prec) :: momentumMax
    INTEGER :: i, nfeatures
    REAL(prec) :: featureParams(1:10)



    ! Get environmental parameters
    CALL config % Get("lsw.environment.g",this % g)

    CALL config % Get("lsw.environment.coriolis",fcori)
    CALL this % SetCoriolis( fcori )

    CALL config % Get("lsw.environment.bathymetry",bathymetry)
    CALL this % SetBathymetry( bathymetry )   

    ! Get initial conditions 
    CALL config % Get("lsw.initial_conditions.u",u)
    CALL config % Get("lsw.initial_conditions.v",v)
    CALL config % Get("lsw.initial_conditions.eta",eta)

    ! Set the initial condition 
    CALL this % SetSolution( (/u,v,eta/) )

    ! Set the time integrator
    CALL config % Get("time_options.integrator",integrator)
    CALL this % SetTimeIntegrator(TRIM(integrator)) ! Set the integrator
    CALL config % Get("time_options.dt",this % dt) ! Set the time step size
    CALL config % Get("time_options.start_time",this % t) ! Set the initial time

  
    CALL this % CalculateEntropy()
    CALL this % ReportEntropy()


  END SUBROUTINE SetInitialConditions_lsw

  SUBROUTINE CalculateEntropy_lsw(this)
  !! Base method for calculating entropy of a model
  !! Calculates the entropy as the integration of the
  !! squared tracer over the domain
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iVar,iEl
    REAL(prec) :: Jacobian,u,v,eta,H
    REAL(prec) :: wi,wj
    REAL(prec) :: entropy

    IF (this % gpuAccel) THEN
      CALL this % solution % interior % UpdateHost()
    END IF

    entropy = 0.0_PREC

    DO iEl = 1,this % geometry % x % nElem
      DO j = 0,this % geometry % x % interp % N
        DO i = 0,this % geometry % x % interp % N

          ! Coordinate mapping Jacobian
          Jacobian = this % geometry % J % interior % hostData(i,j,1,iEl)

          ! Quadrature weights
          wi = this % geometry % x % interp % qWeights % hostData(i)
          wj = this % geometry % x % interp % qWeights % hostData(j)

          ! Solution
          u = this % solution % interior % hostData(i,j,1,iEl)
          v = this % solution % interior % hostData(i,j,2,iEl)
          eta = this % solution % interior % hostData(i,j,3,iEl)

          H = this % H % interior % hostData(i,j,1,iEl)

          entropy = entropy + &
                    0.5_PREC*(H*(u*u + v*v) + this % g*eta*eta)*Jacobian*wi*wj

        END DO
      END DO
    END DO

    CALL this % decomp % GlobalReduce(entropy,this % entropy)

  END SUBROUTINE CalculateEntropy_lsw

  SUBROUTINE SetCoriolisFromEqn_lsw(this,eqn)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn

    ! Copy the equation parser
    CALL this % fCori % SetEquation(1,eqn % equation)

    CALL this % fCori % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % fCori % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % fCori % UpdateDevice()
    END IF

  END SUBROUTINE SetCoriolisFromEqn_lsw

  SUBROUTINE SetCoriolisFromChar_lsw(this,eqnChar)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: eqnChar

    CALL this % fCori % SetEquation(1,TRIM(eqnChar))

    CALL this % fCori % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % fCori % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % fCori % UpdateDevice()
    END IF

  END SUBROUTINE SetCoriolisFromChar_lsw

  SUBROUTINE SetBathymetryFromEqn_lsw(this,eqn)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn

    ! Copy the equation parser
    CALL this % H % SetEquation(1,eqn % equation)

    CALL this % H % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % H % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % H % UpdateDevice()
    END IF

  END SUBROUTINE SetBathymetryFromEqn_lsw

  SUBROUTINE SetBathymetryFromChar_lsw(this,eqnChar)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: eqnChar

    CALL this % H % SetEquation(1,TRIM(eqnChar))

    CALL this % H % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % H % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % H % UpdateDevice()
    END IF

  END SUBROUTINE SetBathymetryFromChar_lsw

  SUBROUTINE SetBathymetryFromConstant_lsw(this,H)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: H

    this % H % interior % hostData = H
    CALL this % H % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % H % UpdateDevice()
    END IF

  END SUBROUTINE SetBathymetryFromConstant_lsw

  SUBROUTINE SetBoundaryCondition_lsw(this)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,i
    INTEGER :: bcid,e2
    REAL(prec) :: u,v,nhat(1:2)

    IF (this % gpuAccel) THEN

      CALL SetBoundaryCondition_lsw_gpu_wrapper(this % solution % boundary % deviceData, &
                                                               this % solution % extBoundary % deviceData, &
                                                               this % geometry % nHat % boundary % deviceData, &
                                                               this % mesh % sideInfo % deviceData, &
                                                               this % solution % interp % N, &
                                                               this % solution % nVar, &
                                                               this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iSide = 1,4
          DO i = 0,this % solution % interp % N

            bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
            e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
            IF (e2 == 0) THEN
              IF (bcid == SELF_BC_RADIATION) THEN

                this % solution % extBoundary % hostData(i,1,iSide,iEl) = 0.0_PREC
                this % solution % extBoundary % hostData(i,2,iSide,iEl) = 0.0_PREC
                this % solution % extBoundary % hostData(i,3,iSide,iEl) = 0.0_PREC

              ELSEIF (bcid == SELF_BC_NONORMALFLOW) THEN

                nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
                u = this % solution % boundary % hostData(i,1,iSide,iEl)
                v = this % solution % boundary % hostData(i,2,iSide,iEl)
                this % solution % extBoundary % hostData(i,1,iSide,iEl) = (nhat(2)**2 - nhat(1)**2)*u - 2.0_PREC*nhat(1)*nhat(2)*v
                this % solution % extBoundary % hostData(i,2,iSide,iEl) = (nhat(1)**2 - nhat(2)**2)*v - 2.0_PREC*nhat(1)*nhat(2)*u
                this % solution % extBoundary % hostData(i,3,iSide,iEl) = this % solution % boundary % hostData(i,3,iSide,iEl)

              ELSE ! Default boundary condition is radiation

                this % solution % extBoundary % hostData(i,1,iSide,iEl) = 0.0_PREC
                this % solution % extBoundary % hostData(i,2,iSide,iEl) = 0.0_PREC
                this % solution % extBoundary % hostData(i,3,iSide,iEl) = 0.0_PREC

              END IF
            END IF

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE SetBoundaryCondition_lsw

  SUBROUTINE DiagnoseGeostrophicVelocity_lsw(this)
  !! Sets the velocity components (solution 1-2) to 0 and then diagnoses the
  !! the velocity field using a balance of the pressure gradient force
  !! and the coriolis force
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    DO iEl = 1,this % source % nElem
      DO j = 0,this % source % interp % N
        DO i = 0,this % source % interp % N

          ! u velocity component
          this % solution % interior % hostData(i,j,1,iEl) = 0.0_PREC

          ! v velocity component
          this % solution % interior % hostData(i,j,2,iEl) = 0.0_PREC

        END DO
      END DO
    END DO

    IF (this % gpuAccel) THEN
      CALL this % solution % interior % UpdateDevice()
    END IF

    ! Calculate tendency
    CALL this % CalculateTendency()

    IF (this % gpuAccel) THEN
      CALL this % dSdt % interior % UpdateHost()
    END IF

    DO iEl = 1,this % source % nElem
      DO j = 0,this % source % interp % N
        DO i = 0,this % source % interp % N

          ! u velocity component = ( - g ta_y ) / f
          this % solution % interior % hostData(i,j,1,iEl) = this % dSdt % interior % hostData(i,j,2,iEl)/ &
                                                             this % fCori % interior % hostData(i,j,1,iEl)

          ! v velocity component = - ( - g ta_x ) / f
          this % solution % interior % hostData(i,j,2,iEl) = -this % dSdt % interior % hostData(i,j,1,iEl)/ &
                                                             this % fCori % interior % hostData(i,j,1,iEl)

        END DO
      END DO
    END DO

    IF (this % gpuAccel) THEN
      CALL this % solution % interior % UpdateDevice()
    END IF

  END SUBROUTINE DiagnoseGeostrophicVelocity_lsw

  SUBROUTINE Source_lsw(this)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    IF (this % gpuAccel) THEN

      CALL Source_lsw_gpu_wrapper(this % source % interior % deviceData, &
                                                 this % solution % interior % deviceData, &
                                                 this % fCori % interior % deviceData, &
                                                 this % source % interp % N, &
                                                 this % solution % nVar, &
                                                 this % solution % nElem)

    ELSE

      DO iEl = 1,this % source % nElem
        DO j = 0,this % source % interp % N
          DO i = 0,this % source % interp % N

            this % source % interior % hostData(i,j,1,iEl) = this % fCori % interior % hostData(i,j,1,iEl)* &
                                                             this % solution % interior % hostData(i,j,2,iEl)
            this % source % interior % hostData(i,j,2,iEl) = -this % fCori % interior % hostData(i,j,1,iEl)* &
                                                             this % solution % interior % hostData(i,j,1,iEl)
            this % source % interior % hostData(i,j,3,iEl) = 0.0_PREC

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE Source_lsw

  SUBROUTINE Flux_lsw(this)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    IF (this % gpuAccel) THEN

      CALL Flux_lsw_gpu_wrapper(this % flux % interior % deviceData, &
                                               this % solution % interior % deviceData, &
                                               this % H % interior % deviceData, &
                                               this % g,this % solution % interp % N, &
                                               this % solution % nVar,this % solution % nElem)

    ELSE
      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO j = 0,this % solution % interp % N
            DO i = 0,this % solution % interp % N

              IF (iVar == 1) THEN ! u-velocity
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                  this % g*this % solution % interior % hostData(i,j,3,iEl)

                this % flux % interior % hostData(2,i,j,iVar,iEl) = 0.0_PREC

              ELSEIF (iVar == 2) THEN ! v-velocity

                this % flux % interior % hostData(1,i,j,iVar,iEl) = 0.0_PREC

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                  this % g*this % solution % interior % hostData(i,j,3,iEl)

              ELSEIF (iVar == 3) THEN ! free surface height
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                  this % H % interior % hostData(i,j,1,iEl)* &
                  this % solution % interior % hostData(i,j,1,iEl)

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                  this % H % interior % hostData(i,j,1,iEl)* &
                  this % solution % interior % hostData(i,j,2,iEl)
              END IF

            END DO
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE Flux_lsw

  SUBROUTINE RiemannSolver_lsw(this)
    IMPLICIT NONE
    CLASS(lsw),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: nhat(1:2),nmag
    REAL(prec) :: c,unL,unR,etaL,etaR,wL,wR

    IF (this % gpuAccel) THEN

      CALL RiemannSolver_lsw_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
                                                        this % solution % boundary % deviceData, &
                                                        this % solution % extBoundary % deviceData, &
                                                        this % H % boundary % deviceData, &
                                                        this % geometry % nHat % boundary % deviceData, &
                                                        this % geometry % nScale % boundary % deviceData, &
                                                        this % g, &
                                                        this % solution % interp % N, &
                                                        this % solution % nVar, &
                                                        this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iSide = 1,4
          DO i = 0,this % solution % interp % N

            ! Get the boundary normals on cell edges from the mesh geometry
            nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
            nmag = this % geometry % nScale % boundary % hostData(i,1,iSide,iEl)
            c = SQRT(this % g*this % H % boundary % hostData(i,1,iSide,iEl))

            ! Calculate the normal velocity at the cell edges
            unL = this % solution % boundary % hostData(i,1,iSide,iEl)*nHat(1) + &
                  this % solution % boundary % hostData(i,2,iSide,iEl)*nHat(2)

            unR = this % solution % extBoundary % hostData(i,1,iSide,iEl)*nHat(1) + &
                  this % solution % extBoundary % hostData(i,2,iSide,iEl)*nHat(2)

            etaL = this % solution % boundary % hostData(i,3,iSide,iEl)
            etaR = this % solution % extBoundary % hostData(i,3,iSide,iEl)

            ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
            wL = 0.5_PREC*(unL/this % g + etaL/c)
            wR = 0.5_PREC*(unR/this % g - etaR/c)

            this % flux % boundaryNormal % hostData(i,1,iSide,iEl) = this % g*c*(wL - wR)*nHat(1)*nmag
            this % flux % boundaryNormal % hostData(i,2,iSide,iEl) = this % g*c*(wL - wR)*nHat(2)*nmag
            this % flux % boundaryNormal % hostData(i,3,iSide,iEl) = c*c*(wL + wR)*nmag

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE RiemannSolver_lsw

END MODULE SELF_lsw
