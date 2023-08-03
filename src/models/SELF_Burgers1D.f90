!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Burgers1D

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model1D
  USE SELF_Config

  IMPLICIT NONE
#include "../SELF_Macros.h"

  TYPE,EXTENDS(Model1D) :: Burgers1D

    REAL(prec) :: viscosity = 0.0_PREC
    TYPE(EquationParser) :: uLEqn
    TYPE(EquationParser) :: uREqn
    REAL(prec) :: uL = 0.0_PREC
    REAL(prec) :: uR = 0.0_PREC

    TYPE(EquationParser) :: uxLEqn
    TYPE(EquationParser) :: uxREqn
    REAL(prec) :: uxL = 0.0_PREC
    REAL(prec) :: uxR = 0.0_PREC

  CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_Burgers1D
    PROCEDURE :: PrintType => PrintType_Burgers1D
    PROCEDURE :: CalculateEntropy => CalculateEntropy_Burgers1D
    PROCEDURE :: SetInitialConditions => SetInitialConditions_Burgers1D

    ! Concretized Methods

    PROCEDURE :: PreTendency => PreTendency_Burgers1D
    PROCEDURE :: SourceMethod => Source_Burgers1D
    PROCEDURE :: FluxMethod => Flux_Burgers1D
    PROCEDURE :: RiemannSolver => RiemannSolver_Burgers1D
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_Burgers1D

  END TYPE Burgers1D

!  INTERFACE
!    SUBROUTINE SetBoundaryCondition_Burgers1D_gpu_wrapper(solution, extBoundary, nHat, sideInfo, N, nVar, nEl) &
!      bind(c,name="SetBoundaryCondition_Burgers1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: solution, extBoundary, nHat, sideInfo
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!    END SUBROUTINE SetBoundaryCondition_Burgers1D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE Source_Burgers1D_gpu_wrapper(source, solution, gradH, g, N, nVar, nEl) &
!      bind(c,name="Source_Burgers1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: source, solution, gradH
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!      REAL(c_prec),VALUE :: g
!    END SUBROUTINE Source_Burgers1D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE Flux_Burgers1D_gpu_wrapper(flux, solution, g, H, N, nVar, nEl) &
!      bind(c,name="Flux_Burgers1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: flux, solution
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!      REAL(c_prec),VALUE :: g, H
!    END SUBROUTINE Flux_Burgers1D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE RiemannSolver_Burgers1D_gpu_wrapper(flux, solution, extBoundary, nHat, nScale, g, H, N, nVar, nEl) &
!      bind(c,name="RiemannSolver_Burgers1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: flux, solution, extBoundary, nHat, nScale
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!      REAL(c_prec),VALUE :: g, H
!    END SUBROUTINE RiemannSolver_Burgers1D_gpu_wrapper
!  END INTERFACE

CONTAINS

  SUBROUTINE Init_Burgers1D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh1D),INTENT(in),TARGET :: mesh
    TYPE(Geometry1D),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp
    ! Local
    INTEGER :: ivar
    CHARACTER(LEN=3) :: ivarChar
    CHARACTER(LEN=25) :: varname
    INTEGER :: nvarloc

    ! Ensure that the number of variables is 3
    ! nvar is unused in this class extension
    nvarloc = 1

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % velocity % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    CALL this % solution % SetName(1,"u")
    CALL this % solution % SetUnits(1,"[]")
    CALL this % solution % SetDescription(1,"null")

  END SUBROUTINE Init_Burgers1D

  SUBROUTINE SetInitialConditions_Burgers1D(this,config)
#undef __FUNC__
#define __FUNC__ "SetInitialConditions"
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    TYPE(SELFConfig),INTENT(inout) :: config
    ! Local
    CHARACTER(LEN=self_EquationLength) ::uLEqn
    CHARACTER(LEN=self_EquationLength) ::uREqn
    CHARACTER(LEN=self_QuadratureTypeCharLength) :: integrator
    CHARACTER(LEN=self_EquationLength) :: u

    INFO("Setting initial conditions")
    ! Set the time integrator
    CALL config % Get("time_options.integrator",integrator)
    CALL this % SetTimeIntegrator(TRIM(integrator)) ! Set the integrator
    CALL config % Get("time_options.dt",this % dt) ! Set the time step size
    CALL config % Get("time_options.start_time",this % t) ! Set the initial time

    ! Get static parameters
    CALL config % Get("brg1d.viscosity",this % viscosity)

    ! Get initial conditions
    CALL config % Get("brg1d.u0",u)
    ! If the character is empty - default the velocity
    ! to zero
    IF (TRIM(u) == "") THEN
      u = "u = 0.0"
    END IF
    INFO(TRIM(u))
    
    ! joe@fluidnumerics.com 7/18/2023
    ! This workaround is needed - feq-parse currently
    ! trips up when we use "t" as an independent variable
    ! with the functions that start with "t". Big problem 
    ! here, for sure.
    !CALL this % solution % SetEquation(1,u)
    this % solution % eqn(1) = EquationParser( TRIM(u), &
                                               (/'x','t'/) )

    CALL this % solution % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % solution % BoundaryInterp(gpuAccel=.FALSE.)

    ! Get boundary conditions for "u"
    CALL config % Get("brg1d.uL",uLEqn)
    CALL config % Get("brg1d.uR",uREqn)
    INFO("uL : "//TRIM(uLEqn))
    INFO("uR : "//TRIM(uREqn))
    this % uLEqn = EquationParser(TRIM(uLEqn), (/'x','t'/))
    this % uREqn = EquationParser(TRIM(uREqn), (/'x','t'/))
 
    this % uL = this % uLEqn % Evaluate((/this % geometry % x % boundary % hostData(1,1,1),this % t/))
    this % uR = this % uREqn % Evaluate((/this % geometry % x % boundary % hostData(1,2,this % geometry % x % nElem), this % t/))

    ! Get boundary conditions for du/dx
    CALL config % Get("brg1d.uxL",uLEqn)
    CALL config % Get("brg1d.uxR",uREqn)
    INFO("uxL : "//TRIM(uLEqn))
    INFO("uxR : "//TRIM(uREqn))
    this % uxLEqn = EquationParser(TRIM(uLEqn), (/'x','t'/))
    this % uxREqn = EquationParser(TRIM(uREqn), (/'x','t'/))
 
    this % uxL = this % uxLEqn % Evaluate((/this % geometry % x % boundary % hostData(1,1,1),this % t/))
    this % uxR = this % uxREqn % Evaluate((/this % geometry % x % boundary % hostData(1,2,this % geometry % x % nElem), this % t/))

    CALL this % CalculateEntropy()
    CALL this % ReportEntropy()

  END SUBROUTINE SetInitialConditions_Burgers1D

  SUBROUTINE PrintType_Burgers1D(this)
#undef __FUNC__
#define __FUNC__ "PrintType"
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(in) :: this

    INFO("Viscous Burgers (1D)")

  END SUBROUTINE PrintType_Burgers1D

  SUBROUTINE CalculateEntropy_Burgers1D(this)
  !! Base method for calculating entropy of a model
  !! Calculates the entropy as the integration of the
  !! squared tracer over the domain
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iVar,iEl
    REAL(prec) :: Jacobian,u
    REAL(prec) :: wi

    this % entropy = 0.0_PREC

    DO iEl = 1,this % geometry % x % nElem
      DO i = 0,this % geometry % x % interp % N

        ! Coordinate mapping Jacobian
        Jacobian = this % geometry % dxds % interior % hostData(i,1,iEl)

        ! Quadrature weights
        wi = this % geometry % x % interp % qWeights % hostData(i)

        ! Solution
        u = this % solution % interior % hostData(i,1,iEl)

        this % entropy = this % entropy + (0.5_PREC*u*u)*wi*Jacobian

      END DO
    END DO

  END SUBROUTINE CalculateEntropy_Burgers1D

  SUBROUTINE PreTendency_Burgers1D(this)
    !! Calculate the velocity and density weighted enthalpy at element interior and element boundaries
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl
    REAL(prec) :: x

    ! Evaluate the boundary conditions at this time level
    x = this % geometry % x % boundary % hostData(1,1,1) ! Left most x-point
    this % uL = this % uLEqn % Evaluate((/x,this % t/))
    this % uxL = this % uxLEqn % Evaluate((/x,this % t/))

    x = this % geometry % x % boundary % hostData(1,2,this % geometry % x % nElem) ! right most x-point
    this % uR = this % uREqn % Evaluate((/x, this % t/))
    this % uxR = this % uxREqn % Evaluate((/x, this % t/))

  END SUBROUTINE PreTendency_Burgers1D

  SUBROUTINE SetBoundaryCondition_Burgers1D(this)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,i
    INTEGER :: bcid,e2
    REAL(prec) :: u,v,nhat(1:2)

    ! Ths uL and uR attributes are set in the pre-tendency
    this % solution % extBoundary % hostData(1,1,1) = this % uL     ! Left most boundary
    this % solution % extBoundary % hostData(1,2,this % solution % nElem) = this % uR     ! Right most boundary

    ! Calculate the average value on the sides
    CALL this % solution % BassiRebaySides(this % gpuAccel)

    ! Calculate the BR gradient
    CALL this % solution % Derivative(this % geometry, &
                                      this % solutionGradient, &
                                      selfWeakBRForm, &
                                      this % gpuAccel)

    ! Boundary interp the BR gradient
    CALL this % solutionGradient % BoundaryInterp(this % gpuAccel)
    CALL this % solutionGradient % SideExchange(this % mesh,this % decomp,this % gpuAccel)

    ! Calculate the average value on the sides
    CALL this % solutionGradient % BassiRebaySides(this % gpuAccel)

    ! Set the external gradient values so that no-flux for viscous fluxes holds 
    this % solutionGradient % avgBoundary % hostData(1,1,1) =  -this % uxL ! 0.0_PREC ! ! left most boundary
    this % solutionGradient % avgBoundary % hostData(1,2,this % solution % nElem) = this % uxR !0.0_PREC !

  END SUBROUTINE SetBoundaryCondition_Burgers1D

  SUBROUTINE Source_Burgers1D(this)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    DO iEl = 1,this % source % nElem
      DO i = 0,this % source % interp % N

        this % source % interior % hostData(i,1,iEl) = 0.0_PREC

      END DO
    END DO

  END SUBROUTINE Source_Burgers1D

  SUBROUTINE Flux_Burgers1D(this)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl

    DO iEl = 1,this % solution % nElem
      DO i = 0,this % solution % interp % N

        ! u^2
        this % flux % interior % hostData(i,1,iEl) = &
          0.5_PREC*this % solution % interior % hostData(i,1,iEl)* &
          this % solution % interior % hostData(i,1,iEl) - &
          this % viscosity*this % solutionGradient % interior % hostData(i,1,iEl)

      END DO
    END DO

  END SUBROUTINE Flux_Burgers1D

  SUBROUTINE RiemannSolver_Burgers1D(this)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: unL,unR
    REAL(prec) :: alpha,nhat
    REAL(prec) :: fluxL
    REAL(prec) :: fluxR
    REAL(prec) :: dun

    DO iEl = 1,this % solution % nElem
      DO iSide = 1,2

        ! Get the boundary normals on cell edges from the mesh geometry
        IF (iSide == 1) THEN
          nhat = -1.0_PREC
        ELSE
          nhat = 1.0_PREC
        END IF

        ! Get the internal and external solution states
        unL = this % solution % boundary % hostData(1,iSide,iEl)
        unR = this % solution % extBoundary % hostData(1,iSide,iEl)

        fluxL = (0.5_PREC*unL*unL)*nHat
        fluxR = (0.5_PREC*unR*unR)*nHat

        alpha = MAX(ABS(unL),ABS(unR))

        ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
        this % flux % boundary % hostData(1,iSide,iEl) = 0.5_PREC*( &
                                                         fluxL + fluxR + alpha*(unL - unR))

        ! Add the viscous flux - nhat is not needed, since avgBoundary accounts for normal direction
        ! Gradient at the boundaries
        dun = this % solutionGradient % avgBoundary % hostData(1,iSide,iEl)
        this % flux % boundary % hostData(1,iSide,iEl) = this % flux % boundary % hostData(1,iSide,iEl) - &
                                                         this % viscosity*dun

      END DO
    END DO

  END SUBROUTINE RiemannSolver_Burgers1D

END MODULE SELF_Burgers1D
