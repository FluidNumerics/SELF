!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_cns2d

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_ECModel2D
  USE SELF_Config

  IMPLICIT NONE

#include "SELF_Macros.h"

  TYPE,EXTENDS(ECModel2D) :: cns2d
    !!
    !! For the solution attribute, we use the following convention
    !! iVar = 1 ~> rho*u (x-momentum)
    !! iVar = 2 ~> rho*v (y-momentum)
    !! iVar = 3 ~> rho (density)
    !! iVar = 4 ~> rho*E (Total Energy)

    !! For the diagnostics attribute, we use the following convention
    !! iVar = 1 ~> Kinetic Energy
    !! iVar = 2 ~> Enthalpy
    !! iVar = 3 ~> Sound Speed
    !! iVar = 4 ~> In-Situ Temperature
    !!
    !! primitive variables
    !!  u = rho*u/rho ~> x-component of velocity
    !!  v = rho*v/rho ~> y-component of velocity
    !!  rho = rho ~> density
    !!  p = p(rho,e) ~> equation of state
    !!
    !! environmentals
    !!  iVar = 1 ~> gravitational potential
    !!  iVar = 2 ~> Linear momentum drag
    !!
    PROCEDURE(Flux_cns2d),POINTER :: HyperbolicBoundaryFlux => LocalLaxFriedrichs
    PROCEDURE(Flux_cns2d),POINTER :: HyperbolicInteriorFlux => ConservativeFlux
    PROCEDURE(Flux_cns2d),POINTER :: ParabolicInteriorFlux => PrimitiveDiffusiveInteriorFlux
    PROCEDURE(Flux_cns2d),POINTER :: ParabolicBoundaryFlux => PrimitiveDiffusiveBoundaryFlux

    TYPE(MappedScalar2D) :: primitive ! Contains primitive variables
    TYPE(MappedScalar2D) :: entropyVars ! Contains entropy variables
    TYPE(MappedScalar2D) :: environmentals ! Functions to describe environmental features (e.g. gravitational potential)

    TYPE(MappedVector2D) :: primitiveGradient ! Gradient of the primitive variables
    TYPE(MappedVector2D) :: environmentalsGradient ! Functions to describe environmental features (e.g. gravitational potential)

    TYPE(MappedScalar2D) :: diagnostics

    TYPE(MappedScalar2D) :: prescribedSolution
    TYPE(MappedScalar2D) :: prescribedPrimitive
    TYPE(MappedScalar2D) :: prescribedDiagnostics

    REAL(prec) :: expansionFactor
    REAL(prec) :: Cp ! Heat capacity at constant pressure ( J/g/K )
    REAL(prec) :: Cv ! Heat capacity at constant volume ( J/g/K )
    REAL(prec) :: R ! Ideal gas constant (J/kg/K)
    REAL(prec) :: viscosity ! Dynamic viscosity (constant; for now) (kg * m^-1 /s^-1 )
    REAL(prec) :: diffusivity ! Thermal diffusivity

  CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_cns2d
    PROCEDURE :: Free => Free_cns2d
    PROCEDURE :: PrintType => PrintType_cns2d
    PROCEDURE :: PreTendency => PreTendency_cns2d
    PROCEDURE :: PreFlux => PreFlux_cns2d
    PROCEDURE :: CalculateEntropy => CalculateEntropy_cns2d
    PROCEDURE :: SetInitialConditions => SetInitialConditions_cns2d

    PROCEDURE :: WriteModel => Write_cns2d
    !PROCEDURE :: WriteTecplot => WriteTecplot_cns2d

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_cns2d
    PROCEDURE :: FluxMethod => FluxMethod_cns2d
    PROCEDURE :: RiemannSolver => RiemannSolver_cns2d
    PROCEDURE :: UpdateBoundary => UpdateBoundary_cns2d

    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_cns2d
    PROCEDURE,PRIVATE :: SetSolutionBoundaryCondition
    PROCEDURE,PRIVATE :: SetPrimitiveBoundaryCondition
    PROCEDURE,PRIVATE :: SetPrimitiveGradientBoundaryCondition
    PROCEDURE,PRIVATE :: SetDiagnosticsBoundaryCondition

    ! New Methods
    PROCEDURE :: CheckMinMax => CheckMinMax_cns2d
    PROCEDURE :: SetMaxCFL => SetMaxCFL_cns2d
    PROCEDURE :: HydrostaticAdjustment => HydrostaticAdjustment_cns2d

    ! Interior Hyperbolic Flux methods
    PROCEDURE,PRIVATE :: ConservativeFlux

    ! Riemann Fluxes for hyperbolic terms
    PROCEDURE,PRIVATE :: LocalLaxFriedrichs

    ! Interior parabolic (diffusive) flux methods
    PROCEDURE,PRIVATE :: PrimitiveDiffusiveInteriorFlux

    ! Boundary parabolic (diffusive) flux methods
    PROCEDURE,PRIVATE :: PrimitiveDiffusiveBoundaryFlux

    ! Methods to set flux models
    GENERIC :: SetFluxMethod => SetFluxMethod_withInt, &
      SetFluxMethod_withChar
    PROCEDURE,PRIVATE :: SetFluxMethod_withInt
    PROCEDURE,PRIVATE :: SetFluxMethod_withChar

    GENERIC :: SetVelocity => SetVelocityFromChar_cns2d
    PROCEDURE,PRIVATE :: SetVelocityFromChar_cns2d

    GENERIC :: SetGravity => SetGravityFromChar_cns2d
    PROCEDURE,PRIVATE :: SetGravityFromChar_cns2d

    GENERIC :: SetDrag => SetDragFromChar_cns2d, &
      SetDragFromConstant_cns2d
    PROCEDURE,PRIVATE :: SetDragFromChar_cns2d
    PROCEDURE,PRIVATE :: SetDragFromConstant_cns2d

    GENERIC :: SetPrescribedSolution => SetPrescribedSolutionFromChar_cns2d, &
      SetPrescribedSolutionFromSolution_cns2d
    PROCEDURE,PRIVATE :: SetPrescribedSolutionFromChar_cns2d
    PROCEDURE,PRIVATE :: SetPrescribedSolutionFromSolution_cns2d

    PROCEDURE :: SetCp => SetCp_cns2d
    PROCEDURE :: SetCv => SetCv_cns2d
    PROCEDURE :: SetGasConstant => SetGasConstant_cns2d
    PROCEDURE :: SetStatic => SetStatic_cns2d

    PROCEDURE,PRIVATE :: AddThermalBubble

    PROCEDURE,PRIVATE :: CalculateDiagnostics
    PROCEDURE,PRIVATE :: ConservativeToPrimitive
    PROCEDURE,PRIVATE :: ConservativeToEntropy

  END TYPE cns2d

  ! ---------------------------------------- !
  ! Riemann Flux integer flags
  !
  INTEGER,PARAMETER :: SELF_NLLF_CIG2D = 500 ! Naive Local Lax Friedrich's Riemann Flux

  ! ---------------------------------------- !
  ! Diagnostics variable indices
  !
  INTEGER,PARAMETER,PRIVATE :: nDiagnostics = 4

  ! ---------------------------------------- !
  ! Variables controlling hydrostatic adjustment
  !
  INTEGER,PARAMETER,PRIVATE :: hydrostaticAdjMaxIters = 1000

  ! ---------------------------------------- !
  !
  REAL(prec),PRIVATE :: Cp_static != 1.005_PREC
  REAL(prec),PRIVATE :: Cv_static != 0.718_PREC
  REAL(prec),PRIVATE :: R_static != 287.0_PREC ! J/kg/K
  REAL(prec),PRIVATE :: rho_static != 1.2754_PREC ! kg/m^3
  REAL(prec),PRIVATE :: T_static != 273.0_PREC ! K
  REAL(prec),PRIVATE :: e_static != 1.5_PREC*R_static*T_static ! J/kg

  INTERFACE
    SUBROUTINE Flux_cns2d(this)
      IMPORT cns2d
      IMPLICIT NONE
      CLASS(cns2d),INTENT(inout) :: this
    END SUBROUTINE Flux_cns2d
  END INTERFACE

  INTERFACE
    SUBROUTINE SetSolutionBoundaryCondition_cns2d_gpu_wrapper(scalar, &
                                                              extScalar,prescribedScalar,nHat,sideInfo,N,nVar,nEl) &
      BIND(c,name="SetSolutionBoundaryCondition_cns2d_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: scalar
      TYPE(C_PTR) :: extScalar
      TYPE(C_PTR) :: prescribedScalar
      TYPE(C_PTR) :: nHat
      TYPE(C_PTR) :: sideInfo
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SetSolutionBoundaryCondition_cns2d_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE SetDiagBoundaryCondition_cns2d_gpu_wrapper(scalar, &
                                                          extScalar,prescribedScalar,nHat,sideInfo,N,nVar,nEl) &
      BIND(c,name="SetDiagBoundaryCondition_cns2d_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: scalar
      TYPE(C_PTR) :: extScalar
      TYPE(C_PTR) :: prescribedScalar
      TYPE(C_PTR) :: nHat
      TYPE(C_PTR) :: sideInfo
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SetDiagBoundaryCondition_cns2d_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Source_cns2d_gpu_wrapper(source,solution,environmentalsGradient,N,nVar,nEl) &
      BIND(c,name="Source_cns2d_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: source,solution,environmentalsGradient
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE Source_cns2d_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ConservativeFlux_gpu_wrapper(flux,solution,primitive,N,nVar,nEl) &
      BIND(c,name="ConservativeFlux_cns2d_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: flux,solution,primitive
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ConservativeFlux_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE LocalLaxFriedrichs_gpu_wrapper(flux, &
                                              solution,extSolution,primitive,extPrimitive,diagnostics, &
                                              extDiagnostics,nHat,nScale,N,nVar,nDiag,nEl) &
      BIND(c,name="LocalLaxFriedrichs_cns2d_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: flux
      TYPE(C_PTR) :: solution
      TYPE(C_PTR) :: extSolution
      TYPE(C_PTR) :: primitive
      TYPE(C_PTR) :: extPrimitive
      TYPE(C_PTR) :: diagnostics
      TYPE(C_PTR) :: extDiagnostics
      TYPE(C_PTR) :: nHat
      TYPE(C_PTR) :: nScale
      INTEGER(C_INT),VALUE :: N,nVar,nDiag,nEl
    END SUBROUTINE LocalLaxFriedrichs_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE CalculateDiagnostics_cns2d_gpu_wrapper(solution, &
                                                      diagnostics,expansionFactor,R,N,nVar, &
                                                      nDiag,nEl) &
      BIND(c,name="CalculateDiagnostics_cns2d_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: solution
      TYPE(C_PTR) :: diagnostics
      REAL(c_prec),VALUE :: expansionFactor,R
      INTEGER(C_INT),VALUE :: N,nVar,nDiag,nEl
    END SUBROUTINE CalculateDiagnostics_cns2d_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ConservativeToPrimitive_cns2d_gpu_wrapper(solution, &
                                                         primitive,expansionFactor,N,nVar,nEl) &
      BIND(c,name="ConservativeToPrimitive_cns2d_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: solution
      TYPE(C_PTR) :: primitive
      REAL(c_prec),VALUE :: expansionFactor
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ConservativeToPrimitive_cns2d_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ConservativeToEntropy_cns2d_gpu_wrapper(solution, &
                                                       entropy,expansionFactor,N,nVar,nEl) &
      BIND(c,name="ConservativeToEntropy_cns2d_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(C_PTR) :: solution
      TYPE(C_PTR) :: entropy
      REAL(c_prec),VALUE :: expansionFactor
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ConservativeToEntropy_cns2d_gpu_wrapper
  END INTERFACE

CONTAINS

  SUBROUTINE Init_cns2d(this,nvar,mesh,geometry,decomp)
#undef __FUNC__
#define __FUNC__ "Init"
    IMPLICIT NONE
    CLASS(cns2d),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh2D),INTENT(in),TARGET :: mesh
    TYPE(SEMQuad),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp
    ! Local
    INTEGER :: ivar
    CHARACTER(LEN=3) :: ivarChar
    CHARACTER(LEN=25) :: varname
    INTEGER :: nvarloc

    ! Ensure that the number of variables is 4
    ! nvar is unused in this class extension
    nvarloc = 4

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % HyperbolicBoundaryFlux => LocalLaxFriedrichs
    this % HyperbolicInteriorFlux => ConservativeFlux
    this % ParabolicInteriorFlux => PrimitiveDiffusiveInteriorFlux
    this % ParabolicBoundaryFlux => PrimitiveDiffusiveBoundaryFlux
    this % gpuAccel = .FALSE.

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % prescribedSolution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    CALL this % diagnostics % Init(geometry % x % interp,nDiagnostics,this % mesh % nElem)
    CALL this % prescribedDiagnostics % Init(geometry % x % interp,nDiagnostics,this % mesh % nElem)

    CALL this % environmentals % Init(geometry % x % interp,2,this % mesh % nElem)

    CALL this % solutionGradient % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % primitiveGradient % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % environmentalsGradient % Init(geometry % x % interp,2,this % mesh % nElem)

    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    CALL this % flux % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    CALL this % primitive % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % prescribedPrimitive % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % entropyVars % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    ! First three variables are treated as u, v, eta
    ! Any additional are treated as passive tracers
    CALL this % solution % SetName(1,"rho*u")
    CALL this % solution % SetUnits(1,"kg/m^2/s")
    CALL this % solution % SetDescription(1,"x-component of the momentum per unit volume")

    CALL this % solution % SetName(2,"rho*v")
    CALL this % solution % SetUnits(2,"kg/m^2/s")
    CALL this % solution % SetDescription(2,"y-component of the momentum per unit volume")

    CALL this % solution % SetName(3,"rho")
    CALL this % solution % SetUnits(3,"kg/m^3")
    CALL this % solution % SetDescription(3,"fluid density; mass per unit volume")

    CALL this % solution % SetName(4,"rho*E")
    CALL this % solution % SetUnits(4,"kg/m/s^2")
    CALL this % solution % SetDescription(4,"Density weighted total (kinetic + internal) energy per unit volume")

    CALL this % diagnostics % SetName(1,"KE")
    CALL this % diagnostics % SetUnits(1,"kg/m/s^2")
    CALL this % diagnostics % SetDescription(1,"Kinetic energy per unit volume")

    CALL this % diagnostics % SetName(2,"H")
    CALL this % diagnostics % SetUnits(2,"kg/m/s^2")
    CALL this % diagnostics % SetDescription(2,"Density weighted enthalpy")

    CALL this % diagnostics % SetName(3,"c")
    CALL this % diagnostics % SetUnits(3,"m/s")
    CALL this % diagnostics % SetDescription(3,"Sound speed")

    CALL this % diagnostics % SetName(4,"T")
    CALL this % diagnostics % SetUnits(4,"K")
    CALL this % diagnostics % SetDescription(4,"In-Situ Temperature")

    CALL this % primitive % SetName(1,"u")
    CALL this % primitive % SetUnits(1,"m/s")
    CALL this % primitive % SetDescription(1,"Fluid velocity x-component")

    CALL this % primitive % SetName(2,"v")
    CALL this % primitive % SetUnits(2,"m/s")
    CALL this % primitive % SetDescription(2,"Fluid velocity x-component")

    CALL this % primitive % SetName(3,"rho")
    CALL this % primitive % SetUnits(3,"kg/m^3")
    CALL this % primitive % SetDescription(3,"Fluid density")

    CALL this % primitive % SetName(4,"P")
    CALL this % primitive % SetUnits(4,"kg/m/s^2")
    CALL this % primitive % SetDescription(4,"Fluid pressure")

    CALL this % environmentals % SetName(1,"gp")
    CALL this % environmentals % SetUnits(1,"m^2/s^2")
    CALL this % environmentals % SetDescription(1,"Gravitational Potential")

    CALL this % environmentals % SetName(2,"Cd")
    CALL this % environmentals % SetUnits(2,"m1/s")
    CALL this % environmentals % SetDescription(2,"Linear momentum drag")

  END SUBROUTINE Init_cns2d

  SUBROUTINE PrintType_cns2d(this)
#undef __FUNC__
#define __FUNC__ "PrintType"
    IMPLICIT NONE
    CLASS(cns2d),INTENT(in) :: this

    INFO("Compressible Ideal Gas (2D)")

  END SUBROUTINE PrintType_cns2d

  SUBROUTINE Free_cns2d(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this

    CALL this % solution % Free()
    CALL this % prescribedSolution % Free()
    CALL this % diagnostics % Free()
    CALL this % prescribedDiagnostics % Free()
    CALL this % environmentals % Free()
    CALL this % environmentalsGradient % Free()
    CALL this % workSol % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()
    CALL this % primitive % Free()
    CALL this % primitiveGradient % Free()
    CALL this % prescribedPrimitive % Free()
    CALL this % entropyVars % Free()

  END SUBROUTINE Free_cns2d

  SUBROUTINE SetInitialConditions_cns2d(this,config)
#undef __FUNC__
#define __FUNC__ "SetInitialConditions"
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    TYPE(SELFConfig),INTENT(inout) :: config
    ! Local
    LOGICAL :: setStaticState
    LOGICAL :: hydrostaticAdjust
    LOGICAL :: found
    CHARACTER(LEN=self_EquationLength) :: u ! x-velocity component
    CHARACTER(LEN=self_EquationLength) :: v ! y-velocity component
    CHARACTER(LEN=self_EquationLength) :: rho ! density
    CHARACTER(LEN=self_EquationLength) :: T  ! in-situ temperature
    CHARACTER(LEN=self_EquationLength) :: gp ! gravitational potential
    CHARACTER(LEN=self_EquationLength) :: featureType
    CHARACTER(LEN=SELF_JSON_DEFAULT_KEY_LENGTH) :: jsonKey
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
    INTEGER :: i,nfeatures
    REAL(prec) :: featureParams(1:10)

    ! Get static parameters
    CALL config % Get("cns2d.initial_conditions.static_state",setStaticState)
    CALL config % Get("cns2d.fluid.Cp",Cp_static)
    CALL config % Get("cns2d.fluid.Cv",Cv_static)
    CALL config % Get("cns2d.fluid.R",R_static)
    CALL config % Get("cns2d.fluid.rho",rho_static)
    CALL config % Get("cns2d.fluid.T",T_static)
    CALL config % Get("cns2d.fluid.energy",e_static)
    CALL config % Get("cns2d.fluid.viscosity",this % viscosity)
    CALL config % Get("cns2d.fluid.diffusivity",this % diffusivity)

    IF (setStaticState) THEN
      ! INFO("Set fluid to static state")
      ! WRITE (tUSC,"(ES16.7E3)") T_static
      ! str = usc2_'T\u2070 = '//TRIM(tUSC)
      ! INFO(str)

      ! WRITE (rhoUSC,"(ES16.7E3)") rho_static
      ! str = usc2_'\u03C1\u2070 = '//TRIM(rhoUSC)
      ! INFO(str)

      ! WRITE (CpUSC,"(ES16.7E3)") Cp_static
      ! str = usc2_'C\u209A = '//TRIM(CpUSC)
      ! INFO(str)

      ! WRITE (CvUSC,"(ES16.7E3)") Cv_static
      ! str = usc2_'C\u1D65 = '//TRIM(CvUSC)
      ! INFO(str)

      CALL this % SetStatic() ! Set field and parameters to STP
    END IF

    ! Get environmental parameters
    CALL config % Get("cns2d.environment.potential",gp)

    ! If the character is empty - default the gravitational
    ! potential to 0
    IF (TRIM(gp) == "") THEN
      gp = "p = 0.0"
    END IF
    ! Configure environmentals
    CALL this % SetGravity(gp)

    CALL config % Get("cns2d.initial_conditions.hydrostatic_adjustment",hydrostaticAdjust)
    IF (hydrostaticAdjust) THEN
      CALL config % Get("cns2d.initial_conditions.hydrostatic_momentum_max",momentumMax)
      CALL this % HydrostaticAdjustment(momentumMax)

    END IF

    ! Get additional initial conditions (add to static state if provided)
    CALL config % Get("cns2d.initial_conditions.u",u)
    CALL config % Get("cns2d.initial_conditions.v",v)
    CALL config % Get("cns2d.initial_conditions.rho",rho)
    CALL config % Get("cns2d.initial_conditions.T",T)

    ! If the character is empty - default the velocity
    ! components to zero
    IF (TRIM(u) == "") THEN
      u = "u = 0.0"
    END IF
    IF (TRIM(u) == "") THEN
      v = "v = 0.0"
    END IF

    velocity = (/u,v/)
    CALL this % SetVelocity(velocity)

    ! Set the time integrator
    CALL config % Get("time_options.integrator",integrator)
    CALL this % SetTimeIntegrator(TRIM(integrator)) ! Set the integrator
    CALL config % Get("time_options.dt",this % dt) ! Set the time step size
    CALL config % Get("time_options.start_time",this % t) ! Set the initial time

    ! Add Features, e.g. thermal bubble
    CALL config % concretization % info('cns2d.initial_conditions.features',n_children=nfeatures)
    DO i = 1,nfeatures
      WRITE (arrayCount,"(I0)") i

      jsonKey = "cns2d.initial_conditions.features["// &
                TRIM(arrayCount)// &
                "].type"
      CALL config % Get(TRIM(jsonKey),featureType)

      IF (UpperCase(TRIM(featureType)) == "THERMAL_BUBBLE") THEN

        jsonKey = "cns2d.initial_conditions.features["// &
                  TRIM(arrayCount)// &
                  "].parameters"
        CALL config % Get(TRIM(jsonKey),featureParams(1:4))

        CALL this % AddThermalBubble((/featureParams(1),featureParams(2)/), &
                                     featureParams(3),featureParams(4))

      END IF

    END DO

    CALL this % SetPrescribedSolution()

    CALL this % CheckMinMax()
    CALL this % CalculateEntropy()
    CALL this % ReportEntropy()

  END SUBROUTINE SetInitialConditions_cns2d

  SUBROUTINE SetFluxMethod_withInt(this,fluxMethod)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    INTEGER,INTENT(in) :: fluxMethod

    SELECT CASE (fluxMethod)

    CASE (SELF_NLLF_CIG2D)
      this % HyperbolicBoundaryFlux => LocalLaxFriedrichs
    CASE DEFAULT
      this % HyperbolicBoundaryFlux => LocalLaxFriedrichs
    END SELECT

  END SUBROUTINE SetFluxMethod_withInt

  SUBROUTINE SetFluxMethod_withChar(this,fluxMethod)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: fluxMethod
    ! Local
    CHARACTER(SELF_INTEGRATOR_LENGTH) :: upperCaseInt

    upperCaseInt = UpperCase(TRIM(fluxMethod))

    SELECT CASE (TRIM(upperCaseInt))

    CASE ("NAIVELLF")
      this % HyperbolicBoundaryFlux => LocalLaxFriedrichs
    CASE DEFAULT
      this % HyperbolicBoundaryFlux => LocalLaxFriedrichs

    END SELECT

  END SUBROUTINE SetFluxMethod_withChar

  SUBROUTINE AddThermalBubble(this,xc,R,Tmax)
    !! Adds a temperature anomaly to the fluid state
    !! The anomaly is a gaussian blob of radius R
    !! centered at xc with an extrema of Tmax.
    !!
    !! The density field is calculated so that the
    !! background pressure field remains undisturbed
    !!
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: xc(1:2)
    REAL(prec),INTENT(in) :: R
    REAL(prec),INTENT(in) :: Tmax
    ! Local
    INTEGER :: i,j,iEl,iVar
    REAL(prec) :: x(1:2)
    REAL(prec) :: Tprime,Rg,T

    Rg = this % R

    DO iEl = 1,this % source % nElem
      DO j = 0,this % source % interp % N
        DO i = 0,this % source % interp % N

          x = this % geometry % x % interior % hostData(1:2,i,j,1,iEl)

          Tprime = Tmax*EXP(-((x(1) - xc(1))**2 + (x(2) - xc(2))**2)/(2.0_PREC*R*R))

          T = this % diagnostics % interior % hostData(i,j,4,iEl) + Tprime

          ! Update the density
          this % solution % interior % hostData(i,j,3,iEl) = &
            this % solution % interior % hostData(i,j,3,iEl)*(1.0_PREC - Tprime/T)

          ! Add internal energy
          this % solution % interior % hostData(i,j,4,iEl) = &
            1.5_PREC*this % solution % interior % hostData(i,j,3,iEl)*Rg*T

        END DO
      END DO
    END DO

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateDevice()
      CALL this % primitive % UpdateDevice()
      CALL this % diagnostics % UpdateDevice()
    END IF

    CALL this % solution % BoundaryInterp(gpuAccel=this % gpuAccel)
    CALL this % CalculateDiagnostics()
    CALL this % ConservativeToPrimitive()
    CALL this % ConservativeToEntropy()

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateHost()
      CALL this % primitive % UpdateHost()
      CALL this % diagnostics % UpdateHost()
    END IF

  END SUBROUTINE AddThermalBubble

  SUBROUTINE SetStatic_cns2d(this)
  !! Sets the default fluid state as "air" at STP with
  !! no motion
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    CALL this % SetCv(Cv_static)
    CALL this % SetCp(Cp_static)
    CALL this % SetGasConstant(R_static)
    DO iEl = 1,this % source % nElem
      DO j = 0,this % source % interp % N
        DO i = 0,this % source % interp % N
          this % solution % interior % hostData(i,j,1,iEl) = 0.0_PREC ! rho*u
          this % solution % interior % hostData(i,j,2,iEl) = 0.0_PREC ! rho*v
          this % solution % interior % hostData(i,j,3,iEl) = rho_static ! rho
          this % solution % interior % hostData(i,j,4,iEl) = rho_static*e_static ! rho*E
        END DO
      END DO
    END DO

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateDevice()
    END IF

    CALL this % solution % BoundaryInterp(gpuAccel=this % gpuAccel)
    CALL this % CalculateDiagnostics()
    CALL this % ConservativeToPrimitive()
    CALL this % ConservativeToEntropy()

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateHost()
      CALL this % primitive % UpdateHost()
      CALL this % diagnostics % UpdateHost()
    END IF

  END SUBROUTINE SetStatic_cns2d

  SUBROUTINE SetVelocityFromChar_cns2d(this,eqnChar)
  !! Sets the fluid velocity field using the provided equation parser
  !! The momentum is then updated using the current fluid density field.
  !! From here, the PreTendency method is called to set other diagnostics
  !!
  !! The total energy field is calculated using the internal energy (diagnosed from the
  !! in-situ temperature) and the new kinetic energy field.
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:2)
    ! Local
    INTEGER :: i,j,iEl,iVar
    REAL(prec) :: rho,u,v,temperature,internalE,KE

    DO iVar = 1,2
      CALL this % primitive % SetEquation(ivar,eqnChar(iVar))
    END DO

    CALL this % primitive % SetInteriorFromEquation(this % geometry,this % t)
    IF (this % gpuAccel) THEN
      CALL this % primitive % UpdateDevice()
    END IF
    CALL this % primitive % BoundaryInterp(this % gpuAccel)

    DO iEl = 1,this % source % nElem
      DO j = 0,this % source % interp % N
        DO i = 0,this % source % interp % N
          rho = this % solution % interior % hostData(i,j,3,iEl)
          u = this % primitive % interior % hostData(i,j,1,iEl)
          v = this % primitive % interior % hostData(i,j,2,iEl)
          temperature = this % diagnostics % interior % hostData(i,j,4,iEl)
          internalE = 1.5_PREC*rho*this % R*temperature ! Internal energy
          KE = rho*(u*u + v*v)*0.5_PREC

          this % solution % interior % hostData(i,j,1,iEl) = rho*u  ! rho*u
          this % solution % interior % hostData(i,j,2,iEl) = rho*v ! rho*v
          this % solution % interior % hostData(i,j,4,iEl) = internalE + KE ! rho*E
        END DO
      END DO
    END DO

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateDevice()
    END IF

    CALL this % solution % BoundaryInterp(gpuAccel=this % gpuAccel)
    CALL this % CalculateDiagnostics()
    CALL this % ConservativeToPrimitive()
    CALL this % ConservativeToEntropy()

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateHost()
      CALL this % primitive % UpdateHost()
      CALL this % diagnostics % UpdateHost()
    END IF

  END SUBROUTINE SetVelocityFromChar_cns2d

  SUBROUTINE SetGravityFromChar_cns2d(this,eqnChar)
    !! Sets the gravitational acceleration from an equation input as a character
    !! The interior points are set and then the strong form of the gradient is used
    !! to calculate the gradient. After, environmentals and environmentalsGradient
    !! fields are copied to device memory (if an accelerator device is present)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    CHARACTER(LEN=*),INTENT(in) :: eqnChar

    CALL this % environmentals % SetEquation(1,eqnChar)

    CALL this % environmentals % SetInteriorFromEquation(this % geometry,this % t)

    CALL this % environmentals % GradientSF(this % geometry, &
                                            this % environmentalsGradient, &
                                            .FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % environmentals % UpdateDevice()
      CALL this % environmentalsGradient % UpdateDevice()
    END IF

  END SUBROUTINE SetGravityFromChar_cns2d

  SUBROUTINE SetDragFromChar_cns2d(this,eqnChar)
    !! Sets the momentum drag from an equation input as a character
    !! The interior points are set and then copied to the device
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    CHARACTER(LEN=*),INTENT(in) :: eqnChar

    CALL this % environmentals % SetEquation(2,eqnChar)
    CALL this % environmentals % SetInteriorFromEquation(this % geometry,this % t)

    IF (this % gpuAccel) THEN
      CALL this % environmentals % UpdateDevice()
    END IF

  END SUBROUTINE SetDragFromChar_cns2d

  SUBROUTINE SetDragFromConstant_cns2d(this,Cd)
    !! Sets the momentum drag from an equation input as a character
    !! The interior points are set and then copied to the device
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: Cd
    ! Local
    INTEGER :: i,j,iEl,iVar

    DO iEl = 1,this % solution % nElem
      DO j = 0,this % solution % interp % N
        DO i = 0,this % solution % interp % N

          this % environmentals % interior % hostData(i,j,2,iEl) = Cd

        END DO
      END DO
    END DO

    IF (this % gpuAccel) THEN
      CALL this % environmentals % UpdateDevice()
    END IF

  END SUBROUTINE SetDragFromConstant_cns2d

  SUBROUTINE SetPrescribedSolutionFromChar_cns2d(this,eqnChar)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % prescribedSolution % nVar)
    ! Local
    INTEGER :: iVar

    DO iVar = 1,this % prescribedSolution % nVar
      CALL this % prescribedSolution % SetEquation(ivar,eqnChar(iVar))
    END DO

    CALL this % prescribedSolution % SetInteriorFromEquation(this % geometry,this % t)
    CALL this % prescribedSolution % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % prescribedSolution % UpdateDevice()
    END IF

  END SUBROUTINE SetPrescribedSolutionFromChar_cns2d

  SUBROUTINE SetPrescribedSolutionFromSolution_cns2d(this)
  !! Sets the prescribed solution using the current solution attribute
  !! This can be useful for situations where you want to set the
  !! boundary conditions and initial conditions to be identical
  !!
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    DO iEl = 1,this % solution % nElem
      DO iVar = 1,this % solution % nvar
        DO j = 0,this % solution % interp % N
          DO i = 0,this % solution % interp % N

            this % prescribedSolution % interior % hostData(i,j,iVar,iEl) = &
              this % solution % interior % hostData(i,j,iVar,iEl)

            this % prescribedPrimitive % interior % hostData(i,j,iVar,iEl) = &
              this % primitive % interior % hostData(i,j,iVar,iEl)

          END DO
        END DO
      END DO
    END DO

    DO iEl = 1,this % diagnostics % nElem
      DO iVar = 1,this % diagnostics % nVar
        DO j = 0,this % diagnostics % interp % N
          DO i = 0,this % diagnostics % interp % N

            this % prescribedDiagnostics % interior % hostData(i,j,iVar,iEl) = &
              this % diagnostics % interior % hostData(i,j,iVar,iEl)

          END DO
        END DO
      END DO
    END DO

    CALL this % prescribedSolution % BoundaryInterp(gpuAccel=.FALSE.)
    CALL this % prescribedPrimitive % BoundaryInterp(gpuAccel=.FALSE.)
    CALL this % prescribedDiagnostics % BoundaryInterp(gpuAccel=.FALSE.)

    IF (this % gpuAccel) THEN
      CALL this % prescribedSolution % UpdateDevice()
      CALL this % prescribedPrimitive % UpdateDevice()
      CALL this % prescribedDiagnostics % UpdateDevice()
    END IF

  END SUBROUTINE SetPrescribedSolutionFromSolution_cns2d

  SUBROUTINE SetCp_cns2d(this,Cp)
  !! Accessor routine to set the heat capacity at constant pressure
  !! Also updates the expansionFactor attribute when called.
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: Cp

    this % Cp = Cp
    IF (this % Cv == 0.0_PREC) THEN
      PRINT *, "Warning : Expansion factor not set; Cv = 0"
    END IF
    this % expansionFactor = this % Cp/this % Cv

  END SUBROUTINE SetCp_cns2d

  SUBROUTINE SetCv_cns2d(this,Cv)
  !! Accessor routine to set the heat capacity at constant volume
  !! Also updates the expansionFactor attribute when called.
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: Cv

    this % Cv = Cv
    this % expansionFactor = this % Cp/this % Cv

  END SUBROUTINE SetCv_cns2d

  SUBROUTINE SetGasConstant_cns2d(this,R)
  !! Accessor routine to set the ideal gas constant
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: R

    this % R = R

  END SUBROUTINE SetGasConstant_cns2d

  SUBROUTINE CalculateEntropy_cns2d(this)
  !!
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,j,i
    REAL(prec) :: Jacobian,wi,wj
    REAL(prec) :: P,rho
    REAL(prec) :: entropy

    IF (this % gpuAccel) THEN
      CALL this % solution % interior % UpdateHost()
      CALL this % primitive % interior % UpdateHost()
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

          rho = this % solution % interior % hostData(i,j,3,iEl)
          P = this % primitive % interior % hostData(i,j,4,iEl)

          entropy = entropy + &
                    rho*(LOG(P) - this % expansionFactor*LOG(rho))/ &
                    (this % expansionFactor - 1.0_PREC)*wi*wj*Jacobian

        END DO
      END DO
    END DO

    CALL this % decomp % GlobalReduce(entropy,this % entropy)

  END SUBROUTINE CalculateEntropy_cns2d

  SUBROUTINE CheckMinMax_cns2d(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iVar

    IF (this % gpuAccel) THEN
      CALL this % solution % UpdateHost()
      CALL this % environmentals % UpdateHost()
      CALL this % primitive % UpdateHost()
      CALL this % diagnostics % UpdateHost()
    END IF

    PRINT *, '---------------------'
    DO iVar = 1,this % solution % nVar
      PRINT *, TRIM(this % solution % meta(iVar) % name)//" (t, min, max) :", &
        this % t, &
        MINVAL(this % solution % interior % hostData(:,:,iVar,:)), &
        MAXVAL(this % solution % interior % hostData(:,:,iVar,:))
    END DO

    DO iVar = 1,this % environmentals % nVar
      PRINT *, TRIM(this % environmentals % meta(iVar) % name)//" (t, min, max) :", &
        this % t, &
        MINVAL(this % environmentals % interior % hostData(:,:,iVar,:)), &
        MAXVAL(this % environmentals % interior % hostData(:,:,iVar,:))
    END DO

    DO iVar = 1,this % primitive % nVar
      PRINT *, TRIM(this % primitive % meta(iVar) % name)//" (t, min, max) :", &
        this % t, &
        MINVAL(this % primitive % interior % hostData(:,:,iVar,:)), &
        MAXVAL(this % primitive % interior % hostData(:,:,iVar,:))
    END DO

    DO iVar = 1,this % diagnostics % nVar
      PRINT *, TRIM(this % diagnostics % meta(iVar) % name)//" (t, min, max) :", &
        this % t, &
        MINVAL(this % diagnostics % interior % hostData(:,:,iVar,:)), &
        MAXVAL(this % diagnostics % interior % hostData(:,:,iVar,:))
    END DO

    PRINT *, '---------------------'

  END SUBROUTINE CheckMinMax_cns2d

  SUBROUTINE SetMaxCFL_cns2d(this,cfl)
    !! This method uses the model grid and sound speed
    !! to set the time step size so that the desired
    !! maximum cfl number fixed.
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: cfl
    REAL(prec) :: dxMin
    REAL(prec) :: diagMax(nDiagnostics)
    REAL(prec) :: currentDt,currentTime,tn
    REAL(prec) :: Cd
    REAL(prec) :: dsdtMax(this % solution % nVar)
    REAL(prec) :: sMax(this % solution % nVar)

    dxMin = this % geometry % CovariantArcMin()
    diagMax = this % diagnostics % AbsMaxInterior() ! Get the absolute max of the diagnostics

    PRINT *, "Min(dx) : ",dxMin
    PRINT *, "Max(c) : ",diagMax(3)

    ! Reassign the time step for the hydrostatic adjustment
    ! so that the max CFL number is 0.5
    currentDt = this % dt
    this % dt = cfl*dxMin/diagMax(3)

    PRINT *, "Adjusted time step size : ",this % dt

  END SUBROUTINE SetMaxCFL_cns2d

  SUBROUTINE HydrostaticAdjustment_cns2d(this,tolerance)
    !! This method can be used to adjust a compressible fluid
    !! to hydrostatic equilibrium. On input, the cns2d
    !! object is expected to have the gravitational potential set to
    !! a non-zero field, the density and temperature fields are
    !! assumed uniform, and the velocity field corresponds to a motionless
    !! fluid.
    !!
    !! To adjust the fluid to hydrostatic equilibrium, an artificial
    !! momentum drag term is introduced to the momentum equations. The model
    !! is stepped forward until the fluid tendency absolute max value is
    !! less than a specified tolerance (optional input).
    !!
    !! The drag coefficient size is chosen to be
    !!
    !!  Cd = rac{max(c)}{min(\Delta x)}
    !!
    !! In this subroutine, the time step is chosen so that the sound-wave
    !! maximum CFL number is 0.75
    !!
    !!   CFL = rac{max(c) \Delta t}{min(\Delta x)} = 0.75
    !!
    !!
    !!
    !! ightarrow\Delta t = 0.75rac{MIN(\Delta x) }{\max{c}}
    !!
    !! After adjustment, the
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: tolerance
    ! Local
    INTEGER :: i
    REAL(prec) :: currentDt,currentTime,tn
    REAL(prec) :: Cd
    REAL(prec) :: dsdtMax(this % solution % nVar)
    REAL(prec) :: sMax(this % solution % nVar)
    REAL(prec) :: dxMin
    REAL(prec) :: diagMax(nDiagnostics)

    dxMin = this % geometry % CovariantArcMin()
    diagMax = this % diagnostics % AbsMaxInterior() ! Get the absolute max of the diagnostics

    ! Reassign the time step for the hydrostatic adjustment
    ! so that the max CFL number is 0.5
    currentDt = this % dt
    CALL this % SetMaxCFL(0.5_PREC)

    ! Calculate the drag coefficient
    Cd = 0.3*diagMax(3)/dxMin

    PRINT *, "Drag coefficient : ",Cd

    CALL this % SetDrag(Cd)

    ! Save the current time
    currentTime = this % t

    tn = 1000.0_PREC*this % dt
    DO i = 1,hydrostaticAdjMaxIters

      this % t = 0.0_PREC
      CALL this % ForwardStep(tn)

      sMax = this % solution % AbsMaxInterior()

      PRINT *, "Momentum Max : ",SQRT(sMax(1)**2 + sMax(2)**2)

      ! Check if momentum max is small in comparison to tolerance
      IF (SQRT(sMax(1)**2 + sMax(2)**2) <= tolerance) THEN
        PRINT *, "Reached hydrostatic equilibrium in ",i," iterations"
        EXIT
      END IF
    END DO

    ! Reset delta t and time
    this % dt = currentDt
    this % t = currentTime

  END SUBROUTINE HydrostaticAdjustment_cns2d

  SUBROUTINE PreTendency_cns2d(this)
    !! Calculate the velocity and density weighted enthalpy at element interior and element boundaries
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this

    CALL this % CalculateDiagnostics()
    CALL this % ConservativeToPrimitive()
    CALL this % ConservativeToEntropy()

  END SUBROUTINE PreTendency_cns2d

  SUBROUTINE CalculateDiagnostics(this)
    !! Calculates
    !!  * kinetic energy
    !!  * speed of sound
    !!  * enthalpy
    !!  * in-situ temperature
    !!
    !! Kinetic Energy
    !!
    !!  We recognize there are two ways we can calculate the kinetic
    !!  energy with the given prognostic and diagnostic variables
    !!  we track.
    !!
    !!  Option 1.
    !!
    !!    Use the velocity field diagnostic with the prognostic density.
    !!    This would result in
    !!
    !!     KE = 0.5_prec*rho*( u*u + v*v )
    !!
    !!    where
    !!
    !!      u = (rho*u)/rho
    !!      v = (rho*v)/rho
    !!
    !!  Option 2.
    !!
    !!    Use the prognostic momentum and density fields. This would
    !!    result in
    !!
    !!     KE = 0.5_prec*( (rho*u)*(rho*u) + (rho*v)*(rho*v) )/rho
    !!
    !!  Analytically, the two options are identical. In floating point
    !!  arithmetic, these are different.
    !!
    !!  It's currently unclear which option is more advantageous (and when),
    !!  and I am arbitrarily implementing Option 2.
    !!
    !!  If you find a good reason to use Option 1, or some other approach to
    !!  calculate kinetic energy that is more advantageous, develop an
    !!  example that highlights the benefits of your approach and open a
    !!  pull request.
    !!
    !! Pressure
    !!
    !!  We use the Ideal Gas Law
    !!
    !!     p = (\gamma-1)*
    !! ho*e
    !!
    !!  where $\gamma = rac{C_p}{C_v}$ is the expansion coefficient,
    !!  $
    !! ho$is the fluid density,and the internal energy.
    !!
    !!  We calculate $
    !! ho*e$as
    !!
    !!    rho*e = (rho*E - 0.5_prec*rho*KE)
    !!
    !!  where rho*E is the total energy, a prognostic variable, modelled
    !!  by this class, and
    !!
    !! Sound Speed
    !!
    !!  The speed of sound is defined through the relation
    !!
    !!    rac{\partial P}{\partial
    !! ho} = c^{2}
    !!
    !!  Then, we have that
    !!
    !!    c = ((\gamma-1)*e)^{1/2}
    !!
    !!  where gamma = Cp/Cv is the expansion coefficient,
    !!  rho is the fluid density, and e is the internal energy.
    !!
    !!  Using the equation of state,
    !!
    !!   (\gamma-1)*e = p/
    !! ho
    !!
    !!  We calculate $ as
    !!
    !!    e = (
    !! ho*E - 0.5_PREC*
    !! ho*KE)/
    !! ho
    !!
    !!
    !!  where rho*E is the total energy, a prognostic variable, modelled
    !!  by this class, $
    !! ho*KE the kinetic energy(a required diagnostic)
    !!  and, $
    !! ho the density(a prognostic variable) .
    !!
    !! In-Situ Temperature
    !!
    !!    T = 2/3*e/R
    !!

    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,j,i
    REAL(prec) :: rho,rhoU,rhoV,u,v,rhoE,rhoKE,p

    IF (this % gpuAccel) THEN

      CALL CalculateDiagnostics_cns2d_gpu_wrapper( &
        this % solution % interior % deviceData, &
        this % diagnostics % interior % deviceData, &
        this % expansionFactor, &
        this % R, &
        this % solution % interp % N, &
        this % solution % nVar, &
        this % diagnostics % nVar, &
        this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO j = 0,this % solution % interp % N
          DO i = 0,this % solution % interp % N

            rhoU = this % solution % interior % hostData(i,j,1,iEl)
            rhoV = this % solution % interior % hostData(i,j,2,iEl)
            rho = this % solution % interior % hostData(i,j,3,iEl)
            rhoE = this % solution % interior % hostData(i,j,4,iEl)
            u = rhoU/rho
            v = rhoV/rho
            rhoKE = 0.5_PREC*(rhoU*u + rhoV*v)
            p = (this % expansionFactor - 1.0_PREC)*(rhoE - rhoKE)

            ! Kinetic Energy
            this % diagnostics % interior % hostData(i,j,1,iEl) = rhoKE

            ! Enthalpy
            this % diagnostics % interior % hostData(i,j,2,iEl) = rhoE + p

            ! Speed of sound
            this % diagnostics % interior % hostData(i,j,3,iEl) = &
              SQRT(this % expansionFactor*p/rho)

            ! In-Situ Temperature
            this % diagnostics % interior % hostData(i,j,4,iEl) = &
              (2.0_PREC/3.0_PREC)*((rhoE - rhoKE)/rho)/this % R

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE CalculateDiagnostics

  SUBROUTINE ConservativeToPrimitive(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,j,i
    REAL(prec) :: rhoU,rhoV,rho,u,v,rhoE,rhoKE,p

    IF (this % gpuAccel) THEN

      CALL ConservativeToPrimitive_cns2d_gpu_wrapper( &
        this % solution % interior % deviceData, &
        this % primitive % interior % deviceData, &
        this % expansionFactor, &
        this % solution % interp % N, &
        this % solution % nVar, &
        this % solution % nElem)

    ELSE
      DO iEl = 1,this % solution % nElem
        DO j = 0,this % solution % interp % N
          DO i = 0,this % solution % interp % N

            rhoU = this % solution % interior % hostData(i,j,1,iEl)
            rhoV = this % solution % interior % hostData(i,j,2,iEl)
            rho = this % solution % interior % hostData(i,j,3,iEl)
            rhoE = this % solution % interior % hostData(i,j,4,iEl)
            u = rhoU/rho
            v = rhoV/rho
            rhoKE = 0.5_PREC*(rhoU*u + rhoV*v)
            p = (this % expansionFactor - 1.0_PREC)*(rhoE - rhoKE)

            this % primitive % interior % hostData(i,j,1,iEl) = u
            this % primitive % interior % hostData(i,j,2,iEl) = v
            this % primitive % interior % hostData(i,j,3,iEl) = rho
            this % primitive % interior % hostData(i,j,4,iEl) = p

          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE ConservativeToPrimitive

  SUBROUTINE ConservativeToEntropy(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,j,i
    REAL(prec) :: rhoU,rhoV,rho,u,v,E,KE,p,s

    IF (this % gpuAccel) THEN

      CALL ConservativeToEntropy_cns2d_gpu_wrapper( &
        this % solution % interior % deviceData, &
        this % entropyVars % interior % deviceData, &
        this % expansionFactor, &
        this % solution % interp % N, &
        this % solution % nVar, &
        this % solution % nElem)

    ELSE
      DO iEl = 1,this % solution % nElem
        DO j = 0,this % solution % interp % N
          DO i = 0,this % solution % interp % N

            rhoU = this % solution % interior % hostData(i,j,1,iEl)
            rhoV = this % solution % interior % hostData(i,j,2,iEl)
            rho = this % solution % interior % hostData(i,j,3,iEl)
            E = this % solution % interior % hostData(i,j,4,iEl)
            u = rhoU/rho
            v = rhoV/rho
            KE = 0.5_PREC*(rhoU*u + rhoV*v)
            p = (this % expansionFactor - 1.0_PREC)*(E - KE)
            s = LOG(p) - this % expansionFactor*LOG(rho)

            this % entropyVars % interior % hostData(i,j,1,iEl) = u*rho/p
            this % entropyVars % interior % hostData(i,j,2,iEl) = v*rho/p
            this % entropyVars % interior % hostData(i,j,3,iEl) = (this % expansionFactor - s)/ &
                                                                  (this % expansionFactor - 1.0_PREC) - KE/p
            this % entropyVars % interior % hostData(i,j,4,iEl) = -rho/p

          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE ConservativeToEntropy

  SUBROUTINE UpdateBoundary_cns2d(this)
  !! Interpolates the solution to element boundaries and performs the sideExchange
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this

    ! Interpolate velocity and required diagnostics to the element boundaries
    CALL this % solution % BoundaryInterp(this % gpuAccel)
    CALL this % primitive % BoundaryInterp(this % gpuAccel)
    CALL this % diagnostics % BoundaryInterp(this % gpuAccel)
    CALL this % entropyVars % BoundaryInterp(this % gpuAccel)

    ! Perform any MPI exchanges for the velocity and the required diagnostics
    ! across shared element faces between neighboring processes.
    CALL this % solution % SideExchange(this % mesh,this % decomp,this % gpuAccel)
    CALL this % primitive % SideExchange(this % mesh,this % decomp,this % gpuAccel)
    CALL this % diagnostics % SideExchange(this % mesh,this % decomp,this % gpuAccel)
    CALL this % entropyVars % SideExchange(this % mesh,this % decomp,this % gpuAccel)

  END SUBROUTINE UpdateBoundary_cns2d

  SUBROUTINE SetSolutionBoundaryCondition(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,i
    INTEGER :: bcid,e2
    REAL(prec) :: u,v,nhat(1:2)

    IF (this % gpuAccel) THEN

      CALL SetSolutionBoundaryCondition_cns2d_gpu_wrapper( &
        this % solution % boundary % deviceData, &
        this % solution % extBoundary % deviceData, &
        this % prescribedSolution % boundary % deviceData, &
        this % geometry % nHat % boundary % deviceData, &
        this % mesh % sideInfo % deviceData, &
        this % solution % interp % N, &
        this % solution % nVar, &
        this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iSide = 1,4
          bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
          e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
          IF (e2 == 0) THEN
            IF (bcid == SELF_BC_NONORMALFLOW) THEN

              DO i = 0,this % solution % interp % N
                nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
                ! Conservative variables
                u = this % solution % boundary % hostData(i,1,iSide,iEl)
                v = this % solution % boundary % hostData(i,2,iSide,iEl)
                this % solution % extBoundary % hostData(i,1,iSide,iEl) = &
                  (nhat(2)**2 - nhat(1)**2)*u - 2.0_PREC*nhat(1)*nhat(2)*v
                this % solution % extBoundary % hostData(i,2,iSide,iEl) = &
                  (nhat(1)**2 - nhat(2)**2)*v - 2.0_PREC*nhat(1)*nhat(2)*u
                this % solution % extBoundary % hostData(i,3,iSide,iEl) = &
                  this % solution % boundary % hostData(i,3,iSide,iEl)
                this % solution % extBoundary % hostData(i,4,iSide,iEl) = &
                  this % solution % boundary % hostData(i,4,iSide,iEl)

              END DO

            ELSEIF (bcid == SELF_BC_PRESCRIBED .OR. bcid == SELF_BC_RADIATION) THEN

              DO i = 0,this % solution % interp % N

                this % solution % extBoundary % hostData(i,1,iSide,iEl) = &
                  this % prescribedSolution % boundary % hostData(i,1,iSide,iEl)
                this % solution % extBoundary % hostData(i,2,iSide,iEl) = &
                  this % prescribedSolution % boundary % hostData(i,2,iSide,iEl)
                this % solution % extBoundary % hostData(i,3,iSide,iEl) = &
                  this % prescribedSolution % boundary % hostData(i,3,iSide,iEl)
                this % solution % extBoundary % hostData(i,4,iSide,iEl) = &
                  this % prescribedSolution % boundary % hostData(i,4,iSide,iEl)

              END DO

            END IF

          END IF

        END DO
      END DO
    END IF

  END SUBROUTINE SetSolutionBoundaryCondition

  SUBROUTINE SetPrimitiveBoundaryCondition(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,i
    INTEGER :: bcid,e2
    REAL(prec) :: u,v,nhat(1:2)

    IF (this % gpuAccel) THEN

      CALL SetSolutionBoundaryCondition_cns2d_gpu_wrapper( &
        this % primitive % boundary % deviceData, &
        this % primitive % extBoundary % deviceData, &
        this % prescribedPrimitive % boundary % deviceData, &
        this % geometry % nHat % boundary % deviceData, &
        this % mesh % sideInfo % deviceData, &
        this % primitive % interp % N, &
        this % primitive % nVar, &
        this % primitive % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iSide = 1,4
          bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
          e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
          IF (e2 == 0) THEN
            IF (bcid == SELF_BC_NONORMALFLOW) THEN

              DO i = 0,this % solution % interp % N
                nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
                u = this % primitive % boundary % hostData(i,1,iSide,iEl)
                v = this % primitive % boundary % hostData(i,2,iSide,iEl)
                this % primitive % extBoundary % hostData(i,1,iSide,iEl) = &
                  (nhat(2)**2 - nhat(1)**2)*u - 2.0_PREC*nhat(1)*nhat(2)*v
                this % primitive % extBoundary % hostData(i,2,iSide,iEl) = &
                  (nhat(1)**2 - nhat(2)**2)*v - 2.0_PREC*nhat(1)*nhat(2)*u
                this % primitive % extBoundary % hostData(i,3,iSide,iEl) = &
                  this % primitive % boundary % hostData(i,3,iSide,iEl)
                this % primitive % extBoundary % hostData(i,4,iSide,iEl) = &
                  this % primitive % boundary % hostData(i,4,iSide,iEl)

              END DO

            ELSEIF (bcid == SELF_BC_PRESCRIBED .OR. bcid == SELF_BC_RADIATION) THEN

              DO i = 0,this % solution % interp % N

                this % primitive % extBoundary % hostData(i,1,iSide,iEl) = &
                  this % prescribedPrimitive % boundary % hostData(i,1,iSide,iEl)
                this % primitive % extBoundary % hostData(i,2,iSide,iEl) = &
                  this % prescribedPrimitive % boundary % hostData(i,2,iSide,iEl)
                this % primitive % extBoundary % hostData(i,3,iSide,iEl) = &
                  this % prescribedPrimitive % boundary % hostData(i,3,iSide,iEl)
                this % primitive % extBoundary % hostData(i,4,iSide,iEl) = &
                  this % prescribedPrimitive % boundary % hostData(i,4,iSide,iEl)

              END DO

            END IF

          END IF

        END DO
      END DO
    END IF

  END SUBROUTINE SetPrimitiveBoundaryCondition

  SUBROUTINE SetDiagnosticsBoundaryCondition(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,i,iVar
    INTEGER :: bcid,e2
    REAL(prec) :: u,v,nhat(1:2)

    IF (this % gpuAccel) THEN

      CALL SetDiagBoundaryCondition_cns2d_gpu_wrapper( &
        this % diagnostics % boundary % deviceData, &
        this % diagnostics % extBoundary % deviceData, &
        this % prescribedDiagnostics % boundary % deviceData, &
        this % geometry % nHat % boundary % deviceData, &
        this % mesh % sideInfo % deviceData, &
        this % diagnostics % interp % N, &
        this % diagnostics % nVar, &
        this % diagnostics % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iSide = 1,4
          bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
          e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
          IF (e2 == 0) THEN
            IF (bcid == SELF_BC_NONORMALFLOW) THEN

              DO iVar = 1,nDiagnostics
                DO i = 0,this % solution % interp % N

                  ! Prolong the diagnostic values to the external state
                  this % diagnostics % extBoundary % hostData(i,iVar,iSide,iEl) = &
                    this % diagnostics % boundary % hostData(i,iVar,iSide,iEl)

                END DO
              END DO

            ELSEIF (bcid == SELF_BC_PRESCRIBED .OR. bcid == SELF_BC_RADIATION) THEN

              DO iVar = 1,nDiagnostics
                DO i = 0,this % solution % interp % N
                  this % diagnostics % extBoundary % hostData(i,iVar,iSide,iEl) = &
                    this % prescribedDiagnostics % boundary % hostData(i,iVar,iSide,iEl)
                END DO
              END DO

            END IF

          END IF

        END DO
      END DO
    END IF

  END SUBROUTINE SetDiagnosticsBoundaryCondition

  SUBROUTINE PreFlux_cns2d(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this


    ! Compute the gradient of the solution using the Bassi-Rebay Method
    ! TO DO : IF (diffusiveVars == primitive)
    CALL this % primitive % GradientBR(this % geometry, &
                                       this % primitiveGradient, &
                                       this % gpuAccel)

    CALL this % primitiveGradient % BoundaryInterp(this % gpuAccel)
    CALL this % primitiveGradient % SideExchange(this % mesh,this % decomp,this % gpuAccel)
    CALL this % SetPrimitiveGradientBoundaryCondition()
    CALL this % primitiveGradient % BassiRebaySides(this % gpuAccel)

  END SUBROUTINE PreFlux_cns2d

  SUBROUTINE SetBoundaryCondition_cns2d(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this

    CALL this % SetSolutionBoundaryCondition()
    CALL this % SetPrimitiveBoundaryCondition()
    CALL this % SetDiagnosticsBoundaryCondition()


  END SUBROUTINE SetBoundaryCondition_cns2d

  SUBROUTINE SetPrimitiveGradientBoundaryCondition(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,i
    INTEGER :: bcid,e2,ivar
    REAL(prec) :: u,v,nhat(1:2)

    ! TO DO : GPU kernel !

    DO iEl = 1,this % solution % nElem
      DO iSide = 1,4
        bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
        e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
        IF (e2 == 0) THEN

          DO ivar = 1,this % solution % nvar
          DO i = 0,this % solution % interp % N

            ! For now, prolong the gradient
            this % solutionGradient % extBoundary % hostData(1:2,i,ivar,iSide,iEl) = &
              this % solutionGradient % boundary % hostData(1:2,i,ivar,iSide,iEl)

          END DO
          END DO

        END IF

      END DO
    END DO

  END SUBROUTINE SetPrimitiveGradientBoundaryCondition

  SUBROUTINE Source_cns2d(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar
    REAL(prec) :: rhou,rhov,rho,gx,gy,Cd

    IF (this % gpuAccel) THEN

      CALL Source_cns2d_gpu_wrapper(this % source % interior % deviceData, &
                                    this % solution % interior % deviceData, &
                                    this % environmentalsGradient % interior % deviceData, &
                                    this % source % interp % N, &
                                    this % source % nVar, &
                                    this % source % nElem)

    ELSE

      DO iEl = 1,this % source % nElem
        DO j = 0,this % source % interp % N
          DO i = 0,this % source % interp % N

            rhou = this % solution % interior % hostData(i,j,1,iEl)
            rhov = this % solution % interior % hostData(i,j,2,iEl)
            rho = this % solution % interior % hostData(i,j,3,iEl)
            gx = this % environmentalsGradient % interior % hostData(1,i,j,1,iEl)
            gy = this % environmentalsGradient % interior % hostData(2,i,j,1,iEl)
            Cd = this % environmentals % interior % hostData(i,j,2,iEl)

            this % source % interior % hostData(i,j,1,iEl) = -rho*gx - Cd*rhou
            this % source % interior % hostData(i,j,2,iEl) = -rho*gy - Cd*rhov
            this % source % interior % hostData(i,j,3,iEl) = 0.0_PREC
            this % source % interior % hostData(i,j,4,iEl) = -rhou*gx - rhou*gy - &
                                                             Cd*(rhou*rhou + rhov*rhov)/rho

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE Source_cns2d

  SUBROUTINE ConservativeFlux(this)
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,n,iEl,iVar

    IF (this % gpuAccel) THEN
      CALL ConservativeFlux_gpu_wrapper(this % flux % physical % deviceData, &
                                        this % solution % interior % deviceData, &
                                        this % primitive % interior % deviceData, &
                                        this % solution % interp % N, &
                                        this % solution % nVar, &
                                        this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1,this % solution % nVar
          DO j = 0,this % solution % interp % N
            DO i = 0,this % solution % interp % N

              IF (iVar == 1) THEN ! rho*u

                DO n = 0,this % solution % interp % N
                  ! rho*u*u + p
                  this % flux % physical % hostData(1,1,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(n,j,1,iEl)* & ! u
                     this % solution % interior % hostData(n,j,1,iEl) + & ! rho*u
                     this % primitive % interior % hostData(n,j,4,iEl))

                  ! rho*u*v
                  this % flux % physical % hostData(2,1,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(n,j,2,iEl)* & ! v
                     this % solution % interior % hostData(n,j,1,iEl)) ! rho*u

                  ! rho*u*u + p
                  this % flux % physical % hostData(1,2,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(i,n,1,iEl)* & ! u
                     this % solution % interior % hostData(i,n,1,iEl) + & ! rho*u
                     this % primitive % interior % hostData(i,n,4,iEl))

                  ! rho*u*v
                  this % flux % physical % hostData(2,2,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(i,n,2,iEl)* & ! v
                     this % solution % interior % hostData(i,n,1,iEl)) ! rho*u
                END DO

              ELSEIF (iVar == 2) THEN ! rho*v

                DO n = 0,this % solution % interp % N
                  ! rho*v*u
                  this % flux % physical % hostData(1,1,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(n,j,1,iEl)* & ! u
                     this % solution % interior % hostData(n,j,2,iEl)) ! rho*v

                  ! rho*v*v + p
                  this % flux % physical % hostData(2,1,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(n,j,2,iEl)* & ! v
                     this % solution % interior % hostData(n,j,2,iEl) + & ! rho*v
                     this % primitive % interior % hostData(n,j,4,iEl)) ! pressure

                  ! rho*v*u
                  this % flux % physical % hostData(1,2,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(i,n,1,iEl)* & ! u
                     this % solution % interior % hostData(i,n,2,iEl)) ! rho*v

                  ! rho*v*v + p
                  this % flux % physical % hostData(2,2,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(i,n,2,iEl)* & ! v
                     this % solution % interior % hostData(i,n,2,iEl) + & ! rho*v
                     this % primitive % interior % hostData(i,n,4,iEl)) ! pressure

                END DO

              ELSEIF (iVar == 3) THEN ! density

                DO n = 0,this % solution % interp % N
                  this % flux % physical % hostData(1,1,n,i,j,iVar,iEl) = &
                    (this % solution % interior % hostData(n,j,1,iEl)) !rho*u

                  this % flux % physical % hostData(2,1,n,i,j,iVar,iEl) = &
                    (this % solution % interior % hostData(n,j,2,iEl)) !rho*v

                  this % flux % physical % hostData(1,2,n,i,j,iVar,iEl) = &
                    (this % solution % interior % hostData(i,n,1,iEl)) !rho*u

                  this % flux % physical % hostData(2,2,n,i,j,iVar,iEl) = &
                    (this % solution % interior % hostData(i,n,2,iEl)) !rho*v
                END DO

              ELSEIF (iVar == 4) THEN ! total energy (rho*u*H)

                DO n = 0,this % solution % interp % N
                  this % flux % physical % hostData(1,1,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(n,j,1,iEl)* &
                     this % diagnostics % interior % hostData(n,j,2,iEl))!rho*u*H

                  this % flux % physical % hostData(2,1,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(n,j,2,iEl)* &
                     this % diagnostics % interior % hostData(n,j,2,iEl)) !rho*v*H

                  this % flux % physical % hostData(1,2,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(i,n,1,iEl)* &
                     this % diagnostics % interior % hostData(i,n,2,iEl))!rho*u*H

                  this % flux % physical % hostData(2,2,n,i,j,iVar,iEl) = &
                    (this % primitive % interior % hostData(i,n,2,iEl)* &
                     this % diagnostics % interior % hostData(i,n,2,iEl)) !rho*v*H
                END DO

              END IF

            END DO
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE ConservativeFlux

  SUBROUTINE FluxMethod_cns2d(this)
    !! This overridden method serves as a wrapper that calls
    !! the HyperbolicInteriorFlux procedure pointer. This design allows
    !! users to dynamically select the type of Riemann solver
    !! to use

    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this

    ! For hyperbolic terms
    CALL this % HyperbolicInteriorFlux()

    ! For parabolic terms
    CALL this % ParabolicInteriorFlux()

  END SUBROUTINE FluxMethod_cns2d

  SUBROUTINE RiemannSolver_cns2d(this)
  !! This overridden method serves as a wrapper that calls
  !! the HyperbolicBoundaryFlux procedure pointer. This design allows
  !! users to dynamically select the type of Riemann solver
  !! to use

    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this

    ! For hyperbolic terms
    CALL this % HyperbolicBoundaryFlux()

    ! For parabolic terms
    CALL this % ParabolicBoundaryFlux() ! Parabolic Boundary Flux

  END SUBROUTINE RiemannSolver_cns2d

  SUBROUTINE LocalLaxFriedrichs(this)
  !! Approximate Riemann Solver for the Compressible Navier-Stokes equations
  !! The Riemann Solver implemented here is the Local Lax-Friedrichs.
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: nhat(1:2),nmag
    REAL(prec) :: cL,cR,unL,unR,HL,HR,HuL,HuR,HvL,HvR
    REAL(prec) :: alpha
    REAL(prec) :: fluxL(1:4)
    REAL(prec) :: fluxR(1:4)
    REAL(prec) :: jump(1:4)

    IF (this % gpuAccel) THEN

      CALL LocalLaxFriedrichs_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
                                          this % solution % boundary % deviceData, &
                                          this % solution % extBoundary % deviceData, &
                                          this % primitive % boundary % deviceData, &
                                          this % primitive % extBoundary % deviceData, &
                                          this % diagnostics % boundary % deviceData, &
                                          this % diagnostics % extBoundary % deviceData, &
                                          this % geometry % nHat % boundary % deviceData, &
                                          this % geometry % nScale % boundary % deviceData, &
                                          this % solution % interp % N, &
                                          this % solution % nVar, &
                                          this % diagnostics % nVar, &
                                          this % solution % nElem)

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iSide = 1,4
          DO i = 0,this % solution % interp % N

            ! Get the boundary normals on cell edges from the mesh geometry
            nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
            nmag = this % geometry % nScale % boundary % hostData(i,1,iSide,iEl)

            ! Calculate the normal velocity at the cell edges
            unL = this % primitive % boundary % hostData(i,1,iSide,iEl)*nHat(1) + &
                  this % primitive % boundary % hostData(i,2,iSide,iEl)*nHat(2)

            unR = this % primitive % extBoundary % hostData(i,1,iSide,iEl)*nHat(1) + &
                  this % primitive % extBoundary % hostData(i,2,iSide,iEl)*nHat(2)

            fluxL(1) = unL*this % solution % boundary % hostData(i,1,iSide,iEl) + &
                       this % primitive % boundary % hostData(i,4,iSide,iEl)*nHat(1)

            fluxL(2) = unL*this % solution % boundary % hostData(i,2,iSide,iEl) + &
                       this % primitive % boundary % hostData(i,4,iSide,iEl)*nHat(2)

            fluxL(3) = this % solution % boundary % hostData(i,1,iSide,iEl)*nHat(1) + &
                       this % solution % boundary % hostData(i,2,iSide,iEl)*nHat(2)

            fluxL(4) = unL*this % diagnostics % boundary % hostData(i,2,iSide,iEl)

            fluxR(1) = unR*this % solution % extBoundary % hostData(i,1,iSide,iEl) + &
                       this % primitive % extBoundary % hostData(i,4,iSide,iEl)*nHat(1)

            fluxR(2) = unR*this % solution % extBoundary % hostData(i,2,iSide,iEl) + &
                       this % primitive % extBoundary % hostData(i,4,iSide,iEl)*nHat(2)

            fluxR(3) = this % solution % extBoundary % hostData(i,1,iSide,iEl)*nHat(1) + &
                       this % solution % extBoundary % hostData(i,2,iSide,iEl)*nHat(2)

            fluxR(4) = unR*this % diagnostics % extBoundary % hostData(i,2,iSide,iEl)

            jump(1:4) = this % solution % boundary % hostData(i,1:4,iSide,iEl) - &
                        this % solution % extBoundary % hostData(i,1:4,iSide,iEl)

            cL = this % diagnostics % boundary % hostData(i,3,iSide,iEl)
            cR = this % diagnostics % extBoundary % hostData(i,3,iSide,iEl)

            alpha = MAX(ABS(unL),ABS(unR)) + MAX(ABS(cL),ABS(cR))

            ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
    this % flux % boundaryNormal % hostData(i,1:4,iSide,iEl) = 0.5_PREC*(fluxL(1:4) + fluxR(1:4) + alpha*jump(1:4))*nmag

          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE LocalLaxFriedrichs

  SUBROUTINE PrimitiveDiffusiveInteriorFlux(this)
    !! This method uses the primitive variables, their gradient,
    !! and the conservative variables (solution) to compute the
    !! diffusive momentum and energy fluxes.
    !!
    !! Thie method assumes that the hyperbolic flux is computed before
    !! calling this method; we add to the boundaryNormal component
    !! of the flux.
    !!
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,i,j,n
    REAL(prec) :: u,v
    REAL(prec) :: dudx,dudy,dvdx,dvdy
    REAL(prec) :: dTdx,dTdy
    REAL(prec) :: tau11,tau12,tau22
    REAL(prec) :: nx,ny,nmag
    REAL(prec) :: viscosity
    REAL(prec) :: diffusivity

    !    IF( this % gpuAccel )THEN

    !    ! TO DO : GPU Accelerated flux
    !    !CALL LocalLaxFriedrichs_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
    !    !       this % solution % boundary % deviceData, &
    !    !       this % solution % extBoundary % deviceData, &
    !    !       this % primitive % boundary % deviceData, &
    !    !       this % primitive % extBoundary % deviceData, &
    !    !       this % diagnostics % boundary % deviceData, &
    !    !       this % diagnostics % extBoundary % deviceData, &
    !    !       this % geometry % nHat % boundary % deviceData, &
    !    !       this % geometry % nScale % boundary % deviceData, &
    !    !       this % solution % interp % N, &
    !    !       this % solution % nVar, &
    !    !       this % diagnostics % nVar, &
    !    !       this % solution % nElem)

    !  ELSE

    ! For now, viscosity and diffusivity are constant values
    ! TO DO : Make spatially varying
    diffusivity = this % diffusivity
    viscosity = this % viscosity
    DO iEl = 1,this % solution % nElem
      DO j = 0,this % solution % interp % N
        DO i = 0,this % solution % interp % N

          DO n = 0,this % solution % interp % N

            u = this % primitive % interior % hostData(n,j,1,iEl)
            v = this % primitive % interior % hostData(n,j,2,iEl)

            ! average x-component of the stress tensor
            dudx = this % primitiveGradient % interior % hostData(1,n,j,1,iEl)
            dudy = this % primitiveGradient % interior % hostData(2,n,j,1,iEl)
            dvdx = this % primitiveGradient % interior % hostData(1,n,j,2,iEl)
            dvdy = this % primitiveGradient % interior % hostData(2,n,j,2,iEl)

            tau11 = 4.0_PREC/3.0_PREC*dudx - 2.0_PREC/3.0_PREC*dvdy
            tau12 = dudy + dvdx ! = tau_21
            tau22 = 4.0_PREC/3.0_PREC*dvdy - 2.0_PREC/3.0_PREC*dudx

            dTdx = this % primitiveGradient % interior % hostData(1,n,j,4,iEl)
            dTdy = this % primitiveGradient % interior % hostData(2,n,j,4,iEl)

            ! x-momentum
            this % flux % physical % hostData(1,1,n,i,j,1,iEl) = &
            this % flux % physical % hostData(1,1,n,i,j,1,iEl) &
            - tau11*viscosity ! x-momentum (x-component)

            this % flux % physical % hostData(2,1,n,i,j,1,iEl) = &
            this % flux % physical % hostData(2,1,n,i,j,1,iEl) &
            - tau12*viscosity ! x-momentum (y-component) 

            ! y-momentum
            this % flux % physical % hostData(1,1,n,i,j,2,iEl) = &
            this % flux % physical % hostData(1,1,n,i,j,2,iEl) &
            - tau12*viscosity ! y-momentum (x-component)

            this % flux % physical % hostData(2,1,n,i,j,2,iEl) = &
            this % flux % physical % hostData(2,1,n,i,j,2,iEl) &
            - tau22*viscosity ! y-momentum (y-component)

            ! Skip density

            ! Total Energy
            this % flux % physical % hostData(1,1,n,i,j,4,iEl) = &
            this % flux % physical % hostData(1,1,n,i,j,4,iEl) &
            - (u*tau11 + v*tau12 + diffusivity*dTdx)*viscosity ! total energy (x-component)

            this % flux % physical % hostData(1,1,n,i,j,4,iEl) = &
            this % flux % physical % hostData(2,1,n,i,j,4,iEl) &
            - (u*tau12 + v*tau22 + diffusivity*dTdy)*viscosity ! total energy (y-component)
              
          END DO

          ! Second part of two-point flux
          DO n = 0,this % solution % interp % N

            u = this % primitive % interior % hostData(i,n,1,iEl)
            v = this % primitive % interior % hostData(i,n,2,iEl)

            ! average x-component of the stress tensor
            dudx = this % primitiveGradient % interior % hostData(1,i,n,1,iEl)
            dudy = this % primitiveGradient % interior % hostData(2,i,n,1,iEl)
            dvdx = this % primitiveGradient % interior % hostData(1,i,n,2,iEl)
            dvdy = this % primitiveGradient % interior % hostData(2,i,n,2,iEl)

            tau11 = 4.0_PREC/3.0_PREC*dudx - 2.0_PREC/3.0_PREC*dvdy
            tau12 = dudy + dvdx ! = tau_21
            tau22 = 4.0_PREC/3.0_PREC*dvdy - 2.0_PREC/3.0_PREC*dudx

            dTdx = this % primitiveGradient % interior % hostData(1,i,n,4,iEl)
            dTdy = this % primitiveGradient % interior % hostData(2,i,n,4,iEl)

            ! x-momentum
            this % flux % physical % hostData(1,2,n,i,j,1,iEl) = &
            this % flux % physical % hostData(1,2,n,i,j,1,iEl) &
            - tau11*viscosity ! x-momentum (x-component)

            this % flux % physical % hostData(2,2,n,i,j,1,iEl) = &
            this % flux % physical % hostData(2,2,n,i,j,1,iEl) &
            - tau12*viscosity ! x-momentum (y-component) 

            ! y-momentum
            this % flux % physical % hostData(1,2,n,i,j,2,iEl) = &
            this % flux % physical % hostData(1,2,n,i,j,2,iEl) &
            - tau12*viscosity ! y-momentum (x-component)

            this % flux % physical % hostData(2,2,n,i,j,2,iEl) = &
            this % flux % physical % hostData(2,2,n,i,j,2,iEl) &
            - tau22*viscosity ! y-momentum (y-component)

            ! Skip density

            ! Total Energy
            this % flux % physical % hostData(1,2,n,i,j,4,iEl) = &
            this % flux % physical % hostData(1,2,n,i,j,4,iEl) &
            - (u*tau11 + v*tau12 + diffusivity*dTdx)*viscosity ! total energy (x-component)

            this % flux % physical % hostData(1,2,n,i,j,4,iEl) = &
            this % flux % physical % hostData(2,2,n,i,j,4,iEl) &
            - (u*tau12 + v*tau22 + diffusivity*dTdy)*viscosity ! total energy (y-component)
              
          END DO


        END DO
      END DO
    END DO

    !   ENDIF

  END SUBROUTINE PrimitiveDiffusiveInteriorFlux

  SUBROUTINE PrimitiveDiffusiveBoundaryFlux(this)
  !! This method uses the primitive variables, their gradient,
  !! and the conservative variables (solution) to compute the
  !! diffusive momentum and energy fluxes.
  !!
  !! Thie method assumes that the hyperbolic flux is computed before
  !! calling this method; we add to the boundaryNormal component
  !! of the flux.
  !!
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl,iSide,i
    REAL(prec) :: u,v
    REAL(prec) :: dudx,dudy,dvdx,dvdy
    REAL(prec) :: dTdx,dTdy
    REAL(prec) :: tau11,tau12,tau22
    REAL(prec) :: nx,ny,nmag
    REAL(prec) :: viscosity
    REAL(prec) :: diffusivity

    !    IF( this % gpuAccel )THEN

    !    ! TO DO : GPU Accelerated flux
    !    !CALL LocalLaxFriedrichs_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
    !    !       this % solution % boundary % deviceData, &
    !    !       this % solution % extBoundary % deviceData, &
    !    !       this % primitive % boundary % deviceData, &
    !    !       this % primitive % extBoundary % deviceData, &
    !    !       this % diagnostics % boundary % deviceData, &
    !    !       this % diagnostics % extBoundary % deviceData, &
    !    !       this % geometry % nHat % boundary % deviceData, &
    !    !       this % geometry % nScale % boundary % deviceData, &
    !    !       this % solution % interp % N, &
    !    !       this % solution % nVar, &
    !    !       this % diagnostics % nVar, &
    !    !       this % solution % nElem)

    !  ELSE

    ! For now, viscosity and diffusivity are constant values
    ! TO DO : Make spatially varying
    diffusivity = this % diffusivity
    viscosity = this % viscosity
    DO iEl = 1,this % solution % nElem
      DO iSide = 1,4
        DO i = 0,this % solution % interp % N

          ! Get the boundary normals on cell edges from the mesh geometry
          nx = this % geometry % nHat % boundary % hostData(1,i,1,iSide,iEl)
          ny = this % geometry % nHat % boundary % hostData(2,i,1,iSide,iEl)
          nmag = this % geometry % nScale % boundary % hostData(i,1,iSide,iEl)

          ! average x-component of the stress tensor
          dudx = 0.5_PREC*(this % primitiveGradient % boundary % hostData(1,i,1,iSide,iEl) + &
                           this % primitiveGradient % extBoundary % hostData(1,i,1,iSide,iEl))

          dudy = 0.5_PREC*(this % primitiveGradient % boundary % hostData(2,i,1,iSide,iEl) + &
                           this % primitiveGradient % boundary % hostData(2,i,1,iSide,iEl))

          dvdx = 0.5_PREC*(this % primitiveGradient % boundary % hostData(1,i,2,iSide,iEl) + &
                           this % primitiveGradient % boundary % hostData(1,i,2,iSide,iEl))

          dvdy = 0.5_PREC*(this % primitiveGradient % boundary % hostData(2,i,2,iSide,iEl) + &
                           this % primitiveGradient % boundary % hostData(2,i,2,iSide,iEl))

          tau11 = 4.0_PREC/3.0_PREC*dudx - 2.0_PREC/3.0_PREC*dvdy
          tau12 = dudy + dvdx ! = tau_21
          tau22 = 4.0_PREC/3.0_PREC*dvdy - 2.0_PREC/3.0_PREC*dudx

          dTdx = 0.5_PREC*(this % primitiveGradient % boundary % hostData(1,i,4,iSide,iEl) + &
                           this % primitiveGradient % boundary % hostData(1,i,4,iSide,iEl))

          dTdy = 0.5_PREC*(this % primitiveGradient % boundary % hostData(2,i,4,iSide,iEl) + &
                           this % primitiveGradient % boundary % hostData(2,i,4,iSide,iEl))

          this % flux % boundaryNormal % hostData(i,1,iSide,iEl) = this % flux % boundaryNormal % hostData(i,1,iSide,iEl) - &
                                                                   viscosity*(tau11*nx + tau12*ny)*nmag ! rhou
          this % flux % boundaryNormal % hostData(i,2,iSide,iEl) = this % flux % boundaryNormal % hostData(i,2,iSide,iEl) - &
                                                                   viscosity*(tau12*nx + tau22*ny)*nmag ! rhov
          !this % flux % boundaryNormal % hostData(i,3,iSide,iEl) =  this % flux % boundaryNormal % hostData(i,3,iSide,iEl)+0.0_prec ! rho
          this % flux % boundaryNormal % hostData(i,4,iSide,iEl) = this % flux % boundaryNormal % hostData(i,4,iSide,iEl) &
                                                                   - (viscosity*((tau11*u + tau12*v)*nx + &
                                                                                 (tau12*u + tau22*v)*ny) + &
                                                                      diffusivity*(dTdx*nx + dTdy*ny))*nmag ! rho*E
        END DO
      END DO
    END DO

!   ENDIF

  END SUBROUTINE PrimitiveDiffusiveBoundaryFlux

  SUBROUTINE Write_cns2d(this,fileName)
#undef __FUNC__
#define __FUNC__ "Write_cns2d"
    IMPLICIT NONE
    CLASS(cns2d),INTENT(inout) :: this
    CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    TYPE(Scalar2D) :: solution
    TYPE(Vector2D) :: x
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

    INFO("Writing pickup file : "//TRIM(pickupFile))

    IF (this % decomp % mpiEnabled) THEN

      CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId,this % decomp % mpiComm)

      ! Write the interpolant to the file
      INFO("Writing interpolant data to file")
      CALL this % solution % interp % WriteHDF5(fileId)

      ! In this section, we write the solution and geometry on the control (quadrature) grid
      ! which can be used for model pickup runs or post-processing
      ! Write the model state to file
      INFO("Writing control grid conservative variables to file")
      CALL CreateGroup_HDF5(fileId,'/controlgrid')
      CALL this % solution % WriteHDF5(fileId,'/controlgrid/solution', &
                                    this % decomp % offsetElem % hostData(this % decomp % rankId),this % decomp % nElem)

      INFO("Writing control grid entropy variables to file")
      CALL this % entropyVars % WriteHDF5(fileId,'/controlgrid/entropy', &
                                    this % decomp % offsetElem % hostData(this % decomp % rankId),this % decomp % nElem)

      INFO("Writing control grid primitive variables to file")
      CALL this % primitive % WriteHDF5(fileId,'/controlgrid/primitive', &
                                    this % decomp % offsetElem % hostData(this % decomp % rankId),this % decomp % nElem)

      ! Write the geometry to file
      INFO("Writing control grid geometry to file")
      CALL CreateGroup_HDF5(fileId,'/controlgrid/geometry')
      CALL this % geometry % x % WriteHDF5(fileId,'/controlgrid/geometry/x', &
                                    this % decomp % offsetElem % hostData(this % decomp % rankId),this % decomp % nElem)

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
      CALL solution % WriteHDF5(fileId,'/targetgrid/solution', &
                                this % decomp % offsetElem % hostData(this % decomp % rankId),this % decomp % nElem)

      ! Map the entropy variables to the target grid
      CALL this % entropyVars % GridInterp(solution,gpuAccel=.FALSE.)
      CALL solution % WriteHDF5(fileId,'/targetgrid/entropy', &
                                this % decomp % offsetElem % hostData(this % decomp % rankId),this % decomp % nElem)

      ! Map the entropy variables to the target grid
      CALL this % primitive % GridInterp(solution,gpuAccel=.FALSE.)
      CALL solution % WriteHDF5(fileId,'/targetgrid/primitive', &
                                this % decomp % offsetElem % hostData(this % decomp % rankId),this % decomp % nElem)

      ! Write the geometry to file
      CALL CreateGroup_HDF5(fileId,'/targetgrid/mesh')
      CALL x % WriteHDF5(fileId,'/targetgrid/mesh/coords', &
                         this % decomp % offsetElem % hostData(this % decomp % rankId),this % decomp % nElem)

      CALL Close_HDF5(fileId)

    ELSE

      CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

      ! Write the interpolant to the file
      INFO("Writing interpolant data to file")
      CALL this % solution % interp % WriteHDF5(fileId)

      ! In this section, we write the solution and geometry on the control (quadrature) grid
      ! which can be used for model pickup runs or post-processing

      ! Write the model state to file
      INFO("Writing control grid conservative to file")
      CALL CreateGroup_HDF5(fileId,'/controlgrid')
      CALL this % solution % WriteHDF5(fileId,'/controlgrid/solution')

      INFO("Writing control grid entropy variables to file")
      CALL this % entropyVars % WriteHDF5(fileId,'/controlgrid/entropy')

      INFO("Writing control grid primitive variables to file")
      CALL this % primitive % WriteHDF5(fileId,'/controlgrid/primitive')

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
      INFO("Writing target grid conservative variables to file")
      CALL CreateGroup_HDF5(fileId,'/targetgrid')
      CALL solution % WriteHDF5(fileId,'/targetgrid/solution')

      INFO("Writing target grid entropy variables to file")
      CALL this % entropyVars % GridInterp(solution,gpuAccel=.FALSE.)
      CALL solution % WriteHDF5(fileId,'/targetgrid/entropy')

      INFO("Writing target grid primitive variables to file")
      CALL this % primitive % GridInterp(solution,gpuAccel=.FALSE.)
      CALL solution % WriteHDF5(fileId,'/targetgrid/primitive')

      ! Write the geometry to file
      INFO("Writing target grid geometry to file")
      CALL CreateGroup_HDF5(fileId,'/targetgrid/geometry')
      CALL x % WriteHDF5(fileId,'/targetgrid/geometry/x')

      CALL Close_HDF5(fileId)

    END IF

    CALL x % Free()
    CALL solution % Free()
    CALL interp % Free()

  END SUBROUTINE Write_cns2d

END MODULE SELF_cns2d
