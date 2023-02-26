!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_CompressibleIdealGas2D

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_ECModel2D

  TYPE,EXTENDS(ECModel2D) :: CompressibleIdealGas2D
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
    PROCEDURE(RiemannFlux_CompressibleIdealGas2D), POINTER :: RiemannFlux => NaiveLLF_CompressibleIdealGas2D

    TYPE(MappedScalar2D) :: primitive ! Contains primitive variables 
    TYPE(MappedScalar2D) :: entropyVars ! Contains entropy variables
    TYPE(MappedScalar2D) :: environmentals ! Functions to describe environmental features (e.g. gravitational potential)
    TYPE(MappedVector2D) :: environmentalsGradient ! Functions to describe environmental features (e.g. gravitational potential)

    TYPE(MappedScalar2D) :: diagnostics

    TYPE(MappedScalar2D) :: prescribedSolution
    TYPE(MappedScalar2D) :: prescribedPrimitive
    TYPE(MappedScalar2D) :: prescribedDiagnostics

    REAL(prec) :: expansionFactor
    REAL(prec) :: Cp ! Heat capacity at constant pressure ( J/g/K )
    REAL(prec) :: Cv ! Heat capacity at constant volume ( J/g/K )
    REAL(prec) :: R ! Ideal gas constant (J/kg/K)

    CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_CompressibleIdealGas2D
    PROCEDURE :: Free => Free_CompressibleIdealGas2D
    PROCEDURE :: PreTendency => PreTendency_CompressibleIdealGas2D
    PROCEDURE :: CalculateEntropy => CalculateEntropy_CompressibleIdealGas2D
    
    PROCEDURE :: WriteTecplot => WriteTecplot_CompressibleIdealGas2D

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_CompressibleIdealGas2D
    PROCEDURE :: FluxMethod => SinglePointFlux_CompressibleIdealGas2D
    PROCEDURE :: RiemannSolver => RiemannSolver_CompressibleIdealGas2D
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_CompressibleIdealGas2D

    ! New Methods
    PROCEDURE :: CheckMinMax => CheckMinMax_CompressibleIdealGas2D
    PROCEDURE :: SetMaxCFL => SetMaxCFL_CompressibleIdealGas2D
    PROCEDURE :: HydrostaticAdjustment => HydrostaticAdjustment_CompressibleIdealGas2D
    
    ! Riemann Fluxes
    PROCEDURE :: NaiveLLF_CompressibleIdealGas2D

    GENERIC :: SetRiemannFlux => SetRiemannFlux_withInt, &
                                    SetRiemannFlux_withChar
    PROCEDURE,PRIVATE :: SetRiemannFlux_withInt
    PROCEDURE,PRIVATE :: SetRiemannFlux_withChar

    GENERIC :: SetVelocity => SetVelocityFromChar_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetVelocityFromChar_CompressibleIdealGas2D

    GENERIC :: SetGravity => SetGravityFromChar_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetGravityFromChar_CompressibleIdealGas2D

    GENERIC :: SetDrag => SetDragFromChar_CompressibleIdealGas2D,&
                              SetDragFromConstant_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetDragFromChar_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetDragFromConstant_CompressibleIdealGas2D

    GENERIC :: SetPrescribedSolution => SetPrescribedSolutionFromChar_CompressibleIdealGas2D,&
                              SetPrescribedSolutionFromSolution_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetPrescribedSolutionFromChar_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetPrescribedSolutionFromSolution_CompressibleIdealGas2D

    PROCEDURE :: SetCp => SetCp_CompressibleIdealGas2D
    PROCEDURE :: SetCv => SetCv_CompressibleIdealGas2D
    PROCEDURE :: SetGasConstant => SetGasConstant_CompressibleIdealGas2D
    PROCEDURE :: SetStatic => SetStatic_CompressibleIdealGas2D
    PROCEDURE :: SetStaticSTP => SetStaticSTP_CompressibleIdealGas2D

    PROCEDURE :: AddThermalBubble

    PROCEDURE,PRIVATE :: CalculateDiagnostics
    PROCEDURE,PRIVATE :: ConservativeToPrimitive
    PROCEDURE,PRIVATE :: ConservativeToEntropy

  END TYPE CompressibleIdealGas2D

  ! ---------------------------------------- !
  ! Riemann Flux integer flags
  !
  INTEGER, PARAMETER :: SELF_NLLF_CIG2D = 500 ! Naive Local Lax Friedrich's Riemann Flux

  ! ---------------------------------------- !
  ! Diagnostics variable indices
  !
  INTEGER, PARAMETER, PRIVATE :: nDiagnostics = 4

  ! ---------------------------------------- !
  ! Variables controlling hydrostatic adjustment
  !
  INTEGER, PARAMETER, PRIVATE :: hydrostaticAdjMaxIters = 1000

  ! ---------------------------------------- ! 
  ! Static fluid state for "air" at stp
  !  > To do : move to json input file
  !
  REAL(prec), PARAMETER :: Cp_stpAir = 1.005_prec
  REAL(prec), PARAMETER :: Cv_stpAir = 0.718_prec
  REAL(prec), PARAMETER :: R_stpAir = 287.0_prec ! J/kg/K
  REAL(prec), PARAMETER :: rho_stpAir = 1.2754_prec ! kg/m^3
  REAL(prec), PARAMETER :: T_stpAir = 273.0_prec ! K
  REAL(prec), PARAMETER :: e_stpAir = 1.5_prec*R_stpAir*T_stpAir ! J/kg

  INTERFACE
    SUBROUTINE RiemannFlux_CompressibleIdealGas2D(this)
      IMPORT CompressibleIdealGas2D
      IMPLICIT NONE
      CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    END SUBROUTINE RiemannFlux_CompressibleIdealGas2D
  END INTERFACE 

  INTERFACE
    SUBROUTINE SetBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper(solution, &
                    extSolution, prescribedSolution, primitive, extPrimitive, &
                    prescribedPrimitive, diagnostics, extDiagnostics, &
                    prescribedDiagnostics, nHat, sideInfo, N, nVar, nDiag, nEl) &
      bind(c,name="SetBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: solution
      TYPE(c_ptr) :: extSolution
      TYPE(c_ptr) :: prescribedSolution
      TYPE(c_ptr) :: primitive
      TYPE(c_ptr) :: extPrimitive
      TYPE(c_ptr) :: prescribedPrimitive
      TYPE(c_ptr) :: diagnostics
      TYPE(c_ptr) :: extDiagnostics
      TYPE(c_ptr) :: prescribedDiagnostics
      TYPE(c_ptr) :: nHat
      TYPE(c_ptr) :: sideInfo
      INTEGER(C_INT),VALUE :: N,nVar,nDiag,nEl
    END SUBROUTINE SetBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Source_CompressibleIdealGas2D_gpu_wrapper(source, solution, environmentalsGradient, N, nVar, nEl) &
      bind(c,name="Source_CompressibleIdealGas2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: source, solution, environmentalsGradient
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE Source_CompressibleIdealGas2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE SinglePointFlux_CompressibleIdealGas2D_gpu_wrapper(flux, solution, primitive, N, nVar, nEl) &
      bind(c,name="SinglePointFlux_CompressibleIdealGas2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, solution, primitive
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SinglePointFlux_CompressibleIdealGas2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE NaiveLLF_CompressibleIdealGas2D_gpu_wrapper(flux, &
                    solution, extSolution, primitive, extPrimitive, diagnostics, &
                    extDiagnostics, nHat, nScale, N, nVar, nDiag, nEl) &
      bind(c,name="NaiveLLF_CompressibleIdealGas2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: flux
      TYPE(c_ptr) :: solution
      TYPE(c_ptr) :: extSolution
      TYPE(c_ptr) :: primitive
      TYPE(c_ptr) :: extPrimitive
      TYPE(c_ptr) :: diagnostics
      TYPE(c_ptr) :: extDiagnostics
      TYPE(c_ptr) :: nHat
      TYPE(c_ptr) :: nScale
      INTEGER(C_INT),VALUE :: N,nVar,nDiag,nEl
    END SUBROUTINE NaiveLLF_CompressibleIdealGas2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE CalculateDiagnostics_CompressibleIdealGas2D_gpu_wrapper(solution, &
                   diagnostics, expansionFactor, R, N, nVar, &
                   nDiag, nEl) &
      bind(c,name="CalculateDiagnostics_CompressibleIdealGas2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: solution
      TYPE(c_ptr) :: diagnostics
      REAL(c_prec),VALUE :: expansionFactor, R
      INTEGER(C_INT),VALUE :: N,nVar,nDiag,nEl
    END SUBROUTINE CalculateDiagnostics_CompressibleIdealGas2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ConservativeToPrimitive_CompressibleIdealGas2D_gpu_wrapper(solution, &
                   primitive, expansionFactor, N, nVar, nEl) &
      bind(c,name="ConservativeToPrimitive_CompressibleIdealGas2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: solution
      TYPE(c_ptr) :: primitive
      REAL(c_prec),VALUE :: expansionFactor
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ConservativeToPrimitive_CompressibleIdealGas2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE ConservativeToEntropy_CompressibleIdealGas2D_gpu_wrapper(solution, &
                   entropy, expansionFactor, N, nVar, nEl) &
      bind(c,name="ConservativeToEntropy_CompressibleIdealGas2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: solution
      TYPE(c_ptr) :: entropy
      REAL(c_prec),VALUE :: expansionFactor
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE ConservativeToEntropy_CompressibleIdealGas2D_gpu_wrapper
  END INTERFACE

CONTAINS

  SUBROUTINE Init_CompressibleIdealGas2D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(out) :: this
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
    this % RiemannFlux => NaiveLLF_CompressibleIdealGas2D
    this % gpuAccel = .FALSE.

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % prescribedSolution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % diagnostics % Init(geometry % x % interp,nDiagnostics,this % mesh % nElem)
    CALL this % prescribedDiagnostics % Init(geometry % x % interp,nDiagnostics,this % mesh % nElem)
    CALL this % environmentals % Init(geometry % x % interp,2,this % mesh % nElem)
    CALL this % environmentalsGradient % Init(geometry % x % interp,2,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
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

  END SUBROUTINE Init_CompressibleIdealGas2D

  SUBROUTINE Free_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this

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
    CALL this % prescribedPrimitive % Free()
    CALL this % entropyVars % Free()

  END SUBROUTINE Free_CompressibleIdealGas2D

  SUBROUTINE SetRiemannFlux_withInt(this,fluxMethod)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    INTEGER, INTENT(in) :: fluxMethod

      SELECT CASE ( fluxMethod )

        CASE ( SELF_NLLF_CIG2D )
          this % RiemannFlux => NaiveLLF_CompressibleIdealGas2D
        CASE DEFAULT
          this % RiemannFlux => NaiveLLF_CompressibleIdealGas2D

      END SELECT

  END SUBROUTINE SetRiemannFlux_withInt

  SUBROUTINE SetRiemannFlux_withChar(this,fluxMethod)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: fluxMethod
    ! Local
    CHARACTER(SELF_INTEGRATOR_LENGTH) :: upperCaseInt

      upperCaseInt = UpperCase(TRIM(fluxMethod))

      SELECT CASE (TRIM(upperCaseInt))

        CASE ( "NAIVELLF" )
          this % RiemannFlux => NaiveLLF_CompressibleIdealGas2D
        CASE DEFAULT
          this % RiemannFlux => NaiveLLF_CompressibleIdealGas2D

      END SELECT

  END SUBROUTINE SetRiemannFlux_withChar

  SUBROUTINE SetStatic_CompressibleIdealGas2D(this)
  !! Sets the default fluid state with uniform 
  !! density and temperature and no motion with
  !! speed of sound as ~ 2 m/s
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      CALL this % SetCv( Cv_stpAir )
      CALL this % SetCp( Cp_stpAir )
      CALL this % SetGasConstant( 1.0_prec )
      
      DO iEl = 1, this % source % nElem
        DO j = 0, this % source % interp % N
          DO i = 0, this % source % interp % N
            this % solution % interior % hostData(i,j,1,iEl) = 0.0_prec ! rho*u
            this % solution % interior % hostData(i,j,2,iEl) = 0.0_prec ! rho*v
            this % solution % interior % hostData(i,j,3,iEl) = rho_stpAir ! rho
            this % solution % interior % hostData(i,j,4,iEl) = rho_stpAir*10.0_prec ! rho*E
          ENDDO
        ENDDO
      ENDDO

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateDevice()
      ENDIF

      CALL this % solution % BoundaryInterp( gpuAccel = this % gpuAccel )
      CALL this % PreTendency( )

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateHost()
        CALL this % primitive % UpdateHost()
        CALL this % diagnostics % UpdateHost()
      ENDIF
      
  END SUBROUTINE SetStatic_CompressibleIdealGas2D

  SUBROUTINE AddThermalBubble(this,xc,R,Tmax)
    !! Adds a temperature anomaly to the fluid state
    !! The anomaly is a gaussian blob of radius R 
    !! centered at xc with an extrema of Tmax.
    !!
    !! The density field is calculated so that the
    !! background pressure field remains undisturbed
    !!
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: xc(1:2)
    REAL(prec), INTENT(in) :: R
    REAL(prec), INTENT(in) :: Tmax
    ! Local
    INTEGER :: i,j,iEl,iVar
    REAL(prec) :: x(1:2)
    REAL(prec) :: Tprime, Rg, T

      Rg = this % R

      DO iEl = 1, this % source % nElem
        DO j = 0, this % source % interp % N
          DO i = 0, this % source % interp % N

            x = this % geometry % x % interior % hostData(1:2,i,j,1,iEl)

            Tprime = Tmax*exp( -( (x(1) - xc(1))**2 + (x(2) - xc(2))**2 )/(2.0_prec*R*R) )

            T = this % diagnostics % interior % hostData(i,j,4,iEl) + Tprime

            ! Update the density
            this % solution % interior % hostData(i,j,3,iEl) = &
               this % solution % interior % hostData(i,j,3,iEl)*(1.0_prec - Tprime/T) 

            ! Add internal energy
            this % solution % interior % hostData(i,j,4,iEl) = &
              1.5_prec*this % solution % interior % hostData(i,j,3,iEl)*Rg*T 

          ENDDO
        ENDDO
      ENDDO

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateDevice()
      ENDIF

      CALL this % solution % BoundaryInterp( gpuAccel = this % gpuAccel )
      CALL this % PreTendency( )

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateHost()
        CALL this % primitive % UpdateHost()
        CALL this % diagnostics % UpdateHost()
      ENDIF
      
  END SUBROUTINE AddThermalBubble

  SUBROUTINE SetStaticSTP_CompressibleIdealGas2D(this)
  !! Sets the default fluid state as "air" at STP with
  !! no motion
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      CALL this % SetCv( Cv_stpAir )
      CALL this % SetCp( Cp_stpAir )
      CALL this % SetGasConstant( R_stpAir )
      
      DO iEl = 1, this % source % nElem
        DO j = 0, this % source % interp % N
          DO i = 0, this % source % interp % N
            this % solution % interior % hostData(i,j,1,iEl) = 0.0_prec ! rho*u
            this % solution % interior % hostData(i,j,2,iEl) = 0.0_prec ! rho*v
            this % solution % interior % hostData(i,j,3,iEl) = rho_stpAir ! rho
            this % solution % interior % hostData(i,j,4,iEl) = rho_stpAir*e_stpAir ! rho*E
          ENDDO
        ENDDO
      ENDDO

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateDevice()
      ENDIF

      CALL this % solution % BoundaryInterp( gpuAccel = this % gpuAccel )
      CALL this % PreTendency( )

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateHost()
        CALL this % primitive % UpdateHost()
        CALL this % diagnostics % UpdateHost()
      ENDIF
      
  END SUBROUTINE SetStaticSTP_CompressibleIdealGas2D

  SUBROUTINE SetVelocityFromChar_CompressibleIdealGas2D(this, eqnChar)
  !! Sets the fluid velocity field using the provided equation parser
  !! The momentum is then updated using the current fluid density field.
  !! From here, the PreTendency method is called to set other diagnostics
  !!
  !! The total energy field is calculated using the internal energy (diagnosed from the
  !! in-situ temperature) and the new kinetic energy field. 
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:2)
    ! Local
    INTEGER :: i,j,iEl,iVar
    REAL(prec) :: rho, u, v, temperature, internalE, KE

      DO iVar = 1, 2
        CALL this % primitive % SetEquation(ivar, eqnChar(iVar))
      ENDDO

      CALL this % primitive % SetInteriorFromEquation( this % geometry, this % t )
      CALL this % primitive % BoundaryInterp( gpuAccel = .FALSE. )

      DO iEl = 1, this % source % nElem
        DO j = 0, this % source % interp % N
          DO i = 0, this % source % interp % N
            rho = this % solution % interior % hostData(i,j,3,iEl)
            u = this % primitive % interior % hostData(i,j,1,iEl)
            v = this % primitive % interior % hostData(i,j,2,iEl)
            temperature = this % diagnostics % interior % hostData(i,j,4,iEl)
            internalE = 1.5_prec*rho*this % R*temperature ! Internal energy
            KE = rho*(u*u+v*v)*0.5_prec

            this % solution % interior % hostData(i,j,1,iEl) = rho*u  ! rho*u
            this % solution % interior % hostData(i,j,2,iEl) = rho*v ! rho*v
            this % solution % interior % hostData(i,j,4,iEl) = internalE + KE ! rho*E
          ENDDO
        ENDDO
      ENDDO

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateDevice()
      ENDIF

      CALL this % solution % BoundaryInterp( gpuAccel = this % gpuAccel )
      CALL this % PreTendency( )

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateHost()
        CALL this % primitive % UpdateHost()
        CALL this % diagnostics % UpdateHost()
      ENDIF
      
  END SUBROUTINE SetVelocityFromChar_CompressibleIdealGas2D

  SUBROUTINE SetGravityFromChar_CompressibleIdealGas2D(this, eqnChar)
    !! Sets the gravitational acceleration from an equation input as a character
    !! The interior points are set and then the strong form of the gradient is used
    !! to calculate the gradient. After, environmentals and environmentalsGradient
    !! fields are copied to device memory (if an accelerator device is present)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    CHARACTER(LEN=*),INTENT(in) :: eqnChar

      CALL this % environmentals % SetEquation(1, eqnChar)

      CALL this % environmentals % SetInteriorFromEquation( this % geometry, this % t )
   
      CALL this % environmentals % GradientSF( this % geometry, &
                                               this % environmentalsGradient, &
                                               .FALSE. )
      
      IF( this % gpuAccel )THEN
        CALL this % environmentals % UpdateDevice()
        CALL this % environmentalsGradient % UpdateDevice()
      ENDIF

  END SUBROUTINE SetGravityFromChar_CompressibleIdealGas2D

  SUBROUTINE SetDragFromChar_CompressibleIdealGas2D(this, eqnChar)
    !! Sets the momentum drag from an equation input as a character
    !! The interior points are set and then copied to the device
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    CHARACTER(LEN=*),INTENT(in) :: eqnChar

      CALL this % environmentals % SetEquation(2, eqnChar)
      CALL this % environmentals % SetInteriorFromEquation( this % geometry, this % t )
      
      IF( this % gpuAccel )THEN
        CALL this % environmentals % UpdateDevice()
      ENDIF
      
  END SUBROUTINE SetDragFromChar_CompressibleIdealGas2D

  SUBROUTINE SetDragFromConstant_CompressibleIdealGas2D(this, Cd)
    !! Sets the momentum drag from an equation input as a character
    !! The interior points are set and then copied to the device
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: Cd
    ! Local
    INTEGER :: i,j,iEl,iVar

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            this % environmentals % interior % hostData(i,j,2,iEl) = Cd

          ENDDO
        ENDDO
      ENDDO

      
      IF( this % gpuAccel )THEN
        CALL this % environmentals % UpdateDevice()
      ENDIF
      
  END SUBROUTINE SetDragFromConstant_CompressibleIdealGas2D

  SUBROUTINE SetPrescribedSolutionFromChar_CompressibleIdealGas2D(this, eqnChar) 
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % prescribedSolution % nVar)
    ! Local
    INTEGER :: iVar

      DO iVar = 1, this % prescribedSolution % nVar
        CALL this % prescribedSolution % SetEquation(ivar, eqnChar(iVar))
      ENDDO

      CALL this % prescribedSolution % SetInteriorFromEquation( this % geometry, this % t )
      CALL this % prescribedSolution % BoundaryInterp( gpuAccel = .FALSE. )

      IF( this % gpuAccel )THEN
        CALL this % prescribedSolution % UpdateDevice()
      ENDIF

  END SUBROUTINE SetPrescribedSolutionFromChar_CompressibleIdealGas2D

  SUBROUTINE SetPrescribedSolutionFromSolution_CompressibleIdealGas2D(this) 
  !! Sets the prescribed solution using the current solution attribute
  !! This can be useful for situations where you want to set the
  !! boundary conditions and initial conditions to be identical
  !!
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nvar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              this % prescribedSolution % interior % hostData(i,j,iVar,iEl) = &
                this % solution % interior % hostData(i,j,iVar,iEl)

              this % prescribedPrimitive % interior % hostData(i,j,iVar,iEl) = &
                this % primitive % interior % hostData(i,j,iVar,iEl)
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      
      DO iEl = 1, this % diagnostics % nElem
        DO iVar = 1, this % diagnostics % nVar
          DO j = 0, this % diagnostics % interp % N
            DO i = 0, this % diagnostics % interp % N

              this % prescribedDiagnostics % interior % hostData(i,j,iVar,iEl) = &
                this % diagnostics % interior % hostData(i,j,iVar,iEl)
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL this % prescribedSolution % BoundaryInterp( gpuAccel = .FALSE. )
      CALL this % prescribedPrimitive % BoundaryInterp( gpuAccel = .FALSE. )
      CALL this % prescribedDiagnostics % BoundaryInterp( gpuAccel = .FALSE. )

      IF( this % gpuAccel )THEN
        CALL this % prescribedSolution % UpdateDevice()
        CALL this % prescribedPrimitive % UpdateDevice()
        CALL this % prescribedDiagnostics % UpdateDevice()
      ENDIF

  END SUBROUTINE SetPrescribedSolutionFromSolution_CompressibleIdealGas2D

  SUBROUTINE SetCp_CompressibleIdealGas2D(this, Cp)
  !! Accessor routine to set the heat capacity at constant pressure
  !! Also updates the expansionFactor attribute when called.
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    REAL(prec), INTENT(in) :: Cp

      this % Cp = Cp
      IF( this % Cv == 0.0_prec )THEN
        PRINT*, "Warning : Expansion factor not set; Cv = 0"
      ENDIF
      this % expansionFactor = this % Cp/this % Cv

  END SUBROUTINE SetCp_CompressibleIdealGas2D

  SUBROUTINE SetCv_CompressibleIdealGas2D(this, Cv)
  !! Accessor routine to set the heat capacity at constant volume
  !! Also updates the expansionFactor attribute when called.
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    REAL(prec), INTENT(in) :: Cv

      this % Cv = Cv
      this % expansionFactor = this % Cp/this % Cv

  END SUBROUTINE SetCv_CompressibleIdealGas2D

  SUBROUTINE SetGasConstant_CompressibleIdealGas2D(this, R)
  !! Accessor routine to set the ideal gas constant
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    REAL(prec), INTENT(in) :: R

      this % R = R

  END SUBROUTINE SetGasConstant_CompressibleIdealGas2D
  
  SUBROUTINE CalculateEntropy_CompressibleIdealGas2D(this)
  !! 
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: Jacobian, wi, wj
    REAL(prec) :: P, rho
    REAL(prec) :: entropy
          
      IF( this % gpuAccel ) THEN
        CALL this % solution % interior % UpdateHost()
        CALL this % diagnostics % interior % UpdateHost()
      ENDIF

      entropy = 0.0_prec

      DO iEl = 1, this % geometry % x % nElem
        DO j = 0, this % geometry % x % interp % N
          DO i = 0, this % geometry % x % interp % N

            ! Coordinate mapping Jacobian
            Jacobian = this % geometry % J % interior % hostData(i,j,1,iEl)

            ! Quadrature weights
            wi = this % geometry % x % interp % qWeights % hostData(i) 
            wj = this % geometry % x % interp % qWeights % hostData(j)
            
            rho = this % solution % interior % hostData(i,j,3,iEl)
            P = this % primitive % interior % hostData(i,j,4,iEl)
            
            entropy = entropy + &
              rho*(log(P) - this % expansionFactor*log(rho))/&
              (this % expansionFactor - 1.0_prec)*wi*wj*Jacobian
          
          ENDDO
        ENDDO
      ENDDO

      CALL this % decomp % GlobalReduce( entropy, this % entropy )

  END SUBROUTINE CalculateEntropy_CompressibleIdealGas2D

  SUBROUTINE CheckMinMax_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    ! Local
    INTEGER :: iVar


    IF( this % gpuAccel )THEN
      CALL this % solution % UpdateHost()
      CALL this % environmentals % UpdateHost()
      CALL this % primitive % UpdateHost()
      CALL this % diagnostics % UpdateHost()
    ENDIF

    PRINT*, '---------------------'
    DO iVar = 1, this % solution % nVar
      PRINT*, TRIM(this % solution % meta(iVar) % name)//" (t, min, max) :", &
              this % t, &
              MINVAL( this % solution % interior % hostData(:,:,iVar,:) ), &
              MAXVAL( this % solution % interior % hostData(:,:,iVar,:) )
    ENDDO

    DO iVar = 1, this % environmentals % nVar
      PRINT*, TRIM(this % environmentals % meta(iVar) % name)//" (t, min, max) :", &
              this % t, &
              MINVAL( this % environmentals % interior % hostData(:,:,iVar,:) ), &
              MAXVAL( this % environmentals % interior % hostData(:,:,iVar,:) )
    ENDDO

    DO iVar = 1, this % primitive % nVar
      PRINT*, TRIM(this % primitive % meta(iVar) % name)//" (t, min, max) :", &
              this % t, &
              MINVAL( this % primitive % interior % hostData(:,:,iVar,:) ), &
              MAXVAL( this % primitive % interior % hostData(:,:,iVar,:) )
    ENDDO

    DO iVar = 1, this % diagnostics % nVar
      PRINT*, TRIM(this % diagnostics % meta(iVar) % name)//" (t, min, max) :", &
              this % t, &
              MINVAL( this % diagnostics % interior % hostData(:,:,iVar,:) ), &
              MAXVAL( this % diagnostics % interior % hostData(:,:,iVar,:) )
    ENDDO

    PRINT*, '---------------------'

  END SUBROUTINE CheckMinMax_CompressibleIdealGas2D

  SUBROUTINE SetMaxCFL_CompressibleIdealGas2D(this, cfl)
    !! This method uses the model grid and sound speed
    !! to set the time step size so that the desired
    !! maximum cfl number fixed.
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    REAL(prec), INTENT(in) :: cfl
    REAL(prec) :: dxMin
    REAL(prec) :: diagMax(nDiagnostics)
    REAL(prec) :: currentDt, currentTime, tn
    REAL(prec) :: Cd
    REAL(prec) :: dsdtMax(this % solution % nVar)
    REAL(prec) :: sMax(this % solution % nVar)

      dxMin = this % geometry % CovariantArcMin()
      diagMax = this % diagnostics % AbsMaxInterior( ) ! Get the absolute max of the diagnostics 

      PRINT*, "Min(dx) : ",dxMin
      PRINT*, "Max(c) : ", diagMax(3)

      ! Reassign the time step for the hydrostatic adjustment
      ! so that the max CFL number is 0.5
      currentDt = this % dt
      this % dt = cfl*dxMin/diagMax(3)
      
      PRINT*, "Adjusted time step size : ", this % dt

  END SUBROUTINE SetMaxCFL_CompressibleIdealGas2D

  SUBROUTINE HydrostaticAdjustment_CompressibleIdealGas2D(this,tolerance)
    !! This method can be used to adjust a compressible fluid
    !! to hydrostatic equilibrium. On input, the CompressibleIdealGas2D
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
    !!  Cd = \frac{max(c)}{min(\Delta x)}
    !!
    !! In this subroutine, the time step is chosen so that the sound-wave
    !! maximum CFL number is 0.75
    !!
    !!   CFL = \frac{max(c) \Delta t}{min(\Delta x)} = 0.75
    !! 
    !!    \rightarrow \Delta t = 0.75\frac{min( \Delta x )}{\max{c}}
    !!
    !! After adjustment, the 
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tolerance
    ! Local
    INTEGER :: i 
    REAL(prec) :: currentDt, currentTime, tn
    REAL(prec) :: Cd
    REAL(prec) :: dsdtMax(this % solution % nVar)
    REAL(prec) :: sMax(this % solution % nVar)
    REAL(prec) :: dxMin
    REAL(prec) :: diagMax(nDiagnostics)

      dxMin = this % geometry % CovariantArcMin()
      diagMax = this % diagnostics % AbsMaxInterior( ) ! Get the absolute max of the diagnostics 

      ! Reassign the time step for the hydrostatic adjustment
      ! so that the max CFL number is 0.5
      currentDt = this % dt
      CALL this % SetMaxCFL( 0.5_prec )
      
      ! Calculate the drag coefficient
      Cd = 0.3*diagMax(3)/dxMin

      PRINT*, "Drag coefficient : ", Cd

      CALL this % SetDrag( Cd )

      ! Save the current time
      currentTime = this % t

      tn = 1000.0_prec*this % dt
      DO i = 1, hydrostaticAdjMaxIters

        this % t = 0.0_prec
        CALL this % ForwardStep( tn )
        
        sMax = this % solution % AbsMaxInterior()

        PRINT*, "Momentum Max : ", SQRT(sMax(1)**2 + sMax(2)**2)
         
        ! Check if momentum max is small in comparison to tolerance
        IF( SQRT(sMax(1)**2 + sMax(2)**2) <= tolerance )THEN
          PRINT*, "Reached hydrostatic equilibrium in ", i, " iterations"
          EXIT
        ENDIF
      ENDDO

      ! Reset delta t and time
      this % dt = currentDt
      this % t = currentTime

  END SUBROUTINE HydrostaticAdjustment_CompressibleIdealGas2D

  SUBROUTINE PreTendency_CompressibleIdealGas2D(this)
    !! Calculate the velocity and density weighted enthalpy at element interior and element boundaries
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this

      CALL this % CalculateDiagnostics()
      CALL this % ConservativeToPrimitive()
      CALL this % ConservativeToEntropy()

      ! Interpolate velocity and required diagnostics to the element boundaries
      CALL this % primitive % BoundaryInterp(this % gpuAccel)
      CALL this % diagnostics % BoundaryInterp(this % gpuAccel)
      CALL this % entropyVars % BoundaryInterp(this % gpuAccel)

      ! Perform any MPI exchanges for the velocity and the required diagnostics
      ! across shared element faces between neighboring processes.
      CALL this % primitive % SideExchange(this % mesh, this % decomp, this % gpuAccel)
      CALL this % diagnostics % SideExchange(this % mesh, this % decomp, this % gpuAccel)
      CALL this % entropyVars % SideExchange(this % mesh, this % decomp, this % gpuAccel)

  END SUBROUTINE PreTendency_CompressibleIdealGas2D

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
    !!     p = (\gamma-1)*\rho*e
    !! 
    !!  where $\gamma = \frac{C_p}{C_v}$ is the expansion coefficient,
    !!  $\rho$ is the fluid density, and $e$ is the internal energy.
    !! 
    !!  We calculate $rho*e$ as
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
    !!    \frac{\partial P}{\partial \rho} = c^{2}
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
    !!   (\gamma-1)*e = p/\rho
    !! 
    !!  We calculate $e$ as
    !! 
    !!    e = (\rho*E - 0.5_prec*\rho*KE)/\rho
    !!    
    !! 
    !!  where rho*E is the total energy, a prognostic variable, modelled
    !!  by this class, $\rho*KE$ is the kinetic energy (a required diagnostic)
    !!  and, $\rho$ is the density ( a prognostic variable ).
    !!
    !! In-Situ Temperature
    !!
    !!    T = 2/3*e/R
    !!

    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rho, rhoU, rhoV, u, v, rhoE, rhoKE, p

      IF( this % gpuAccel )THEN

        CALL CalculateDiagnostics_CompressibleIdealGas2D_gpu_wrapper( & 
          this % solution % interior % deviceData, &
          this % diagnostics % interior % deviceData, &
          this % expansionFactor, &
          this % R, &
          this % solution % interp % N, &
          this % solution % nVar, &
          this % diagnostics % nVar, &
          this % solution % nElem )

      ELSE

        DO iEl = 1, this % solution % nElem
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              rhoU = this % solution % interior % hostData(i,j,1,iEl)
              rhoV = this % solution % interior % hostData(i,j,2,iEl)
              rho = this % solution % interior % hostData(i,j,3,iEl)
              rhoE = this % solution % interior % hostData(i,j,4,iEl)
              u = rhoU/rho
              v = rhoV/rho
              rhoKE = 0.5_prec*(rhoU*u+rhoV*v)
              p = (this % expansionFactor - 1.0_prec)*(rhoE - rhoKE)

              ! Kinetic Energy
              this % diagnostics % interior % hostData(i,j,1,iEl) = rhoKE

              ! Enthalpy
              this % diagnostics % interior % hostData(i,j,2,iEl) = rhoE + p

              ! Speed of sound
              this % diagnostics % interior % hostData(i,j,3,iEl) = &
                      sqrt(this % expansionFactor*p/rho)

              ! In-Situ Temperature
              this % diagnostics % interior % hostData(i,j,4,iEl) = &
                      (2.0_prec/3.0_prec)*((rhoE - rhoKE)/rho)/this % R

            ENDDO
          ENDDO
        ENDDO

      ENDIF

  END SUBROUTINE CalculateDiagnostics

  SUBROUTINE ConservativeToPrimitive(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rhoU, rhoV, rho, u, v, rhoE, rhoKE, p

      IF( this % gpuAccel )THEN

        CALL ConservativeToPrimitive_CompressibleIdealGas2D_gpu_wrapper( & 
          this % solution % interior % deviceData, &
          this % primitive % interior % deviceData, &
          this % expansionFactor, &
          this % solution % interp % N, &
          this % solution % nVar, &
          this % solution % nElem )

      ELSE
        DO iEl = 1, this % solution % nElem
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              rhoU = this % solution % interior % hostData(i,j,1,iEl)
              rhoV = this % solution % interior % hostData(i,j,2,iEl)
              rho = this % solution % interior % hostData(i,j,3,iEl)
              rhoE = this % solution % interior % hostData(i,j,4,iEl)
              u = rhoU/rho
              v = rhoV/rho
              rhoKE = 0.5_prec*(rhoU*u+rhoV*v)
              p = (this % expansionFactor - 1.0_prec)*(rhoE - rhoKE)

              this % primitive % interior % hostData(i,j,1,iEl) = u 
              this % primitive % interior % hostData(i,j,2,iEl) = v
              this % primitive % interior % hostData(i,j,3,iEl) = rho 
              this % primitive % interior % hostData(i,j,4,iEl) = p

            ENDDO
          ENDDO
        ENDDO
      ENDIF

  END SUBROUTINE ConservativeToPrimitive

  SUBROUTINE ConservativeToEntropy(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rhoU, rhoV, rho, u, v, E, KE, p, s

      IF( this % gpuAccel )THEN

        CALL ConservativeToEntropy_CompressibleIdealGas2D_gpu_wrapper( & 
          this % solution % interior % deviceData, &
          this % entropyVars % interior % deviceData, &
          this % expansionFactor, &
          this % solution % interp % N, &
          this % solution % nVar, &
          this % solution % nElem )

      ELSE
        DO iEl = 1, this % solution % nElem
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              rhoU = this % solution % interior % hostData(i,j,1,iEl)
              rhoV = this % solution % interior % hostData(i,j,2,iEl)
              rho = this % solution % interior % hostData(i,j,3,iEl)
              E = this % solution % interior % hostData(i,j,4,iEl)
              u = rhoU/rho
              v = rhoV/rho
              KE = 0.5_prec*(rhoU*u+rhoV*v)
              p = (this % expansionFactor - 1.0_prec)*(E - KE)
              s = log(p) - this % expansionFactor*log(rho)

              this % entropyVars % interior % hostData(i,j,1,iEl) = u*rho/p 
              this % entropyVars % interior % hostData(i,j,2,iEl) = v*rho/p
              this % entropyVars % interior % hostData(i,j,3,iEl) = (this % expansionFactor - s)/&
                     (this % expansionFactor - 1.0_prec) - KE/p
              this % entropyVars % interior % hostData(i,j,4,iEl) = -rho/p

            ENDDO
          ENDDO
        ENDDO
      ENDIF

  END SUBROUTINE ConservativeToEntropy

  SUBROUTINE SetBoundaryCondition_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, i
    INTEGER :: bcid, e2
    REAL(prec) :: u, v, nhat(1:2)


      IF( this % gpuAccel )THEN

        CALL SetBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper(&
          this % solution % boundary % deviceData, &
          this % solution % extBoundary % deviceData, &
          this % prescribedSolution % boundary % deviceData, &
          this % primitive % boundary % deviceData, &
          this % primitive % extBoundary % deviceData, &
          this % prescribedPrimitive % boundary % deviceData, &
          this % diagnostics % boundary % deviceData, &
          this % diagnostics % extBoundary % deviceData, &
          this % prescribedDiagnostics % boundary % deviceData, &
          this % geometry % nHat % boundary % deviceData, &
          this % mesh % sideInfo % deviceData, &
          this % solution % interp % N, &
          this % solution % nVar, &
          this % diagnostics % nVar, &
          this % solution % nElem )

      ELSE

        DO iEl = 1, this % solution % nElem
          DO iSide = 1, 4
            bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
            e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
            IF( e2 == 0 )THEN
              IF( bcid == SELF_BC_NONORMALFLOW )THEN

                DO i = 0, this % solution % interp % N
                  nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
                  ! Conservative variables
                  u = this % solution % boundary % hostData(i,1,iSide,iEl) 
                  v = this % solution % boundary % hostData(i,2,iSide,iEl) 
                  this % solution % extBoundary % hostData(i,1,iSide,iEl) = &
                          (nhat(2)**2 - nhat(1)**2)*u - 2.0_prec*nhat(1)*nhat(2)*v
                  this % solution % extBoundary % hostData(i,2,iSide,iEl) = &
                          (nhat(1)**2 - nhat(2)**2)*v - 2.0_prec*nhat(1)*nhat(2)*u
                  this % solution % extBoundary % hostData(i,3,iSide,iEl) = &
                          this % solution % boundary % hostData(i,3,iSide,iEl)
                  this % solution % extBoundary % hostData(i,4,iSide,iEl) = &
                          this % solution % boundary % hostData(i,4,iSide,iEl)
                  
                  ! Primitive variables
                  u = this % primitive % boundary % hostData(i,1,iSide,iEl) 
                  v = this % primitive % boundary % hostData(i,2,iSide,iEl) 
                  this % primitive % extBoundary % hostData(i,1,iSide,iEl) = &
                          (nhat(2)**2 - nhat(1)**2)*u - 2.0_prec*nhat(1)*nhat(2)*v
                  this % primitive % extBoundary % hostData(i,2,iSide,iEl) = &
                          (nhat(1)**2 - nhat(2)**2)*v - 2.0_prec*nhat(1)*nhat(2)*u
                  this % primitive % extBoundary % hostData(i,3,iSide,iEl) = &
                          this % primitive % boundary % hostData(i,3,iSide,iEl)
                  this % primitive % extBoundary % hostData(i,4,iSide,iEl) = &
                          this % primitive % boundary % hostData(i,4,iSide,iEl)
                  
                  ! Prolong the diagnostic values to the external state
                  this % diagnostics % extBoundary % hostData(i,1:nDiagnostics,iSide,iEl) = &
                    this % diagnostics % boundary % hostData(i,1:nDiagnostics,iSide,iEl)

                ENDDO
             
              ELSEIF( bcid == SELF_BC_PRESCRIBED .OR. bcid == SELF_BC_RADIATION )THEN

                DO i = 0, this % solution % interp % N

                  this % solution % extBoundary % hostData(i,1,iSide,iEl) = &
                          this % prescribedSolution % boundary % hostData(i,1,iSide,iEl)
                  this % solution % extBoundary % hostData(i,2,iSide,iEl) = &
                          this % prescribedSolution % boundary % hostData(i,2,iSide,iEl)
                  this % solution % extBoundary % hostData(i,3,iSide,iEl) = &
                          this % prescribedSolution % boundary % hostData(i,3,iSide,iEl)
                  this % solution % extBoundary % hostData(i,4,iSide,iEl) = &
                          this % prescribedSolution % boundary % hostData(i,4,iSide,iEl)
                      

                  this % primitive % extBoundary % hostData(i,1,iSide,iEl) = &
                          this % prescribedPrimitive % boundary % hostData(i,1,iSide,iEl)
                  this % primitive % extBoundary % hostData(i,2,iSide,iEl) = &
                          this % prescribedPrimitive % boundary % hostData(i,2,iSide,iEl)
                  this % primitive % extBoundary % hostData(i,3,iSide,iEl) = &
                          this % prescribedPrimitive % boundary % hostData(i,3,iSide,iEl)
                  this % primitive % extBoundary % hostData(i,4,iSide,iEl) = &
                          this % prescribedPrimitive % boundary % hostData(i,4,iSide,iEl)

                  this % diagnostics % extBoundary % hostData(i,1:nDiagnostics,iSide,iEl) = &
                    this % prescribedDiagnostics % boundary % hostData(i,1:nDiagnostics,iSide,iEl)

                ENDDO

              ENDIF

            ENDIF

          ENDDO
        ENDDO
      ENDIF


  END SUBROUTINE SetBoundaryCondition_CompressibleIdealGas2D 

  SUBROUTINE Source_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar
    REAL(prec) :: rhou, rhov, rho, gx, gy, Cd

    IF( this % gpuAccel )THEN

      CALL Source_CompressibleIdealGas2D_gpu_wrapper( this % source % interior % deviceData, &
              this % solution % interior % deviceData, &
              this % environmentalsGradient % interior % deviceData, &
              this % source % interp % N, &
              this % source % nVar, &
              this % source % nElem )

    ELSE

      DO iEl = 1, this % source % nElem
          DO j = 0, this % source % interp % N
            DO i = 0, this % source % interp % N


              rhou = this % solution % interior % hostData(i,j,1,iEl)
              rhov = this % solution % interior % hostData(i,j,2,iEl)
              rho = this % solution % interior % hostData(i,j,3,iEl)
              gx = this % environmentalsGradient % interior % hostData(1,i,j,1,iEl)
              gy = this % environmentalsGradient % interior % hostData(2,i,j,1,iEl) 
              Cd = this % environmentals % interior % hostData(i,j,2,iEl)

              this % source % interior % hostData(i,j,1,iEl) = -rho*gx -Cd*rhou ! (\rho u)_t = -\rho gx 
              this % source % interior % hostData(i,j,2,iEl) = -rho*gy -Cd*rhov! (\rho v)_t = -\rho gy 
              this % source % interior % hostData(i,j,3,iEl) = 0.0_prec ! (\rho )_t = 0 
              this % source % interior % hostData(i,j,4,iEl) = -rhou*gx-rhou*gy ! (\rho E )_t = -\rho u g_x - \rho u g_

            ENDDO
          ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE Source_CompressibleIdealGas2D

  SUBROUTINE SinglePointFlux_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,n,iEl,iVar

      IF (this % gpuAccel) THEN
        CALL SinglePointFlux_CompressibleIdealGas2D_gpu_wrapper( this % flux % interior % deviceData,&
                this % solution % interior % deviceData, &
                this % primitive % interior % deviceData, &
                this % solution % interp % N, &
                this % solution % nVar, &
                this % solution % nElem )

      ELSE

        DO iEl = 1, this % solution % nElem
          DO iVar = 1, this % solution % nVar
            DO j = 0, this % solution % interp % N
              DO i = 0, this % solution % interp % N

                IF ( iVar == 1 )THEN ! rho*u

                  DO n = 0, this % solution % interp % N
                    ! rho*u*u + p 
                    this % flux % physical % hostData(1,1,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(n,j,1,iEl)*& ! u
                          this % solution % interior % hostData(n,j,1,iEl)+& ! rho*u
                          this % primitive % interior % hostData(n,j,4,iEl) )

                    ! rho*u*v
                    this % flux % physical % hostData(2,1,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(n,j,2,iEl)*& ! v
                          this % solution % interior % hostData(n,j,1,iEl) ) ! rho*u

                    ! rho*u*u + p 
                    this % flux % physical % hostData(1,2,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(i,n,1,iEl)*& ! u
                          this % solution % interior % hostData(i,n,1,iEl)+& ! rho*u
                          this % primitive % interior % hostData(i,n,4,iEl) )

                    ! rho*u*v
                    this % flux % physical % hostData(2,2,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(i,n,2,iEl)*& ! v
                          this % solution % interior % hostData(i,n,1,iEl) ) ! rho*u
                  ENDDO

                ELSEIF ( iVar == 2 )THEN ! rho*v

                  DO n = 0, this % solution % interp % N
                    ! rho*v*u
                    this % flux % physical % hostData(1,1,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(n,j,1,iEl)*& ! u
                          this % solution % interior % hostData(n,j,2,iEl) ) ! rho*v

                    ! rho*v*v + p
                    this % flux % physical % hostData(2,1,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(n,j,2,iEl)*& ! v
                          this % solution % interior % hostData(n,j,2,iEl)+& ! rho*v
                          this % primitive % interior % hostData(n,j,4,iEl) ) ! pressure

                    ! rho*v*u
                    this % flux % physical % hostData(1,2,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(i,n,1,iEl)*& ! u
                          this % solution % interior % hostData(i,n,2,iEl) ) ! rho*v

                    ! rho*v*v + p
                    this % flux % physical % hostData(2,2,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(i,n,2,iEl)*& ! v
                          this % solution % interior % hostData(i,n,2,iEl)+& ! rho*v
                          this % primitive % interior % hostData(i,n,4,iEl) ) ! pressure

                  ENDDO


                ELSEIF ( iVar == 3 )THEN ! density
                        
                  DO n = 0, this % solution % interp % N
                    this % flux % physical % hostData(1,1,n,i,j,iVar,iEl) = &
                          ( this % solution % interior % hostData(n,j,1,iEl) ) !rho*u

                    this % flux % physical % hostData(2,1,n,i,j,iVar,iEl) = &
                          ( this % solution % interior % hostData(n,j,2,iEl) ) !rho*v

                    this % flux % physical % hostData(1,2,n,i,j,iVar,iEl) = &
                          ( this % solution % interior % hostData(i,n,1,iEl) ) !rho*u

                    this % flux % physical % hostData(2,2,n,i,j,iVar,iEl) = &
                          ( this % solution % interior % hostData(i,n,2,iEl) ) !rho*v
                  ENDDO

                ELSEIF ( iVar == 4 )THEN ! total energy (rho*u*H)

                  DO n = 0, this % solution % interp % N
                    this % flux % physical % hostData(1,1,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(n,j,1,iEl)*&
                          this % diagnostics % interior % hostData(n,j,2,iEl) )!rho*u*H

                    this % flux % physical % hostData(2,1,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(n,j,2,iEl)*&
                          this % diagnostics % interior % hostData(n,j,2,iEl) ) !rho*v*H

                    this % flux % physical % hostData(1,2,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(i,n,1,iEl)*&
                          this % diagnostics % interior % hostData(i,n,2,iEl) )!rho*u*H

                    this % flux % physical % hostData(2,2,n,i,j,iVar,iEl) = &
                          ( this % primitive % interior % hostData(i,n,2,iEl)*&
                          this % diagnostics % interior % hostData(i,n,2,iEl) ) !rho*v*H
                  ENDDO

                ENDIF

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

  END SUBROUTINE SinglePointFlux_CompressibleIdealGas2D

  SUBROUTINE RiemannSolver_CompressibleIdealGas2D(this)
  !! This overridden method serves as a wrapper that calls
  !! the RiemannFlux procedure pointer. This design allows
  !! users to dynamically select the type of Riemann solver
  !! to use

    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this

      CALL this % RiemannFlux()

  END SUBROUTINE RiemannSolver_CompressibleIdealGas2D

  SUBROUTINE NaiveLLF_CompressibleIdealGas2D(this)
  !! Approximate Riemann Solver for the Compressible Navier-Stokes equations
  !! The Riemann Solver implemented here is the Local Lax-Friedrichs.
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: nhat(1:2), nmag
    REAL(prec) :: cL, cR, unL, unR, HL, HR, HuL, HuR, HvL, HvR
    REAL(prec) :: alpha
    REAL(prec) :: fluxL(1:4)
    REAL(prec) :: fluxR(1:4)
    REAL(prec) :: jump(1:4)


    IF( this % gpuAccel )THEN

      CALL NaiveLLF_CompressibleIdealGas2D_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
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

      DO iEl = 1, this % solution % nElem
        DO iSide = 1, 4
          DO i = 0, this % solution % interp % N

             ! Get the boundary normals on cell edges from the mesh geometry
             nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
             nmag = this % geometry % nScale % boundary % hostData(i,1,iSide,iEl)

             ! Calculate the normal velocity at the cell edges
             unL = this % primitive % boundary % hostData(i,1,iSide,iEl)*nHat(1)+&
                   this % primitive % boundary % hostData(i,2,iSide,iEl)*nHat(2)

             unR = this % primitive % extBoundary % hostData(i,1,iSide,iEl)*nHat(1)+&
                   this % primitive % extBoundary % hostData(i,2,iSide,iEl)*nHat(2)

             fluxL(1) = unL*this % solution % boundary % hostData(i,1,iSide,iEl) +&
                        this % primitive % boundary % hostData(i,4,iSide,iEl)*nHat(1)

             fluxL(2) = unL*this % solution % boundary % hostData(i,2,iSide,iEl) +&
                        this % primitive % boundary % hostData(i,4,iSide,iEl)*nHat(2)

             fluxL(3) = this % solution % boundary % hostData(i,1,iSide,iEl)*nHat(1)+&
                        this % solution % boundary % hostData(i,2,iSide,iEl)*nHat(2)
                        
             fluxL(4) = unL*this % diagnostics % boundary % hostData(i,2,iSide,iEl)

             fluxR(1) = unR*this % solution % extBoundary % hostData(i,1,iSide,iEl) +&
                       this % primitive % extBoundary % hostData(i,4,iSide,iEl)*nHat(1)

             fluxR(2) = unR*this % solution % extBoundary % hostData(i,2,iSide,iEl) +&
                        this % primitive % extBoundary % hostData(i,4,iSide,iEl)*nHat(2)

             fluxR(3) = this % solution % extBoundary % hostData(i,1,iSide,iEl)*nHat(1)+&
                        this % solution % extBoundary % hostData(i,2,iSide,iEl)*nHat(2)
                        
             fluxR(4) = unR*this % diagnostics % extBoundary % hostData(i,2,iSide,iEl)

             jump(1:4) = this % solution % boundary % hostData(i,1:4,iSide,iEl)-&
                         this % solution % extBoundary % hostData(i,1:4,iSide,iEl)

             cL = this % diagnostics % boundary % hostData(i,3,iSide,iEl)
             cR = this % diagnostics % extBoundary % hostData(i,3,iSide,iEl)

             alpha = MAX( ABS(unL), ABS(unR) ) + MAX( ABS(cL), ABS(cR) )

             ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
             this % flux % boundaryNormal % hostData(i,1:4,iSide,iEl) =  0.5_prec*( fluxL(1:4) + fluxR(1:4) + alpha*jump(1:4) )*nmag

          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE NaiveLLF_CompressibleIdealGas2D

  SUBROUTINE WriteTecplot_CompressibleIdealGas2D(this, filename)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    CHARACTER(*), INTENT(in), OPTIONAL :: filename
    ! Local
    CHARACTER(8) :: zoneID
    INTEGER :: fUnit
    INTEGER :: iEl, i, j, iVar 
    CHARACTER(LEN=self_FileNameLength) :: tecFile
    CHARACTER(LEN=self_TecplotHeaderLength) :: tecHeader
    CHARACTER(LEN=self_FormatLength) :: fmat
    CHARACTER(13) :: timeStampString
    CHARACTER(5) :: rankString

    IF( PRESENT(filename) )THEN
      tecFile = filename
    ELSE
      timeStampString = TimeStamp(this % t, 's')

      IF( this % decomp % mpiEnabled )THEN
        WRITE(rankString,'(I5.5)') this % decomp % rankId 
        tecFile = 'solution.'//rankString//'.'//timeStampString//'.tec'
      ELSE
        tecFile = 'solution.'//timeStampString//'.tec'
      ENDIF

    ENDIF
                      
    IF( this % gpuAccel )THEN
      CALL this % solution % interior % UpdateHost()
    ENDIF

    CALL this % solution % WriteTecplot( this % geometry, &
                                         this % decomp, &
                                         tecFile )
 
    IF( this % decomp % mpiEnabled )THEN
      WRITE(rankString,'(I5.5)') this % decomp % rankId 
      tecFile = 'diagnostics.'//rankString//'.'//timeStampString//'.tec'
    ELSE
      tecFile = 'diagnostics.'//timeStampString//'.tec'
    ENDIF
                      
    IF( this % gpuAccel )THEN
      CALL this % diagnostics % interior % UpdateHost()
    ENDIF

    CALL this % diagnostics % WriteTecplot( this % geometry, &
                                         this % decomp, &
                                         tecFile )

    IF( this % decomp % mpiEnabled )THEN
      WRITE(rankString,'(I5.5)') this % decomp % rankId 
      tecFile = 'primitive.'//rankString//'.'//timeStampString//'.tec'
    ELSE
      tecFile = 'primitive.'//timeStampString//'.tec'
    ENDIF
                      
    IF( this % gpuAccel )THEN
      CALL this % primitive % interior % UpdateHost()
    ENDIF

    CALL this % primitive % WriteTecplot( this % geometry, &
                                         this % decomp, &
                                         tecFile )

  END SUBROUTINE WriteTecplot_CompressibleIdealGas2D

END MODULE SELF_CompressibleIdealGas2D
