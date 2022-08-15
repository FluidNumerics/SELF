!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Model

  USE SELF_SupportRoutines
  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_HDF5
  USE HDF5
  USE FEQParse

! //////////////////////////////////////////////// !
!   Time integration parameters

!     Runge-Kutta 3rd Order, low storage constants
  REAL(prec),PARAMETER,PRIVATE :: rk3_a(1:3) = (/0.0_prec,-5.0_prec/9.0_prec,-153.0_prec/128.0_prec/)
  REAL(prec),PARAMETER,PRIVATE :: rk3_b(1:3) = (/0.0_prec,1.0_prec/3.0_prec,3.0_prec/4.0_prec/)
  REAL(prec),PARAMETER,PRIVATE :: rk3_g(1:3) = (/1.0_prec/3.0_prec,15.0_prec/16.0_prec,8.0_prec/15.0_prec/)

! 
  INTEGER, PARAMETER :: SELF_EULER = 100
  INTEGER, PARAMETER :: SELF_RK2 = 200
  INTEGER, PARAMETER :: SELF_RK3 = 300
  INTEGER, PARAMETER :: SELF_RK4 = 400
  INTEGER, PARAMETER :: SELF_AB2 = 201
  INTEGER, PARAMETER :: SELF_AB3 = 301
  INTEGER, PARAMETER :: SELF_AB4 = 401

  INTEGER, PARAMETER :: SELF_INTEGRATOR_LENGTH = 10 ! max length of integrator methods when specified as char
  INTEGER, PARAMETER :: SELF_EQUATION_LENGTH = 500

! //////////////////////////////////////////////// !
!   Boundary Condition parameters
!

  ! Conditions on the solution
  INTEGER, PARAMETER :: SELF_BC_PRESCRIBED = 100
  INTEGER, PARAMETER :: SELF_BC_RADIATION = 101
  INTEGER, PARAMETER :: SELF_BC_NONORMALFLOW = 102

  ! Conditions on the solution gradients
  INTEGER, PARAMETER :: SELF_BC_PRESCRIBED_STRESS = 200
  INTEGER, PARAMETER :: SELF_BC_NOSTRESS = 201

! //////////////////////////////////////////////// !
!   Model Formulations
!
  INTEGER, PARAMETER :: SELF_FORMULATION_LENGTH = 30 ! max length of integrator methods when specified as char
  INTEGER, PARAMETER :: SELF_CONSERVATIVE_FLUX = 0
  INTEGER, PARAMETER :: SELF_SPLITFORM_FLUX = 1


  TYPE,ABSTRACT :: Model
    LOGICAL :: gpuAccel
    INTEGER :: fluxDivMethod

    ! Time integration attributes
    PROCEDURE(SELF_timeIntegrator), POINTER :: timeStepper => Euler_timeStepper 
    REAL(prec) :: dt
    REAL(prec) :: t

    ! Standard Diagnostics
    REAL(prec) :: entropy ! Mathematical entropy function for the model

    ! Domain Decomposition
    TYPE(MPILayer),POINTER :: decomp

    CONTAINS

    PROCEDURE :: ForwardStep => ForwardStep_Model

!    PROCEDURE :: Null_timeStepper 

    PROCEDURE :: Euler_timeStepper 

    ! Runge-Kutta methods
    PROCEDURE :: LowStorageRK3_timeStepper 
    PROCEDURE(UpdateGRK3),DEFERRED :: UpdateGRK3

    PROCEDURE :: PreTendency => PreTendency_Model
    PROCEDURE :: SourceMethod => Source_Model
    PROCEDURE :: FluxMethod => Flux_Model
    PROCEDURE :: RiemannSolver => RiemannSolver_Model
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_Model

    PROCEDURE :: ReportEntropy => ReportEntropy_Model
    PROCEDURE :: CalculateEntropy => CalculateEntropy_Model

    PROCEDURE(UpdateSolution),DEFERRED :: UpdateSolution
    PROCEDURE(CalculateTendency),DEFERRED :: CalculateTendency
    PROCEDURE(ReadModel),DEFERRED :: ReadModel
    PROCEDURE(WriteModel),DEFERRED :: WriteModel
    PROCEDURE(WriteTecplot),DEFERRED :: WriteTecplot

    GENERIC :: SetTimeIntegrator => SetTimeIntegrator_withInt, &
                                    SetTimeIntegrator_withChar
    PROCEDURE,PRIVATE :: SetTimeIntegrator_withInt
    PROCEDURE,PRIVATE :: SetTimeIntegrator_withChar

    PROCEDURE :: SetSimulationTime
    PROCEDURE :: GetSimulationTime
!    PROCEDURE :: SetTimeStep
!    PROCEDURE :: GetTimeStep

    GENERIC :: SetFluxMethod => SetFluxMethod_withInt, &
                                    SetFluxMethod_withChar
    PROCEDURE,PRIVATE :: SetFluxMethod_withInt
    PROCEDURE,PRIVATE :: SetFluxMethod_withChar

    PROCEDURE :: EnableGPUAccel => EnableGPUAccel_Model
    PROCEDURE :: DisableGPUAccel => DisableGPUAccel_Model

  END TYPE Model

  TYPE,EXTENDS(Model),ABSTRACT :: Model1D
    TYPE(MappedScalar1D) :: solution
    TYPE(MappedScalar1D) :: solutionGradient
    TYPE(MappedScalar1D) :: velocity
    TYPE(MappedScalar1D) :: flux
    TYPE(MappedScalar1D) :: source
    TYPE(MappedScalar1D) :: fluxDivergence
    TYPE(MappedScalar1D) :: dSdt
    TYPE(MappedScalar1D) :: workSol
    TYPE(Mesh1D),POINTER :: mesh
    TYPE(Geometry1D),POINTER :: geometry

    CONTAINS

    PROCEDURE :: Init => Init_Model1D
    PROCEDURE :: Free => Free_Model1D

    PROCEDURE :: UpdateHost => UpdateHost_Model1D
    PROCEDURE :: UpdateDevice => UpdateDevice_Model1D

    PROCEDURE :: UpdateSolution => UpdateSolution_Model1D
    PROCEDURE :: UpdateGRK3 => UpdateGRK3_Model1D
    PROCEDURE :: CalculateTendency => CalculateTendency_Model1D
    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_Model1D

    GENERIC :: SetSolution => SetSolutionFromChar_Model1D,&
                              SetSolutionFromEqn_Model1D
    PROCEDURE,PRIVATE :: SetSolutionFromChar_Model1D
    PROCEDURE,PRIVATE :: SetSolutionFromEqn_Model1D

    GENERIC :: SetVelocityField => SetVelocityFieldFromChar_Model1D,&
                              SetVelocityFieldFromEqn_Model1D
    PROCEDURE,PRIVATE :: SetVelocityFieldFromChar_Model1D
    PROCEDURE,PRIVATE :: SetVelocityFieldFromEqn_Model1D

!    PROCEDURE :: ReprojectFlux => ReprojectFlux_Model1D

    PROCEDURE :: ReadModel => Read_Model1D
    PROCEDURE :: WriteModel => Write_Model1D
    PROCEDURE :: WriteTecplot => WriteTecplot_Model1D

  END TYPE Model1D

  TYPE,EXTENDS(Model) :: Model2D
    TYPE(MappedScalar2D) :: solution
    TYPE(MappedVector2D) :: solutionGradient
    TYPE(MappedVector2D) :: velocity
    TYPE(MappedVector2D) :: compVelocity
    TYPE(MappedVector2D) :: flux
    TYPE(MappedScalar2D) :: source
    TYPE(MappedScalar2D) :: fluxDivergence
    TYPE(MappedScalar2D) :: dSdt
    TYPE(MappedScalar2D) :: workSol
    TYPE(Mesh2D),POINTER :: mesh
    TYPE(SEMQuad),POINTER :: geometry

    CONTAINS

    PROCEDURE :: Init => Init_Model2D
    PROCEDURE :: Free => Free_Model2D

    PROCEDURE :: UpdateHost => UpdateHost_Model2D
    PROCEDURE :: UpdateDevice => UpdateDevice_Model2D

    PROCEDURE :: UpdateSolution => UpdateSolution_Model2D
    PROCEDURE :: UpdateGRK3 => UpdateGRK3_Model2D
    PROCEDURE :: CalculateTendency => CalculateTendency_Model2D
    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_Model2D

    GENERIC :: SetSolution => SetSolutionFromChar_Model2D,&
                              SetSolutionFromEqn_Model2D
    PROCEDURE,PRIVATE :: SetSolutionFromChar_Model2D
    PROCEDURE,PRIVATE :: SetSolutionFromEqn_Model2D

    GENERIC :: SetVelocityField => SetVelocityFieldFromChar_Model2D,&
                              SetVelocityFieldFromEqn_Model2D
    PROCEDURE,PRIVATE :: SetVelocityFieldFromChar_Model2D
    PROCEDURE,PRIVATE :: SetVelocityFieldFromEqn_Model2D

    PROCEDURE :: ReprojectFlux => ReprojectFlux_Model2D

    PROCEDURE :: ReadModel => Read_Model2D
    PROCEDURE :: WriteModel => Write_Model2D
    PROCEDURE :: WriteTecplot => WriteTecplot_Model2D

  END TYPE Model2D

  INTERFACE
    SUBROUTINE SELF_timeIntegrator(this,tn)
      USE SELF_Constants, ONLY : prec
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model),INTENT(inout) :: this
      REAL(prec), INTENT(in) :: tn
    END SUBROUTINE SELF_timeIntegrator
  END INTERFACE 

  INTERFACE 
    SUBROUTINE UpdateGRK3( this, m )
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
    END SUBROUTINE UpdateGRK3
  END INTERFACE

  INTERFACE 
    SUBROUTINE UpdateSolution( this, dt )
      USE SELF_Constants, ONLY : prec
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model),INTENT(inout) :: this
      REAL(prec),OPTIONAL,INTENT(in) :: dt
    END SUBROUTINE UpdateSolution
  END INTERFACE

  INTERFACE 
    SUBROUTINE CalculateTendency( this )
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model),INTENT(inout) :: this
    END SUBROUTINE CalculateTendency
  END INTERFACE

  INTERFACE
    SUBROUTINE WriteModel(this, filename)
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model), INTENT(inout) :: this
      CHARACTER(*), INTENT(in), OPTIONAL :: filename
    END SUBROUTINE WriteModel
  END INTERFACE

  INTERFACE
    SUBROUTINE ReadModel(this, filename)
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model), INTENT(inout) :: this
      CHARACTER(*), INTENT(in) :: filename
    END SUBROUTINE ReadModel
  END INTERFACE


  INTERFACE
    SUBROUTINE WriteTecplot(this, filename)
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model), INTENT(inout) :: this
      CHARACTER(*), INTENT(in), OPTIONAL :: filename
    END SUBROUTINE WriteTecplot
  END INTERFACE


  INTERFACE
    SUBROUTINE UpdateSolution_Model1D_gpu_wrapper(solution, dSdt, dt, N, nVar, nEl) &
      bind(c,name="UpdateSolution_Model1D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: solution, dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: dt
    END SUBROUTINE UpdateSolution_Model1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE UpdateGRK3_Model1D_gpu_wrapper(grk3, solution, dSdt, rk3_a, rk3_g, dt, N, nVar, nEl) &
      bind(c,name="UpdateGRK3_Model1D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: grk3, solution, dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: rk3_a, rk3_g, dt
    END SUBROUTINE UpdateGRK3_Model1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE UpdateSolution_Model2D_gpu_wrapper(solution, dSdt, dt, N, nVar, nEl) &
      bind(c,name="UpdateSolution_Model2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: solution, dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: dt
    END SUBROUTINE UpdateSolution_Model2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE UpdateGRK3_Model2D_gpu_wrapper(grk3, solution, dSdt, rk3_a, rk3_g, dt, N, nVar, nEl) &
      bind(c,name="UpdateGRK3_Model2D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: grk3, solution, dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: rk3_a, rk3_g, dt
    END SUBROUTINE UpdateGRK3_Model2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE CalculateDSDt_Model1D_gpu_wrapper(fluxDivergence, source, dSdt, N, nVar, nEl) &
      bind(c,name="CalculateDSDt_Model1D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: fluxDivergence, source, dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE CalculateDSDt_Model1D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE CalculateDSDt_Model2D_gpu_wrapper(fluxDivergence, source, dSdt, N, nVar, nEl) &
      bind(c,name="CalculateDSDt_Model2D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: fluxDivergence, source, dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE CalculateDSDt_Model2D_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE CalculateDSDt_Model3D_gpu_wrapper(fluxDivergence, source, dSdt, N, nVar, nEl) &
      bind(c,name="CalculateDSDt_Model3D_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: fluxDivergence, source, dSdt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE CalculateDSDt_Model3D_gpu_wrapper
  END INTERFACE

CONTAINS

  SUBROUTINE PreTendency_Model(this)
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this

      RETURN

  END SUBROUTINE PreTendency_Model

  SUBROUTINE Source_Model(this)
    !!
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this

      RETURN

  END SUBROUTINE Source_Model

  SUBROUTINE RiemannSolver_Model(this)
    !!
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this

      RETURN

  END SUBROUTINE RiemannSolver_Model

  SUBROUTINE Flux_Model(this)
    !!
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this

      RETURN

  END SUBROUTINE Flux_Model

  SUBROUTINE SetBoundaryCondition_Model(this)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this

      RETURN

  END SUBROUTINE SetBoundaryCondition_Model 
  
  SUBROUTINE SetTimeIntegrator_withInt(this,integrator)
    !! Sets the time integrator method, using an integer flag
    !!
    !! Valid options for `integrator` are
    !!
    !!    SELF_EULER
    !!    SELF_RK3
    !!    SELF_RK4
    !!
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    INTEGER, INTENT(in) :: integrator

      SELECT CASE ( integrator )

        CASE ( SELF_EULER )
          this % timeStepper => Euler_timeStepper
        CASE ( SELF_RK3 )
          this % timeStepper => LowStorageRK3_timeStepper
        CASE DEFAULT
          this % timeStepper => LowStorageRK3_timeStepper

      END SELECT


  END SUBROUTINE SetTimeIntegrator_withInt

  SUBROUTINE SetTimeIntegrator_withChar(this,integrator)
    !! Sets the time integrator method, using a character input
    !!
    !! Valid options for integrator are
    !!
    !!   "euler"
    !!   "rk3"
    !!   "rk4"
    !!
    !! Note that the character provided is not case-sensitive
    !!
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: integrator
    ! Local
    CHARACTER(SELF_INTEGRATOR_LENGTH) :: upperCaseInt

      upperCaseInt = UpperCase(TRIM(integrator))

      SELECT CASE (TRIM(upperCaseInt))

        CASE ("EULER")
          this % timeStepper => Euler_timeStepper

        CASE ("RK3")
          this % timeStepper => LowStorageRK3_timeStepper

        CASE DEFAULT
          this % timeStepper => LowStorageRK3_timeStepper


      END SELECT

  END SUBROUTINE SetTimeIntegrator_withChar

  SUBROUTINE SetFluxMethod_withInt(this,fluxDivMethod)
    !! Sets the method for calculating the flux divergence, using an integer flag
    !!
    !! Valid options for `fluxDivMethod` are
    !!
    !!    SELF_CONSERVATIVE_FLUX
    !!    SELF_SPLITFORM_FLUX
    !!
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    INTEGER, INTENT(in) :: fluxDivMethod

      this % fluxDivMethod = fluxDivMethod

  END SUBROUTINE SetFluxMethod_withInt

  SUBROUTINE SetFluxMethod_withChar(this,fluxDivMethod)
    !! Sets the method for calculating the flux divergence, using a character input
    !!
    !! Valid options for flux method are
    !!
    !!   "conservative"
    !!   "split" or "splitform" or "split form" or "split-form"
    !!
    !! Note that the character provided is not case-sensitive
    !!
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: fluxDivMethod
    ! Local
    CHARACTER(SELF_FORMULATION_LENGTH) :: upperCaseInt

      upperCaseInt = UpperCase(TRIM(fluxDivMethod))

      SELECT CASE (TRIM(upperCaseInt))

        CASE ("CONSERVATIVE")
          this % fluxDivMethod = SELF_CONSERVATIVE_FLUX

        CASE ("SPLIT")
          this % fluxDivMethod = SELF_SPLITFORM_FLUX

        CASE ("SPLITFORM")
          this % fluxDivMethod = SELF_SPLITFORM_FLUX

        CASE ("SPLIT FORM")
          this % fluxDivMethod = SELF_SPLITFORM_FLUX

        CASE ("SPLIT-FORM")
          this % fluxDivMethod = SELF_SPLITFORM_FLUX

        CASE DEFAULT
          this % fluxDivMethod = SELF_CONSERVATIVE_FLUX

      END SELECT

  END SUBROUTINE SetFluxMethod_withChar

  SUBROUTINE GetSimulationTime(this,t)
    !! Returns the current simulation time stored in the model % t attribute
    IMPLICIT NONE
    CLASS(Model),INTENT(in) :: this
    REAL(prec),INTENT(out) :: t

      t = this % t

  END SUBROUTINE GetSimulationTime

  SUBROUTINE SetSimulationTime(this,t)
    !! Sets the model % t attribute with the provided simulation time
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec),INTENT(in) :: t

      this % t = t

  END SUBROUTINE SetSimulationTime

  SUBROUTINE EnableGPUAccel_Model(this)
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this

    IF (GPUAvailable()) THEN
      this % gpuAccel = .TRUE.
    ELSE
      this % gpuAccel = .FALSE.
      PRINT*, 'Warning : GPU acceleration requested, but no GPU is available'
    ENDIF

  END SUBROUTINE EnableGPUAccel_Model

  SUBROUTINE DisableGPUAccel_Model(this)
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this

    this % gpuAccel = .FALSE.

  END SUBROUTINE DisableGPUAccel_Model

  SUBROUTINE CalculateEntropy_Model(this)
  !! Base method for calculating entropy of a model
  !! When this method is not overridden, the entropy
  !! is simply set to 0.0. When you develop a model
  !! built on top of this abstract class or one of its
  !! children, it is recommended that you define a
  !! convex mathematical entropy function that is used
  !! as a measure of the model stability.
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this

      this % entropy = 0.0_prec

  END SUBROUTINE CalculateEntropy_Model

  SUBROUTINE ReportEntropy_Model(this)
  !! Base method for reporting the entropy of a model
  !! to stdout. Only override this procedure if additional
  !! reporting is needed. Alternatively, if you think
  !! additional reporting would be valuable for all models,
  !! open a pull request with modifications to this base 
  !! method.
    USE, INTRINSIC :: ISO_FORTRAN_ENV
    IMPLICIT NONE
    CLASS(Model), INTENT(in) :: this
    ! Local
    INTEGER, PARAMETER :: ucs2 = selected_char_kind('ISO_10646')
    CHARACTER(KIND=ucs2, len=20) :: modelTime
    CHARACTER(KIND=ucs2, len=20) :: entropy
    CHARACTER(KIND=ucs2, len=:), ALLOCATABLE :: str

    IF( this % decomp % rankId == 0 )THEN
      ! Copy the time and entropy to a string
      WRITE(modelTime,"(ES16.7E3)") this % t
      WRITE(entropy,"(ES16.7E3)") this % entropy
  
      ! Write the output to STDOUT 
      OPEN(output_unit, ENCODING='utf-8')
      str = ucs2_'t\u1D62 ='//TRIM(modelTime)
      WRITE(output_unit,'(A)',ADVANCE='no') str
      str = ucs2_'  |  e\u1D62 ='//TRIM(entropy)
      WRITE(output_unit,'(A)',ADVANCE='yes') str
    ENDIF

  END SUBROUTINE ReportEntropy_Model

  ! ////////////////////////////////////// !
  !       Time Integrators                 !

  SUBROUTINE ForwardStep_Model(this,tn,dt,ioInterval)
  !!  Forward steps the model using the associated tendency procedure and time integrator
  !!
  !!  If the final time `tn` is provided, the model is forward stepped to that final time,
  !!  otherwise, the model is forward stepped only a single time step
  !!  
  !!  If a time step is provided through the interface, the model time step size is updated
  !!  and that time step is used to update the model
  !!
  !! If ioInterval is provided, file IO will be conducted every ioInterval seconds until tn 
  !! is reached
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec),OPTIONAL,INTENT(in) :: tn
    REAL(prec),OPTIONAL,INTENT(in) :: dt
    REAL(prec),OPTIONAL,INTENT(in) :: ioInterval
    ! Local
    REAL(prec) :: targetTime, tNext
    INTEGER :: i, nIO

    IF (PRESENT(dt)) THEN
      this % dt = dt
    ENDIF

    IF (PRESENT(tn)) THEN
      targetTime = tn
    ELSE
      targetTime = this % t + this % dt
    ENDIF

    IF (PRESENT(ioInterval)) THEN
      nIO = INT( (targetTime - this % t)/ioInterval )
      DO i = 1, nIO
        tNext = this % t + ioInterval
        CALL this % timeStepper(tNext)
        this % t = tNext
        CALL this % WriteModel()
        CALL this % WriteTecplot()
        CALL this % CalculateEntropy()
        CALL this % ReportEntropy()
      ENDDO

    ELSE
      CALL this % timeStepper(targetTime)
      this % t = targetTime
      CALL this % CalculateEntropy()
      CALL this % ReportEntropy()
    ENDIF

  END SUBROUTINE ForwardStep_Model

!  SUBROUTINE Null_timeStepper(this,tn)
!    IMPLICIT NONE
!    CLASS(Model),INTENT(inout) :: this
!    REAL(prec), INTENT(in) :: tn
!
!      this % t = tn
!
!  END SUBROUTINE Null_timeStepper

  SUBROUTINE Euler_timeStepper(this,tn)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tn
    ! Local
    REAL(prec) :: tRemain
    REAL(prec) :: dtLim

    dtLim = this % dt ! Get the max time step size from the dt attribute
    DO WHILE (this % t < tn)

      tRemain = tn - this % t
      this % dt = MIN( dtLim, tRemain )
      CALL this % CalculateTendency()
      CALL this % UpdateSolution()
      this % t = this % t + this % dt

    ENDDO 

    this % dt = dtLim

  END SUBROUTINE Euler_timeStepper

  SUBROUTINE LowStorageRK3_timeStepper(this,tn)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tn
    ! Local
    INTEGER :: m
    REAL(prec) :: tRemain
    REAL(prec) :: dtLim

    dtLim = this % dt ! Get the max time step size from the dt attribute
    DO WHILE (this % t < tn)

      tRemain = tn - this % t
      this % dt = MIN( dtLim, tRemain )
      DO m = 1, 3
        CALL this % CalculateTendency()
        CALL this % UpdateGRK3(m)
        this % t = this % t + rk3_b(m)*this % dt
      ENDDO

    ENDDO 

    this % dt = dtLim

  END SUBROUTINE LowStorageRK3_timeStepper

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
    CALL this % velocity % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_Model1D

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

  SUBROUTINE SetSolutionFromEqn_Model1D(this, eqn) 
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

      ! Copy the equation parser
      DO iVar = 1, this % solution % nVar
        CALL this % solution % SetEquation(ivar, eqn(iVar) % equation)
      ENDDO

      CALL this % solution % SetInteriorFromEquation( this % geometry, this % t )
      CALL this % solution % BoundaryInterp( gpuAccel = .FALSE. )


      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateDevice()
      ENDIF



  END SUBROUTINE SetSolutionFromEqn_Model1D 

  SUBROUTINE SetSolutionFromChar_Model1D(this, eqnChar) 
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

      DO iVar = 1, this % solution % nVar
        CALL this % solution % SetEquation(ivar, eqnChar(iVar))
      ENDDO

      CALL this % solution % SetInteriorFromEquation( this % geometry, this % t )
      CALL this % solution % BoundaryInterp( gpuAccel = .FALSE. )


      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateDevice()
      ENDIF

  END SUBROUTINE SetSolutionFromChar_Model1D

  SUBROUTINE SetVelocityFieldFromEqn_Model1D(this, eqn) 
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn

      ! Copy the equation parser
      ! Set the x-component of the velocity
      CALL this % velocity % SetEquation(1,eqn % equation)

      ! Set the velocity values using the equation parser
      CALL this % velocity % SetInteriorFromEquation( this % geometry, this % t )

      CALL this % velocity % BoundaryInterp( gpuAccel = .FALSE. )

      IF( this % gpuAccel )THEN
        CALL this % velocity % UpdateDevice()
      ENDIF

  END SUBROUTINE SetVelocityFieldFromEqn_Model1D 

  SUBROUTINE SetVelocityFieldFromChar_Model1D(this, eqnChar) 
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar

      ! Set the x-component of the velocity
      CALL this % velocity % SetEquation(1,eqnChar)

      ! Set the velocity values using the equation parser
      CALL this % velocity % SetInteriorFromEquation( this % geometry, this % t )

      CALL this % velocity % BoundaryInterp( gpuAccel = .FALSE. )

      IF( this % gpuAccel )THEN
        CALL this % velocity % UpdateDevice()
      ENDIF

  END SUBROUTINE SetVelocityFieldFromChar_Model1D

  SUBROUTINE UpdateSolution_Model1D(this,dt)
    !! Computes a solution update as `s=s+dt*dsdt`, where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    REAL(prec),OPTIONAL,INTENT(in) :: dt
    ! Local
    REAL(prec) :: dtLoc
    INTEGER :: i, iVar, iEl

    IF (PRESENT(dt)) THEN
      dtLoc = dt
    ELSE 
      dtLoc = this % dt
    ENDIF

    IF (this % gpuAccel) THEN

      CALL UpdateSolution_Model1D_gpu_wrapper( this % solution % interior % deviceData, &
                                      this % dSdt % interior % deviceData, &
                                      dtLoc, &
                                      this % solution % interp % N, &
                                      this % solution % nVar, &
                                      this % solution % nElem ) 
                                      

    ELSE

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO i = 0, this % solution % interp % N

            this % solution % interior % hostData(i,iVar,iEl) = &
                this % solution % interior % hostData(i,iVar,iEl) +&
                dtLoc*this % dSdt % interior % hostData(i,iVar,iEl)

          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE UpdateSolution_Model1D

  SUBROUTINE UpdateGRK3_Model1D(this,m)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    INTEGER, INTENT(in) :: m
    ! Local
    INTEGER :: i, iVar, iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGRK3_Model1D_gpu_wrapper( this % workSol % interior % deviceData, &
                                           this % solution % interior % deviceData, &
                                           this % dSdt % interior % deviceData, &
                                           rk3_a(m),rk3_g(m),this % dt, &
                                           this % solution % interp % N, &
                                           this % solution % nVar, &
                                           this % solution % nElem ) 
                                      

    ELSE

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO i = 0, this % solution % interp % N

            this % workSol % interior % hostData(i,iVar,iEl) = rk3_a(m)*&
                   this % workSol % interior % hostData(i,iVar,iEl) + &
                   this % dSdt % interior % hostData(i,iVar,iEl)


            this % solution % interior % hostData(i,iVar,iEl) = &
                    this % solution % interior % hostData(i,iVar,iEl) + &
                    rk3_g(m)*this % dt*this % workSol % interior % hostData(i,iVar,iEl)

          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE UpdateGRK3_Model1D

  SUBROUTINE CalculateFluxDivergence_Model1D(this)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this

    CALL this % flux % Derivative(this % geometry, &
                                  this % fluxDivergence, &
                                  selfWeakDGForm,&
                                  this % gpuAccel)


  END SUBROUTINE CalculateFluxDivergence_Model1D

  SUBROUTINE CalculateTendency_Model1D(this)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i, iVar, iEl

    CALL this % PreTendency()
    CALL this % solution % BoundaryInterp(this % gpuAccel)
    CALL this % solution % SideExchange(this % mesh, this % decomp, this % gpuAccel)
    CALL this % SetBoundaryCondition()
    CALL this % SourceMethod()
    CALL this % RiemannSolver()
    CALL this % FluxMethod()
    CALL this % CalculateFluxDivergence()

    IF( this % gpuAccel )THEN

      CALL CalculateDSDt_Model1D_gpu_wrapper( this % fluxDivergence % interior % deviceData, &
                                      this % source % interior % deviceData, &
                                      this % dSdt % interior % deviceData, &
                                      this % solution % interp % N, &
                                      this % solution % nVar, &
                                      this % solution % nElem ) 
                                      
    ELSE

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO i = 0, this % solution % interp % N

            this % dSdt % interior % hostData(i,iVar,iEl) = &
                    this % source % interior % hostData(i,iVar,iEl) -&
                    this % fluxDivergence % interior % hostData(i,iVar,iEl)

          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE CalculateTendency_Model1D

  SUBROUTINE Write_Model1D(this,fileName)
    IMPLICIT NONE
    CLASS(Model1D),INTENT(inout) :: this
    CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:3)
    INTEGER(HID_T) :: xOffset(1:3)
    INTEGER(HID_T) :: bOffset(1:3)
    INTEGER(HID_T) :: bxOffset(1:3)
    INTEGER(HID_T) :: solGlobalDims(1:3)
    INTEGER(HID_T) :: xGlobalDims(1:3)
    INTEGER(HID_T) :: bGlobalDims(1:3)
    INTEGER(HID_T) :: bxGlobalDims(1:3)
    INTEGER :: firstElem
    ! Local
    CHARACTER(LEN=self_FileNameLength) :: pickupFile
    CHARACTER(13) :: timeStampString

    IF( PRESENT(filename) )THEN
      pickupFile = filename
    ELSE
      timeStampString = TimeStamp(this % t, 's')
      pickupFile = 'solution.'//timeStampString//'.h5'
    ENDIF

    IF( this % gpuAccel ) THEN
      CALL this % solution % UpdateHost()
      CALL this % solutionGradient % UpdateHost()
    ENDIF

    IF (this % decomp % mpiEnabled) THEN

      CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId,this % decomp % mpiComm)

      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)
      solOffset(1:3) = (/0,1,firstElem/)
      solGlobalDims(1:3) = (/this % solution % interp % N, &
                             this % solution % nVar, &
                             this % decomp % nElem/)


      xOffset(1:3) = (/0,1,firstElem/)
      xGlobalDims(1:3) = (/this % solution % interp % N, &
                           this % solution % nVar, &
                           this % decomp % nElem/)

      ! Offsets and dimensions for element boundary data
      bOffset(1:3) = (/1,1,firstElem/)
      bGlobalDims(1:3) = (/this % solution % nVar, &
                           2,&
                           this % decomp % nElem/)

      bxOffset(1:3) = (/1,1,firstElem/)
      bxGlobalDims(1:3) = (/this % solution % nVar, &
                           2,&
                           this % decomp % nElem/)

      
      CALL CreateGroup_HDF5(fileId,'/quadrature')

      IF( this % decomp % rankId == 0 )THEN
        CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
                             this % solution % interp % controlPoints)

        CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
                             this % solution % interp % qWeights)

        CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
                             this % solution % interp % dgMatrix)

        CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
                             this % solution % interp % dMatrix)
      ENDIF

      CALL CreateGroup_HDF5(fileId,'/state')

      CALL CreateGroup_HDF5(fileId,'/state/interior')

      CALL CreateGroup_HDF5(fileId,'/state/boundary')

      CALL CreateGroup_HDF5(fileId,'/mesh')

      CALL CreateGroup_HDF5(fileId,'/mesh/interior')

      CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

      CALL WriteArray_HDF5(fileId,'/state/interior/solution', &
                           this % solution % interior,solOffset,solGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solution', &
                           this % solution % boundary,bOffset,bGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/fluxDivergence', &
                           this % fluxDivergence % interior,solOffset,solGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/flux', &
                           this % flux % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/flux', &
                           this % flux % boundary,bxOffset,bxGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/solutionGradient', &
                           this % solutionGradient % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solutionGradient', &
                           this % solutionGradient % boundary,bxOffset,bxGlobalDims)

      CALL WriteArray_HDF5(fileId,'/mesh/interior/x', &
                           this % geometry % x % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/mesh/boundary/x', &
                           this % geometry % x % boundary,bxOffset,bxGlobalDims)

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

      CALL CreateGroup_HDF5(fileId,'/state')

      CALL CreateGroup_HDF5(fileId,'/state/interior')

      CALL CreateGroup_HDF5(fileId,'/state/boundary')

      CALL CreateGroup_HDF5(fileId,'/mesh')

      CALL CreateGroup_HDF5(fileId,'/mesh/interior')

      CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

      CALL WriteArray_HDF5(fileId,'/state/interior/solution',this % solution % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solution',this % solution % boundary)

      CALL WriteArray_HDF5(fileId,'/state/interior/fluxDivergence',this % fluxDivergence % interior)

      CALL WriteArray_HDF5(fileId,'/state/interior/flux',this % flux % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/flux',this % flux % boundary)

      CALL WriteArray_HDF5(fileId,'/state/interior/solutionGradient',this % solutionGradient % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solutionGradient',this % solutionGradient % boundary)

      CALL WriteArray_HDF5(fileId,'/mesh/interior/x',this % geometry % x % interior)

      CALL WriteArray_HDF5(fileId,'/mesh/boundary/x',this % geometry % x % boundary)

      CALL Close_HDF5(fileId)

    END IF

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

    CALL ReadAttribute_HDF5(fileId,'N',N)

    IF (this % solution % interp % N /= N) THEN
      STOP 'Error : Solution polynomial degree does not match input file'
    END IF

    IF (this % decomp % mpiEnabled) THEN
      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId) + 1
      solOffset(1:3) = (/0,1,firstElem/)
      CALL ReadArray_HDF5(fileId,'/state/interior/solution', &
                          this % solution % interior,solOffset)
    ELSE
      CALL ReadArray_HDF5(fileId,'/state/interior/solution',this % solution % interior)
    END IF

    CALL Close_HDF5(fileId)
    IF( this % gpuAccel )THEN
      CALL this % solution % interior % UpdateDevice()
    ENDIF

  END SUBROUTINE Read_Model1D

  SUBROUTINE WriteTecplot_Model1D(this, filename)
    IMPLICIT NONE
    CLASS(Model1D), INTENT(inout) :: this
    CHARACTER(*), INTENT(in), OPTIONAL :: filename
    ! Local
    CHARACTER(8) :: zoneID
    INTEGER :: fUnit
    INTEGER :: iEl, i, iVar
    CHARACTER(LEN=self_FileNameLength) :: tecFile
    CHARACTER(LEN=self_TecplotHeaderLength) :: tecHeader
    CHARACTER(LEN=self_FormatLength) :: fmat
    CHARACTER(13) :: timeStampString
    CHARACTER(5) :: rankString
    TYPE(Scalar1D) :: solution
    TYPE(Scalar1D) :: x
    TYPE(Lagrange),TARGET :: interp

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
                      
    ! Create an interpolant for the uniform grid
    CALL interp % Init(this % solution % interp % M,&
            this % solution % interp % targetNodeType,&
            this % solution % interp % N, &
            this % solution % interp % controlNodeType)

    CALL solution % Init( interp, &
            this % solution % nVar, this % solution % nElem )

    CALL x % Init( interp, 1, this % solution % nElem )

    ! Map the mesh positions to the target grid
    CALL this % geometry % x % GridInterp(x, gpuAccel=.FALSE.)

    ! Map the solution to the target grid
    CALL this % solution % GridInterp(solution,gpuAccel=.FALSE.)
   
    ! Let's write some tecplot!! 
     OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(tecFile), &
      FORM='formatted', &
      STATUS='replace')

    tecHeader = 'VARIABLES = "X"'
    DO iVar = 1, this % solution % nVar
      tecHeader = TRIM(tecHeader)//', "'//TRIM(this % solution % meta(iVar) % name)//'"'
    ENDDO

    WRITE(fUnit,*) TRIM(tecHeader) 

    ! Create format statement
    WRITE(fmat,*) this % solution % nvar+1
    fmat = '('//TRIM(fmat)//'(ES16.7E3,1x))'

    WRITE(fUnit,*) TRIM(tecHeader) 

    DO iEl = 1, this % solution % nElem

      ! TO DO :: Get the global element ID 
      WRITE(zoneID,'(I8.8)') iEl
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',this % solution % interp % M+1

      DO i = 0, this % solution % interp % M

        WRITE(fUnit,fmat) x % interior % hostData(i,1,iEl), &
                          solution % interior % hostData(i,1:this % solution % nvar,iEl)

      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)

    CALL x % Free()
    CALL solution % Free() 
    CALL interp % Free()

  END SUBROUTINE WriteTecplot_Model1D

  SUBROUTINE Init_Model2D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(out) :: this
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
    this % fluxDivMethod = SELF_CONSERVATIVE_FLUX 

    CALL this % solution % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % velocity % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % compVelocity % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nVar,this % mesh % nElem)

    ! set default metadata
    DO ivar = 1,nvar
      WRITE(ivarChar,'(I3.3)') ivar
      varname="solution"//TRIM(ivarChar)
      CALL this % solution % SetName(ivar,varname)
      CALL this % solution % SetUnits(ivar,"[null]")
    ENDDO

  END SUBROUTINE Init_Model2D

  SUBROUTINE Free_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % solution % Free()
    CALL this % workSol % Free()
    CALL this % velocity % Free()
    CALL this % compVelocity % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_Model2D

  SUBROUTINE UpdateHost_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % mesh % UpdateHost()
    CALL this % geometry % UpdateHost()
    CALL this % solution % UpdateHost()
    CALL this % dSdt % UpdateHost()
    CALL this % solution % UpdateHost()
    CALL this % velocity % UpdateHost()
    CALL this % solutionGradient % UpdateHost()
    CALL this % flux % UpdateHost()
    CALL this % source % UpdateHost()
    CALL this % fluxDivergence % UpdateHost()

  END SUBROUTINE UpdateHost_Model2D

  SUBROUTINE UpdateDevice_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % mesh % UpdateDevice()
    CALL this % geometry % UpdateDevice()
    CALL this % dSdt % UpdateDevice()
    CALL this % solution % UpdateDevice()
    CALL this % velocity % UpdateDevice()
    CALL this % solutionGradient % UpdateDevice()
    CALL this % flux % UpdateDevice()
    CALL this % source % UpdateDevice()
    CALL this % fluxDivergence % UpdateDevice()

  END SUBROUTINE UpdateDevice_Model2D

  SUBROUTINE SetSolutionFromEqn_Model2D(this, eqn) 
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

      ! Copy the equation parser
      DO iVar = 1, this % solution % nVar
        CALL this % solution % SetEquation(ivar, eqn(iVar) % equation)
      ENDDO

      CALL this % solution % SetInteriorFromEquation( this % geometry, this % t )
      CALL this % solution % BoundaryInterp( gpuAccel = .FALSE. )

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateDevice()
      ENDIF

  END SUBROUTINE SetSolutionFromEqn_Model2D 

  SUBROUTINE SetVelocityFieldFromEqn_Model2D(this, eqn) 
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn(1:2)

      ! Copy the equation parser
      ! Set the x-component of the velocity
      CALL this % velocity % SetEquation(1,1,eqn(1) % equation)

      ! Set the y-component of the velocity
      CALL this % velocity % SetEquation(2,1,eqn(2) % equation)

      ! Set the velocity values using the equation parser
      CALL this % velocity % SetInteriorFromEquation( this % geometry, this % t )

      CALL this % velocity % BoundaryInterp( gpuAccel = .FALSE. )

      IF( this % gpuAccel )THEN
        CALL this % velocity % UpdateDevice()
      ENDIF

  END SUBROUTINE SetVelocityFieldFromEqn_Model2D 

  SUBROUTINE SetVelocityFieldFromChar_Model2D(this, eqnChar) 
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:2)

      ! Set the x-component of the velocity
      CALL this % velocity % SetEquation(1,1,eqnChar(1))

      ! Set the y-component of the velocity
      CALL this % velocity % SetEquation(2,1,eqnChar(2))

      ! Set the velocity values using the equation parser
      CALL this % velocity % SetInteriorFromEquation( this % geometry, this % t )

      CALL this % velocity % BoundaryInterp( gpuAccel = .FALSE. )

      IF( this % gpuAccel )THEN
        CALL this % velocity % UpdateDevice()
      ENDIF

  END SUBROUTINE SetVelocityFieldFromChar_Model2D

  SUBROUTINE SetSolutionFromChar_Model2D(this, eqnChar) 
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % solution % nVar)
    ! Local
    INTEGER :: iVar

      DO iVar = 1, this % solution % nVar
        CALL this % solution % SetEquation(ivar, eqnChar(iVar))
      ENDDO

      CALL this % solution % SetInteriorFromEquation( this % geometry, this % t )
      CALL this % solution % BoundaryInterp( gpuAccel = .FALSE. )

      IF( this % gpuAccel )THEN
        CALL this % solution % UpdateDevice()
      ENDIF


  END SUBROUTINE SetSolutionFromChar_Model2D

  SUBROUTINE UpdateSolution_Model2D(this,dt)
    !! Computes a solution update as `s=s+dt*dsdt`, where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    REAL(prec),OPTIONAL,INTENT(in) :: dt
    ! Local
    REAL(prec) :: dtLoc
    INTEGER :: i, j, iVar, iEl

    IF (PRESENT(dt)) THEN
      dtLoc = dt
    ELSE 
      dtLoc = this % dt
    ENDIF

    IF (this % gpuAccel) THEN

      CALL UpdateSolution_Model2D_gpu_wrapper( this % solution % interior % deviceData, &
                                      this % dSdt % interior % deviceData, &
                                      dtLoc, &
                                      this % solution % interp % N, &
                                      this % solution % nVar, &
                                      this % solution % nElem ) 
                                      

    ELSE

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              this % solution % interior % hostData(i,j,iVar,iEl) = &
                  this % solution % interior % hostData(i,j,iVar,iEl) +&
                  dtLoc*this % dSdt % interior % hostData(i,j,iVar,iEl)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE UpdateSolution_Model2D

  SUBROUTINE UpdateGRK3_Model2D(this,m)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    INTEGER, INTENT(in) :: m
    ! Local
    INTEGER :: i, j, iVar, iEl

    IF (this % gpuAccel) THEN

      CALL UpdateGRK3_Model2D_gpu_wrapper( this % workSol % interior % deviceData, &
                                           this % solution % interior % deviceData, &
                                           this % dSdt % interior % deviceData, &
                                           rk3_a(m),rk3_g(m),this % dt, &
                                           this % solution % interp % N, &
                                           this % solution % nVar, &
                                           this % solution % nElem ) 
                                      

    ELSE

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              this % workSol % interior % hostData(i,j,iVar,iEl) = rk3_a(m)*&
                     this % workSol % interior % hostData(i,j,iVar,iEl) + &
                     this % dSdt % interior % hostData(i,j,iVar,iEl)


              this % solution % interior % hostData(i,j,iVar,iEl) = &
                      this % solution % interior % hostData(i,j,iVar,iEl) + &
                      rk3_g(m)*this % dt*this % workSol % interior % hostData(i,j,iVar,iEl)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE UpdateGRK3_Model2D

  SUBROUTINE ReprojectFlux_Model2D(this) 
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

      CALL this % flux % ContravariantProjection(this % geometry, this % gpuAccel)

  END SUBROUTINE ReprojectFlux_Model2D

  SUBROUTINE CalculateFluxDivergence_Model2D(this)
    !! Calculates the divergence of the flux vector using either the split-form or conservative formulation.
    !! If the split-form is used, you need to set the velocity field
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

      IF (this % fluxDivMethod == SELF_SPLITFORM_FLUX) THEN
        CALL this % velocity % ContravariantProjection(this % geometry, this % gpuAccel)

        IF (this % gpuAccel) THEN
          CALL this % flux % interp % VectorDGDivergence_2D(this % flux % interior % deviceData, &
                                                           this % solution % interior % deviceData, &
                                                           this % compVelocity % interior % deviceData, &
                                                           this % flux % boundaryNormal % deviceData, &
                                                           this % fluxDivergence % interior % deviceData, &
                                                           this % flux % nvar, &
                                                           this % flux % nelem)
        ELSE
          CALL this % flux % interp % VectorDGDivergence_2D(this % flux % interior % hostData, &
                                                           this % solution % interior % hostData, &
                                                           this % compVelocity % interior % hostData, &
                                                           this % flux % boundaryNormal % hostData, &
                                                           this % fluxDivergence % interior % hostData, &
                                                           this % flux % nvar, &
                                                           this % flux % nelem)
        END IF

      ELSE ! Conservative Form

        CALL this % flux % Divergence(this % geometry, &
                                      this % fluxDivergence, &
                                      selfWeakDGForm,&
                                      this % gpuAccel)
      ENDIF

  END SUBROUTINE CalculateFluxDivergence_Model2D

  SUBROUTINE CalculateTendency_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, iVar, iEl

    CALL this % PreTendency()
    CALL this % solution % BoundaryInterp(this % gpuAccel)
    CALL this % solution % SideExchange(this % mesh, this % decomp, this % gpuAccel)
    CALL this % SetBoundaryCondition()
    CALL this % SourceMethod()
    CALL this % RiemannSolver()
    CALL this % FluxMethod()
    CALL this % flux % ContravariantProjection(this % geometry, this % gpuAccel)
    CALL this % CalculateFluxDivergence()

    IF( this % gpuAccel )THEN

      CALL CalculateDSDt_Model2D_gpu_wrapper( this % fluxDivergence % interior % deviceData, &
                                      this % source % interior % deviceData, &
                                      this % dSdt % interior % deviceData, &
                                      this % solution % interp % N, &
                                      this % solution % nVar, &
                                      this % solution % nElem ) 
                                      
    ELSE

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              this % dSdt % interior % hostData(i,j,iVar,iEl) = &
                      this % source % interior % hostData(i,j,iVar,iEl) -&
                      this % fluxDivergence % interior % hostData(i,j,iVar,iEl)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE CalculateTendency_Model2D

  SUBROUTINE Write_Model2D(this,fileName)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:4)
    INTEGER(HID_T) :: xOffset(1:5)
    INTEGER(HID_T) :: bOffset(1:4)
    INTEGER(HID_T) :: bxOffset(1:5)
    INTEGER(HID_T) :: solGlobalDims(1:4)
    INTEGER(HID_T) :: xGlobalDims(1:5)
    INTEGER(HID_T) :: bGlobalDims(1:4)
    INTEGER(HID_T) :: bxGlobalDims(1:5)
    INTEGER :: firstElem
    ! Local
    CHARACTER(LEN=self_FileNameLength) :: pickupFile
    CHARACTER(13) :: timeStampString

    IF( PRESENT(filename) )THEN
      pickupFile = filename
    ELSE
      timeStampString = TimeStamp(this % t, 's')
      pickupFile = 'solution.'//timeStampString//'.h5'
    ENDIF

    IF( this % gpuAccel ) THEN
      CALL this % solution % UpdateHost()
      CALL this % solutionGradient % UpdateHost()
    ENDIF

    IF (this % decomp % mpiEnabled) THEN

      CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId,this % decomp % mpiComm)

      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)
      solOffset(1:4) = (/0,0,0,firstElem/)
      solGlobalDims(1:4) = (/this % solution % interp % N+1, &
                             this % solution % interp % N+1, &
                             this % solution % nVar, &
                             this % decomp % nElem/)



      xOffset(1:5) = (/0,0,0,0,firstElem/)
      xGlobalDims(1:5) = (/2, &
                           this % solution % interp % N+1, &
                           this % solution % interp % N+1, &
                           this % solution % nVar, &
                           this % decomp % nElem/)

      ! Offsets and dimensions for element boundary data
      bOffset(1:4) = (/0,0,0,firstElem/)
      bGlobalDims(1:4) = (/this % solution % interp % N+1, &
                           this % solution % nVar, &
                           4,&
                           this % decomp % nElem/)

      bxOffset(1:5) = (/0,0,0,0,firstElem/)
      bxGlobalDims(1:5) = (/2,&
                           this % solution % interp % N+1, &
                           this % solution % nVar, &
                           4,&
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


      CALL CreateGroup_HDF5(fileId,'/state')

      CALL CreateGroup_HDF5(fileId,'/state/interior')

      CALL CreateGroup_HDF5(fileId,'/state/boundary')

      CALL CreateGroup_HDF5(fileId,'/mesh')

      CALL CreateGroup_HDF5(fileId,'/mesh/interior')

      CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

      CALL WriteArray_HDF5(fileId,'/state/interior/solution', &
                           this % solution % interior,solOffset,solGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solution', &
                           this % solution % boundary,bOffset,bGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/fluxDivergence', &
                           this % fluxDivergence % interior,solOffset,solGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/flux', &
                           this % flux % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/flux', &
                           this % flux % boundary,bxOffset,bxGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/solutionGradient', &
                           this % solutionGradient % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solutionGradient', &
                           this % solutionGradient % boundary,bxOffset,bxGlobalDims)

      CALL WriteArray_HDF5(fileId,'/mesh/interior/x', &
                           this % geometry % x % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/mesh/boundary/x', &
                           this % geometry % x % boundary,bxOffset,bxGlobalDims)

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

      CALL CreateGroup_HDF5(fileId,'/state')

      CALL CreateGroup_HDF5(fileId,'/state/interior')

      CALL CreateGroup_HDF5(fileId,'/state/boundary')

      CALL CreateGroup_HDF5(fileId,'/mesh')

      CALL CreateGroup_HDF5(fileId,'/mesh/interior')

      CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

      CALL WriteArray_HDF5(fileId,'/state/interior/solution',this % solution % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solution',this % solution % boundary)

      CALL WriteArray_HDF5(fileId,'/state/interior/fluxDivergence',this % fluxDivergence % interior)

      CALL WriteArray_HDF5(fileId,'/state/interior/flux',this % flux % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/flux',this % flux % boundary)

      CALL WriteArray_HDF5(fileId,'/state/interior/solutionGradient',this % solutionGradient % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solutionGradient',this % solutionGradient % boundary)

      CALL WriteArray_HDF5(fileId,'/mesh/interior/x',this % geometry % x % interior)

      CALL WriteArray_HDF5(fileId,'/mesh/boundary/x',this % geometry % x % boundary)

      CALL Close_HDF5(fileId)

    END IF

  END SUBROUTINE Write_Model2D

  SUBROUTINE Read_Model2D(this,fileName)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
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
      CALL ReadArray_HDF5(fileId,'/state/interior/solution', &
                          this % solution % interior,solOffset)
    ELSE
      CALL ReadArray_HDF5(fileId,'/state/interior/solution',this % solution % interior)
    END IF

    CALL Close_HDF5(fileId)

    IF( this % gpuAccel )THEN
      CALL this % solution % interior % UpdateDevice()
    ENDIF

  END SUBROUTINE Read_Model2D

  SUBROUTINE WriteTecplot_Model2D(this, filename)
    IMPLICIT NONE
    CLASS(Model2D), INTENT(inout) :: this
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
    TYPE(Scalar2D) :: solution
    TYPE(Vector2D) :: x
    TYPE(Lagrange),TARGET :: interp

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
                      
    ! Create an interpolant for the uniform grid
    CALL interp % Init(this % solution % interp % M,&
            this % solution % interp % targetNodeType,&
            this % solution % interp % N, &
            this % solution % interp % controlNodeType)

    CALL solution % Init( interp, &
            this % solution % nVar, this % solution % nElem )

    CALL x % Init( interp, 1, this % solution % nElem )

    ! Map the mesh positions to the target grid
    CALL this % geometry % x % GridInterp(x, gpuAccel=.FALSE.)

    ! Map the solution to the target grid
    CALL this % solution % GridInterp(solution,gpuAccel=.FALSE.)
   
     OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(tecFile), &
      FORM='formatted', &
      STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y"'
    DO iVar = 1, this % solution % nVar
      tecHeader = TRIM(tecHeader)//', "'//TRIM(this % solution % meta(iVar) % name)//'"'
    ENDDO

    WRITE(fUnit,*) TRIM(tecHeader) 

    ! Create format statement
    WRITE(fmat,*) this % solution % nvar+2
    fmat = '('//TRIM(fmat)//'(ES16.7E3,1x))'

    DO iEl = 1, this % solution % nElem

      ! TO DO :: Get the global element ID 
      WRITE(zoneID,'(I8.8)') iEl
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',this % solution % interp % M+1,&
                                                 ', J=',this % solution % interp % M+1

      DO j = 0, this % solution % interp % M
        DO i = 0, this % solution % interp % M

          WRITE(fUnit,fmat) x % interior % hostData(1,i,j,1,iEl), &
                            x % interior % hostData(2,i,j,1,iEl), &
                            solution % interior % hostData(i,j,1:this % solution % nvar,iEl)

        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)

    CALL x % Free()
    CALL solution % Free() 
    CALL interp % Free()

  END SUBROUTINE WriteTecplot_Model2D

END MODULE SELF_Model
