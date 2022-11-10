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


  ! Runge-Kutta 2nd Order (Low Storage)
  REAL(prec),PARAMETER :: rk2_a(1:2) = (/0.0_prec,-0.5_prec/)
  REAL(prec),PARAMETER :: rk2_b(1:2) = (/0.5_prec,0.5_prec/)
  REAL(prec),PARAMETER :: rk2_g(1:2) = (/0.5_prec,1.0_prec/)

  ! Williamson's Runge-Kutta 3rd Order (Low Storage)
  REAL(prec),PARAMETER :: rk3_a(1:3) = (/0.0_prec,-5.0_prec/9.0_prec,-153.0_prec/128.0_prec/)
  REAL(prec),PARAMETER :: rk3_b(1:3) = (/0.0_prec,1.0_prec/3.0_prec,3.0_prec/4.0_prec/)
  REAL(prec),PARAMETER :: rk3_g(1:3) = (/1.0_prec/3.0_prec,15.0_prec/16.0_prec,8.0_prec/15.0_prec/)

  ! Carpenter-Kennedy Runge-Kuttta 4th Order (Low Storage)
  REAL(prec),PARAMETER :: rk4_a(1:5) = (/0.0_prec, &
          -1.0_prec, &
          -1.0_prec/3.0_prec+2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
          -2.0_prec**(1.0_prec/3.0_prec)-2.0_prec**(2.0_prec/3.0_prec)-2.0_prec, &
          -1.0_prec+2.0_prec**(1.0_prec/3.0_prec) /)

  REAL(prec),PARAMETER :: rk4_b(1:5) = (/ 0.0_prec, &
          2.0_prec/3.0_prec+2.0_prec**(1.0_prec/3.0_prec)/3.0_prec+2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
          2.0_prec/3.0_prec+2.0_prec**(1.0_prec/3.0_prec)/3.0_prec+2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
          1.0_prec/3.0_prec-2.0_prec**(1.0_prec/3.0_prec)/3.0_prec-2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
          1.0_prec/)
  
  REAL(prec),PARAMETER :: rk4_g(1:5) = (/&
          2.0_prec/3.0_prec+2.0_prec**(1.0_prec/3.0_prec)/3.0_prec+2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
          -2.0_prec**(2.0_prec/3.0_prec)/6.0_prec+1.0_prec/6.0_prec, &
          -1.0_prec/3.0_prec-2.0_prec*2.0_prec**(1.0_prec/3.0_prec)/3.0_prec-2.0_prec**(2.0_prec/3.0_prec)/3.0_prec, &
          1.0_prec/3.0_prec-2.0_prec**(1.0_prec/3.0_prec)/3.0_prec-2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
          1.0_prec/3.0_prec+2.0_prec**(1.0_prec/3.0_prec)/6.0_prec+2.0_prec**(2.0_prec/3.0_prec)/12.0_prec /)



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


  TYPE,ABSTRACT :: Model
    LOGICAL :: gpuAccel

    ! Time integration attributes
    PROCEDURE(SELF_timeIntegrator), POINTER :: timeIntegrator => Euler_timeIntegrator 
    REAL(prec) :: dt
    REAL(prec) :: t

    ! Standard Diagnostics
    REAL(prec) :: entropy ! Mathematical entropy function for the model

    ! Domain Decomposition
    TYPE(MPILayer),POINTER :: decomp

    CONTAINS

    PROCEDURE :: ForwardStep => ForwardStep_Model

    PROCEDURE :: Euler_timeIntegrator 

    ! Adams-Bashforth Methods
    PROCEDURE(ResizePrevSol),DEFERRED :: ResizePrevSol

    PROCEDURE :: AdamsBashforth2_timeIntegrator
    PROCEDURE(UpdateGAB),DEFERRED :: UpdateGAB2

    PROCEDURE :: AdamsBashforth3_timeIntegrator
    PROCEDURE(UpdateGAB),DEFERRED :: UpdateGAB3

    PROCEDURE :: AdamsBashforth4_timeIntegrator
    PROCEDURE(UpdateGAB),DEFERRED :: UpdateGAB4

    ! Runge-Kutta methods
    PROCEDURE :: LowStorageRK2_timeIntegrator 
    PROCEDURE(UpdateGRK),DEFERRED :: UpdateGRK2

    PROCEDURE :: LowStorageRK3_timeIntegrator 
    PROCEDURE(UpdateGRK),DEFERRED :: UpdateGRK3

    PROCEDURE :: LowStorageRK4_timeIntegrator 
    PROCEDURE(UpdateGRK),DEFERRED :: UpdateGRK4

!    PROCEDURE :: CrankNicholson_timeIntegrator

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

    PROCEDURE :: EnableGPUAccel => EnableGPUAccel_Model
    PROCEDURE :: DisableGPUAccel => DisableGPUAccel_Model

  END TYPE Model

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
    SUBROUTINE ResizePrevSol( this, m )
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
    END SUBROUTINE ResizePrevSol
  END INTERFACE

  INTERFACE 
    SUBROUTINE UpdateGAB( this, m )
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
    END SUBROUTINE UpdateGAB
  END INTERFACE

  INTERFACE 
    SUBROUTINE UpdateGRK( this, m )
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model),INTENT(inout) :: this
      INTEGER,INTENT(in) :: m
    END SUBROUTINE UpdateGRK
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
          this % timeIntegrator => Euler_timeIntegrator
        CASE ( SELF_AB2 )
          this % timeIntegrator => AdamsBashforth2_timeIntegrator
          CALL this % ResizePrevSol(2)
        CASE ( SELF_AB3 )
          this % timeIntegrator => AdamsBashforth3_timeIntegrator
          CALL this % ResizePrevSol(3)
        CASE ( SELF_AB4 )
          this % timeIntegrator => AdamsBashforth4_timeIntegrator
          CALL this % ResizePrevSol(4)
        CASE ( SELF_RK2 )
          this % timeIntegrator => LowStorageRK2_timeIntegrator
        CASE ( SELF_RK3 )
          this % timeIntegrator => LowStorageRK3_timeIntegrator
        CASE ( SELF_RK4 )
          this % timeIntegrator => LowStorageRK4_timeIntegrator
        CASE DEFAULT
          this % timeIntegrator => LowStorageRK3_timeIntegrator

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
          this % timeIntegrator => Euler_timeIntegrator

        CASE ("AB2")
          this % timeIntegrator => AdamsBashforth2_timeIntegrator
          CALL this % ResizePrevSol(2)

        CASE ("AB3")
          this % timeIntegrator => AdamsBashforth3_timeIntegrator
          CALL this % ResizePrevSol(3)

        CASE ("AB4")
          this % timeIntegrator => AdamsBashforth4_timeIntegrator
          CALL this % ResizePrevSol(4)

        CASE ("RK2")
          this % timeIntegrator => LowStorageRK2_timeIntegrator

        CASE ("RK3")
          this % timeIntegrator => LowStorageRK3_timeIntegrator

        CASE ("RK4")
          this % timeIntegrator => LowStorageRK4_timeIntegrator

        CASE DEFAULT
          this % timeIntegrator => LowStorageRK3_timeIntegrator

      END SELECT

  END SUBROUTINE SetTimeIntegrator_withChar

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
        CALL this % timeIntegrator(tNext)
        this % t = tNext
        CALL this % WriteModel()
        CALL this % WriteTecplot()
        CALL this % CalculateEntropy()
        CALL this % ReportEntropy()
      ENDDO

    ELSE
      CALL this % timeIntegrator(targetTime)
      this % t = targetTime
      CALL this % CalculateEntropy()
      CALL this % ReportEntropy()
    ENDIF

  END SUBROUTINE ForwardStep_Model

  SUBROUTINE Euler_timeIntegrator(this,tn)
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

  END SUBROUTINE Euler_timeIntegrator

  SUBROUTINE AdamsBashforth2_timeIntegrator(this,tn)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tn
    ! Local
    INTEGER :: m
    REAL(prec) :: tRemain
    REAL(prec) :: dtLim
    REAL(prec) :: t0

    dtLim = this % dt ! Get the max time step size from the dt attribute
    t0 = this % t

    ! Do a single step with RK2
    ! Initialize the PrevSol attribute
    CALL this % UpdateGAB2(0)
    CALL this % LowStorageRK2_timeIntegrator(t0+this%dt)

    DO WHILE (this % t < tn)

      t0 = this % t
      tRemain = tn - this % t
      this % dt = MIN( dtLim, tRemain )

      CALL this % UpdateGAB2(2) ! Store the solution in PrevSol and store the interpolated
                                ! solution in the solution attribute for tendency calculation
      CALL this % CalculateTendency()
      CALL this % UpdateGAB2(1) ! Reset the solution from the PrevSol
      CALL this % UpdateSolution()

      this % t = t0 + this % dt

    ENDDO 

    this % dt = dtLim

  END SUBROUTINE AdamsBashforth2_timeIntegrator

  SUBROUTINE AdamsBashforth3_timeIntegrator(this,tn)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tn
    ! Local
    INTEGER :: m
    REAL(prec) :: tRemain
    REAL(prec) :: dtLim
    REAL(prec) :: t0

    dtLim = this % dt ! Get the max time step size from the dt attribute

    ! Do two time steps with RK3
    ! Initialize the PrevSol attribute
    t0 = this % t
    CALL this % UpdateGAB3(0)
    CALL this % LowStorageRK3_timeIntegrator(t0+this%dt)

    t0 = this % t
    CALL this % UpdateGAB3(1)
    CALL this % LowStorageRK3_timeIntegrator(t0+this%dt)

    DO WHILE (this % t < tn)

      t0 = this % t
      tRemain = tn - this % t
      this % dt = MIN( dtLim, tRemain )

      CALL this % UpdateGAB3(3) ! Store the solution in PrevSol and store the interpolated
                                ! solution in the solution attribute for tendency calculation
      CALL this % CalculateTendency()
      CALL this % UpdateGAB3(2) ! Reset the solution from the PrevSol
      CALL this % UpdateSolution()

      this % t = t0 + this % dt

    ENDDO 

    this % dt = dtLim

  END SUBROUTINE AdamsBashforth3_timeIntegrator

  SUBROUTINE AdamsBashforth4_timeIntegrator(this,tn)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tn
    ! Local
    INTEGER :: m
    REAL(prec) :: tRemain
    REAL(prec) :: dtLim
    REAL(prec) :: t0

    dtLim = this % dt ! Get the max time step size from the dt attribute

    ! Do three time steps with RK4
    ! Initialize the PrevSol attribute
    t0 = this % t
    CALL this % UpdateGAB4(0)
    CALL this % LowStorageRK4_timeIntegrator(t0+this%dt)

    t0 = this % t
    CALL this % UpdateGAB4(1)
    CALL this % LowStorageRK4_timeIntegrator(t0+this%dt)

    t0 = this % t
    CALL this % UpdateGAB4(2)
    CALL this % LowStorageRK4_timeIntegrator(t0+this%dt)

    DO WHILE (this % t < tn)

      t0 = this % t
      tRemain = tn - this % t
      this % dt = MIN( dtLim, tRemain )

      CALL this % UpdateGAB4(4) ! Store the solution in PrevSol and store the interpolated
                                ! solution in the solution attribute for tendency calculation
      CALL this % CalculateTendency()
      CALL this % UpdateGAB4(3) ! Reset the solution from the PrevSol
      CALL this % UpdateSolution()

      this % t = t0 + this % dt

    ENDDO 

    this % dt = dtLim

  END SUBROUTINE AdamsBashforth4_timeIntegrator

  SUBROUTINE LowStorageRK2_timeIntegrator(this,tn)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tn
    ! Local
    INTEGER :: m
    REAL(prec) :: tRemain
    REAL(prec) :: dtLim
    REAL(prec) :: t0

    dtLim = this % dt ! Get the max time step size from the dt attribute
    DO WHILE (this % t < tn)

      t0 = this % t
      tRemain = tn - this % t
      this % dt = MIN( dtLim, tRemain )
      DO m = 1, 2
        CALL this % CalculateTendency()
        CALL this % UpdateGRK2(m)
        this % t = t0 + rk2_b(m)*this % dt
      ENDDO

      this % t = t0 + this % dt

    ENDDO 

    this % dt = dtLim

  END SUBROUTINE LowStorageRK2_timeIntegrator

  SUBROUTINE LowStorageRK3_timeIntegrator(this,tn)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tn
    ! Local
    INTEGER :: m
    REAL(prec) :: tRemain
    REAL(prec) :: dtLim
    REAL(prec) :: t0

    dtLim = this % dt ! Get the max time step size from the dt attribute
    DO WHILE (this % t < tn)

      t0 = this % t
      tRemain = tn - this % t
      this % dt = MIN( dtLim, tRemain )
      DO m = 1, 3
        CALL this % CalculateTendency()
        CALL this % UpdateGRK3(m)
        this % t = t0 + rk3_b(m)*this % dt
      ENDDO

      this % t = t0 + this % dt

    ENDDO 

    this % dt = dtLim

  END SUBROUTINE LowStorageRK3_timeIntegrator

  SUBROUTINE LowStorageRK4_timeIntegrator(this,tn)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec), INTENT(in) :: tn
    ! Local
    INTEGER :: m
    REAL(prec) :: tRemain
    REAL(prec) :: dtLim
    REAL(prec) :: t0

    dtLim = this % dt ! Get the max time step size from the dt attribute
    DO WHILE (this % t < tn)

      t0 = this % t
      tRemain = tn - this % t
      this % dt = MIN( dtLim, tRemain )
      DO m = 1, 5
        CALL this % CalculateTendency()
        CALL this % UpdateGRK4(m)
        this % t = t0 + rk4_b(m)*this % dt
      ENDDO

      this % t = t0 + this % dt

    ENDDO 

    this % dt = dtLim

  END SUBROUTINE LowStorageRK4_timeIntegrator

!  SUBROUTINE CrankNicholson_timeIntegrator(this,tn)
!    !! Solves the equation formed by the Crank Nicholson method
!    !! using JFNK, where the Krylov solver is chosen to be
!    !! BiCG-Stabilized.
!    IMPLICIT NONE
!    CLASS(Model),INTENT(inout) :: this
!    REAL(prec), INTENT(in) :: tn
!    ! Local
!    INTEGER :: m
!    REAL(prec) :: tRemain
!    REAL(prec) :: dtLim
!    REAL(prec) :: t0
!
!    dtLim = this % dt ! Get the max time step size from the dt attribute
!    DO WHILE (this % t < tn)
!
!      t0 = this % t
!      tRemain = tn - this % t
!      this % dt = MIN( dtLim, tRemain )
!      ! Copy existing solution to old solution
!
!      ! Evaluate tendency with old solution
!
!      ! Calculate r_k (fixed) -> Store in PrevSol
!
!      ! Calculate Fk(m-1)
!
!      DO m = 1, SELF_maxJFNKiterations
!
!        
!        ! Linear iterations on Jk(m-1) dS(m) = -Fk(m-1)
!        CALL this % JFNKLinearSolver(t0+dt) ! Use PrevSol to store Fk(m-1), sk(m-1), dSm, rk
!                                            ! Linear solver updates
!                                            ! dSm, sk(m), and Fk(m)
!
!        ! Check for convergence
!        !  >> global reduction on dSm (either l2 or lmax)
!        !  >> global reduction on Fkm1 (either l2 or lmax)
!        !
!        !    -> Both done as reduction on PrevSol attribute - reduction over grid, not variables
!        !
!
!      ENDDO
!
!      this % t = t0 + this % dt
!
!    ENDDO 
!
!    this % dt = dtLim
!
!  END SUBROUTINE CrankNicholson_timeIntegrator

END MODULE SELF_Model
