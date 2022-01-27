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
  REAL(prec),PARAMETER :: rk3_a(1:3) = (/0.0_prec,-5.0_prec/9.0_prec,-153.0_prec/128.0_prec/)
  REAL(prec),PARAMETER :: rk3_b(1:3) = (/0.0_prec,1.0_prec/3.0_prec,3.0_prec/4.0_prec/)
  REAL(prec),PARAMETER :: rk3_g(1:3) = (/1.0_prec/3.0_prec,15.0_prec/16.0_prec,8.0_prec/15.0_prec/)

! 
  INTEGER, PARAMETER :: SELF_EULER = 100
  INTEGER, PARAMETER :: SELF_RK3 = 300
  INTEGER, PARAMETER :: SELF_RK4 = 400

  INTEGER, PARAMETER :: SELF_INTEGRATOR_LENGTH = 10 ! max length of integrator methods when specified as char
  INTEGER, PARAMETER :: SELF_EQUATION_LENGTH = 500

! //////////////////////////////////////////////// !

  TYPE,ABSTRACT :: Model
    LOGICAL :: gpuAccel

    ! Time integration attributes
    INTEGER :: timeIntegrator
    REAL(prec) :: dt
    REAL(prec) :: t

    CONTAINS

    PROCEDURE :: ForwardStep => ForwardStep_Model
    PROCEDURE :: ForwardStepEuler => ForwardStepEuler_Model 

    PROCEDURE(UpdateSolution),DEFERRED :: UpdateSolution
    PROCEDURE(CalculateTendency),DEFERRED :: CalculateTendency

    GENERIC :: SetTimeIntegrator => SetTimeIntegrator_withInt, &
                                    SetTimeIntegrator_withChar
    PROCEDURE,PRIVATE :: SetTimeIntegrator_withInt
    PROCEDURE,PRIVATE :: SetTimeIntegrator_withChar

    PROCEDURE :: SetSimulationTime
    PROCEDURE :: GetSimulationTime
!    PROCEDURE :: SetTimeStep
!    PROCEDURE :: GetTimeStep


    PROCEDURE :: EnableGPUAccel => EnableGPUAccel_Model
    PROCEDURE :: DisableGPUAccel => DisableGPUAccel_Model

  END TYPE Model


  TYPE,EXTENDS(Model),ABSTRACT :: Model2D
    TYPE(MappedScalar2D) :: solution
    TYPE(MappedVector2D) :: solutionGradient
    TYPE(MappedVector2D) :: flux
    TYPE(MappedScalar2D) :: source
    TYPE(MappedScalar2D) :: fluxDivergence
    TYPE(MappedScalar2D) :: dSdt
    TYPE(MPILayer),POINTER :: decomp
    TYPE(Mesh2D),POINTER :: mesh
    TYPE(SEMQuad),POINTER :: geometry

    CONTAINS

    PROCEDURE :: Init => Init_Model2D
    PROCEDURE :: Free => Free_Model2D

    PROCEDURE :: UpdateHost => UpdateHost_Model2D
    PROCEDURE :: UpdateDevice => UpdateDevice_Model2D

    PROCEDURE :: UpdateSolution => UpdateSolution_Model2D
    PROCEDURE :: CalculateTendency => CalculateTendency_Model2D
    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_Model2D

    GENERIC :: SetSolution => SetSolutionFromChar_Model2D,&
                              SetSolutionFromEqn_Model2D
    PROCEDURE,PRIVATE :: SetSolutionFromChar_Model2D
    PROCEDURE,PRIVATE :: SetSolutionFromEqn_Model2D

    PROCEDURE(Source2D),DEFERRED :: Source2D
    PROCEDURE(Flux2D),DEFERRED :: Flux2D
    PROCEDURE(RiemannSolver2D),DEFERRED :: RiemannSolver2D

    PROCEDURE :: ReprojectFlux => ReprojectFlux_Model2D

    PROCEDURE :: Read => Read_Model2D
    PROCEDURE :: Write => Write_Model2D

  END TYPE Model2D

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
    SUBROUTINE Flux2D( this )
      IMPORT Model2D
      IMPLICIT NONE
      CLASS(Model2D),INTENT(inout) :: this
    END SUBROUTINE Flux2D
  END INTERFACE

  INTERFACE 
    SUBROUTINE Source2D( this )
      IMPORT Model2D
      IMPLICIT NONE
      CLASS(Model2D),INTENT(inout) :: this
    END SUBROUTINE Source2D
  END INTERFACE

  INTERFACE 
    SUBROUTINE RiemannSolver2D( this )
      IMPORT Model2D
      IMPLICIT NONE
      CLASS(Model2D),INTENT(inout) :: this
    END SUBROUTINE RiemannSolver2D
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

      this % timeIntegrator = integrator

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

      SELECT CASE (upperCaseInt)

        CASE ("EULER")
          this % timeIntegrator = SELF_EULER

        CASE ("RK3")
          this % timeIntegrator = SELF_RK3

        CASE ("RK4")
          this % timeIntegrator = SELF_RK4

        CASE DEFAULT
          this % timeIntegrator = SELF_EULER

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
      ! TO DO : Warning to user that no GPU is available
    ENDIF

  END SUBROUTINE EnableGPUAccel_Model

  SUBROUTINE DisableGPUAccel_Model(this)
    IMPLICIT NONE
    CLASS(Model), INTENT(inout) :: this

    this % gpuAccel = .FALSE.

  END SUBROUTINE DisableGPUAccel_Model

  ! ////////////////////////////////////// !
  !       Time Integrators                 !

  SUBROUTINE ForwardStep_Model(this,tn,dt)
  !!  Forward steps the model using the associated tendency procedure and time integrator
  !!
  !!  If the final time `tn` is provided, the model is forward stepped to that final time,
  !!  otherwise, the model is forward stepped only a single time step
  !!  
  !!  If a time step is provided through the interface, the model time step size is updated
  !!  and that time step is used to update the model
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    REAL(prec),OPTIONAL,INTENT(in) :: tn
    REAL(prec),OPTIONAL,INTENT(in) :: dt
    ! Local
    INTEGER :: nSteps
    

    IF (PRESENT(dt)) THEN
      this % dt = dt
    ENDIF

    IF (PRESENT(tn)) THEN
      nSteps = INT( (tn - this % t)/(this % dt) )
    ELSE
      nSteps = 1
    ENDIF

    SELECT CASE (this % timeIntegrator)

      CASE (SELF_EULER)

        CALL this % ForwardStepEuler(nSteps)

!      CASE RK3
!
!        CALL this % ForwardStepRK3(nSteps)

      CASE DEFAULT
        ! TODO : Warn user that time integrator not valid, default to Euler
        CALL this % ForwardStepEuler(nSteps)

    END SELECT

  END SUBROUTINE ForwardStep_Model

  SUBROUTINE ForwardStepEuler_Model(this,nSteps)
    IMPLICIT NONE
    CLASS(Model),INTENT(inout) :: this
    INTEGER,INTENT(in) :: nSteps
    ! Local
    INTEGER :: i

    DO i = 1, nSteps

      CALL this % CalculateTendency()
      CALL this % UpdateSolution()

    ENDDO 

  END SUBROUTINE ForwardStepEuler_Model


  SUBROUTINE Init_Model2D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh2D),INTENT(in),TARGET :: mesh
    TYPE(SEMQuad),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry

    CALL this % solution % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nVar,this % mesh % nElem)

  END SUBROUTINE Init_Model2D

  SUBROUTINE Free_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % solution % Free()
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

  END SUBROUTINE SetSolutionFromEqn_Model2D 

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

  SUBROUTINE ReprojectFlux_Model2D(this) 
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iVar, j, i
    REAL(prec) :: Fx, Fy

      DO iEl = 1,this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              Fx = this % flux % interior % hostData(1,i,j,iVar,iEl)
              Fy = this % flux % interior % hostData(2,i,j,iVar,iEl)

              this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                this % geometry % dsdx % interior % hostData(1,1,i,j,1,iel)*Fx + &
                this % geometry % dsdx % interior % hostData(2,1,i,j,1,iel)*Fy 

              this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                this % geometry % dsdx % interior % hostData(1,2,i,j,1,iel)*Fx + &
                this % geometry % dsdx % interior % hostData(2,2,i,j,1,iel)*Fy 


            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE ReprojectFlux_Model2D

  SUBROUTINE CalculateFluxDivergence_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % flux % Divergence(this % geometry, &
                                  this % fluxDivergence, &
                                  selfWeakDGForm,&
                                  this % gpuAccel)

  END SUBROUTINE CalculateFluxDivergence_Model2D

  SUBROUTINE CalculateTendency_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, iVar, iEl

!      CALL this % solution % AverageSides()
!      CALL this % solution % DiffSides()
!      CALL this % SetBoundaryCondition()
      CALL this % Source2D()
      CALL this % RiemannSolver2D()
      CALL this % Flux2D()
      CALL this % ReprojectFlux()
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
    CLASS(Model2D),INTENT(in) :: this
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

    IF (this % decomp % mpiEnabled) THEN

      CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId,this % decomp % mpiComm)

      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)
      solOffset(1:4) = (/0,0,1,firstElem/)
      solGlobalDims(1:4) = (/this % solution % interp % N, &
                             this % solution % interp % N, &
                             this % solution % nVar, &
                             this % decomp % nElem/)


      xOffset(1:5) = (/1,0,0,1,firstElem/)
      xGlobalDims(1:5) = (/2, &
                           this % solution % interp % N, &
                           this % solution % interp % N, &
                           this % solution % nVar, &
                           this % decomp % nElem/)

      ! Offsets and dimensions for element boundary data
      bOffset(1:4) = (/0,1,1,firstElem/)
      bGlobalDims(1:4) = (/this % solution % interp % N, &
                           this % solution % nVar, &
                           4,&
                           this % decomp % nElem/)

      bxOffset(1:5) = (/1,0,1,1,firstElem/)
      bxGlobalDims(1:5) = (/2,&
                           this % solution % interp % N, &
                           this % solution % nVar, &
                           4,&
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

  END SUBROUTINE Read_Model2D

END MODULE SELF_Model
