! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Model

  use SELF_SupportRoutines
  use SELF_Metadata
  use SELF_MPI
  use SELF_HDF5
  use HDF5
  use FEQParse

  implicit none

#include "SELF_Macros.h"

! //////////////////////////////////////////////// !
!   Time integration parameters

  ! Runge-Kutta 2nd Order (Low Storage)
  real(prec),parameter :: rk2_a(1:2) = (/0.0_prec,-0.5_prec/)
  real(prec),parameter :: rk2_b(1:2) = (/0.5_prec,0.5_prec/)
  real(prec),parameter :: rk2_g(1:2) = (/0.5_prec,1.0_prec/)

  ! Williamson's Runge-Kutta 3rd Order (Low Storage)
  real(prec),parameter :: rk3_a(1:3) = (/0.0_prec,-5.0_prec/9.0_prec,-153.0_prec/128.0_prec/)
  real(prec),parameter :: rk3_b(1:3) = (/0.0_prec,1.0_prec/3.0_prec,3.0_prec/4.0_prec/)
  real(prec),parameter :: rk3_g(1:3) = (/1.0_prec/3.0_prec,15.0_prec/16.0_prec,8.0_prec/15.0_prec/)

  ! Carpenter-Kennedy Runge-Kuttta 4th Order (Low Storage)
  real(prec),parameter :: rk4_a(1:5) = (/0.0_prec, &
                                         -1.0_prec, &
                                         -1.0_prec/3.0_prec+ &
                                         2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
                                         -2.0_prec**(1.0_prec/3.0_prec)- &
                                         2.0_prec**(2.0_prec/3.0_prec)-2.0_prec, &
                                         -1.0_prec+2.0_prec**(1.0_prec/3.0_prec)/)

  real(prec),parameter :: rk4_b(1:5) = (/ &
                          0.0_prec, &
                          2.0_prec/3.0_prec+2.0_prec**(1.0_prec/3.0_prec)/3.0_prec+ &
                          2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
                          2.0_prec/3.0_prec+2.0_prec**(1.0_prec/3.0_prec)/3.0_prec+ &
                          2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
                          1.0_prec/3.0_prec-2.0_prec**(1.0_prec/3.0_prec)/3.0_prec- &
                          2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
                          1.0_prec/)

  real(prec),parameter :: rk4_g(1:5) = (/ &
                          2.0_prec/3.0_prec+2.0_prec**(1.0_prec/3.0_prec)/3.0_prec+ &
                          2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
                          -2.0_prec**(2.0_prec/3.0_prec)/6.0_prec+1.0_prec/6.0_prec, &
                          -1.0_prec/3.0_prec-2.0_prec*2.0_prec**(1.0_prec/3.0_prec)/3.0_prec- &
                          2.0_prec**(2.0_prec/3.0_prec)/3.0_prec, &
                          1.0_prec/3.0_prec-2.0_prec**(1.0_prec/3.0_prec)/3.0_prec- &
                          2.0_prec**(2.0_prec/3.0_prec)/6.0_prec, &
                          1.0_prec/3.0_prec+2.0_prec**(1.0_prec/3.0_prec)/6.0_prec+ &
                          2.0_prec**(2.0_prec/3.0_prec)/12.0_prec/)

!
  integer,parameter :: SELF_EULER = 100
  integer,parameter :: SELF_RK2 = 200
  integer,parameter :: SELF_RK3 = 300
  integer,parameter :: SELF_RK4 = 400
  ! integer,parameter :: SELF_AB2 = 201
  ! integer,parameter :: SELF_AB3 = 301
  ! integer,parameter :: SELF_AB4 = 401

  integer,parameter :: SELF_INTEGRATOR_LENGTH = 10 ! max length of integrator methods when specified as char
  integer,parameter :: SELF_EQUATION_LENGTH = 500

! //////////////////////////////////////////////// !
!   Boundary Condition parameters
!

  ! Conditions on the solution
  integer,parameter :: SELF_BC_PRESCRIBED = 100
  integer,parameter :: SELF_BC_RADIATION = 101
  integer,parameter :: SELF_BC_NONORMALFLOW = 102

  ! Conditions on the solution gradients
  integer,parameter :: SELF_BC_PRESCRIBED_STRESS = 200
  integer,parameter :: SELF_BC_NOSTRESS = 201

! //////////////////////////////////////////////// !
!   Model Formulations
!
  integer,parameter :: SELF_FORMULATION_LENGTH = 30 ! max length of integrator methods when specified as char

  type,abstract :: Model

    ! Time integration attributes
    procedure(SELF_timeIntegrator),pointer :: timeIntegrator => Euler_timeIntegrator
    real(prec) :: dt
    real(prec) :: t
    integer :: ioIterate = 0
    logical :: gradient_enabled = .false.

    ! Standard Diagnostics
    real(prec) :: entropy ! Mathematical entropy function for the model

    ! Domain Decomposition
    type(MPILayer),pointer :: decomp

  contains

    procedure :: IncrementIOCounter

    procedure :: PrintType => PrintType_Model

    procedure :: SetInitialConditions => SetInitialConditions_Model

    procedure :: ForwardStep => ForwardStep_Model

    procedure :: Euler_timeIntegrator

    ! Adams-Bashforth Methods
    procedure(ResizePrevSol),deferred :: ResizePrevSol

    ! Runge-Kutta methods
    procedure :: LowStorageRK2_timeIntegrator
    procedure(UpdateGRK),deferred :: UpdateGRK2

    procedure :: LowStorageRK3_timeIntegrator
    procedure(UpdateGRK),deferred :: UpdateGRK3

    procedure :: LowStorageRK4_timeIntegrator
    procedure(UpdateGRK),deferred :: UpdateGRK4

    procedure :: PreTendency => PreTendency_Model
    procedure :: PreFlux => PreFlux_Model
    procedure :: SourceMethod => Source_Model
    procedure :: FluxMethod => Flux_Model
    procedure :: RiemannSolver => RiemannSolver_Model
    procedure :: UpdateBoundary => UpdateBoundary_Model
    procedure :: SetBoundaryCondition => SetBoundaryCondition_Model
    procedure :: SetGradientBoundaryCondition => SetGradientBoundaryCondition_Model

    procedure :: ReportEntropy => ReportEntropy_Model
    procedure :: CalculateEntropy => CalculateEntropy_Model

    procedure(UpdateSolution),deferred :: UpdateSolution
    procedure(CalculateTendency),deferred :: CalculateTendency
    procedure(ReadModel),deferred :: ReadModel
    procedure(WriteModel),deferred :: WriteModel
    procedure(WriteTecplot),deferred :: WriteTecplot

    generic :: SetTimeIntegrator => SetTimeIntegrator_withInt, &
      SetTimeIntegrator_withChar
    procedure,private :: SetTimeIntegrator_withInt
    procedure,private :: SetTimeIntegrator_withChar

    procedure :: SetSimulationTime
    procedure :: GetSimulationTime

  endtype Model

  interface
    subroutine SELF_timeIntegrator(this,tn)
      use SELF_Constants,only:prec
      import Model
      implicit none
      class(Model),intent(inout) :: this
      real(prec),intent(in) :: tn
    endsubroutine SELF_timeIntegrator
  endinterface

  interface
    subroutine ResizePrevSol(this,m)
      import Model
      implicit none
      class(Model),intent(inout) :: this
      integer,intent(in) :: m
    endsubroutine ResizePrevSol
  endinterface

  interface
    subroutine UpdateGRK(this,m)
      import Model
      implicit none
      class(Model),intent(inout) :: this
      integer,intent(in) :: m
    endsubroutine UpdateGRK
  endinterface

  interface
    subroutine UpdateSolution(this,dt)
      use SELF_Constants,only:prec
      import Model
      implicit none
      class(Model),intent(inout) :: this
      real(prec),optional,intent(in) :: dt
    endsubroutine UpdateSolution
  endinterface

  interface
    subroutine CalculateTendency(this)
      import Model
      implicit none
      class(Model),intent(inout) :: this
    endsubroutine CalculateTendency
  endinterface

  interface
    subroutine WriteModel(this,filename)
      import Model
      implicit none
      class(Model),intent(inout) :: this
      character(*),intent(in),optional :: filename
    endsubroutine WriteModel
  endinterface

  interface
    subroutine ReadModel(this,filename)
      import Model
      implicit none
      class(Model),intent(inout) :: this
      character(*),intent(in) :: filename
    endsubroutine ReadModel
  endinterface

  interface
    subroutine WriteTecplot(this,filename)
      import Model
      implicit none
      class(Model),intent(inout) :: this
      character(*),intent(in),optional :: filename
    endsubroutine WriteTecplot
  endinterface

contains

  subroutine IncrementIOCounter(this)
    implicit none
    class(Model),intent(inout) :: this

    ! Increment the ioIterate
    this%ioIterate = this%ioIterate+1

  endsubroutine IncrementIOCounter

  function GetBCFlagForChar(charFlag) result(intFlag)
  !! This method is used to return the integer flag from a char for boundary conditions
  !!
    implicit none
    character(*),intent(in) :: charFlag
    integer :: intFlag

    select case(UpperCase(trim(charFlag)))

    case("PRESCRIBED")
      intFlag = SELF_BC_PRESCRIBED

    case("RADIATION")
      intFlag = SELF_BC_RADIATION

    case("NO_NORMAL_FLOW")
      intFlag = SELF_BC_NONORMALFLOW

    case("PRESCRIBED_STRESS")
      intFlag = SELF_BC_PRESCRIBED_STRESS

    case("NO_STRESS")
      intFlag = SELF_BC_NOSTRESS

    case DEFAULT
      intFlag = 0

    endselect

  endfunction GetBCFlagForChar

  subroutine PrintType_Model(this)
#undef __FUNC__
#define __FUNC__ "PrintType"
    implicit none
    class(Model),intent(in) :: this

    INFO("None")

  endsubroutine PrintType_Model

  subroutine PreTendency_Model(this)
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    implicit none
    class(Model),intent(inout) :: this

    return

  endsubroutine PreTendency_Model

  subroutine PreFlux_Model(this)
    !! PreFlux is a template routine that is used to house any additional calculations
    !! that you want to execute just before the calculation of flux terms.
    !! This default PreFlux simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    implicit none
    class(Model),intent(inout) :: this

    return

  endsubroutine PreFlux_Model

  subroutine Source_Model(this)
    !!
    implicit none
    class(Model),intent(inout) :: this

    return

  endsubroutine Source_Model

  subroutine RiemannSolver_Model(this)
    !!
    implicit none
    class(Model),intent(inout) :: this

    return

  endsubroutine RiemannSolver_Model

  subroutine Flux_Model(this)
    !!
    implicit none
    class(Model),intent(inout) :: this

    return

  endsubroutine Flux_Model

  subroutine UpdateBoundary_Model(this)
  !!
    implicit none
    class(Model),intent(inout) :: this

    return
  endsubroutine UpdateBoundary_Model

  subroutine SetBoundaryCondition_Model(this)
    implicit none
    class(Model),intent(inout) :: this

    return

  endsubroutine SetBoundaryCondition_Model

  subroutine SetGradientBoundaryCondition_Model(this)
    implicit none
    class(Model),intent(inout) :: this

    return

  endsubroutine SetGradientBoundaryCondition_Model

  subroutine SetTimeIntegrator_withInt(this,integrator)
    !! Sets the time integrator method, using an integer flag
    !!
    !! Valid options for  are
    !!
    !!    SELF_EULER
    !!    SELF_RK3
    !!    SELF_RK4
    !!
    implicit none
    class(Model),intent(inout) :: this
    integer,intent(in) :: integrator

    select case(integrator)

    case(SELF_EULER)
      this%timeIntegrator => Euler_timeIntegrator
    case(SELF_RK2)
      this%timeIntegrator => LowStorageRK2_timeIntegrator
    case(SELF_RK3)
      this%timeIntegrator => LowStorageRK3_timeIntegrator
    case(SELF_RK4)
      this%timeIntegrator => LowStorageRK4_timeIntegrator
    case DEFAULT
      this%timeIntegrator => LowStorageRK3_timeIntegrator

    endselect

  endsubroutine SetTimeIntegrator_withInt

  subroutine SetTimeIntegrator_withChar(this,integrator)
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
    implicit none
    class(Model),intent(inout) :: this
    character(*),intent(in) :: integrator
    ! Local
    character(SELF_INTEGRATOR_LENGTH) :: upperCaseInt

    upperCaseInt = UpperCase(trim(integrator))

    select case(trim(upperCaseInt))

    case("EULER")
      this%timeIntegrator => Euler_timeIntegrator

    case("RK2")
      this%timeIntegrator => LowStorageRK2_timeIntegrator

    case("RK3")
      this%timeIntegrator => LowStorageRK3_timeIntegrator

    case("RK4")
      this%timeIntegrator => LowStorageRK4_timeIntegrator

    case DEFAULT
      this%timeIntegrator => LowStorageRK3_timeIntegrator

    endselect

  endsubroutine SetTimeIntegrator_withChar

  subroutine GetSimulationTime(this,t)
    !! Returns the current simulation time stored in the model % t attribute
    implicit none
    class(Model),intent(in) :: this
    real(prec),intent(out) :: t

    t = this%t

  endsubroutine GetSimulationTime

  subroutine SetSimulationTime(this,t)
    !! Sets the model % t attribute with the provided simulation time
    implicit none
    class(Model),intent(inout) :: this
    real(prec),intent(in) :: t

    this%t = t

  endsubroutine SetSimulationTime

  subroutine SetInitialConditions_Model(this)
#undef __FUNC__
#define __FUNC__ "SetInitialConditions"
    implicit none
    class(Model),intent(inout) :: this

    INFO("No model, so nothing to set")

  endsubroutine SetInitialConditions_Model

  subroutine CalculateEntropy_Model(this)
  !! Base method for calculating entropy of a model
  !! When this method is not overridden, the entropy
  !! is simply set to 0.0. When you develop a model
  !! built on top of this abstract class or one of its
  !! children, it is recommended that you define a
  !! convex mathematical entropy function that is used
  !! as a measure of the model stability.
    implicit none
    class(Model),intent(inout) :: this

    this%entropy = 0.0_prec

  endsubroutine CalculateEntropy_Model

  subroutine ReportEntropy_Model(this)
#undef __FUNC__
#define __FUNC__ "ReportEntropy"
  !! Base method for reporting the entropy of a model
  !! to stdout. Only override this procedure if additional
  !! reporting is needed. Alternatively, if you think
  !! additional reporting would be valuable for all models,
  !! open a pull request with modifications to this base
  !! method.
    implicit none
    class(Model),intent(in) :: this
    ! Local
    character(len=20) :: modelTime
    character(len=20) :: entropy
    character(len=:),allocatable :: str

    if(this%decomp%rankId == 0) then
      ! Copy the time and entropy to a string
      write(modelTime,"(ES16.7E3)") this%t
      write(entropy,"(ES16.7E3)") this%entropy

      ! Write the output to STDOUT
      open(output_unit,ENCODING='utf-8')
      write(output_unit,'("INFO : [",A,"] : ")',ADVANCE='no') __FUNC__
      str = 'tᵢ ='//trim(modelTime)
      write(output_unit,'(A)',ADVANCE='no') str
      str = '  |  eᵢ ='//trim(entropy)
      write(output_unit,'(A)',ADVANCE='yes') str
    endif

  endsubroutine ReportEntropy_Model

  ! ////////////////////////////////////// !
  !       Time Integrators                 !

  subroutine ForwardStep_Model(this,tn,dt,ioInterval)
  !!  Forward steps the model using the associated tendency procedure and time integrator
  !!
  !!  If the final time  is provided, the model is forward stepped to that final time,
  !!  otherwise, the model is forward stepped only a single time step
  !!
  !!  If a time step is provided through the interface, the model time step size is updated
  !!  and that time step is used to update the model
  !!
  !! If ioInterval is provided, file IO will be conducted every ioInterval seconds until tn
  !! is reached
    implicit none
    class(Model),intent(inout) :: this
    real(prec),optional,intent(in) :: tn
    real(prec),optional,intent(in) :: dt
    real(prec),optional,intent(in) :: ioInterval
    ! Local
    real(prec) :: targetTime,tNext
    integer :: i,nIO

    if(present(dt)) then
      this%dt = dt
    endif

    if(present(tn)) then
      targetTime = tn
    else
      targetTime = this%t+this%dt
    endif

    if(present(ioInterval)) then
      nIO = int((targetTime-this%t)/ioInterval)
      do i = 1,nIO
        tNext = this%t+ioInterval
        call this%timeIntegrator(tNext)
        this%t = tNext
        call this%WriteModel()
        call this%WriteTecplot()
        call this%IncrementIOCounter()
        call this%CalculateEntropy()
        call this%ReportEntropy()
      enddo

    else
      call this%timeIntegrator(targetTime)
      this%t = targetTime
      call this%CalculateEntropy()
      call this%ReportEntropy()
    endif

  endsubroutine ForwardStep_Model

  subroutine Euler_timeIntegrator(this,tn)
    implicit none
    class(Model),intent(inout) :: this
    real(prec),intent(in) :: tn
    ! Local
    real(prec) :: tRemain
    real(prec) :: dtLim

    dtLim = this%dt ! Get the max time step size from the dt attribute
    do while(this%t < tn)

      tRemain = tn-this%t
      this%dt = min(dtLim,tRemain)
      call this%CalculateTendency()
      call this%UpdateSolution()
      this%t = this%t+this%dt

    enddo

    this%dt = dtLim

  endsubroutine Euler_timeIntegrator

  subroutine LowStorageRK2_timeIntegrator(this,tn)
    implicit none
    class(Model),intent(inout) :: this
    real(prec),intent(in) :: tn
    ! Local
    integer :: m
    real(prec) :: tRemain
    real(prec) :: dtLim
    real(prec) :: t0

    dtLim = this%dt ! Get the max time step size from the dt attribute
    do while(this%t < tn)

      t0 = this%t
      tRemain = tn-this%t
      this%dt = min(dtLim,tRemain)
      do m = 1,2
        call this%CalculateTendency()
        call this%UpdateGRK2(m)
        this%t = t0+rk2_b(m)*this%dt
      enddo

      this%t = t0+this%dt

    enddo

    this%dt = dtLim

  endsubroutine LowStorageRK2_timeIntegrator

  subroutine LowStorageRK3_timeIntegrator(this,tn)
    implicit none
    class(Model),intent(inout) :: this
    real(prec),intent(in) :: tn
    ! Local
    integer :: m
    real(prec) :: tRemain
    real(prec) :: dtLim
    real(prec) :: t0

    dtLim = this%dt ! Get the max time step size from the dt attribute
    do while(this%t < tn)

      t0 = this%t
      tRemain = tn-this%t
      this%dt = min(dtLim,tRemain)
      do m = 1,3
        call this%CalculateTendency()
        call this%UpdateGRK3(m)
        this%t = t0+rk3_b(m)*this%dt
      enddo

      this%t = t0+this%dt

    enddo

    this%dt = dtLim

  endsubroutine LowStorageRK3_timeIntegrator

  subroutine LowStorageRK4_timeIntegrator(this,tn)
    implicit none
    class(Model),intent(inout) :: this
    real(prec),intent(in) :: tn
    ! Local
    integer :: m
    real(prec) :: tRemain
    real(prec) :: dtLim
    real(prec) :: t0

    dtLim = this%dt ! Get the max time step size from the dt attribute
    do while(this%t < tn)

      t0 = this%t
      tRemain = tn-this%t
      this%dt = min(dtLim,tRemain)
      do m = 1,5
        call this%CalculateTendency()
        call this%UpdateGRK4(m)
        this%t = t0+rk4_b(m)*this%dt
      enddo

      this%t = t0+this%dt

    enddo

    this%dt = dtLim

  endsubroutine LowStorageRK4_timeIntegrator

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

endmodule SELF_Model
