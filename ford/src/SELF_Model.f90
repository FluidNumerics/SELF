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
    logical :: prescribed_bcs_enabled = .true.
    logical :: tecplot_enabled = .true.
    integer :: nvar
    ! Standard Diagnostics
    real(prec) :: entropy ! Mathematical entropy function for the model

  contains

    procedure :: IncrementIOCounter

    procedure :: PrintType => PrintType_Model

    procedure :: SetNumberOfVariables => SetNumberOfVariables_Model

    procedure :: AdditionalInit => AdditionalInit_Model
    procedure :: AdditionalFree => AdditionalFree_Model
    procedure :: AdditionalOutput => AdditionalOutput_Model

    procedure :: ForwardStep => ForwardStep_Model

    procedure :: Euler_timeIntegrator

    ! Runge-Kutta methods
    procedure :: LowStorageRK2_timeIntegrator
    procedure(UpdateGRK),deferred :: UpdateGRK2

    procedure :: LowStorageRK3_timeIntegrator
    procedure(UpdateGRK),deferred :: UpdateGRK3

    procedure :: LowStorageRK4_timeIntegrator
    procedure(UpdateGRK),deferred :: UpdateGRK4

    procedure :: PreTendency => PreTendency_Model
    procedure :: entropy_func => entropy_func_Model

    procedure :: flux1D => flux1d_Model
    procedure :: flux2D => flux2d_Model
    procedure :: flux3D => flux3d_Model

    procedure :: riemannflux1d => riemannflux1d_Model
    procedure :: riemannflux2d => riemannflux2d_Model
    procedure :: riemannflux3d => riemannflux3d_Model

    procedure :: source1d => source1d_Model
    procedure :: source2d => source2d_Model
    procedure :: source3d => source3d_Model

    ! Boundary condition functions (hyperbolic)
    procedure :: hbc1d_Prescribed => hbc1d_Prescribed_Model
    procedure :: hbc1d_Radiation => hbc1d_Generic_Model
    procedure :: hbc1d_NoNormalFlow => hbc1d_Generic_Model
    procedure :: hbc2d_Prescribed => hbc2d_Prescribed_Model
    procedure :: hbc2d_Radiation => hbc2d_Generic_Model
    procedure :: hbc2d_NoNormalFlow => hbc2d_Generic_Model
    procedure :: hbc3d_Prescribed => hbc3d_Prescribed_Model
    procedure :: hbc3d_Radiation => hbc3d_Generic_Model
    procedure :: hbc3d_NoNormalFlow => hbc3d_Generic_Model

    ! Boundary condition functions (parabolic)
    procedure :: pbc1d_Prescribed => pbc1d_Prescribed_Model
    procedure :: pbc1d_Radiation => pbc1d_Generic_Model
    procedure :: pbc1d_NoNormalFlow => pbc1d_Generic_Model
    procedure :: pbc2d_Prescribed => pbc2d_Prescribed_Model
    procedure :: pbc2d_Radiation => pbc2d_Generic_Model
    procedure :: pbc2d_NoNormalFlow => pbc2d_Generic_Model
    procedure :: pbc3d_Prescribed => pbc3d_Prescribed_Model
    procedure :: pbc3d_Radiation => pbc3d_Generic_Model
    procedure :: pbc3d_NoNormalFlow => pbc3d_Generic_Model

    procedure :: ReportEntropy => ReportEntropy_Model
    procedure :: ReportMetrics => ReportMetrics_Model
    procedure :: ReportUserMetrics => ReportUserMetrics_Model
    procedure :: CalculateEntropy => CalculateEntropy_Model

    procedure(UpdateSolution),deferred :: UpdateSolution
    procedure(CalculateTendency),deferred :: CalculateTendency
    procedure(ReadModel),deferred :: ReadModel
    procedure(WriteModel),deferred :: WriteModel
    procedure(WriteTecplot),deferred :: WriteTecplot

    generic :: SetTimeIntegrator => SetTimeIntegrator_withChar
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

  subroutine SetNumberOfVariables_Model(this)
    implicit none
    class(Model),intent(inout) :: this

    this%nvar = 1

  endsubroutine SetNumberOfVariables_Model

  subroutine AdditionalInit_Model(this)
    implicit none
    class(Model),intent(inout) :: this
    return
  endsubroutine AdditionalInit_Model

  subroutine AdditionalFree_Model(this)
    implicit none
    class(Model),intent(inout) :: this
    return
  endsubroutine AdditionalFree_Model

  subroutine AdditionalOutput_Model(this,fileid)
    implicit none
    class(Model),intent(inout) :: this
    integer(HID_T),intent(in) :: fileid
    return
  endsubroutine AdditionalOutput_Model

  subroutine PrintType_Model(this)
    implicit none
    class(Model),intent(in) :: this

    print*,__FILE__//" : Model : No model type"

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

  pure function entropy_func_Model(this,s) result(e)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec) :: e

    e = 0.0_prec

  endfunction entropy_func_Model

  pure function riemannflux1d_Model(this,sL,sR,dsdx,nhat) result(flux)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar)
    real(prec),intent(in) :: nhat
    real(prec) :: flux(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      flux(ivar) = 0.0_prec
    enddo

  endfunction riemannflux1d_Model

  pure function riemannflux2d_Model(this,sL,sR,dsdx,nhat) result(flux)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: flux(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      flux(ivar) = 0.0_prec
    enddo

  endfunction riemannflux2d_Model

  pure function riemannflux3d_Model(this,sL,sR,dsdx,nhat) result(flux)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: sL(1:this%nvar)
    real(prec),intent(in) :: sR(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: flux(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      flux(ivar) = 0.0_prec
    enddo

  endfunction riemannflux3d_Model

  pure function flux1d_Model(this,s,dsdx) result(flux)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar)
    real(prec) :: flux(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      flux(ivar) = 0.0_prec
    enddo

  endfunction flux1d_Model

  pure function flux2d_Model(this,s,dsdx) result(flux)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec) :: flux(1:this%nvar,1:2)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      flux(ivar,1:2) = 0.0_prec
    enddo

  endfunction flux2d_Model

  pure function flux3d_Model(this,s,dsdx) result(flux)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec) :: flux(1:this%nvar,1:3)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      flux(ivar,1:3) = 0.0_prec
    enddo

  endfunction flux3d_Model

  pure function source1d_Model(this,s,dsdx) result(source)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar)
    real(prec) :: source(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      source(ivar) = 0.0_prec
    enddo

  endfunction source1d_Model

  pure function source2d_Model(this,s,dsdx) result(source)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec) :: source(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      source(ivar) = 0.0_prec
    enddo

  endfunction source2d_Model

  pure function source3d_Model(this,s,dsdx) result(source)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec) :: source(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      source(ivar) = 0.0_prec
    enddo

  endfunction source3d_Model

  pure function hbc1d_Generic_Model(this,s,nhat) result(exts)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: nhat
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      exts(ivar) = 0.0_prec
    enddo

  endfunction hbc1d_Generic_Model

  pure function hbc1d_Prescribed_Model(this,x,t) result(exts)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: x
    real(prec),intent(in) :: t
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      exts(ivar) = 0.0_prec
    enddo

  endfunction hbc1d_Prescribed_Model

  pure function hbc2d_Generic_Model(this,s,nhat) result(exts)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      exts(ivar) = 0.0_prec
    enddo

  endfunction hbc2d_Generic_Model

  pure function hbc2d_Prescribed_Model(this,x,t) result(exts)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: x(1:2)
    real(prec),intent(in) :: t
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      exts(ivar) = 0.0_prec
    enddo

  endfunction hbc2d_Prescribed_Model

  pure function hbc3d_Generic_Model(this,s,nhat) result(exts)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: s(1:this%nvar)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      exts(ivar) = 0.0_prec
    enddo

  endfunction hbc3d_Generic_Model

  pure function hbc3d_Prescribed_Model(this,x,t) result(exts)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: x(1:3)
    real(prec),intent(in) :: t
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      exts(ivar) = 0.0_prec
    enddo

  endfunction hbc3d_Prescribed_Model

  pure function pbc1d_Generic_Model(this,dsdx,nhat) result(extDsdx)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: dsdx(1:this%nvar)
    real(prec),intent(in) :: nhat
    real(prec) :: extDsdx(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      extDsdx(ivar) = dsdx(ivar)
    enddo

  endfunction pbc1d_Generic_Model

  pure function pbc1d_Prescribed_Model(this,x,t) result(extDsdx)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: x
    real(prec),intent(in) :: t
    real(prec) :: extDsdx(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      extDsdx(ivar) = 0.0_prec
    enddo

  endfunction pbc1d_Prescribed_Model

  pure function pbc2d_Generic_Model(this,dsdx,nhat) result(extDsdx)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: dsdx(1:this%nvar,1:2)
    real(prec),intent(in) :: nhat(1:2)
    real(prec) :: extDsdx(1:this%nvar,1:2)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      extDsdx(ivar,1:2) = dsdx(ivar,1:2)
    enddo

  endfunction pbc2d_Generic_Model

  pure function pbc2d_Prescribed_Model(this,x,t) result(extDsdx)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: x(1:2)
    real(prec),intent(in) :: t
    real(prec) :: extDsdx(1:this%nvar,1:2)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      extDsdx(ivar,1:2) = 0.0_prec
    enddo

  endfunction pbc2d_Prescribed_Model

  pure function pbc3d_Generic_Model(this,dsdx,nhat) result(extDsdx)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: dsdx(1:this%nvar,1:3)
    real(prec),intent(in) :: nhat(1:3)
    real(prec) :: extDsdx(1:this%nvar,1:3)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      extDsdx(ivar,1:3) = dsdx(ivar,1:3)
    enddo

  endfunction pbc3d_Generic_Model

  pure function pbc3d_Prescribed_Model(this,x,t) result(extDsdx)
    class(Model),intent(in) :: this
    real(prec),intent(in) :: x(1:3)
    real(prec),intent(in) :: t
    real(prec) :: extDsdx(1:this%nvar,1:3)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      extDsdx(ivar,1:3) = 0.0_prec
    enddo

  endfunction pbc3d_Prescribed_Model

  subroutine SetTimeIntegrator_withChar(this,integrator)
    !! Sets the time integrator method, using a character input
    !!
    !! Valid options for integrator are
    !!
    !!   "euler"
    !!   "rk2"
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

    ! Copy the time and entropy to a string
    write(modelTime,"(ES16.7E3)") this%t
    write(entropy,"(ES16.7E3)") this%entropy

    ! Write the output to STDOUT
    open(output_unit,ENCODING='utf-8')
    write(output_unit,'(1x,A," : ")',ADVANCE='no') __FILE__
    str = 'tᵢ ='//trim(modelTime)
    write(output_unit,'(A)',ADVANCE='no') str
    str = '  |  eᵢ ='//trim(entropy)
    write(output_unit,'(A)',ADVANCE='yes') str

  endsubroutine ReportEntropy_Model

  subroutine ReportMetrics_Model(this)
      !! Method that can be overridden by users to
      !! report their own custom metrics after file io
    implicit none
    class(Model),intent(inout) :: this
    return
  endsubroutine ReportMetrics_Model

  subroutine ReportUserMetrics_Model(this)
    !! Method that can be overridden by users to
    !! report their own custom metrics after file io
    implicit none
    class(Model),intent(inout) :: this
    return
  endsubroutine ReportUserMetrics_Model

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
    real(prec),intent(in) :: tn
    real(prec),intent(in) :: dt
    real(prec),intent(in) :: ioInterval
    ! Local
    real(prec) :: targetTime,tNext
    integer :: i,nIO
    character(10) :: ntimesteps
    real(prec) :: t1,t2
    character(len=:),allocatable :: str
    character(len=20) :: modelTime

    this%dt = dt
    targetTime = tn

    write(ntimesteps,"(I10)") int(ioInterval/this%dt)
    nIO = int((targetTime-this%t)/ioInterval)
    do i = 1,nIO

      tNext = this%t+ioInterval
      call cpu_time(t1)
      call this%timeIntegrator(tNext)
      call cpu_time(t2)

      open(output_unit,ENCODING='utf-8')
      write(output_unit,'(A)',ADVANCE='no') ' --------------------------------------------------------'
      write(output_unit,'(A)',ADVANCE='yes') '----------------------------------------------------------------------- '
      write(output_unit,'(A)',ADVANCE='no') ' <><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
      write(output_unit,'(A)',ADVANCE='yes') '<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> '
      write(output_unit,'(A)',ADVANCE='no') ' --------------------------------------------------------'
      write(output_unit,'(A)',ADVANCE='yes') '------------------------------------------------------------------------ '
      write(modelTime,"(ES16.7E3)") this%t
      ! Write the wall-time
      write(output_unit,'(1x, A," : ")',ADVANCE='no') __FILE__
      str = 'tᵢ ='//trim(modelTime)
      write(output_unit,'(A)',ADVANCE='no') str
      str = '  | Time to complete '//trim(ntimesteps)//' time steps (s) : '
      write(output_unit,'(A)',ADVANCE='no') str
      write(str,"(ES16.7E3)") t2-t1
      write(output_unit,'(A)',ADVANCE='yes') str

      ! Wall-time per time step
      write(output_unit,'(1x, A," : ")',ADVANCE='no') __FILE__
      str = 'tᵢ ='//trim(modelTime)
      write(output_unit,'(A)',ADVANCE='no') str
      str = '  | Wall-time per time step : '
      write(output_unit,'(A)',ADVANCE='no') str
      write(str,"(ES16.7E3)")(t2-t1)/(floor(ioInterval/this%dt))
      write(output_unit,'(A)',ADVANCE='yes') str

      ! Wall-time per simulation time
      write(output_unit,'(1x, A," : ")',ADVANCE='no') __FILE__
      str = 'tᵢ ='//trim(modelTime)
      write(output_unit,'(A)',ADVANCE='no') str
      str = '  | Wall-time per simulation time : '
      write(output_unit,'(A)',ADVANCE='no') str
      write(str,"(ES16.7E3)")(t2-t1)/(ioInterval)
      write(output_unit,'(A)',ADVANCE='yes') str

      this%t = tNext

      call this%CalculateEntropy()
      call this%ReportEntropy()
      call this%ReportMetrics()

      call this%WriteModel()
      if(this%tecplot_enabled) then
        call this%WriteTecplot()
      endif
      call this%IncrementIOCounter()

    enddo

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

endmodule SELF_Model
