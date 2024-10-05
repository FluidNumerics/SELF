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

module burgers1d_shock_model
!! This module can be used for simulating stationary and traveling shocks in 1-D
!! that obey burger's equation. In this module, we provide a class that is a type
!! extension of the self_Burgers1D class to provide prescribed boundary conditions
!! to produce a traveling or stationary shock.
!!
!! The initial and boundary conditions are set using an exact solution,
!!
!! \begin{equation}
!! u(x,t) = s - 0.5(u_l-u_r) \tanh\left( \frac{(x - st - x_0)(u_l-u_r)}{4\nu} \right)
!! \end{equation}
!!
!! where $s = \frac{u_l + u_r}{2}$ is the shock speed.

use self_Burgers1D

implicit none

  type, extends(burgers1D) :: burgers1d_shock
    real(prec) :: ul = 1.0_prec
    real(prec) :: ur = 0.0_prec
    real(prec) :: x0 = 0.1_prec

    contains

    procedure :: hbc1d_Prescribed => hbc1d_Prescribed_burgers1d_shock ! override for the hyperbolic boundary conditions
    procedure :: pbc1d_Prescribed => pbc1d_Prescribed_burgers1d_shock ! override for the parabolic boundary conditions

  endtype burgers1d_shock

  contains

  pure function hbc1d_Prescribed_burgers1d_shock(this,x,t) result(exts)
  class(burgers1d_shock),intent(in) :: this
  real(prec),intent(in) :: x
  real(prec),intent(in) :: t
  real(prec) :: exts(1:this%nvar)
  ! Local
  real(prec) :: jump, s

    jump = this%ul - this%ur
    s = 0.5_prec*(this%ul + this%ur)
    exts(1) = s - 0.5_prec*tanh( (x-s*t-this%x0)*jump/(4.0_prec*this%nu) )

  endfunction hbc1d_Prescribed_burgers1d_shock

  pure function pbc1d_Prescribed_burgers1d_shock(this,x,t) result(extDsdx)
  class(burgers1d_shock),intent(in) :: this
  real(prec),intent(in) :: x
  real(prec),intent(in) :: t
  real(prec) :: extDsdx(1:this%nvar)
  ! Local
  real(prec) :: jump, s, r, drdx

    jump = this%ul - this%ur
    s = 0.5_prec*(this%ul + this%ur)
    r = (x-s*t-this%x0)*jump/(4.0_prec*this%nu)
    drdx = 1.0_prec/(4.0_prec*this%nu)
    extDsdx(1) = -0.5_prec*drdx*( sech( r ) )**2

  endfunction pbc1d_Prescribed_burgers1d_shock

  pure real(prec) function sech(x) result(fx)
  real(prec),intent(in) :: x
  fx = 2.0_prec/(exp(x)+exp(-x))
  endfunction
 
endmodule burgers1d_shock_model

program traveling_shock

    use self_data
    use burgers1d_shock_model
  
    implicit none
    character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
    integer,parameter :: nvar = 1
    integer,parameter :: nelem = 10
    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 10
    real(prec),parameter :: nu = 0.01_prec ! diffusivity
    real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-5) ! time-step size
    real(prec),parameter :: endtime = 2.0_prec
    real(prec),parameter :: iointerval = 0.05_prec
    type(burgers1d_shock) :: modelobj
    type(Lagrange),target :: interp
    type(Mesh1D),target :: mesh
    type(Geometry1D),target :: geometry
    real(prec) :: jump,s
  
    ! Create a mesh using the built-in
    ! uniform mesh generator.
    ! The domain is set to x in [0,1]
    ! We use `nelem` elements
    call mesh%UniformBlockMesh(nGeo=1, &
                               nElem=nelem, &
                               x=(/0.0_prec,1.0_prec/))

    ! Set the left and right boundary conditions to prescribed                               
    call mesh%ResetBoundaryConditionType(SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED)

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)
  
    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)
  
    ! Initialize the model
    call modelobj%Init(nvar,mesh,geometry)
    modelobj%gradient_enabled = .true.
    !Set the diffusivity
    modelobj%nu = nu
  
    ! Set the initial condition
    jump = modelobj%ul - modelobj%ur
    s = 0.5_prec*(modelobj%ul + modelobj%ur)
    modelobj%solution%interior(:,:,1) = s - 0.5_prec*tanh( &
      (geometry%x%interior(:,:,1)-modelobj%x0)*jump/(4.0_prec*modelobj%nu) )
  
    print*,"min, max (interior)", &
      minval(modelobj%solution%interior), &
      maxval(modelobj%solution%interior)
  
    call modelobj%CalculateEntropy()
    call modelobj%ReportEntropy()

    !Write the initial condition
    call modelobj%WriteModel()
    call modelobj%WriteTecplot()
    call modelobj%IncrementIOCounter()
  
    ! Set the model's time integration method
    call modelobj%SetTimeIntegrator(integrator)
  
    ! forward step the model to `endtime` using a time step
    ! of `dt` and outputing model data every `iointerval`
    call modelobj%ForwardStep(endtime,dt,iointerval)
    
    print*,"min, max (interior)", &
      minval(modelobj%solution%interior), &
      maxval(modelobj%solution%interior)
  
    ! ef = modelobj%entropy
  
    ! if(ef > e0) then
    !   print*,"Error: Final entropy greater than initial entropy! ",e0,ef
    !   stop 1
    ! endif
  
    ! Clean up
    call modelobj%free()
    call mesh%free()
    call geometry%free()
    call interp%free()
  
  endprogram traveling_shock
  