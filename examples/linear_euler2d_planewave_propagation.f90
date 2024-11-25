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
module lineareuler2d_planewave_prop_model
!! This module can be used for simulating plane wave propagation in a 2-D domain
!! We use a type extension of the linearEuler2D class to add parameters for a simple
!! plane wave solution, including the x and y wave numbers and the peak pressure.
!! This information, combined with the sound speed and reference density, is used
!! to compute initial and boundary conditions for a plane wave.
!!
!! Any boundary that is set as a prescribed boundary condition will be prescribed
!! the exact plane wave solution at each model time step.
!!
!! See section 5.4.3 of D. A. Kopriva, "Implementing Spectral Methods for Partial Differential Equations" (2009)

  use self_lineareuler2d

  implicit none

  type,extends(lineareuler2d) :: lineareuler2d_planewave
    real(prec) :: wx = 0.0_prec ! Wave number in the x-direction
    real(prec) :: wy = 0.0_prec ! Wave number in the y-direction
    real(prec) :: p = 0.0_prec ! Peak pressure amplitude
    real(prec) :: x0 = 0.0_prec ! x component of the wave center position
    real(prec) :: y0 = 0.0_prec ! y component of the wave center position
    real(prec) :: L = 1.0_prec ! Halfwidth of the plane wave

  contains

    procedure :: setInitialCondition
    procedure :: hbc2d_Prescribed => hbc2d_Prescribed_lineareuler2d_planewave ! override for the hyperbolic boundary conditions

  endtype lineareuler2d_planewave

contains

  subroutine setInitialCondition(this)
    implicit none
    class(lineareuler2d_planewave),intent(inout) :: this
    ! Local
    integer :: i,j,iel
    real(prec) :: p,rho,u,v,x,y,phase,shape

    p = this%p
    rho = this%p/this%c/this%c
    u = this%p*this%wx/this%c
    v = this%p*this%wy/this%c

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)

      x = this%geometry%x%interior(i,j,iel,1,1)
      y = this%geometry%x%interior(i,j,iel,1,2)
      phase = this%wx*(x-this%x0)+this%wy*(y-this%y0)-this%c*this%t
      shape = exp(-phase*phase/(this%L*this%L))

      this%solution%interior(i,j,iel,1) = rho*shape ! density
      this%solution%interior(i,j,iel,2) = u*shape ! u
      this%solution%interior(i,j,iel,3) = v*shape ! v
      this%solution%interior(i,j,iel,4) = p*shape ! pressure

    enddo

    call this%solution%UpdateDevice()

  endsubroutine setInitialCondition

  pure function hbc2d_Prescribed_lineareuler2d_planewave(this,x,t) result(exts)
    class(lineareuler2d_planewave),intent(in) :: this
    real(prec),intent(in) :: x(1:2)
    real(prec),intent(in) :: t
    real(prec) :: exts(1:this%nvar)
    ! Local
    real(prec) :: p,rho,u,v,phase,shape

    p = this%p
    rho = this%p/this%c/this%c
    u = this%p*this%wx/this%c
    v = this%p*this%wy/this%c

    phase = this%wx*(x(1)-this%x0)+this%wy*(x(2)-this%y0)-this%c*this%t
    shape = exp(-phase*phase/(this%L*this%L))

    exts(1) = rho*shape ! density
    exts(2) = u*shape ! u
    exts(3) = v*shape ! v
    exts(4) = p*shape ! pressure

  endfunction hbc2d_Prescribed_lineareuler2d_planewave

endmodule lineareuler2d_planewave_prop_model

program LinearEuler_Example

  use self_data
  use lineareuler2d_planewave_prop_model

  implicit none
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 15
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-4) ! time-step size
  real(prec),parameter :: endtime = 0.05_prec
  real(prec),parameter :: iointerval = 0.05_prec
  real(prec) :: e0,ef ! Initial and final entropy
  type(lineareuler2d_planewave) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  ! Create a structured mesh
  bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                SELF_BC_PRESCRIBED, & ! East
                SELF_BC_PRESCRIBED, & ! North
                SELF_BC_PRESCRIBED] ! West

  call mesh%StructuredMesh(10,10,2,2,0.05_prec,0.05_prec,bcids)

  ! Create an interpolant
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize the model
  call modelobj%Init(mesh,geometry)

  ! Set the plane wave parameters
  modelobj%x0 = 0.2_prec
  modelobj%y0 = 0.2_prec
  modelobj%L = 0.2_prec/(2.0_prec*sqrt(log(2.0_prec)))
  modelobj%wx = sqrt(2.0_prec)/2.0_prec
  modelobj%wy = sqrt(2.0_prec)/2.0_prec
  modelobj%p = 10.0_prec**(-4)

  ! Set the initial condition
  call modelobj%setInitialCondition()

  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  e0 = modelobj%entropy
  ! Set the model's time integration method
  call modelobj%SetTimeIntegrator(integrator)

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy

  if(ef /= ef) then
    print*,"Error: Final entropy is not finite or not a number",ef
    stop 1
  endif

  ! Clean up
  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram LinearEuler_Example
