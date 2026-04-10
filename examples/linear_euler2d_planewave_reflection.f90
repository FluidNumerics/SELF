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
module lineareuler2d_planewave_model
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
!!
!! The exact solution is derived using method of images, assuming we have a no-normal-flow boundary condition
!! on the "eastern" boundary of the domain.

  use self_lineareuler2d
  use SELF_BoundaryConditions

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
    procedure :: AdditionalInit => AdditionalInit_lineareuler2d_planewave

  endtype lineareuler2d_planewave

contains

  subroutine AdditionalInit_lineareuler2d_planewave(this)
    implicit none
    class(lineareuler2d_planewave),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Register the parent class NoNormalFlow BC first
    call AdditionalInit_LinearEuler2D_t(this)

    ! Register the prescribed BC
    bcfunc => hbc2d_Prescribed_lineareuler2d_planewave
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_PRESCRIBED,"prescribed",bcfunc)

  endsubroutine AdditionalInit_lineareuler2d_planewave

  subroutine setInitialCondition(this)
    implicit none
    class(lineareuler2d_planewave),intent(inout) :: this
    ! Local
    integer :: i,j,iel
    real(prec) :: p,rho,u,v,x,y
    real(prec) :: phi,phr,shi,shr

    p = this%p
    rho = this%p/this%c/this%c
    u = this%p*this%wx/this%c
    v = this%p*this%wy/this%c

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)

      x = this%geometry%x%interior(i,j,iel,1,1)
      y = this%geometry%x%interior(i,j,iel,1,2)

      ! Incident wabe
      phi = this%wx*(x-this%x0)+this%wy*(y-this%y0)-this%c*this%t
      shi = exp(-phi*phi/(this%L*this%L))

      ! Reflected wave
      phr = -this%wx*(x-(2.0_prec-this%x0))+this%wy*(y-this%y0)-this%c*this%t
      shr = exp(-phr*phr/(this%L*this%L))

      this%solution%interior(i,j,iel,1) = rho*(shi+shr) ! density
      this%solution%interior(i,j,iel,2) = u*(shi-shr) ! u
      this%solution%interior(i,j,iel,3) = v*(shi+shr) ! v
      this%solution%interior(i,j,iel,4) = p*(shi+shr) ! pressure

    enddo

    call this%solution%UpdateDevice()

  endsubroutine setInitialCondition

  subroutine hbc2d_Prescribed_lineareuler2d_planewave(bc,mymodel)
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,i,iEl,j
    real(prec) :: x(1:2)
    real(prec) :: p,rho,u,v,phase,shi,shr

    select type(m => mymodel)
    class is(lineareuler2d_planewave)
      p = m%p
      rho = m%p/m%c/m%c
      u = m%p*m%wx/m%c
      v = m%p*m%wy/m%c

      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        j = bc%sides(n)
        do i = 1,m%solution%interp%N+1
          x = m%geometry%x%boundary(i,j,iEl,1,1:2)

          ! Incident wave
          phase = m%wx*(x(1)-m%x0)+m%wy*(x(2)-m%y0)-m%c*m%t
          shi = exp(-phase*phase/(m%L*m%L))

          ! Reflected wave
          phase = -m%wx*(x(1)+m%x0-2.0_prec)+m%wy*(x(2)-m%y0)-m%c*m%t
          shr = exp(-phase*phase/(m%L*m%L))

          m%solution%extBoundary(i,j,iEl,1) = rho*(shi+shr) ! density
          m%solution%extBoundary(i,j,iEl,2) = u*(shi-shr) ! u
          m%solution%extBoundary(i,j,iEl,3) = v*(shi+shr) ! v
          m%solution%extBoundary(i,j,iEl,4) = p*(shi+shr) ! pressure
        enddo
      enddo
    endselect

  endsubroutine hbc2d_Prescribed_lineareuler2d_planewave

endmodule lineareuler2d_planewave_model

program LinearEuler_Example

  use self_data
  use lineareuler2d_planewave_model

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
                SELF_BC_NONORMALFLOW, & ! East
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
