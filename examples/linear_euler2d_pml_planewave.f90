! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2026 Fluid Numerics LLC
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
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
module lineareuler2d_pml_planewave_model
!! Right-going Gaussian acoustic pulse on a 2D domain with a perfectly
!! matched layer (Hu 2001, unsplit form) on the east side. The pulse
!! propagates into the PML region and is absorbed without producing a
!! visible reflection back into the interior.

  use self_lineareuler2d_pml

  implicit none

  type,extends(LinearEuler2D_PML) :: lineareuler2d_pml_planewave
    real(prec) :: c = 1.0_prec ! Reference sound speed
    real(prec) :: amp = 1.0e-3_prec ! Peak pressure amplitude
    real(prec) :: x0 = 1.0_prec ! Pulse centre (x)
    real(prec) :: y0 = 1.0_prec ! Pulse centre (y)
    real(prec) :: L = 0.25_prec ! Gaussian halfwidth
  contains
    procedure :: setInitialCondition
  endtype lineareuler2d_pml_planewave

contains

  subroutine setInitialCondition(this)
    !! Right-going acoustic wave packet:
    !!   p(x,y,0) = amp * exp(-((x-x0)^2 + (y-y0)^2)/L^2)
    !!   u        = p/(rho0*c)
    !!   v        = 0
    !! consistent with a planewave travelling in +x at sound speed c.
    implicit none
    class(lineareuler2d_pml_planewave),intent(inout) :: this
    ! Local
    integer :: i,j,iel
    real(prec) :: x,y,r2,p_loc,u_loc

    do concurrent(i=1:this%solution%N+1,j=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem)
      x = this%geometry%x%interior(i,j,iel,1,1)
      y = this%geometry%x%interior(i,j,iel,1,2)
      r2 = (x-this%x0)**2+(y-this%y0)**2
      p_loc = this%amp*exp(-r2/(this%L*this%L))
      u_loc = p_loc/(this%rho0*this%c)

      this%solution%interior(i,j,iel,1) = u_loc
      this%solution%interior(i,j,iel,2) = 0.0_prec
      this%solution%interior(i,j,iel,3) = p_loc
      this%solution%interior(i,j,iel,4) = this%c
      ! PML auxiliaries start at zero
      this%solution%interior(i,j,iel,5:7) = 0.0_prec
    enddo

    call this%solution%UpdateDevice()

  endsubroutine setInitialCondition

  subroutine tagPMLElements(mesh,x_interior_max)
    !! Tag mesh elements whose centroid lies east of x_interior_max
    !! with material name "pml". All other elements remain "interior".
    !! Centroid is the corner-node average (linear elements: nGeo=1).
    use SELF_Mesh_2D
    implicit none
    class(Mesh2D),intent(inout) :: mesh
    real(prec),intent(in) :: x_interior_max
    ! Local
    integer :: iel
    real(prec) :: xc

    if(allocated(mesh%materialNames)) deallocate(mesh%materialNames)
    mesh%nMaterials = 2
    allocate(mesh%materialNames(1:2))
    mesh%materialNames(1) = "interior"
    mesh%materialNames(2) = "pml"

    do iel = 1,mesh%nElem
      xc = sum(mesh%nodeCoords(1,:,:,iel))/real(size(mesh%nodeCoords,2)*size(mesh%nodeCoords,3),prec)
      if(xc > x_interior_max) then
        mesh%elemMaterial(iel) = 2 ! pml
      else
        mesh%elemMaterial(iel) = 1 ! interior
      endif
    enddo

  endsubroutine tagPMLElements

endmodule lineareuler2d_pml_planewave_model

program LinearEuler_PML_Example

  use self_data
  use lineareuler2d_pml_planewave_model

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 10
  real(prec),parameter :: dt = 5.0e-3_prec
  real(prec),parameter :: endtime = 4.0_prec
  real(prec),parameter :: iointerval = 0.5_prec

  ! Interior + PML layout (one element = 0.5 wide, total 12 elements in x)
  !  Interior: x in [0, 4]   (8 elements)
  !  PML east: x in [4, 6]   (4 elements, width 2)
  real(prec),parameter :: x_interior_max = 4.0_prec
  real(prec),parameter :: pml_width = 2.0_prec
  real(prec),parameter :: sigma_max = 6.0_prec

  type(lineareuler2d_pml_planewave) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)
  real(prec) :: e0,ef

  ! Structured mesh: 12 elements in x, 4 in y, all 0.5 wide.
  bcids(1:4) = [SELF_BC_NONORMALFLOW, & ! south
                SELF_BC_NONORMALFLOW, & ! east  (behind the PML)
                SELF_BC_NONORMALFLOW, & ! north
                SELF_BC_NONORMALFLOW] ! west
  call mesh%StructuredMesh(12,4,1,1,0.5_prec,0.5_prec,bcids)

  ! Tag the eastern strip as the "pml" material.
  call tagPMLElements(mesh,x_interior_max)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%rho0 = 1.0_prec

  ! Configure the PML profile: damping ramps from 0 at x=x_interior_max
  ! up to sigma_max at x=x_interior_max + pml_width with a cubic ramp.
  call modelobj%SetPMLProfile(x_interior_min=0.0_prec, &
                              x_interior_max=x_interior_max, &
                              y_interior_min=0.0_prec, &
                              y_interior_max=2.0_prec, &
                              pml_width=pml_width, &
                              sigma_max=sigma_max, &
                              ramp_exponent=3)

  call modelobj%setInitialCondition()

  call modelobj%CalculateEntropy()
  call modelobj%ReportEntropy()
  e0 = modelobj%entropy

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(endtime,dt,iointerval)
  ef = modelobj%entropy

  if(ef /= ef) then
    print*,"Error: Final entropy is inf or nan",ef
    stop 1
  endif

  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram LinearEuler_PML_Example
