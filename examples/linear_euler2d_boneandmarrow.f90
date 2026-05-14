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
module lineareuler2d_boneandmarrow_model
!! Heterogeneous-medium linear-acoustics test on the
!! share/mesh/MultiMaterial2D/BoneAndMarrow.mesh ISM-MM mesh.
!!
!! The mesh tags every element with one of three materials:
!! "Muscle" (background annulus, r in roughly [6.5, 14]),
!! "Bone" (disk of radius 6.5 about (-0.5, 0)), and
!! "Marrow" (smaller disk of radius 1.5 about (-1.5, 0.5))
!! nested inside the bone region. We map each material to a
!! representative sound speed and write that into solution(...,5),
!! which is held fixed in time by the LinearEuler2D model.
!!
!! Initial condition: a small Gaussian pressure / density bump is
!! placed in the Muscle region, well outside the bone region. The
!! transient acoustic pulse propagates outward, refracts/reflects
!! at the material interfaces (because c is discontinuous across
!! them), and radiates out through the outer boundary.

  use self_lineareuler2d
  use SELF_BoundaryConditions

  implicit none

  type,extends(lineareuler2d) :: lineareuler2d_boneandmarrow
    ! Sound speeds in normalized units, preserving real-tissue ratios
    ! (bone is ~2.3x faster than muscle, marrow is ~0.92x).
    real(prec) :: c_muscle = 1.0_prec
    real(prec) :: c_bone = 2.3_prec
    real(prec) :: c_marrow = 0.92_prec
    real(prec) :: bump_x0 = 10.0_prec ! Pulse center, well into the muscle annulus
    real(prec) :: bump_y0 = 0.0_prec
    real(prec) :: bump_L = 0.6_prec ! Halfwidth (e-folding length) of the bump
    real(prec) :: bump_amp = 1.0e-3_prec ! Density-perturbation amplitude

  contains

    procedure :: setInitialCondition
    procedure :: AdditionalInit => AdditionalInit_lineareuler2d_boneandmarrow

  endtype lineareuler2d_boneandmarrow

contains

  subroutine AdditionalInit_lineareuler2d_boneandmarrow(this)
    !! Registers the radiation BC for the CPU path. The base class
    !! `_t` only wires up no_normal_flow; we want all four sides of
    !! this disk-shaped domain to behave as transparent outflow.
    implicit none
    class(lineareuler2d_boneandmarrow),intent(inout) :: this
    procedure(SELF_bcMethod),pointer :: bcfunc

    call AdditionalInit_LinearEuler2D_t(this)

    bcfunc => hbc2d_Radiation_lineareuler2d_boneandmarrow
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      SELF_BC_RADIATION,"radiation",bcfunc)

  endsubroutine AdditionalInit_lineareuler2d_boneandmarrow

  subroutine hbc2d_Radiation_lineareuler2d_boneandmarrow(bc,mymodel)
    !! Radiation (zero-state) BC for the four prognostic variables.
    !! The sound speed (variable 5) is copied through from the
    !! interior so that the face Riemann flux sees a consistent c.
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    integer :: n,i,iEl,j

    select type(m => mymodel)
    class is(lineareuler2d_boneandmarrow)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        j = bc%sides(n)
        do i = 1,m%solution%interp%N+1
          m%solution%extBoundary(i,j,iEl,1) = 0.0_prec ! rho
          m%solution%extBoundary(i,j,iEl,2) = 0.0_prec ! u
          m%solution%extBoundary(i,j,iEl,3) = 0.0_prec ! v
          m%solution%extBoundary(i,j,iEl,4) = 0.0_prec ! p
          m%solution%extBoundary(i,j,iEl,5) = m%solution%boundary(i,j,iEl,5) ! c
        enddo
      enddo
    endselect

  endsubroutine hbc2d_Radiation_lineareuler2d_boneandmarrow

  subroutine setInitialCondition(this)
    !! Material-aware initial condition: stamp `c` from the per-
    !! element material id, then add a Gaussian rho/p bump only in
    !! Muscle elements (so the pulse starts cleanly in a uniform
    !! medium and the bone/marrow inclusions are seen as scatterers).
    implicit none
    class(lineareuler2d_boneandmarrow),intent(inout) :: this
    integer :: i,j,iel,matid
    real(prec) :: c_mat,x,y,r2,shape
    character(LEN=64) :: matname

    do iel = 1,this%mesh%nElem
      matid = this%mesh%elemMaterial(iel)
      matname = this%mesh%materialNames(matid)
      select case(trim(matname))
      case("Muscle"); c_mat = this%c_muscle
      case("Bone"); c_mat = this%c_bone
      case("Marrow"); c_mat = this%c_marrow
      case default; c_mat = this%c_muscle
      endselect

      do j = 1,this%solution%N+1
        do i = 1,this%solution%N+1
          x = this%geometry%x%interior(i,j,iel,1,1)
          y = this%geometry%x%interior(i,j,iel,1,2)

          if(trim(matname) == "Muscle") then
            r2 = (x-this%bump_x0)**2+(y-this%bump_y0)**2
            shape = this%bump_amp*exp(-r2/(this%bump_L*this%bump_L))
          else
            shape = 0.0_prec
          endif

          this%solution%interior(i,j,iel,1) = shape ! density perturbation
          this%solution%interior(i,j,iel,2) = 0.0_prec ! u
          this%solution%interior(i,j,iel,3) = 0.0_prec ! v
          this%solution%interior(i,j,iel,4) = shape*c_mat*c_mat ! p = rho * c^2 (acoustic)
          this%solution%interior(i,j,iel,5) = c_mat ! sound speed for this material
        enddo
      enddo
    enddo

    call this%solution%UpdateDevice()

  endsubroutine setInitialCondition

endmodule lineareuler2d_boneandmarrow_model

program LinearEuler_BoneAndMarrow

  use self_data
  use lineareuler2d_boneandmarrow_model

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  ! CFL for DGSEM: dt < dx_min / (c_max * (N+1)^2). With curved elements
  ! near the material interfaces dx_min can be much smaller than the
  ! background 1.0, so we run conservatively.
  real(prec),parameter :: dt = 1.0e-3_prec
  real(prec),parameter :: endtime = 0.5_prec
  real(prec),parameter :: iointerval = 0.1_prec
  real(prec) :: e0,ef
  type(lineareuler2d_boneandmarrow) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  character(LEN=255) :: WORKSPACE

  ! Read the multi-material ISM-MM mesh. Boundary names from the
  ! .mesh file (e.g. "outer") become bc indices 1..nBCs; remap them
  ! all to SELF_BC_RADIATION since the disk has one transparent
  ! outer boundary.
  call get_environment_variable("WORKSPACE",WORKSPACE)
  call mesh%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/MultiMaterial2D/BoneAndMarrow.mesh")
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%tecplot_enabled = .false.
  modelobj%rho0 = 1.0_prec

  call modelobj%setInitialCondition()

  call modelobj%WriteModel()
  call modelobj%IncrementIOCounter()

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
  if(ef > 10.0_prec*e0) then
    print*,"Error: Final entropy grew unphysically (e0=",e0," ef=",ef,")"
    stop 1
  endif

  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram LinearEuler_BoneAndMarrow
