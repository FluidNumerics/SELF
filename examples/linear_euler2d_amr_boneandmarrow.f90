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

module lineareuler2d_amr_boneandmarrow_model
!! Adaptive mesh refinement on a curved, unstructured, multi-material mesh: the BoneAndMarrow
!! HOHQMesh disk (radius 14, 973 elements, nGeo = 4) with three media of distinct sound
!! speeds and densities:
!!
!!   Muscle (background annulus)             : c = 1.00, rho0 = 1.00
!!   Bone   (disk, radius 6.5 at (-0.5, 0))  : c = 2.30, rho0 = 1.80
!!   Marrow (disk, radius 1.5 at (-1.5,0.5)) : c = 0.92, rho0 = 0.93
!!
!! (normalized units preserving real-tissue ratios). The acoustic impedance Z = rho0*c jumps
!! by ~4.1x at the muscle/bone interface, so the incident wave splits into a strong reflection
!! and a fast transmitted front.
!!
!! A Gaussian pressure pulse is released at the 3 o'clock position (x,y) = (10,0) in the
!! muscle annulus. Travel times from the source (leading edge): first bone contact at
!! t ~ 3.6-4.0, the transmitted front reaches the marrow at t ~ 6.6 and exits the far side of
!! the bone at t ~ 9.7. The AMR controller tracks all of it with the modal-decay indicator on
!! pressure: the refined band follows the incident annulus, then splits to follow the
!! reflected front in muscle and the transmitted front in bone - where the same wave is 2.3x
!! longer and the indicator refines correspondingly less.
!!
!! Beyond the physics, this exercises the AMR machinery on a base mesh that structured tests
!! cannot: exact isoparametric subdivision of curved (nGeo = 4) elements, and rotated face
!! adjacencies (base flip = 1) in the forest's neighbour search and mortar emission.

  use self_lineareuler2d

  implicit none

  type,extends(lineareuler2d) :: lineareuler2d_bone_amr
    ! Distinct sound speeds per material, normalized to muscle (real-tissue ratios: cortical
    ! bone is ~2.3x faster than muscle, marrow ~0.92x).
    real(prec) :: c_muscle = 1.0_prec
    real(prec) :: c_bone = 2.3_prec
    real(prec) :: c_marrow = 0.92_prec
    ! Distinct background densities (bone ~1.8x muscle, marrow ~0.93x).
    real(prec) :: rho0_muscle = 1.0_prec
    real(prec) :: rho0_bone = 1.8_prec
    real(prec) :: rho0_marrow = 0.93_prec
    real(prec) :: bump_x0 = 10.0_prec ! 3 o'clock position in the muscle annulus
    real(prec) :: bump_y0 = 0.0_prec
    real(prec) :: bump_L = 0.4_prec ! pulse half-width (e-folding length)
    real(prec) :: bump_amp = 1.0e-3_prec ! pressure amplitude in units of rho'*c^2

  contains
    procedure :: setInitialCondition

  endtype lineareuler2d_bone_amr

contains

  subroutine setInitialCondition(this)
    !! Material-aware initial condition: stamp c (solution variable 4) and rho0 (variable 5)
    !! from each element's material id and add the Gaussian pressure bump in Muscle elements
    !! only. Works unchanged on any emitted AMR mesh because leaf elements inherit their root
    !! element's material id, and children of a common parent never straddle a material
    !! interface (materials are element-aligned on the base mesh).
    implicit none
    class(lineareuler2d_bone_amr),intent(inout) :: this
    integer :: i,j,iel,matid
    real(prec) :: c_mat,rho0_mat,x,y,r2,shape
    character(LEN=64) :: matname

    do iel = 1,this%mesh%nElem
      matid = this%mesh%elemMaterial(iel)
      matname = this%mesh%materialNames(matid)
      select case(trim(matname))
      case("Muscle"); c_mat = this%c_muscle; rho0_mat = this%rho0_muscle
      case("Bone"); c_mat = this%c_bone; rho0_mat = this%rho0_bone
      case("Marrow"); c_mat = this%c_marrow; rho0_mat = this%rho0_marrow
      case default; c_mat = this%c_muscle; rho0_mat = this%rho0_muscle
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
          this%solution%interior(i,j,iel,1) = 0.0_prec ! u
          this%solution%interior(i,j,iel,2) = 0.0_prec ! v
          this%solution%interior(i,j,iel,3) = shape*c_mat*c_mat ! p = rho' * c^2
          this%solution%interior(i,j,iel,4) = c_mat ! material sound speed
          this%solution%interior(i,j,iel,5) = rho0_mat ! material background density
        enddo
      enddo
    enddo

    call this%solution%UpdateDevice()

  endsubroutine setInitialCondition

endmodule lineareuler2d_amr_boneandmarrow_model

program LinearEuler2D_AMR_BoneAndMarrow

  use iso_fortran_env,only:int64,real64,output_unit
  use self_data
  use lineareuler2d_amr_boneandmarrow_model
  use SELF_Geometry_2D
  use SELF_AMRController_2D

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 4
  integer,parameter :: targetDegree = 4
  integer,parameter :: maxLevel = 2
  real(prec),parameter :: dtBase = 5.0e-4_prec ! stable on the (curved) base mesh
  real(prec),parameter :: epochLength = 0.05_prec ! adaptation cadence
  ! CI-sized default (debug builds run every example ~5-10x slower than release, and the
  ! curved 973-element base mesh makes geometry generation itself substantial; see issue
  ! #157). One epoch still exercises initial static refinement to the level cap plus one
  ! dynamic adaptation on the unstructured mesh. See SELF_AMR_BONE_EPOCHS below.
  integer,parameter :: defaultEpochs = 1

  type(lineareuler2d_bone_amr) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  type(AMRController2D) :: controller
  integer :: nEpochs,epoch,i,envstat,nElemBase
  logical :: adapted
  real(prec) :: e0,ef,dt
  character(32) :: envstr
  character(LEN=255) :: WORKSPACE
  ! Wall-clock instrumentation of the epoch loop. system_clock with an integer(int64) count
  ! is a monotonic wall clock; ForwardStep's own TIMER macro (src/SELF_Macros.h) expands to
  ! cpu_time in non-multithreaded builds and so cannot be used for GPU timing.
  integer(int64) :: c0clk,c1clk,c2clk,crate
  integer(int64) :: nStepsTotal
  real(real64) :: tFwd,tAdapt,tSim

  ! Number of adaptation epochs (each epochLength = 0.05 time units). The default is
  ! CI-sized: the pulse refines and begins to propagate, exercising the full AMR loop on the
  ! curved unstructured mesh. For the full reflection/transmission movie through the bone and
  ! marrow, set SELF_AMR_BONE_EPOCHS=200 (t = 10).
  nEpochs = defaultEpochs
  call get_environment_variable("SELF_AMR_BONE_EPOCHS",envstr,status=envstat)
  if(envstat == 0) then
    read(envstr,*) nEpochs
  endif

  ! Multi-material ISM-MM mesh; all outer boundaries are transparent (radiation).
  call get_environment_variable("WORKSPACE",WORKSPACE)
  call mesh%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/MultiMaterial2D/BoneAndMarrow.mesh")
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)
  nElemBase = mesh%nElem

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
  call modelobj%SetTimeIntegrator(integrator)

  call controller%Init(modelobj,refineThreshold=-3.0_prec,coarsenThreshold=-8.0_prec, &
                       ivar=3,maxLevel=maxLevel,nHalo=2)

  ! ---- Initial condition + static refinement to the pulse ----
  ! Re-evaluated analytically after each mesh change so the indicator sees the true pulse.
  call modelobj%setInitialCondition()
  do i = 1,maxLevel+2
    call controller%Adapt(modelobj,adapted)
    if(.not. adapted) exit
    call modelobj%setInitialCondition()
  enddo
  print*,"initial adaptation: nElem =",nElemBase," -> ",modelobj%mesh%nElem, &
    ", max level =",controller%forest%MaxLevel()
  if(controller%forest%MaxLevel() < 1 .or. modelobj%mesh%nElem <= nElemBase) then
    print*,"ERROR: the indicator did not refine around the pulse."
    stop 1
  endif

  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy

  ! Snapshot the initial condition (epoch snapshots follow from ForwardStep's file IO).
  call modelobj%WriteModel()
  call modelobj%IncrementIOCounter()

  ! ---- Time-step through adaptation epochs ----
  ! The epoch loop is timed in two buckets: time integration (ForwardStep, which also writes
  ! one snapshot per epoch) and adaptation (Adapt). Adapt opens with solution%UpdateHost(), a
  ! blocking device-to-host copy, so the c1clk boundary lands after the GPU has drained.
  call system_clock(count_rate=crate)
  tFwd = 0.0_real64
  tAdapt = 0.0_real64
  nStepsTotal = 0
  do epoch = 1,nEpochs
    dt = controller%RecommendedTimeStep(dtBase)
    call system_clock(c0clk)
    ! ioInterval = (tn - t) exactly, so ForwardStep's io count is exactly 1 per epoch
    ! regardless of accumulated floating-point drift in t.
    call modelobj%ForwardStep(tn=modelobj%t+epochLength,dt=dt, &
                              ioInterval=(modelobj%t+epochLength)-modelobj%t)
    call system_clock(c1clk)
    call controller%Adapt(modelobj,adapted)
    call system_clock(c2clk)
    tFwd = tFwd+real(c1clk-c0clk,real64)/real(crate,real64)
    tAdapt = tAdapt+real(c2clk-c1clk,real64)/real(crate,real64)
    nStepsTotal = nStepsTotal+int(epochLength/dt,int64)
    print*,"epoch",epoch,": t =",modelobj%t,", dt =",dt, &
      ", nElem =",modelobj%mesh%nElem,", adapted =",adapted
  enddo

  ! ---- Wall-clock summary (fixed BENCH_ keys for machine parsing) ----
  tSim = real(nEpochs,real64)*real(epochLength,real64)
  write(output_unit,'(A)') "BENCH_case = boneandmarrow"
  write(output_unit,'(A,I0)') "BENCH_nEpochs = ",nEpochs
  write(output_unit,'(A,I0)') "BENCH_nSteps = ",nStepsTotal
  write(output_unit,'(A,I0)') "BENCH_nElemBase = ",nElemBase
  write(output_unit,'(A,I0)') "BENCH_nElemFinal = ",modelobj%mesh%nElem
  write(output_unit,'(A,ES16.7E3)') "BENCH_simTime_s = ",tSim
  write(output_unit,'(A,ES16.7E3)') "BENCH_tForwardStep_s = ",tFwd
  write(output_unit,'(A,ES16.7E3)') "BENCH_tAdapt_s = ",tAdapt
  write(output_unit,'(A,ES16.7E3)') "BENCH_tEpochLoop_s = ",tFwd+tAdapt
  write(output_unit,'(A,ES16.7E3)') "BENCH_wallPerStep_s = ",tFwd/real(nStepsTotal,real64)
  write(output_unit,'(A,ES16.7E3)') "BENCH_wallPerSimTime = ",tFwd/tSim
  write(output_unit,'(A,ES16.7E3)') "BENCH_wallPerSimTimeIncAMR = ",(tFwd+tAdapt)/tSim
  write(output_unit,'(A,ES16.7E3)') "BENCH_amrFraction = ",tAdapt/(tFwd+tAdapt)

  ! ---- Integration checks ----
  call modelobj%CalculateEntropy()
  ef = modelobj%entropy
  print*,"initial, final acoustic energy :",e0,ef
  if(ef /= ef) then
    print*,"ERROR: acoustic energy is NaN."
    stop 1
  endif
  if(ef > e0*(1.0_prec+100.0_prec*epsilon(1.0_prec))) then
    print*,"ERROR: acoustic energy grew over the adaptive run."
    stop 1
  endif

  call modelobj%Free()
  call controller%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram LinearEuler2D_AMR_BoneAndMarrow
