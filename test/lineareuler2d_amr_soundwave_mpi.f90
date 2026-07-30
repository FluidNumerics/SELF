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

program LinearEuler2D_AMR_SoundWave
!! Regression test for the 2-D AMR controller driving a live LinearEuler2D model: a spherical
!! acoustic pulse, deliberately under-resolved on a coarse base mesh, is refined by the
!! modal-decay indicator and time-stepped through several adaptation epochs on the resulting
!! nonconforming (mortar) meshes.
!!
!! Assertions:
!!  (a) the initial adaptation refines the mesh around the pulse (element count grows and the
!!      forest reaches a level > 0), and the level cap is respected;
!!  (b) the solution transfer at an adaptation epoch conserves the Jacobian-weighted global
!!      integral of every prognostic variable to roundoff;
!!  (c) the model's time and model-specific parameters (rho0) survive regridding, and the
!!      recommended time step follows the level-based bound dtBase/2**MaxLevel;
!!  (d) at least one mid-run epoch adapts the mesh (the refined band follows the wave);
!!  (e) the acoustic energy (entropy) stays finite, NaN-free, and does not exceed its initial
!!      value: the interior upwind flux and radiation boundaries are dissipative, prolongation
!!      preserves the discrete L2 energy, and restriction (an L2 projection) cannot increase it.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Geometry_2D
  use SELF_AMRController_2D
  use SELF_RefinementIndicator_2D
  use mpi

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 7
  integer,parameter :: maxLevel = 2
  integer,parameter :: nEpochs = 5
  real(prec),parameter :: dx = 0.25_prec ! 4 x 4 base mesh spanning 1 m x 1 m
  real(prec),parameter :: dtBase = 1.0e-3_prec ! stable on the level-0 mesh
  real(prec),parameter :: epochLength = 0.02_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-4_prec
  real(prec),parameter :: Lr = 0.05_prec ! pulse half-width << dx: under-resolved at level 0
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 1.0e-10_prec
#else
  real(prec),parameter :: tolerance = 1.0e-3_prec
#endif

  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  type(AMRController2D) :: controller
  integer :: bcids(1:4)
  integer :: epoch,i,nElemPrev,nAdapted
  logical :: adapted
  real(prec) :: e0,ef,dt,tExpected
  real(prec) :: intBefore(1:3),intAfter(1:3)
  integer :: r

  r = 0

  ! Radiation (outflow) conditions on all domain boundaries
  bcids(1:4) = [SELF_BC_RADIATION, & ! south
                SELF_BC_RADIATION, & ! east
                SELF_BC_RADIATION, & ! north
                SELF_BC_RADIATION] ! west

  call mesh%StructuredMesh(4,4,1,1,dx,dx,bcids)
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%tecplot_enabled = .false.
  modelobj%rho0 = rho0
  call modelobj%SetTimeIntegrator(integrator)

  ! Recommended starting thresholds from the AMR documentation; pressure (variable 3) drives
  ! the indicator.
  call controller%Init(modelobj,refineThreshold=-3.0_prec,coarsenThreshold=-8.0_prec, &
                       ivar=3,maxLevel=maxLevel,nHalo=1)

  ! ---- (a) Initial adaptation: refine until the (re-evaluated) pulse is resolved ----
  ! After each mesh change the initial condition is re-evaluated analytically on the new mesh,
  ! so the indicator sees the true pulse rather than its coarse-mesh interpolant.
  call modelobj%SphericalSoundWave(amp,Lr,0.5_prec,0.5_prec,c0)
  nElemPrev = modelobj%mesh%nElem
  do i = 1,maxLevel+2
    call controller%Adapt(modelobj,adapted)
    if(.not. adapted) exit
    call modelobj%SphericalSoundWave(amp,Lr,0.5_prec,0.5_prec,c0)
  enddo
  print*,"initial adaptation: nElem ",nElemPrev," -> ",modelobj%mesh%nElem, &
    ", max level ",controller%forest%MaxLevel()
  if(modelobj%mesh%nElem <= nElemPrev .or. controller%forest%MaxLevel() < 1) then
    print*,"FAIL: initial adaptation did not refine around the pulse"
    r = 1
  endif
  if(controller%forest%MaxLevel() > maxLevel) then
    print*,"FAIL: refinement exceeded the level cap"
    r = 1
  endif

  ! ---- (c) Level-based time step ----
  dt = controller%RecommendedTimeStep(dtBase)
  if(abs(dt-dtBase/real(2**controller%forest%MaxLevel(),prec)) > 0.0_prec) then
    print*,"FAIL: RecommendedTimeStep does not follow dtBase/2**MaxLevel"
    r = 1
  endif

  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy

  ! ---- Time-step through adaptation epochs ----
  nAdapted = 0
  do epoch = 1,nEpochs
    dt = controller%RecommendedTimeStep(dtBase)
    ! Drive the epoch with ioInterval = (tn - t) exactly, so ForwardStep's io count
    ! int((tn-t)/ioInterval) is exactly 1 regardless of accumulated floating-point drift in t.
    tExpected = modelobj%t
    call modelobj%ForwardStep(tn=modelobj%t+epochLength,dt=dt, &
                              ioInterval=(modelobj%t+epochLength)-modelobj%t)
    if(modelobj%t <= tExpected) then
      print*,"FAIL: epoch",epoch,"did not advance the model time",tExpected,modelobj%t
      r = 1
    endif

    tExpected = modelobj%t

    ! (b) Conservation of every prognostic variable across the epoch's transfer.
    do i = 1,3
      intBefore(i) = GlobalIntegral(i)
    enddo
    call controller%Adapt(modelobj,adapted)
    if(adapted) then
      nAdapted = nAdapted+1
      do i = 1,3
        intAfter(i) = GlobalIntegral(i)
        if(abs(intAfter(i)-intBefore(i)) > tolerance*max(abs(intBefore(i)),1.0_prec)) then
          print*,"FAIL: adaptation did not conserve variable ",i,intBefore(i),intAfter(i)
          r = 1
        endif
      enddo
      ! (c) Time state and model parameters survive the regrid.
      if(modelobj%t /= tExpected) then
        print*,"FAIL: model time changed across a regrid",tExpected,modelobj%t
        r = 1
      endif
      if(modelobj%rho0 /= rho0) then
        print*,"FAIL: model parameter rho0 changed across a regrid"
        r = 1
      endif
    endif

    call modelobj%CalculateEntropy()
    print*,"epoch",epoch,": t =",modelobj%t,", nElem =",modelobj%mesh%nElem, &
      ", adapted =",adapted,", entropy =",modelobj%entropy
  enddo

  ! ---- (d) The mesh evolved during the run ----
  if(nAdapted < 1) then
    print*,"FAIL: no mid-run epoch adapted the mesh"
    r = 1
  endif

  ! ---- Forest replication across ranks (no-op on one rank) ----
  ! Every rank must hold the identical forest: compare a leaf-list checksum via min/max
  ! reductions.
  if(modelobj%mesh%decomp%mpiEnabled) then
    block
      integer :: chk,chkmin,chkmax,ierror
      chk = controller%forest%nLeaves+sum(controller%forest%leaf(1:controller%forest%nLeaves))
      call mpi_allreduce(chk,chkmin,1,MPI_INTEGER,MPI_MIN,modelobj%mesh%decomp%mpiComm,ierror)
      call mpi_allreduce(chk,chkmax,1,MPI_INTEGER,MPI_MAX,modelobj%mesh%decomp%mpiComm,ierror)
      if(chkmin /= chkmax) then
        print*,"FAIL: forest is not identical across ranks"
        r = 1
      endif
    endblock
  endif

  ! ---- (e) Entropy is finite, NaN-free, and non-increasing overall ----
  ef = modelobj%entropy
  if(ef /= ef) then
    print*,"FAIL: entropy is NaN"
    r = 1
  endif
  if(ef > e0+tolerance*max(abs(e0),1.0_prec)) then
    print*,"FAIL: entropy grew across the AMR run :",e0,ef
    r = 1
  endif
  print*,"initial, final entropy :",e0,ef

  if(r == 0) print*,"LINEAR EULER 2D AMR SOUNDWAVE CHECKS PASSED"

  call modelobj%Free()
  call controller%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

  if(r /= 0) stop 1

contains

  function GlobalIntegral(ivar) result(intU)
    !! Jacobian-weighted global integral of solution variable ivar on the model's current
    !! mesh; reduced over all ranks when the mesh is decomposed (mirroring CalculateEntropy).
    implicit none
    integer,intent(in) :: ivar
    real(prec) :: intU
    ! Local
    integer :: iEl,ii,jj,ierror
    real(prec) :: intLocal

    ! Synchronize the host mirror before reading it. The multi-rank path still transfers on the
    ! host (Stage-5 v1 migration allgathers there), so this is currently redundant, but a caller
    ! reading solution%interior must not depend on where the last write happened to land.
    call modelobj%solution%UpdateHost()

    intLocal = 0.0_prec
    do iEl = 1,modelobj%mesh%nElem
      do jj = 1,controlDegree+1
        do ii = 1,controlDegree+1
          intLocal = intLocal+modelobj%solution%interior(ii,jj,iEl,ivar)* &
                     modelobj%geometry%J%interior(ii,jj,iEl,1)* &
                     interp%qWeights(ii)*interp%qWeights(jj)
        enddo
      enddo
    enddo

    if(modelobj%mesh%decomp%mpiEnabled) then
      call mpi_allreduce(intLocal,intU,1,modelobj%mesh%decomp%mpiPrec,MPI_SUM, &
                         modelobj%mesh%decomp%mpiComm,ierror)
    else
      intU = intLocal
    endif

  endfunction GlobalIntegral

endprogram LinearEuler2D_AMR_SoundWave
