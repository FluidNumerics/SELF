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

program LTS_2D_Conservation
!! Local time stepping (AMR Stage 7) on a 2:1 nonconforming mesh: does subcycling the fine
!! level stay conservative, and does the flux-register reflux stay quiet when it should?
!!
!! Assertions:
!!
!!  (a) The RK3 stage weights returned by rk3StageWeight really do reproduce a full
!!      low-storage RK3 step. Everything about the flux registers rests on the telescoping
!!      identity U^{n+1} = U^n + dt * sum_m c_m R_m, so it is checked directly against the
!!      integrator's own rk3_a/rk3_b/rk3_g on a scalar linear ODE, to round-off.
!!
!!  (b) On a conforming (single-level) mesh, ForwardStepLTS reproduces the ordinary RK3
!!      ForwardStep BITWISE. With one level there is no recursion, no mortar and no reflux,
!!      so the two paths must be the same arithmetic in the same order - this is the check
!!      that the LTS driver itself introduces no drift.
!!
!!  (c) On the three-element SimpleMortarMesh (big element at level 0, two small elements at
!!      level 1, so the fine level really does subcycle) with reflecting boundaries, the
!!      domain integral of the acoustic pressure is conserved across many macro steps. The
!!      mortar projection makes each interface conservative in space; refluxing is what
!!      makes it conservative in time when the two sides step at different rates.
!!
!!  (d) The same conservation statement holds on a THREE-level mesh emitted from a real
!!      quadtree forest, where the recursion nests twice (level 2 takes four substeps per
!!      macro step) and level 1 is simultaneously a small side of one interface and the big
!!      side of another - the configuration in which the register bookkeeping and the
!!      parent-window interpolation are easiest to get backwards.
!!
!!  (e) The solution stays finite and NaN-free.
!!
!! (c) is the assertion that actually exercises the flux registers: without refluxing, the
!! coarse side would keep the lagged interface flux it used during its own stages and the
!! integral would drift at the level of the lag.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Model
  use SELF_QuadTreeMesh_2D
  use SELF_AdaptiveMesh_2D
  use SELF_LocalTimeStepping_2D

  implicit none

  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 10
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-4_prec
  real(prec),parameter :: Lr = 0.05_prec
  !! Stable on the level-0 (2*dx) element; level 1 subcycles at dtBase/2.
  real(prec),parameter :: dtBase = 5.0e-4_prec
  real(prec),parameter :: endtime = 0.05_prec
  real(prec),parameter :: iointerval = 0.05_prec

  type(Mesh2D),target :: keepalive
  integer :: keepaliveBCs(1:4)

  call CheckStageWeights()

  keepaliveBCs(1:4) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                       SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]
  ! SELF releases MPI with the last live mesh; keep one alive across both cases.
  call keepalive%StructuredMesh(2,2,1,1,dx,dx,keepaliveBCs)

  call CheckConformingEquivalence()
  call CheckMortarConservation()
  call CheckThreeLevelConservation()

  call keepalive%Free()

  print*,"LTS 2D: stage weights, conforming equivalence, and two- and three-level ", &
    "mortar conservation all pass"

contains

  subroutine CheckStageWeights()
    !! Integrate y' = lambda*y over one step, once with the integrator's own 2N recursion and
    !! once with the telescoped weights, and require the two to agree to round-off.
    implicit none
    ! Local
    integer :: m
    real(prec) :: lam,dt,y2N,yWeighted,g,r,tol

    lam = -0.7_prec
    dt = 0.31_prec

    ! Williamson 2N : g is the low-storage register, r the stage tendency
    y2N = 1.0_prec
    g = 0.0_prec
    do m = 1,3
      r = lam*y2N
      g = rk3_a(m)*g+r
      y2N = y2N+rk3_g(m)*dt*g
    enddo

    ! Telescoped form, re-running the same stage states so that the tendencies match
    yWeighted = 1.0_prec
    block
      real(prec) :: ys,gs,rr(1:3)
      ys = 1.0_prec
      gs = 0.0_prec
      do m = 1,3
        rr(m) = lam*ys
        gs = rk3_a(m)*gs+rr(m)
        ys = ys+rk3_g(m)*dt*gs
      enddo
      do m = 1,3
        yWeighted = yWeighted+dt*rk3StageWeight(m)*rr(m)
      enddo
    endblock

    tol = 8.0_prec*epsilon(1.0_prec)
    if(abs(y2N-yWeighted) > tol) then
      print*,"Error: rk3StageWeight does not reproduce a full low-storage RK3 step."
      print*,"  2N form, weighted form, difference :",y2N,yWeighted,y2N-yWeighted
      stop 1
    endif

    ! Consistency : the weights of any single step must sum to one.
    if(abs(rk3StageWeight(1)+rk3StageWeight(2)+rk3StageWeight(3)-1.0_prec) > tol) then
      print*,"Error: rk3 stage weights do not sum to one :", &
        rk3StageWeight(1),rk3StageWeight(2),rk3StageWeight(3)
      stop 1
    endif

  endsubroutine CheckStageWeights

  subroutine CheckConformingEquivalence()
    !! One refinement level everywhere : ForwardStepLTS must be the ordinary RK3 step,
    !! bitwise. Any difference here is the LTS driver perturbing the single-rate path.
    implicit none
    ! Local
    type(LinearEuler2D) :: mLTS,mRef
    type(Lagrange),target :: interp
    type(Mesh2D),target :: meshLTS,meshRef
    type(SEMQuad),target :: geomLTS,geomRef
    type(LTSSchedule2D) :: sched
    integer :: bcids(1:4)
    integer :: Np,nEl,nVar

    bcids(1:4) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                  SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]

    call interp%Init(N=controlDegree,controlNodeType=GAUSS, &
                     M=targetDegree,targetNodeType=UNIFORM)

    call meshLTS%StructuredMesh(4,4,1,1,dx,dx,bcids)
    call geomLTS%Init(interp,meshLTS%nElem)
    call geomLTS%GenerateFromMesh(meshLTS)
    call mLTS%Init(meshLTS,geomLTS)
    mLTS%prescribed_bcs_enabled = .false.
    mLTS%tecplot_enabled = .false.
    mLTS%rho0 = rho0
    call mLTS%SphericalSoundWave(amp,Lr,2.0_prec*dx,2.0_prec*dx,c0)

    call meshRef%StructuredMesh(4,4,1,1,dx,dx,bcids)
    call geomRef%Init(interp,meshRef%nElem)
    call geomRef%GenerateFromMesh(meshRef)
    call mRef%Init(meshRef,geomRef)
    mRef%prescribed_bcs_enabled = .false.
    mRef%tecplot_enabled = .false.
    mRef%rho0 = rho0
    call mRef%SphericalSoundWave(amp,Lr,2.0_prec*dx,2.0_prec*dx,c0)

    call ForwardStepLTS(mLTS,sched,tn=endtime,dtBase=dtBase,ioInterval=iointerval)

    call mRef%SetTimeIntegrator('rk3')
    call mRef%ForwardStep(tn=endtime,dt=dtBase,ioInterval=iointerval)

    Np = controlDegree+1
    nEl = meshLTS%nElem
    nVar = mLTS%solution%nVar

    if(any(mLTS%solution%interior(1:Np,1:Np,1:nEl,1:nVar) /= &
           mRef%solution%interior(1:Np,1:Np,1:nEl,1:nVar))) then
      print*,"Error: on a single-level mesh, ForwardStepLTS did not reproduce ForwardStep."
      print*,"  max |difference| : ", &
        maxval(abs(mLTS%solution%interior(1:Np,1:Np,1:nEl,1:nVar)- &
                   mRef%solution%interior(1:Np,1:Np,1:nEl,1:nVar)))
      stop 1
    endif

    call sched%Free()
    call mLTS%Free()
    call mRef%Free()
    call geomLTS%Free()
    call geomRef%Free()
    call meshLTS%Free()
    call meshRef%Free()
    call interp%Free()

  endsubroutine CheckConformingEquivalence

  subroutine CheckMortarConservation()
    !! Two levels across a 2:1 mortar, reflecting boundaries : the domain integral of the
    !! pressure must not drift while the fine level subcycles.
    implicit none
    ! Local
    type(LinearEuler2D) :: modelobj
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(LTSSchedule2D) :: sched
    integer :: bcids(1:4)
    real(prec) :: m0,mf,e0,ef,scale,defect

    ! No-normal-flow on every side, so the acoustic mass has nowhere to leave and any drift
    ! in its integral is a defect of the scheme rather than of the physics.
    bcids(1:4) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                  SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]

    call mesh%SimpleMortarMesh(dx,bcids)

    if(mesh%maxElemLevel /= 1) then
      print*,"Error: SimpleMortarMesh should report a level cap of 1, got ",mesh%maxElemLevel
      stop 1
    endif

    call interp%Init(N=controlDegree,controlNodeType=GAUSS, &
                     M=targetDegree,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call modelobj%Init(mesh,geometry)
    modelobj%prescribed_bcs_enabled = .false.
    modelobj%tecplot_enabled = .false.
    modelobj%rho0 = rho0

    ! Pulse centred on the mortar, so the interface carries a real signal throughout.
    call modelobj%SphericalSoundWave(amp,Lr,2.0_prec*dx,dx,c0)

    m0 = PressureIntegral(modelobj,geometry)
    call modelobj%CalculateEntropy()
    e0 = modelobj%entropy

    call ForwardStepLTS(modelobj,sched,tn=endtime,dtBase=dtBase,ioInterval=iointerval)

    mf = PressureIntegral(modelobj,geometry)
    call modelobj%CalculateEntropy()
    ef = modelobj%entropy

    if(ef /= ef) then
      print*,"Error: entropy is NaN after local time stepping across the mortar."
      stop 1
    endif
    if(ef > e0) then
      print*,"Error: entropy grew under local time stepping :",e0,ef
      stop 1
    endif

    ! Scale the defect by the pulse amplitude times the domain measure, so the tolerance is
    ! a relative one rather than a magic absolute number.
    scale = amp*c0*c0*(6.0_prec*dx*dx)
    defect = abs(mf-m0)/scale
    print*,"LTS mortar conservation : integral before, after, relative defect :",m0,mf,defect
    if(defect > 1.0e-10_prec) then
      print*,"Error: the pressure integral drifted across the level interface. Refluxing ", &
        "should hold this at round-off."
      stop 1
    endif

    call sched%Free()
    call modelobj%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

  endsubroutine CheckMortarConservation

  subroutine CheckThreeLevelConservation()
    !! Three refinement levels emitted from a real quadtree forest, so the recursion nests
    !! twice (level 2 subcycles four times per macro step) and level 1 acts as both a big
    !! side and a small side at once - the configuration where the register bookkeeping and
    !! the parent-window interpolation could most easily be wired up back to front.
    !!
    !! The same conservation statement as the two-level case must still hold at round-off.
    implicit none
    ! Local
    type(LinearEuler2D) :: modelobj
    type(Lagrange),target :: interp
    type(Mesh2D),target :: baseMesh,mesh
    type(SEMQuad),target :: geometry
    type(QuadTreeMesh2D) :: forest
    type(LTSSchedule2D) :: sched
    integer :: bcids(1:4)
    integer,allocatable :: flag(:)
    real(prec) :: m0,mf,e0,ef,scale,defect

    bcids(1:4) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                  SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]

    ! 4x4 base mesh; refine one corner twice so the forest carries levels 0, 1 and 2.
    call baseMesh%StructuredMesh(4,4,1,1,dx,dx,bcids)
    call forest%Init(baseMesh)
    allocate(flag(1:forest%nLeaves)); flag = QUADTREE_KEEP; flag(1) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag); deallocate(flag)
    allocate(flag(1:forest%nLeaves)); flag = QUADTREE_KEEP; flag(1) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag); deallocate(flag)
    call forest%Balance2to1()
    call EmitMesh(forest,baseMesh,mesh)

    if(mesh%maxElemLevel /= 2) then
      print*,"Error: the emitted forest should span three levels, got a cap of ", &
        mesh%maxElemLevel
      stop 1
    endif
    if(mesh%nMortars <= 0) then
      print*,"Error: a 2:1-balanced three-level forest must emit mortars."
      stop 1
    endif

    call interp%Init(N=controlDegree,controlNodeType=GAUSS, &
                     M=targetDegree,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call modelobj%Init(mesh,geometry)
    modelobj%prescribed_bcs_enabled = .false.
    modelobj%tecplot_enabled = .false.
    modelobj%rho0 = rho0

    ! Centre the pulse on the refined corner so the signal actually crosses both level jumps.
    call modelobj%SphericalSoundWave(amp,Lr,dx,dx,c0)

    m0 = PressureIntegral(modelobj,geometry)
    call modelobj%CalculateEntropy()
    e0 = modelobj%entropy

    call ForwardStepLTS(modelobj,sched,tn=endtime,dtBase=dtBase,ioInterval=iointerval)

    mf = PressureIntegral(modelobj,geometry)
    call modelobj%CalculateEntropy()
    ef = modelobj%entropy

    if(ef /= ef) then
      print*,"Error: entropy is NaN after three-level local time stepping."
      stop 1
    endif
    if(ef > e0) then
      print*,"Error: entropy grew under three-level local time stepping :",e0,ef
      stop 1
    endif

    scale = amp*c0*c0*(16.0_prec*dx*dx)
    defect = abs(mf-m0)/scale
    print*,"LTS three-level conservation : integral before, after, relative defect :", &
      m0,mf,defect
    if(defect > 1.0e-10_prec) then
      print*,"Error: the pressure integral drifted across the three-level interfaces."
      stop 1
    endif

    call sched%Free()
    call modelobj%Free()
    call geometry%Free()
    call mesh%Free()
    call forest%Free()
    call baseMesh%Free()
    call interp%Free()

  endsubroutine CheckThreeLevelConservation

  real(prec) function PressureIntegral(modelobj,geometry) result(total)
    !! Jacobian-and-quadrature-weighted integral of the pressure (variable 3) over the whole
    !! rank-local domain. With reflecting boundaries this is a conserved quantity of the
    !! linearized acoustics, so it is the sharpest available probe of interface conservation.
    implicit none
    type(LinearEuler2D),intent(in) :: modelobj
    type(SEMQuad),intent(in) :: geometry
    ! Local
    integer :: iel,i,j

    total = 0.0_prec
    do iel = 1,geometry%nelem
      do j = 1,modelobj%solution%interp%N+1
        do i = 1,modelobj%solution%interp%N+1
          total = total+modelobj%solution%interior(i,j,iel,3)* &
                  abs(geometry%J%interior(i,j,iel,1))* &
                  modelobj%solution%interp%qWeights(i)* &
                  modelobj%solution%interp%qWeights(j)
        enddo
      enddo
    enddo

  endfunction PressureIntegral

endprogram LTS_2D_Conservation
