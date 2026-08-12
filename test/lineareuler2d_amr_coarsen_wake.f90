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

program LinearEuler2D_AMR_CoarsenWake
!! Regression test for AMR coarsening BEHIND a propagating front (issue #162).
!!
!! The modal-decay indicator's smoothness ratio is scale-free, so before the amplitude gate existed
!! an element carrying low-amplitude grid-scale residue - what a wave leaves behind once it has
!! passed - looked exactly as under-resolved as the front itself, and was flagged REFINE forever.
!! The measured symptom was an element count that is monotone non-decreasing for a whole run:
!! refinement tracked the wave outward and the refined region grew to the entire swept area, so the
!! saving AMR delivered was set by how little of the domain the wave had reached rather than by
!! adaptivity.
!!
!! The field here is prescribed rather than time-integrated: an expanding ring pulse whose radius is
!! stepped outward, on top of a fixed grid-scale residue at 1e-8 of the ring's amplitude. That is a
!! deliberate choice, and it isolates what the issue is about:
!!
!!   - the adaptation loop under test is indicator -> flags -> forest -> emitted mesh, which the
!!     prescribed field drives exactly as a solver would, at a fraction of the cost;
!!   - the residue is what actually pins the mesh, and prescribing it makes its amplitude an input
!!     rather than an artefact of how long the run happens to be;
!!   - the front stays inside the domain throughout. That matters: the gate normalizes against the
!!     field's peak, so once a wave has radiated away entirely the scale collapses onto the
!!     residue's own peak and no relative gate can discriminate. A propagating-wave problem - the
!!     case the issue is about, and the case AMR is for - always has its front in the domain.
!!
!! Assertions:
!!  (a) the ring is refined (an over-aggressive gate would suppress the front itself);
!!  (b) the mesh is released WHERE THE FRONT HAS LEFT: at the final radius every element well
!!      inside the ring is back at the base (level-0) size;
!!  (c) the refined region follows the ring instead of filling the domain: the element count stays
!!      far below the saturated (uniformly refined) mesh the ungated run degenerates to;
!!  (d) control - with the gate disabled (relativeEnergyFloor = 0, i.e. the pre-fix behaviour) the
!!      same sweep leaves the wake refined, so the assertions above have teeth.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Geometry_2D
  use SELF_AMRController_2D
  use SELF_RefinementIndicator_2D

  implicit none

  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 7
  integer,parameter :: maxLevel = 1
  integer,parameter :: nElemX = 16 ! base-mesh elements per direction
  real(prec),parameter :: Lx = 1.0_prec ! domain extent (m)
  real(prec),parameter :: dx = Lx/real(nElemX,prec)
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: rho0 = 1.0_prec
  ! Ring amplitude and half-width. Lr is below the base-mesh nodal spacing (dx/controlDegree =
  ! 12.5 mm), so the ring is under-resolved at level 0 and the indicator must refine it.
  real(prec),parameter :: amp = 1.0e6_prec
  real(prec),parameter :: Lr = 0.01_prec
  ! Grid-scale residue left everywhere, at 1e-8 of the ring's amplitude: the highest tensor Legendre
  ! mode on each element, for which the smoothness ratio is exactly 1 whatever the mesh. Its modal
  ! energy is 1e-4 - twelve orders of magnitude above machine epsilon in real64, and still ~1e3
  ! above it in real32, so the absolute floor alone never touches it and the ungated control really
  ! does refine on it - against a ring-element energy of order 1e11. That energy ratio, ~1e-15, is
  ! what the relative floor gates, with three orders of margin below the 1e-12 default.
  real(prec),parameter :: wakeAmp = 1.0e-2_prec
  ! Radii the ring is stepped through (m). The last one leaves a wake several base elements deep.
  real(prec),parameter :: ringRadius(1:5) = &
                          [0.05_prec,0.10_prec,0.15_prec,0.20_prec,0.25_prec]
  ! Elements whose centroid is inside this radius are behind the final ring by at least three base
  ! elements, so nothing there can be held refined by the front or its halo.
  real(prec),parameter :: wakeRadius = 0.05_prec

  integer :: nGated,nUngated
  real(prec) :: hGated,hUngated
  integer :: r
  ! SELF owns the MPI lifecycle and keeps it alive only while a mesh is live, so a run that tears
  ! its mesh down between the two sweeps would finalize MPI and leave the second sweep running
  ! against a finalized library (see test/domaindecomposition_two_meshes.f90). This mesh exists
  ! solely to hold that reference open across both sweeps.
  type(Mesh2D),target :: keepAlive
  integer :: keepAliveBcids(1:4)

  r = 0

  keepAliveBcids(1:4) = SELF_BC_NONORMALFLOW
  call keepAlive%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,keepAliveBcids)

  ! ---- The shipped amplitude gate ----
  call RunSweep(SELF_AMR_DEFAULT_RELFLOOR,nGated,hGated)
  print*,"with the amplitude gate    : nElem =",nGated, &
    ", smallest element near the centre / base =",hGated

  ! ---- (d) Control: the pre-fix behaviour ----
  call RunSweep(0.0_prec,nUngated,hUngated)
  print*,"with the gate disabled     : nElem =",nUngated, &
    ", smallest element near the centre / base =",hUngated

  ! ---- (b) The wake is released ----
  if(hGated < 0.75_prec) then
    print*,"FAIL: the wake behind the ring was not coarsened back to the base mesh"
    r = 1
  endif

  ! ---- (d) ... and would not have been before the gate existed ----
  if(hUngated >= 0.75_prec) then
    print*,"FAIL: the control run also released the wake, so this test proves nothing"
    r = 1
  endif
  if(nGated >= nUngated) then
    print*,"FAIL: the gate did not reduce the element count",nGated,nUngated
    r = 1
  endif
  ! ---- (c) The refined region follows the ring instead of filling the domain ----
  ! Without the gate every element is eventually flagged REFINE and the mesh saturates at the
  ! level cap (4*nElemX**2 elements), so AMR degenerates into a uniformly fine mesh. The gated
  ! run must stay far below that; the margin here is large (roughly 2.3x), so half is a
  ! comfortable bound that still fails loudly if the gate stops working.
  block
    integer :: nSaturated
    nSaturated = 4*nElemX*nElemX
    print*,"element count : gated",nGated,", ungated",nUngated, &
      ", saturated",nSaturated,", base",nElemX*nElemX
    if(nGated > nSaturated/2) then
      print*,"FAIL: the refined region is filling the domain rather than following the ring"
      r = 1
    endif
  endblock

  if(r == 0) print*,"LINEAR EULER 2D AMR COARSEN WAKE CHECKS PASSED"

  call keepAlive%Free()

  if(r /= 0) stop 1

contains

  subroutine RunSweep(relativeEnergyFloor,nElemFinal,hRelCentre)
    !! Step an expanding ring through the domain, adapting to it at each radius, and report the
    !! final element count together with the size of the smallest element near the domain centre
    !! relative to a base element.
    implicit none
    real(prec),intent(in) :: relativeEnergyFloor
    integer,intent(out) :: nElemFinal
    real(prec),intent(out) :: hRelCentre
    ! Local
    type(LinearEuler2D) :: modelobj
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(AMRController2D) :: controller
    integer :: bcids(1:4)
    integer :: stage,i
    logical :: adapted

    bcids(1:4) = [SELF_BC_RADIATION, & ! south
                  SELF_BC_RADIATION, & ! east
                  SELF_BC_RADIATION, & ! north
                  SELF_BC_RADIATION] ! west

    call mesh%StructuredMesh(nElemX,nElemX,1,1,dx,dx,bcids)
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

    ! Recommended starting thresholds from the AMR documentation; pressure (variable 3) drives the
    ! indicator, and nHalo = 1 keeps the refined band one element wider than the ring.
    call controller%Init(modelobj,refineThreshold=-3.0_prec,coarsenThreshold=-8.0_prec, &
                         ivar=3,maxLevel=maxLevel,nHalo=1, &
                         relativeEnergyFloor=relativeEnergyFloor)

    do stage = 1,size(ringRadius)
      ! Adapt to convergence at this radius, re-evaluating the field on each new mesh so the
      ! indicator sees the true ring rather than its coarse-mesh interpolant.
      call SetRingField(modelobj,interp,ringRadius(stage))
      do i = 1,maxLevel+2
        call controller%Adapt(modelobj,adapted)
        if(.not. adapted) exit
        call SetRingField(modelobj,interp,ringRadius(stage))
      enddo
      print*,"  ring radius",ringRadius(stage),": nElem =",modelobj%mesh%nElem, &
        ", max level =",controller%forest%MaxLevel()

      ! (a) The ring itself must be refined - a gate that suppressed the front would show up here.
      if(controller%forest%MaxLevel() < 1) then
        print*,"FAIL: the ring was never refined at radius",ringRadius(stage)
        r = 1
      endif
    enddo

    nElemFinal = modelobj%mesh%nElem
    hRelCentre = RelativeSizeNearCentre(modelobj,interp)

    call modelobj%Free()
    call controller%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

  endsubroutine RunSweep

  subroutine SetRingField(modelobj,interp,radius)
    !! Prescribe an expanding ring pulse of the given radius on the model's current mesh, plus a
    !! fixed grid-scale residue. The residue is written in REFERENCE coordinates (the highest
    !! tensor Legendre mode of the element), so it stays grid-scale on whatever mesh the previous
    !! adaptation produced - which is precisely the residue that a scale-free smoothness ratio
    !! cannot distinguish from a genuine front.
    implicit none
    type(LinearEuler2D),intent(inout) :: modelobj
    type(Lagrange),intent(in) :: interp
    real(prec),intent(in) :: radius
    ! Local
    integer :: i,j,iEl
    real(prec) :: x,y,rr,ring,residue

    do iEl = 1,modelobj%mesh%nElem
      do j = 1,controlDegree+1
        do i = 1,controlDegree+1
          x = modelobj%geometry%x%interior(i,j,iEl,1,1)-0.5_prec*Lx
          y = modelobj%geometry%x%interior(i,j,iEl,1,2)-0.5_prec*Lx
          rr = sqrt(x*x+y*y)
          ring = amp*exp(-log(2.0_prec)*(rr-radius)**2/Lr**2)
          residue = wakeAmp*nlegendre(controlDegree,interp%controlPoints(i))* &
                    nlegendre(controlDegree,interp%controlPoints(j))
          modelobj%solution%interior(i,j,iEl,1) = 0.0_prec
          modelobj%solution%interior(i,j,iEl,2) = 0.0_prec
          modelobj%solution%interior(i,j,iEl,3) = ring+residue
          modelobj%solution%interior(i,j,iEl,4) = c0
          modelobj%solution%interior(i,j,iEl,5) = rho0
        enddo
      enddo
    enddo
    call modelobj%solution%UpdateDevice() ! no-op on CPU; pushes interior to device on GPU

  endsubroutine SetRingField

  function RelativeSizeNearCentre(modelobj,interp) result(hRel)
    !! Size of the smallest element whose centroid lies within wakeRadius of the domain centre, as
    !! a fraction of a base (level-0) element. Element size is read back from the geometry the
    !! model is actually running on rather than inferred from a refinement level. The span is
    !! measured across control points, which on Gauss nodes stop short of the element edges; the
    !! reference span is measured the same way, so that factor cancels.
    implicit none
    type(LinearEuler2D),intent(in) :: modelobj
    type(Lagrange),intent(in) :: interp
    real(prec) :: hRel
    ! Local
    integer :: iEl,ii,jj,np
    real(prec) :: xc,yc,xmin,xmax,rc,hNear,hBase

    np = controlDegree+1
    hNear = huge(1.0_prec)
    hBase = dx*0.5_prec*(interp%controlPoints(np)-interp%controlPoints(1))
    do iEl = 1,modelobj%mesh%nElem
      xc = 0.0_prec
      yc = 0.0_prec
      xmin = huge(1.0_prec)
      xmax = -huge(1.0_prec)
      do jj = 1,np
        do ii = 1,np
          xc = xc+modelobj%geometry%x%interior(ii,jj,iEl,1,1)
          yc = yc+modelobj%geometry%x%interior(ii,jj,iEl,1,2)
          xmin = min(xmin,modelobj%geometry%x%interior(ii,jj,iEl,1,1))
          xmax = max(xmax,modelobj%geometry%x%interior(ii,jj,iEl,1,1))
        enddo
      enddo
      xc = xc/real(np*np,prec)
      yc = yc/real(np*np,prec)
      rc = sqrt((xc-0.5_prec*Lx)**2+(yc-0.5_prec*Lx)**2)
      if(rc <= wakeRadius) hNear = min(hNear,xmax-xmin)
    enddo

    hRel = hNear/hBase

  endfunction RelativeSizeNearCentre

  pure function nlegendre(p,x) result(Lp)
    !! L2-normalized Legendre polynomial on [-1,1], matching the basis used by the indicator.
    implicit none
    integer,intent(in) :: p
    real(prec),intent(in) :: x
    real(prec) :: Lp
    integer :: k
    real(prec) :: lkm1,lk,lkp1

    if(p == 0) then
      Lp = 1.0_prec
    elseif(p == 1) then
      Lp = x
    else
      lkm1 = 1.0_prec
      lk = x
      do k = 1,p-1
        lkp1 = (real(2*k+1,prec)*x*lk-real(k,prec)*lkm1)/real(k+1,prec)
        lkm1 = lk
        lk = lkp1
      enddo
      Lp = lk
    endif
    Lp = Lp*sqrt((2.0_prec*real(p,prec)+1.0_prec)/2.0_prec)

  endfunction nlegendre

endprogram LinearEuler2D_AMR_CoarsenWake
