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

program GeometryReuse2D_MPI
!! Multi-rank validation of the AMR Stage 6c geometry reuse.
!!
!! On one rank, reuse is unconditional for every element the transfer plan marks
!! SELF_TRANSFER_COPY. On several ranks it is not: sourceElem indexes the GLOBAL old element list
!! while a rank holds only its own contiguous slice, and EmitMesh re-decomposes the leaf list
!! every epoch, so refinement anywhere shifts the rank boundaries and a kept element's source can
!! end up on another rank. The controller therefore additionally requires the source to be
!! locally owned - a range test against the previous decomposition - and regenerates the rest.
!!
!! The single-rank tests cannot reach that branch, because there the range test is always true.
!! This test exists to exercise it and to prove the resulting MIXED assembly (some elements
!! copied forward, some regenerated, in one geometry) is exact.
!!
!! Assertions:
!!  (a) both branches of the predicate actually fire: across the run at least one element was
!!      reused and at least one was regenerated, on at least one rank. Without this the rest of
!!      the test could pass while silently covering nothing.
!!  (b) the rank-local geometry the controller assembled is EXACTLY equal, element by element and
!!      bit for bit, to geometry generated from scratch for the same mesh. Exactness (rather than
!!      a tolerance) is the right bar: an unchanged leaf's node coordinates are reproducible
!!      bit-for-bit from the forest, and per-element generation is element-local, so any
!!      discrepancy is a bug rather than round-off.
!!  (c) Jacobians stay strictly positive on the assembled geometry.

  use iso_fortran_env,only:int64
  use self_data
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Geometry_2D
  use SELF_AMRController_2D
  use mpi

  implicit none

  integer :: exit_code

  exit_code = geometry_2d_reuse_mpi()
  if(exit_code /= 0) stop 1

contains

  integer function geometry_2d_reuse_mpi() result(r)
    implicit none

    character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
    integer,parameter :: controlDegree = 4
    integer,parameter :: targetDegree = 4
    integer,parameter :: maxLevel = 2
    integer,parameter :: nEpochs = 4
    integer,parameter :: nElemX = 8 ! 8 x 8 base mesh over 1 m x 1 m
    real(prec),parameter :: dx = 0.125_prec
    real(prec),parameter :: dtBase = 1.0e-3_prec
    real(prec),parameter :: epochLength = 0.02_prec
    real(prec),parameter :: c0 = 1.0_prec
    real(prec),parameter :: rho0 = 1.0_prec
    real(prec),parameter :: amp = 1.0e-4_prec
    real(prec),parameter :: Lr = 0.04_prec ! under-resolved at level 0, so the indicator refines

    type(LinearEuler2D) :: modelobj
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(AMRController2D) :: controller
    type(SEMQuad) :: ref
    type(Lagrange),pointer :: interpPtr
    integer :: bcids(1:4)
    integer :: epoch,iel,nBad,ierror,nRanks,rankId
    integer(int64) :: reusedTot,genTot
    logical :: adapted
    real(prec) :: dt,minJ

    r = 0

    bcids(1:4) = [SELF_BC_RADIATION,SELF_BC_RADIATION,SELF_BC_RADIATION,SELF_BC_RADIATION]
    call mesh%StructuredMesh(nElemX,nElemX,1,1,dx,dx,bcids)
    call interp%Init(N=controlDegree,controlNodeType=GAUSS, &
                     M=targetDegree,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    nRanks = mesh%decomp%nRanks
    rankId = mesh%decomp%rankId

    call modelobj%Init(mesh,geometry)
    modelobj%prescribed_bcs_enabled = .false.
    modelobj%tecplot_enabled = .false.
    modelobj%rho0 = rho0
    call modelobj%SetTimeIntegrator(integrator)

    call controller%Init(modelobj,refineThreshold=-3.0_prec,coarsenThreshold=-8.0_prec, &
                         ivar=3,maxLevel=maxLevel,nHalo=1)

    call modelobj%SphericalSoundWave(amp,Lr,0.5_prec,0.5_prec,c0)
    call controller%Adapt(modelobj,adapted)
    call modelobj%SphericalSoundWave(amp,Lr,0.5_prec,0.5_prec,c0)

    do epoch = 1,nEpochs
      dt = controller%RecommendedTimeStep(dtBase)
      call modelobj%ForwardStep(tn=modelobj%t+epochLength,dt=dt, &
                                ioInterval=(modelobj%t+epochLength)-modelobj%t)
      call controller%Adapt(modelobj,adapted)
    enddo

    ! ---- (a) both branches of the reuse predicate must have fired ----
    ! Reduced over ranks: on a small mesh a given rank may see only one of the two, but the run
    ! as a whole must exercise both or this test proves nothing.
    call mpi_allreduce(controller%nGeomReused,reusedTot,1,MPI_INTEGER8,MPI_SUM, &
                       mesh%decomp%mpiComm,ierror)
    call mpi_allreduce(controller%nGeomGenerated,genTot,1,MPI_INTEGER8,MPI_SUM, &
                       mesh%decomp%mpiComm,ierror)

    if(rankId == 0) then
      print*,"ranks =",nRanks,": elements reused =",reusedTot,", regenerated =",genTot
    endif
    if(reusedTot == 0) then
      print*,"FAIL: no element was ever reused; the reuse path was not exercised"
      r = 1
    endif
    if(nRanks > 1 .and. genTot == 0) then
      print*,"FAIL: no element was ever regenerated on ",nRanks," ranks, so the", &
        " migrated-source fallback was not exercised"
      r = 1
    endif

    ! ---- (b) the assembled rank-local geometry must equal a from-scratch generation exactly ----
    interpPtr => interp
    call ref%Init(interpPtr,controller%activeMesh%nElem)
    call ref%GenerateFromMesh(controller%activeMesh)

    nBad = 0
    do iel = 1,controller%activeMesh%nElem
      if(.not. SameElement(controller%activeGeom,iel,ref,iel)) then
        nBad = nBad+1
        if(nBad <= 5) print*,"FAIL: rank",rankId,"element",iel, &
          "differs from a from-scratch generation"
      endif
    enddo
    if(nBad /= 0) then
      print*,"FAIL: rank",rankId,"assembled geometry differs on",nBad,"elements"
      r = 1
    endif

    ! ---- (c) positive Jacobians ----
    minJ = minval(controller%activeGeom%J%interior)
    if(minJ <= 0.0_prec) then
      print*,"FAIL: rank",rankId,"non-positive Jacobian",minJ
      r = 1
    endif

    if(r == 0 .and. rankId == 0) then
      print*,"PASS: multi-rank geometry assembly is exact; both reuse and regeneration exercised"
    endif

    ! The controller owns the mesh the model is currently bound to, so the model must be freed
    ! first; after controller%Free the model's mesh/geometry pointers dangle.
    call ref%Free()
    call modelobj%Free()
    call controller%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

  endfunction geometry_2d_reuse_mpi

  logical function SameElement(a,ia,b,ib) result(same)
    !! Exact comparison of every geometry quantity of element ia of a against element ib of b.
    implicit none
    type(SEMQuad),intent(in) :: a
    integer,intent(in) :: ia
    type(SEMQuad),intent(in) :: b
    integer,intent(in) :: ib

    same = all(a%x%interior(:,:,ia,:,:) == b%x%interior(:,:,ib,:,:)) .and. &
           all(a%x%boundary(:,:,ia,:,:) == b%x%boundary(:,:,ib,:,:)) .and. &
           all(a%dsdx%interior(:,:,ia,:,:,:) == b%dsdx%interior(:,:,ib,:,:,:)) .and. &
           all(a%dsdx%boundary(:,:,ia,:,:,:) == b%dsdx%boundary(:,:,ib,:,:,:)) .and. &
           all(a%nHat%boundary(:,:,ia,:,:) == b%nHat%boundary(:,:,ib,:,:)) .and. &
           all(a%nScale%boundary(:,:,ia,:) == b%nScale%boundary(:,:,ib,:)) .and. &
           all(a%J%interior(:,:,ia,:) == b%J%interior(:,:,ib,:)) .and. &
           all(a%J%boundary(:,:,ia,:) == b%J%boundary(:,:,ib,:))

  endfunction SameElement

endprogram GeometryReuse2D_MPI
