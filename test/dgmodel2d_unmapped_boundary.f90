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

program DGModel2D_Unmapped_Boundary
!! Regression test for issue #180 : a mesh face whose bcid has no registered boundary
!! condition used to be skipped silently, leaving its exterior state unwritten.
!!
!! MapBoundaryConditions now also scans the mesh in the reverse direction and counts the
!! boundary edges no registration covers, so ForwardStep can warn about them once.
!!
!! Assertions:
!!  (a) tagging every domain edge with a bcid no model registers counts every one of them,
!!      and reports one of the offending ids;
!!  (b) ForwardStep consumes the warning exactly once;
!!  (c) re-tagging to a registered bcid and re-mapping clears the count, and re-arms the
!!      one-shot warning (this is the path Regrid takes).

  use self_data
  use self_lineareuler2d
  use self_mesh_2d

  implicit none

  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  integer,parameter :: nxPerTile = 5
  integer,parameter :: nTileX = 2
  !! Every element on the perimeter contributes one boundary edge per domain side.
  integer,parameter :: nBoundaryEdges = 4*nxPerTile*nTileX
  !! A bcid no SELF model registers. The point of the test: any integer is a legal bcid, and
  !! one that nothing registers must be reported rather than silently ignored.
  integer,parameter :: bogusBCID = 999
  real(prec),parameter :: dt = 1.0e-4_prec
  real(prec),parameter :: endtime = 2.5e-3_prec
  real(prec),parameter :: iointerval = 1.0e-3_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: rho0 = 1.0_prec

  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  bcids(1:4) = [bogusBCID,bogusBCID,bogusBCID,bogusBCID] ! south, east, north, west

  call mesh%StructuredMesh(nxPerTile,nxPerTile,nTileX,nTileX,0.1_prec,0.1_prec,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Init runs AdditionalInit (which registers no_normal_flow and radiation) and then
  ! MapBoundaryConditions, so the count is available immediately.
  call modelobj%Init(mesh,geometry)
  modelobj%tecplot_enabled = .false.
  modelobj%rho0 = rho0

  ! (a) every domain edge is unmapped, and the offending bcid is named
  if(modelobj%nUnmappedBoundaries /= nBoundaryEdges) then
    print*,"Error: expected ",nBoundaryEdges," unmapped boundary edges, counted ", &
      modelobj%nUnmappedBoundaries
    stop 1
  endif
  if(modelobj%unmappedBoundaryID /= bogusBCID) then
    print*,"Error: expected the unregistered bcid ",bogusBCID," to be reported, got ", &
      modelobj%unmappedBoundaryID
    stop 1
  endif
  if(modelobj%unmappedBoundariesReported) then
    print*,"Error: the warning was marked reported before ForwardStep ran."
    stop 1
  endif

  call modelobj%SphericalSoundWave(1.0e-4_prec,0.02_prec,0.5_prec,0.5_prec,c0)
  call modelobj%SetTimeIntegrator(integrator)

  ! (b) ForwardStep emits the warning, once
  call modelobj%ForwardStep(tn=endtime,dt=dt,ioInterval=iointerval)
  if(.not. modelobj%unmappedBoundariesReported) then
    print*,"Error: ForwardStep did not report the unmapped boundary edges."
    stop 1
  endif

  ! (c) re-tagging to a registered bcid and re-mapping clears the count. This is the path
  ! Regrid takes: Free the lists, AdditionalInit, MapBoundaryConditions.
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)
  call modelobj%MapBoundaryConditions()
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: re-tagging to a registered bcid left ",modelobj%nUnmappedBoundaries, &
      " unmapped boundary edges."
    stop 1
  endif
  if(modelobj%unmappedBoundariesReported) then
    print*,"Error: MapBoundaryConditions did not re-arm the one-shot warning."
    stop 1
  endif

  call modelobj%ForwardStep(tn=2.0_prec*endtime,dt=dt,ioInterval=iointerval)
  if(modelobj%unmappedBoundariesReported) then
    print*,"Error: the warning fired on a mesh with no unmapped boundary edges."
    stop 1
  endif

  print*,"dgmodel2d_unmapped_boundary : all assertions passed"

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram DGModel2D_Unmapped_Boundary
