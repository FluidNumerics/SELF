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

program DGModel2D_Unmapped_Boundary_MPI
!! Multi-rank companion to dgmodel2d_unmapped_boundary (issue #180).
!!
!! Each rank owns a slice of the mesh, so a bcid that appears on no local face may still be
!! present on another rank. MapBoundaryConditions therefore reduces the count across the
!! communicator. This test asserts on EVERY rank that the reported count is the GLOBAL
!! number of domain boundary edges - a per-rank tally would be smaller on every rank.
!!
!! Assertions:
!!  (a) every rank sees the global unmapped-edge count and the offending bcid;
!!  (b) re-tagging to a registered bcid clears the count on every rank.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d

  implicit none

  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  integer,parameter :: nxPerTile = 5
  integer,parameter :: nTileX = 2
  !! Global count: one boundary edge per perimeter element per domain side.
  integer,parameter :: nBoundaryEdges = 4*nxPerTile*nTileX
  integer,parameter :: bogusBCID = 999

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

  call modelobj%Init(mesh,geometry)
  modelobj%tecplot_enabled = .false.

  ! (a) the count is global, and identical on every rank
  if(modelobj%nUnmappedBoundaries /= nBoundaryEdges) then
    print*,"Error: rank ",mesh%decomp%rankId," expected the global count ",nBoundaryEdges, &
      " unmapped boundary edges, saw ",modelobj%nUnmappedBoundaries
    stop 1
  endif
  if(modelobj%unmappedBoundaryID /= bogusBCID) then
    print*,"Error: rank ",mesh%decomp%rankId," expected the unregistered bcid ",bogusBCID, &
      ", got ",modelobj%unmappedBoundaryID
    stop 1
  endif

  ! (b) re-tagging to a registered bcid clears the count everywhere
  ! As in the serial twin: MapBoundaryConditions re-maps on the host only, so on a GPU
  ! build the device BC arrays are left stale. Nothing steps after this point.
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)
  call modelobj%MapBoundaryConditions()
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: rank ",mesh%decomp%rankId," still reports ", &
      modelobj%nUnmappedBoundaries," unmapped boundary edges after re-tagging."
    stop 1
  endif

  if(mesh%decomp%rankId == 0) then
    print*,"dgmodel2d_unmapped_boundary_mpi : all assertions passed"
  endif

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram DGModel2D_Unmapped_Boundary_MPI
