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

program DGModel3D_Unmapped_Boundary_MPI
!! Multi-rank companion to dgmodel3d_unmapped_boundary. MapBoundaryConditions is a separate
!! routine per dimension, so the 2-D MPI test does not cover these collectives: the sum over
!! ranks, and the presence-aware selection of the representative bcid.
!!
!! The id is deliberately NOT reduced over its value. A bcid is any integer, so no value can
!! mark "absent"; ranks agree by MPI_MIN on the lowest-numbered rank that actually holds an
!! offender and that rank broadcasts its first one. A negative id is the case that a max()
!! over values would get wrong, so this test uses one.
!!
!! Assertions, on EVERY rank:
!!  (a) the count is the GLOBAL number of unmapped faces, not the rank-local tally;
!!  (b) the representative id is the offending one, including when it is negative;
!!  (c) re-tagging to a registered bcid clears the count everywhere.

  use self_data
  use self_lineareuler3d
  use self_mesh_3d

  implicit none

  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  integer,parameter :: nPerTile = 2
  integer,parameter :: nTile = 2
  !! Global: n^3 elements, 6 faces of n^2 each on the domain boundary.
  integer,parameter :: nSide = nPerTile*nTile
  integer,parameter :: nBoundaryFaces = 6*nSide*nSide
  integer,parameter :: negativeBCID = -7
  real(prec),parameter :: dx = 0.1_prec

  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)

  bcids(1:6) = negativeBCID
  call mesh%StructuredMesh(nPerTile,nPerTile,nPerTile,nTile,nTile,nTile,dx,dx,dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%tecplot_enabled = .false.

  ! (a) the count is global and identical on every rank
  if(modelobj%nUnmappedBoundaries /= nBoundaryFaces) then
    print*,"Error: rank ",mesh%decomp%rankId," expected the global count ",nBoundaryFaces, &
      " unmapped boundary faces, saw ",modelobj%nUnmappedBoundaries
    stop 1
  endif

  ! (b) every rank agrees on the representative id, and it survives being negative
  if(modelobj%unmappedBoundaryID /= negativeBCID) then
    print*,"Error: rank ",mesh%decomp%rankId," expected the unregistered bcid ",negativeBCID, &
      ", got ",modelobj%unmappedBoundaryID
    stop 1
  endif

  ! (c) re-tagging to a bcid LinearEuler3D registers clears the count everywhere
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)
  call modelobj%MapBoundaryConditions()
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: rank ",mesh%decomp%rankId," still reports ", &
      modelobj%nUnmappedBoundaries," unmapped boundary faces after re-tagging."
    stop 1
  endif
  if(modelobj%unmappedBoundaryID /= -1) then
    print*,"Error: rank ",mesh%decomp%rankId," reports bcid ",modelobj%unmappedBoundaryID, &
      " with a zero count; it should be reset to -1 when no rank holds an offender."
    stop 1
  endif

  if(mesh%decomp%rankId == 0) then
    print*,"dgmodel3d_unmapped_boundary_mpi : all assertions passed"
  endif

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram DGModel3D_Unmapped_Boundary_MPI
