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

program DGModel3D_Unmapped_Boundary
!! The 3-D half of the issue #180 regression. MapBoundaryConditions is a separate routine per
!! dimension - the 3-D one walks six faces per element rather than four - so the 2-D test does
!! not cover it.
!!
!! Assertions:
!!  (a) tagging every domain face with a bcid no model registers counts every one of them,
!!      and reports that id;
!!  (b) the id reported is the offending one even when it is negative, which a max() against
!!      a sentinel would get wrong;
!!  (c) re-tagging to a registered bcid and re-mapping clears the count and re-arms the
!!      one-shot warning.

  use self_data
  use self_lineareuler3d
  use self_mesh_3d

  implicit none

  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  integer,parameter :: nPerTile = 2
  integer,parameter :: nTile = 1
  !! n^3 elements, so 6 faces of n^2 each on the domain boundary.
  integer,parameter :: nSide = nPerTile*nTile
  integer,parameter :: nBoundaryFaces = 6*nSide*nSide
  integer,parameter :: bogusBCID = 999
  integer,parameter :: negativeBCID = -7 !! legal, and below the -1 "absent" marker
  real(prec),parameter :: dx = 0.1_prec

  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)

  bcids(1:6) = bogusBCID
  call mesh%StructuredMesh(nPerTile,nPerTile,nPerTile,nTile,nTile,nTile,dx,dx,dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Init runs AdditionalInit (LinearEuler3D registers radiation only) and then
  ! MapBoundaryConditions, so the count is available immediately.
  call modelobj%Init(mesh,geometry)
  modelobj%tecplot_enabled = .false.

  ! (a) every domain face is unmapped, and the offending bcid is named
  if(modelobj%nUnmappedBoundaries /= nBoundaryFaces) then
    print*,"Error: expected ",nBoundaryFaces," unmapped boundary faces, counted ", &
      modelobj%nUnmappedBoundaries
    stop 1
  endif
  if(modelobj%unmappedBoundaryID /= bogusBCID) then
    print*,"Error: expected the unregistered bcid ",bogusBCID," to be reported, got ", &
      modelobj%unmappedBoundaryID
    stop 1
  endif
  if(modelobj%unmappedBoundariesReported) then
    print*,"Error: the warning was marked reported before anything reported it."
    stop 1
  endif

  ! (b) a negative bcid is legal and must be reported as itself. Taking a maximum against
  ! the -1 initializer would report -1 here, which reads as "no offender".
  call mesh%ResetBoundaryConditionType(negativeBCID)
  call modelobj%MapBoundaryConditions()
  if(modelobj%nUnmappedBoundaries /= nBoundaryFaces) then
    print*,"Error: expected ",nBoundaryFaces," unmapped boundary faces for a negative bcid, ", &
      "counted ",modelobj%nUnmappedBoundaries
    stop 1
  endif
  if(modelobj%unmappedBoundaryID /= negativeBCID) then
    print*,"Error: expected the negative bcid ",negativeBCID," to be reported, got ", &
      modelobj%unmappedBoundaryID
    stop 1
  endif

  ! The warning path itself. ReportUnmappedBoundaries reads the counters and prints, so it
  ! exercises the one-shot latch without a time step - and on a GPU build the device copies
  ! of bc%elements/bc%sides are refreshed only by Init and Regrid, never by a bare re-map.
  call modelobj%ReportUnmappedBoundaries()
  if(.not. modelobj%unmappedBoundariesReported) then
    print*,"Error: ReportUnmappedBoundaries did not consume the warning."
    stop 1
  endif

  ! (c) re-tagging to a bcid LinearEuler3D does register clears the count
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)
  call modelobj%MapBoundaryConditions()
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: re-tagging to a registered bcid left ",modelobj%nUnmappedBoundaries, &
      " unmapped boundary faces."
    stop 1
  endif
  if(modelobj%unmappedBoundariesReported) then
    print*,"Error: MapBoundaryConditions did not re-arm the one-shot warning."
    stop 1
  endif

  call modelobj%ReportUnmappedBoundaries()
  if(modelobj%unmappedBoundariesReported) then
    print*,"Error: the warning fired on a mesh with no unmapped boundary faces."
    stop 1
  endif

  print*,"dgmodel3d_unmapped_boundary : all assertions passed"

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram DGModel3D_Unmapped_Boundary
