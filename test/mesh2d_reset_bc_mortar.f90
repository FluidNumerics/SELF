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

program Mesh2D_Reset_BC_Mortar
!! ResetBoundaryConditionType must tag physical boundaries and nothing else.
!!
!! Mortar sides carry sideInfo(3) = 0 exactly like a physical boundary - sideInfo(1) holds
!! the mortar index - so a routine that keys only on sideInfo(3) writes a bcid onto the
!! nonconforming interface, and MapBoundaryConditions then hands those sides to a boundary
!! condition that overwrites what the mortar exchange just wrote.
!!
!! The exclusion cannot key on sideInfo(1) alone. That row is the HOPr "side type", copied
!! verbatim by Read_HOPr, so it can be nonzero on an ordinary face; an ungated test would
!! silently skip real boundaries on a HOPr mesh. It is a mortar marker only when the mesh
!! actually has mortars.
!!
!! Assertions:
!!  (a) on a mortar mesh, every mortar side keeps sideInfo(5) = 0 after a reset;
!!  (b) on that same mesh, every non-mortar physical boundary does get the new bcid;
!!  (c) on a mesh with no mortars, a nonzero sideInfo(1) does NOT suppress the tagging -
!!      this is the assertion that fails if the nMortars gate is ever dropped.

  use SELF_Constants
  use SELF_Mesh_2D

  implicit none

  real(prec),parameter :: dx = 0.1_prec
  integer,parameter :: initialBCID = 101
  integer,parameter :: resetBCID = 102
  integer,parameter :: hoprSideType = 7 ! stand-in for a HOPr "side type" on an ordinary face

  ! Two mesh objects rather than one reused object: SELF finalizes MPI when the last live
  ! domain decomposition is freed, so a Free followed by another Init in the same program
  ! calls into MPI after MPI_Finalize. Both are freed at the end, as in
  ! domaindecomposition_two_meshes.
  type(Mesh2D),target :: mortarMesh,structMesh
  integer :: bcids(1:4)
  integer :: m,k,e,s,iEl,iSide
  integer :: nMortarSides,nBoundarySides
  logical,allocatable :: isMortar(:,:)

  ! ---------------------------------------------------------------------------------------
  ! (a) and (b) : a mesh that really has mortars
  ! ---------------------------------------------------------------------------------------
  bcids(1:4) = [initialBCID,initialBCID,initialBCID,initialBCID] ! south, east, north, west
  call mortarMesh%SimpleMortarMesh(dx,bcids)

  if(mortarMesh%nMortars <= 0) then
    print*,"Error: SimpleMortarMesh produced no mortars; the guard is not exercised."
    stop 1
  endif

  ! Mark the local sides that take part in a mortar, from mortarInfo: the big side is
  ! (1,2) and the two small sides are (2*k+1, 2*k+2/10) for k = 1,2. Element ids in
  ! mortarInfo are global; this test is serial, so they are already local.
  allocate(isMortar(1:4,1:mortarMesh%nElem))
  isMortar = .false.
  do m = 1,mortarMesh%nMortars
    e = mortarMesh%mortarInfo(1,m)
    s = mortarMesh%mortarInfo(2,m)
    isMortar(s,e) = .true.
    do k = 1,2
      e = mortarMesh%mortarInfo(2*k+1,m)
      s = mortarMesh%mortarInfo(2*k+2,m)/10
      isMortar(s,e) = .true.
    enddo
  enddo

  call mortarMesh%ResetBoundaryConditionType(resetBCID)

  nMortarSides = 0
  nBoundarySides = 0
  do iEl = 1,mortarMesh%nElem
    do iSide = 1,4

      if(isMortar(iSide,iEl)) then
        nMortarSides = nMortarSides+1
        ! (a) a mortar side is an interior face and must never be tagged
        if(mortarMesh%sideInfo(5,iSide,iEl) /= 0) then
          print*,"Error: mortar side (",iSide,",",iEl,") was tagged with bcid ", &
            mortarMesh%sideInfo(5,iSide,iEl)
          stop 1
        endif

      elseif(mortarMesh%sideInfo(3,iSide,iEl) == 0) then
        nBoundarySides = nBoundarySides+1
        ! (b) a genuine physical boundary must be tagged; the guard must not over-reach
        if(mortarMesh%sideInfo(5,iSide,iEl) /= resetBCID) then
          print*,"Error: boundary side (",iSide,",",iEl,") holds bcid ", &
            mortarMesh%sideInfo(5,iSide,iEl)," but should hold ",resetBCID
          stop 1
        endif
      endif

    enddo
  enddo

  if(nMortarSides /= 3*mortarMesh%nMortars) then
    print*,"Error: expected ",3*mortarMesh%nMortars," mortar sides, marked ",nMortarSides
    stop 1
  endif
  if(nBoundarySides <= 0) then
    print*,"Error: no physical boundaries found; assertion (b) is vacuous."
    stop 1
  endif

  deallocate(isMortar)

  ! ---------------------------------------------------------------------------------------
  ! (c) : no mortars, but a nonzero sideInfo(1), as a HOPr mesh carries
  ! ---------------------------------------------------------------------------------------
  bcids(1:4) = [initialBCID,initialBCID,initialBCID,initialBCID]
  call structMesh%StructuredMesh(2,2,1,1,dx,dx,bcids)

  if(structMesh%nMortars /= 0) then
    print*,"Error: the structured mesh reports ",structMesh%nMortars," mortars; (c) is not the "// &
      "no-mortar case it is meant to be."
    stop 1
  endif

  ! A structured mesh leaves sideInfo(1) at 0. Imitate what Read_HOPr would have left
  ! there, on every side, so the reset has to rely on the nMortars gate to ignore it.
  structMesh%sideInfo(1,1:4,1:structMesh%nElem) = hoprSideType

  call structMesh%ResetBoundaryConditionType(resetBCID)

  nBoundarySides = 0
  do iEl = 1,structMesh%nElem
    do iSide = 1,4
      if(structMesh%sideInfo(3,iSide,iEl) == 0) then
        nBoundarySides = nBoundarySides+1
        if(structMesh%sideInfo(5,iSide,iEl) /= resetBCID) then
          print*,"Error: a nonzero sideInfo(1) suppressed the tagging of boundary side (", &
            iSide,",",iEl,") on a mesh with no mortars. The nMortars gate is missing."
          stop 1
        endif
      endif
    enddo
  enddo

  if(nBoundarySides /= 8) then
    print*,"Error: expected 8 boundary sides on a 2x2 structured mesh, found ",nBoundarySides
    stop 1
  endif

  ! Freed together at the end: the last Free is what finalizes MPI, so neither mesh may be
  ! released while the other still has to be built.
  call mortarMesh%Free()
  call structMesh%Free()

  print*,"mesh2d_reset_bc_mortar : all assertions passed"

endprogram Mesh2D_Reset_BC_Mortar
