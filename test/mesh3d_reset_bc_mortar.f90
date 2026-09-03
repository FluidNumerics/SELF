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

program Mesh3D_Reset_BC_Mortar
!! The 3-D twin of mesh2d_reset_bc_mortar, and the regression for the behaviour change this
!! branch makes in 3-D.
!!
!! ResetBoundaryConditionType_Mesh3D_t used to exclude mortar faces with a bare
!! sideInfo(1) == 0. That row is the HOPr "side type", copied verbatim by Read_HOPr, so it can
!! be nonzero on an ordinary face - which made the whole routine a silent no-op on every HOPr
!! 3-D mesh. The exclusion is now gated on the mesh actually having mortars.
!!
!! Assertions:
!!  (a) on a mortar mesh, every mortar face keeps sideInfo(5) = 0 after a reset;
!!  (b) on that same mesh, every non-mortar physical boundary does get the new bcid;
!!  (c) on a mesh with no mortars, a nonzero sideInfo(1) does NOT suppress the tagging -
!!      this is the assertion that fails if the nMortars gate is ever dropped.

  use SELF_Constants
  use SELF_Mesh_3D

  implicit none

  real(prec),parameter :: dx = 0.1_prec
  integer,parameter :: initialBCID = 101
  integer,parameter :: resetBCID = 102
  integer,parameter :: hoprSideType = 7 ! stand-in for a HOPr "side type" on an ordinary face

  ! Two mesh objects rather than one reused object: SELF finalizes MPI when the last live
  ! domain decomposition is freed, so a Free followed by another Init in the same program
  ! calls into MPI after MPI_Finalize.
  type(Mesh3D),target :: mortarMesh,structMesh
  integer :: bcids(1:6)
  integer :: m,q,e,s,iEl,iSide
  integer :: nMortarFaces,nBoundaryFaces
  logical,allocatable :: isMortar(:,:)

  ! ---------------------------------------------------------------------------------------
  ! (a) and (b) : a mesh that really has mortars
  ! ---------------------------------------------------------------------------------------
  bcids(1:6) = initialBCID
  call mortarMesh%SimpleMortarMesh(dx,bcids)

  if(mortarMesh%nMortars <= 0) then
    print*,"Error: SimpleMortarMesh produced no mortars; the guard is not exercised."
    stop 1
  endif

  ! Mark the local faces taking part in a mortar, from mortarInfo: the big face is (1,2) and
  ! the four small faces are (2*q+1, 2*q+2/10) for q = 1..4. Element ids there are global;
  ! this test is serial, so they are already local.
  allocate(isMortar(1:6,1:mortarMesh%nElem))
  isMortar = .false.
  do m = 1,mortarMesh%nMortars
    e = mortarMesh%mortarInfo(1,m)
    s = mortarMesh%mortarInfo(2,m)
    isMortar(s,e) = .true.
    do q = 1,4
      e = mortarMesh%mortarInfo(2*q+1,m)
      s = mortarMesh%mortarInfo(2*q+2,m)/10
      isMortar(s,e) = .true.
    enddo
  enddo

  call mortarMesh%ResetBoundaryConditionType(resetBCID)

  nMortarFaces = 0
  nBoundaryFaces = 0
  do iEl = 1,mortarMesh%nElem
    do iSide = 1,6

      if(isMortar(iSide,iEl)) then
        nMortarFaces = nMortarFaces+1
        ! (a) a mortar face is an interior face and must never be tagged
        if(mortarMesh%sideInfo(5,iSide,iEl) /= 0) then
          print*,"Error: mortar face (",iSide,",",iEl,") was tagged with bcid ", &
            mortarMesh%sideInfo(5,iSide,iEl)
          stop 1
        endif

      elseif(mortarMesh%sideInfo(3,iSide,iEl) == 0) then
        nBoundaryFaces = nBoundaryFaces+1
        ! (b) a genuine physical boundary must be tagged; the guard must not over-reach
        if(mortarMesh%sideInfo(5,iSide,iEl) /= resetBCID) then
          print*,"Error: boundary face (",iSide,",",iEl,") holds bcid ", &
            mortarMesh%sideInfo(5,iSide,iEl)," but should hold ",resetBCID
          stop 1
        endif
      endif

    enddo
  enddo

  if(nMortarFaces /= 5*mortarMesh%nMortars) then
    print*,"Error: expected ",5*mortarMesh%nMortars," mortar faces, marked ",nMortarFaces
    stop 1
  endif
  if(nBoundaryFaces <= 0) then
    print*,"Error: no physical boundaries found; assertion (b) is vacuous."
    stop 1
  endif

  deallocate(isMortar)

  ! ---------------------------------------------------------------------------------------
  ! (c) : no mortars, but a nonzero sideInfo(1), as a HOPr mesh carries
  ! ---------------------------------------------------------------------------------------
  bcids(1:6) = initialBCID
  call structMesh%StructuredMesh(2,2,2,1,1,1,dx,dx,dx,bcids)

  if(structMesh%nMortars /= 0) then
    print*,"Error: the structured mesh reports ",structMesh%nMortars," mortars; (c) is not "// &
      "the no-mortar case it is meant to be."
    stop 1
  endif

  ! A structured mesh leaves sideInfo(1) at 0. Imitate what Read_HOPr would have left there,
  ! on every face, so the reset has to rely on the nMortars gate to ignore it.
  structMesh%sideInfo(1,1:6,1:structMesh%nElem) = hoprSideType

  call structMesh%ResetBoundaryConditionType(resetBCID)

  nBoundaryFaces = 0
  do iEl = 1,structMesh%nElem
    do iSide = 1,6
      if(structMesh%sideInfo(3,iSide,iEl) == 0) then
        nBoundaryFaces = nBoundaryFaces+1
        if(structMesh%sideInfo(5,iSide,iEl) /= resetBCID) then
          print*,"Error: a nonzero sideInfo(1) suppressed the tagging of boundary face (", &
            iSide,",",iEl,") on a mesh with no mortars. The nMortars gate is missing."
          stop 1
        endif
      endif
    enddo
  enddo

  if(nBoundaryFaces /= 24) then
    print*,"Error: expected 24 boundary faces on a 2x2x2 structured mesh, found ", &
      nBoundaryFaces
    stop 1
  endif

  ! Freed together at the end: the last Free is what finalizes MPI.
  call mortarMesh%Free()
  call structMesh%Free()

  print*,"mesh3d_reset_bc_mortar : all assertions passed"

endprogram Mesh3D_Reset_BC_Mortar
