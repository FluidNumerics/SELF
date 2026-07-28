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

module SELF_AdaptiveMesh_2D
!! Emit a solver-ready Mesh2D_t from an adaptively refined, 2:1-balanced quad-forest (AMR
!! Stage 4b). This closes the adaptive-refinement loop: with the Stage-1 indicator flagging
!! elements, Stage 2b mutating the forest, and Stage 4a balancing it, EmitMesh produces the
!! nonconforming mesh - leaf geometry, conforming-side connectivity, and a mortar table for the
!! 2:1 hanging faces - that the existing mortar solver machinery already handles.
!!
!! Each leaf of the forest becomes an element (in leaf-list order). For every leaf face the
!! Stage-4a FaceNeighbor query classifies the face and drives the emitted connectivity:
!!
!!   - domain boundary      -> sideInfo(3)=0, sideInfo(5)=base BC id
!!   - same-level leaf       -> conforming interior side (sideInfo(3)=neighbour, (4)=10*side+flip)
!!   - one-level-finer face  -> this leaf is the BIG side of a 2:1 mortar; the two small elements
!!                              are the finer neighbour node's children on the shared face
!!   - one-level-coarser face-> this leaf is a SMALL side; filled when its big side is processed
!!
!! Mortar sides carry sideInfo(1)=mortar index and sideInfo(3)=sideInfo(5)=0 so the conforming
!! side-exchange machinery skips them, exactly as in the hand-built SimpleMortarMesh. The mortar
!! table follows the same 8-row layout (big elem/side; small elem + 10*side+flip per sub-edge;
!! two sub-edge global side ids), with sub-edge 1 covering the big edge coordinate [-1,0] and
!! sub-edge 2 covering [0,1], and the small-side flips inherited from the base face flip.
!!
!! Serial only, matching the rest of the AMR stack (dynamic MPI re-partitioning is Stage 5): the
!! output replicates the input mesh's single-rank decomposition rather than re-decomposing.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_QuadTreeMesh_2D
  use SELF_RefinementPrimitives_2D,only:childOfSide

  implicit none

contains

  subroutine EmitMesh(forest,baseMesh,outMesh)
    !! Build outMesh (a conforming-or-mortar Mesh2D_t) from a 2:1-balanced forest. baseMesh is the
    !! mesh the forest was initialised from (supplies BC/material metadata and the serial MPI
    !! state). The forest must already be balanced (MaxLevelJump <= 1); EmitMesh does not mutate it.
    implicit none
    type(QuadTreeMesh2D),intent(in) :: forest
    type(Mesh2D),intent(in) :: baseMesh
    type(Mesh2D),intent(out) :: outMesh
    ! Local
    integer :: nEl,nGeo,nBCs,li,s,node,nbr,ns,nf,e,ne,k
    integer :: m,nMortar,gid,gidA,gidB,t1,t2,c1,c2,es1,es2
    type(Lagrange) :: geomInterp
    integer,allocatable :: leafIdx(:)
    integer,allocatable :: si(:,:,:)
    integer,allocatable :: minfo(:,:)
    real(prec),allocatable :: coords(:,:,:)

    if(forest%MaxLevelJump() > 1) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : EmitMesh requires a 2:1-balanced forest; call Balance2to1 first.'
      stop 1
    endif
    if(baseMesh%decomp%nRanks > 1) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : EmitMesh is serial-only pending AMR Stage 5 (MPI repartitioning).'
      stop 1
    endif

    nEl = forest%nLeaves
    nGeo = forest%nGeo
    nBCs = baseMesh%nBCs

    call geomInterp%Init(nGeo,forest%quadrature,nGeo,forest%quadrature)

    ! node id -> emitted element id (leaf-list order); 0 for non-leaf nodes.
    allocate(leafIdx(1:forest%nNodes))
    leafIdx = 0
    do li = 1,nEl
      leafIdx(forest%leaf(li)) = li
    enddo

    ! ---- Classify every leaf face; build sideInfo and the mortar table ----
    allocate(si(1:5,1:4,1:nEl))
    si = 0
    allocate(minfo(1:8,1:4*nEl)) ! upper bound: at most one mortar per leaf face
    nMortar = 0
    gid = 0

    do li = 1,nEl
      node = forest%leaf(li)
      e = li
      do s = 1,4
        if(si(2,s,e) /= 0) cycle ! already filled (conforming partner, or small side of a mortar)

        call forest%FaceNeighbor(node,s,nbr,ns,nf)

        if(nbr == 0) then
          ! Physical domain boundary.
          gid = gid+1
          si(2,s,e) = gid
          si(5,s,e) = forest%rootBC(s,forest%rootElem(node))

        elseif(forest%child(1,nbr) == 0) then
          ! Neighbour is a leaf.
          ne = leafIdx(nbr)
          if(forest%level(nbr) == forest%level(node)) then
            ! Conforming same-level interior side; assign a shared global id to both.
            gid = gid+1
            si(2,s,e) = gid
            si(3,s,e) = ne
            si(4,s,e) = 10*ns+nf
            si(2,ns,ne) = gid
            si(3,ns,ne) = e
            si(4,ns,ne) = 10*s+nf
          else
            ! Neighbour is one level coarser: this is a SMALL side; its big side fills it later.
            cycle
          endif

        else
          ! Neighbour node is internal (finer) -> this leaf is the BIG side of a 2:1 mortar.
          nMortar = nMortar+1
          m = nMortar
          ! The two small elements are the finer neighbour's children on its side ns.
          ! Big edge coordinate [-1,0] (sub-edge 1) maps to neighbour sub-position t1, and
          ! [0,1] (sub-edge 2) to t2, reversed when the shared face has flip 1.
          if(nf == 0) then
            t1 = 1; t2 = 2
          else
            t1 = 2; t2 = 1
          endif
          c1 = forest%child(childOfSide(t1,ns),nbr)
          c2 = forest%child(childOfSide(t2,ns),nbr)
          es1 = leafIdx(c1)
          es2 = leafIdx(c2)

          gid = gid+1; gidA = gid ! sub-edge 1 (shared by the big side and small 1)
          gid = gid+1; gidB = gid ! sub-edge 2

          ! Big side.
          si(1,s,e) = m
          si(2,s,e) = gidA
          ! Small sides (both on neighbour local side ns).
          si(1,ns,es1) = m
          si(2,ns,es1) = gidA
          si(1,ns,es2) = m
          si(2,ns,es2) = gidB

          minfo(1,m) = e
          minfo(2,m) = s
          minfo(3,m) = es1
          minfo(4,m) = 10*ns+nf
          minfo(5,m) = es2
          minfo(6,m) = 10*ns+nf
          minfo(7,m) = gidA
          minfo(8,m) = gidB
        endif
      enddo
    enddo

    ! ---- Allocate and populate the output mesh (serial decomposition, as UniformRefineMesh) ----
    outMesh%decomp%mpiComm = baseMesh%decomp%mpiComm
    outMesh%decomp%mpiPrec = baseMesh%decomp%mpiPrec
    outMesh%decomp%rankId = baseMesh%decomp%rankId
    outMesh%decomp%nRanks = baseMesh%decomp%nRanks
    outMesh%decomp%mpiEnabled = baseMesh%decomp%mpiEnabled
    outMesh%decomp%initialized = .true.
    allocate(outMesh%decomp%offsetElem(1:outMesh%decomp%nRanks+1))
    call outMesh%decomp%GenerateDecomposition(nEl,64*max(gid,1))

    call outMesh%Init(nGeo,nEl,4*nEl,4*nEl,nBCs)
    outMesh%nUniqueSides = gid
    outMesh%quadrature = forest%quadrature

    ! Leaf geometry.
    allocate(coords(1:2,1:nGeo+1,1:nGeo+1))
    do li = 1,nEl
      call forest%LeafCoords(li,geomInterp,coords)
      outMesh%nodeCoords(1:2,1:nGeo+1,1:nGeo+1,li) = coords(1:2,1:nGeo+1,1:nGeo+1)
    enddo

    outMesh%sideInfo(1:5,1:4,1:nEl) = si(1:5,1:4,1:nEl)
    outMesh%globalNodeIDs = 0 ! node ids are unused by the solver (flips are set directly)

    ! Boundary-condition metadata.
    if(nBCs > 0) then
      outMesh%BCType(1:4,1:nBCs) = baseMesh%BCType(1:4,1:nBCs)
      do k = 1,nBCs
        outMesh%BCNames(k) = baseMesh%BCNames(k)
      enddo
    endif

    ! Material table: each leaf inherits its root element's material.
    outMesh%nMaterials = baseMesh%nMaterials
    if(allocated(outMesh%materialNames)) deallocate(outMesh%materialNames)
    allocate(outMesh%materialNames(1:baseMesh%nMaterials))
    outMesh%materialNames(1:baseMesh%nMaterials) = baseMesh%materialNames(1:baseMesh%nMaterials)
    do li = 1,nEl
      outMesh%elemMaterial(li) = baseMesh%elemMaterial(forest%rootElem(forest%leaf(li)))
    enddo

    ! Mortar table.
    outMesh%nMortars = nMortar
    if(associated(outMesh%mortarInfo)) deallocate(outMesh%mortarInfo)
    if(nMortar > 0) then
      allocate(outMesh%mortarInfo(1:8,1:nMortar))
      outMesh%mortarInfo(1:8,1:nMortar) = minfo(1:8,1:nMortar)
    else
      outMesh%mortarInfo => null()
    endif

    deallocate(leafIdx,si,minfo,coords)
    call geomInterp%Free()

    call outMesh%UpdateDevice()

  endsubroutine EmitMesh

endmodule SELF_AdaptiveMesh_2D
