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

module SELF_AdaptiveMesh_3D
!! Emit a solver-ready Mesh3D_t from an adaptively refined, 2:1-balanced octree forest,
!! the direct analogue of SELF_AdaptiveMesh_2D. This closes the 3-D adaptive-refinement
!! loop: with the indicator flagging elements, the forest mutating, and Balance2to1
!! balancing it, EmitMesh produces the nonconforming mesh - leaf geometry,
!! conforming-face connectivity, and a mortar table for the 2:1 hanging faces - that the
!! 3-D mortar solver machinery already handles.
!!
!! Each leaf of the forest becomes an element (in leaf-list order). For every leaf face
!! the FaceNeighbor query classifies the face and drives the emitted connectivity:
!!
!!   - domain boundary       -> sideInfo(3)=0, sideInfo(5)=base BC id
!!   - same-level leaf        -> conforming interior face (sideInfo(3)=nbr, (4)=10*face+flip)
!!   - one-level-finer face   -> this leaf is the BIG face of a 2:1 mortar; the four small
!!                               elements are the finer neighbour node's children on the
!!                               shared face
!!   - one-level-coarser face -> this leaf is a SMALL face; filled when its big face is
!!                               processed
!!
!! Mortar faces carry sideInfo(1)=mortar index and sideInfo(3)=sideInfo(5)=0 so the
!! conforming side-exchange machinery skips them, exactly as in the hand-built
!! SimpleMortarMesh. The mortar table follows the 14-row layout documented on Mesh3D_t
!! (big elem/face; small elem + 10*face+flip per sub-face in big-face quadrant order;
!! four sub-face global side ids), with the small-face flips inherited from the shared
!! face's flip and the quadrant-to-child pairing given by faceQuadPerm/childOfFace.
!!
!! Decomposition: every rank builds the same GLOBAL connectivity and mortar tables
!! deterministically from the (rank-replicated) forest, generates a fresh contiguous
!! decomposition over the leaf list - leaf-list order is Morton order within each root
!! tree, so contiguous ranges are space-filling-curve partitions - and stores only its
!! local slice of the element-sized arrays, exactly as the built-in mesh constructors do.
!! sideInfo(3) carries global element ids, nUniqueSides is the global side count, and
!! mortarInfo/nMortars are replicated in full with global ids on every rank, which is
!! what SideExchange/MortarExchange require. Repartitioning is implicit: each epoch's
!! emitted mesh is re-decomposed over the new leaf list, so equal-count partitions move
!! with the refinement.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_OctreeMesh_3D
  use SELF_RefinementPrimitives_3D,only:childOfFace,faceQuadPerm

  implicit none

contains

  subroutine EmitMesh(forest,baseMesh,outMesh)
    !! Build outMesh (a conforming-or-mortar Mesh3D_t) from a 2:1-balanced forest.
    !! baseMesh is the mesh the forest was initialised from (supplies BC metadata and the
    !! communicator; on nRanks > 1 the forest must be rank-replicated so every rank emits
    !! identical global tables). The forest must already be balanced (MaxLevelJump <= 1);
    !! EmitMesh does not mutate it.
    implicit none
    type(OctreeMesh3D),intent(in) :: forest
    type(Mesh3D),intent(in) :: baseMesh
    type(Mesh3D),intent(out) :: outMesh
    ! Local
    integer :: nEl,nGeo,nBCs,li,s,node,nbr,ns,nf,e,ne,k
    integer :: m,nMortar,gid,q,tq,cq,esq
    integer :: eFirst,eLast,nLocal
    type(Lagrange) :: geomInterp
    integer,allocatable :: leafIdx(:)
    integer,allocatable :: si(:,:,:)
    integer,allocatable :: minfo(:,:)
    integer :: gidq(1:4)
    real(prec),allocatable :: coords(:,:,:,:)

    if(forest%MaxLevelJump() > 1) then
      print*,__FILE__,':',__LINE__, &
        ' : Error : EmitMesh requires a 2:1-balanced forest; call Balance2to1 first.'
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
    allocate(si(1:5,1:6,1:nEl))
    si = 0
    allocate(minfo(1:14,1:6*nEl)) ! upper bound: at most one mortar per leaf face
    nMortar = 0
    gid = 0

    do li = 1,nEl
      node = forest%leaf(li)
      e = li
      do s = 1,6
        if(si(2,s,e) /= 0) cycle ! already filled (conforming partner, or small face of a mortar)

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
            ! Conforming same-level interior face; assign a shared global id to both. The
            ! reverse flip is the inverse permutation of nf (flips 5 and 7 are mutually
            ! inverse; all others are involutions).
            gid = gid+1
            si(2,s,e) = gid
            si(3,s,e) = ne
            si(4,s,e) = 10*ns+nf
            si(2,ns,ne) = gid
            si(3,ns,ne) = e
            si(4,ns,ne) = 10*s+invFlip3D(nf)
          else
            ! Neighbour is one level coarser: this is a SMALL face; its big face fills it later.
            cycle
          endif

        else
          ! Neighbour node is internal (finer) -> this leaf is the BIG face of a 2:1 mortar.
          nMortar = nMortar+1
          m = nMortar

          minfo(1,m) = e
          minfo(2,m) = s
          si(1,s,e) = m

          ! The four small elements are the finer neighbour's children on its face ns.
          ! Big-face quadrant q pairs with the neighbour-face quadrant faceQuadPerm(q,nf),
          ! and the flip stored with each small face is nf (big-face coordinates to
          ! small-face coordinates, the receiver-to-donor convention).
          do q = 1,4
            tq = faceQuadPerm(q,nf)
            cq = forest%child(childOfFace(tq,ns),nbr)
            esq = leafIdx(cq)
            gid = gid+1
            gidq(q) = gid
            minfo(2*q+1,m) = esq
            minfo(2*q+2,m) = 10*ns+nf
            minfo(10+q,m) = gid
            si(1,ns,esq) = m
            si(2,ns,esq) = gid
          enddo
          ! The big face shares sub-face 1's global side id.
          si(2,s,e) = gidq(1)
        endif
      enddo
    enddo

    ! ---- Allocate and populate the output mesh (fresh contiguous decomposition) ----
    ! Initialize on the base mesh's communicator so MPI is reused (not re-initialized) and
    ! the process-wide live-decomposition count stays correct across mesh lifetimes. The
    ! decomposition is regenerated over the (global) leaf list, and this rank stores only
    ! its contiguous slice eFirst:eLast, exactly as the built-in mesh constructors do.
    call outMesh%decomp%Init(comm=baseMesh%decomp%mpiComm)
    call outMesh%decomp%GenerateDecomposition(nEl,64*max(gid,1))
    eFirst = outMesh%decomp%offsetElem(outMesh%decomp%rankId+1)+1
    eLast = outMesh%decomp%offsetElem(outMesh%decomp%rankId+2)
    nLocal = eLast-eFirst+1

    call outMesh%Init(nGeo,nLocal,6*nLocal,8*nLocal,nBCs)
    outMesh%nGlobalElem = nEl
    outMesh%nUniqueSides = gid ! GLOBAL side count on every rank (the MPI tag stride)
    outMesh%quadrature = forest%quadrature

    ! Leaf geometry (rank-local leaves only).
    allocate(coords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1))
    do li = eFirst,eLast
      call forest%LeafCoords(li,geomInterp,coords)
      outMesh%nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,li-eFirst+1) = &
        coords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1)
    enddo

    ! Local slice of the global side table; sideInfo(3) keeps GLOBAL neighbour element
    ! ids, which is what SideExchange consumes (locality decided through
    ! decomp%elemToRank).
    outMesh%sideInfo(1:5,1:6,1:nLocal) = si(1:5,1:6,eFirst:eLast)
    outMesh%globalNodeIDs = 0 ! node ids are unused by the solver (flips are set directly)
    outMesh%elemInfo = 0

    ! Boundary-condition metadata (replicated on every rank).
    if(nBCs > 0) then
      outMesh%BCType(1:4,1:nBCs) = baseMesh%BCType(1:4,1:nBCs)
      do k = 1,nBCs
        outMesh%BCNames(k) = baseMesh%BCNames(k)
      enddo
    endif

    ! Material table: each leaf inherits its root element's material (rootMaterial is
    ! global on the forest, so this works for any decomposition of the emitted mesh).
    outMesh%nMaterials = baseMesh%nMaterials
    if(allocated(outMesh%materialNames)) deallocate(outMesh%materialNames)
    allocate(outMesh%materialNames(1:baseMesh%nMaterials))
    outMesh%materialNames(1:baseMesh%nMaterials) = baseMesh%materialNames(1:baseMesh%nMaterials)
    do li = eFirst,eLast
      outMesh%elemMaterial(li-eFirst+1) = forest%rootMaterial(forest%rootElem(forest%leaf(li)))
    enddo

    ! Mortar table.
    outMesh%nMortars = nMortar
    if(associated(outMesh%mortarInfo)) deallocate(outMesh%mortarInfo)
    if(nMortar > 0) then
      allocate(outMesh%mortarInfo(1:14,1:nMortar))
      outMesh%mortarInfo(1:14,1:nMortar) = minfo(1:14,1:nMortar)
    else
      outMesh%mortarInfo => null()
    endif

    deallocate(leafIdx,si,minfo,coords)
    call geomInterp%Free()

    call outMesh%UpdateDevice()

  endsubroutine EmitMesh

  pure function invFlip3D(f) result(fi)
    !! Inverse of a SELF face flip: the flip that maps donor-face indices back to
    !! receiver-face indices. Flips 0..4 and 6 are involutions; 5 and 7 invert each other.
    implicit none
    integer,intent(in) :: f
    integer :: fi
    integer,parameter :: inv(0:7) = [0,1,2,3,4,7,6,5]
    fi = inv(f)
  endfunction invFlip3D

endmodule SELF_AdaptiveMesh_3D
