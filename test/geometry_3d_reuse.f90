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

program test

  integer :: exit_code

  exit_code = geometry_3d_reuse()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function geometry_3d_reuse() result(r)
    !! Validates the premise of the AMR Stage 6c geometry reuse: that an element the adaptation did
    !! not change has bit-for-bit identical geometry in the new epoch, so it may be copied forward
    !! instead of regenerated.
    !!
    !! The claim rests on two properties of the implementation, and this test is the direct check
    !! of their combination:
    !!
    !!   1. A leaf's mesh node coordinates come from LeafCoords, a pure function of the root
    !!      element's coordinates, the leaf's level and its octant path, evaluated in a fixed
    !!      order. Root coordinates are never mutated and node ids, levels and octants are stable
    !!      across forest mutations.
    !!   2. Per-element geometry generation is strictly element-local - GenerateFromMesh,
    !!      CalculateMetricTerms and CalculateContravariantBasis contain no neighbour coupling, no
    !!      side pairing and no reduction, and the normal-orientation convention is element-local.
    !!
    !! Together these say a kept element's whole geometry block is reproducible EXACTLY, not merely
    !! to round-off. The test therefore asserts exact equality: anything less means one of the two
    !! properties has been broken, and a tolerance would hide it.
    !!
    !! Unlike the 2-D test, which drives the reuse path through the transfer plan, no 3-D transfer
    !! plan (or AMR controller) exists yet, so the reuse composition is exercised directly the way
    !! the 2-D controller composes it: kept elements are identified by matching leaf node ids
    !! between epochs and copied forward with CopyElements, while the changed elements are
    !! generated COMPACTED through GenerateFromNodeCoords and scattered into place.
    !!
    !! What is checked, on a two-epoch adaptation of a structured base mesh:
    !!
    !!   (a) For every kept element, the geometry assembled by CopyElements from the previous
    !!       epoch equals - bit for bit - a full GenerateFromMesh reference built for the SAME
    !!       new mesh. This is the reuse predicate itself.
    !!   (b) Jacobians are strictly positive on the assembled geometry.
    !!   (c) The discrete metric identity sum_i d(J dsdx(i,j))/dxi_i = 0 holds per element for each
    !!       j, which the rest of the suite never checks. A geometry assembled from mismatched
    !!       pieces would violate it even where volumes and global integrals still balanced.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_OctreeMesh_3D
    use SELF_AdaptiveMesh_3D
    implicit none

    integer,parameter :: N = 3 ! polynomial degree
    integer,parameter :: nElemX = 2 ! base mesh is nElemX x nElemX x nElemX
    real(prec),parameter :: metricTol = 1.0e-9_prec

    type(Lagrange),target :: interp
    type(Lagrange),pointer :: interpPtr
    type(Mesh3D),target :: baseMesh
    type(Mesh3D),target :: meshOld,meshNew
    type(SEMHex),target :: geomOld,geomNew,geomRef,genGeom
    type(OctreeMesh3D) :: forest
    integer :: bcids(1:6)
    integer,allocatable :: flag(:),oldLeaf(:)
    integer,allocatable :: srcIdx(:),dstIdx(:),genIdx(:)
    real(prec),allocatable :: genCoords(:,:,:,:,:)
    integer :: li,lo,nOld,nCopy,nGen,nBad,iel,j,k,nGeo
    real(prec) :: minJ

    r = 0

    bcids(1:6) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED, &
                  SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]
    call baseMesh%StructuredMesh(nElemX,nElemX,nElemX,1,1,1, &
                                 1.0_prec/real(nElemX,prec),1.0_prec/real(nElemX,prec), &
                                 1.0_prec/real(nElemX,prec),bcids)
    call interp%Init(N=N,controlNodeType=GAUSS,M=N,targetNodeType=UNIFORM)
    interpPtr => interp

    call forest%Init(baseMesh)

    ! ---- Epoch 1: refine a couple of elements, emit, generate geometry ----
    allocate(flag(1:forest%nLeaves))
    flag = OCTREE_KEEP
    flag(1) = OCTREE_REFINE
    flag(min(6,forest%nLeaves)) = OCTREE_REFINE
    call forest%AdaptFromFlags(flag)
    call forest%Balance2to1()
    deallocate(flag)

    call EmitMesh(forest,baseMesh,meshOld)
    call geomOld%Init(interpPtr,meshOld%nElem)
    call geomOld%GenerateFromMesh(meshOld)

    ! ---- Epoch 2: refine elsewhere, so most leaves are carried over unchanged ----
    nOld = forest%nLeaves
    allocate(oldLeaf(1:nOld))
    oldLeaf(1:nOld) = forest%leaf(1:nOld)

    allocate(flag(1:nOld))
    flag = OCTREE_KEEP
    flag(nOld) = OCTREE_REFINE
    call forest%AdaptFromFlags(flag)
    call forest%Balance2to1()
    deallocate(flag)

    call EmitMesh(forest,baseMesh,meshNew)
    nGeo = meshNew%nGeo

    ! Full reference for the SAME new mesh: what the reuse path must reproduce exactly.
    call geomRef%Init(interpPtr,meshNew%nElem)
    call geomRef%GenerateFromMesh(meshNew)

    print*,"epoch 1 elements =",meshOld%nElem,", epoch 2 elements =",meshNew%nElem

    ! ---- Assemble the epoch-2 geometry incrementally ----
    ! Elements are emitted in leaf-list order, and node ids are stable across forest
    ! mutations, so a new leaf whose node id appears in the previous epoch's leaf list is
    ! the same element and its geometry may be copied forward; everything else is generated.
    call geomNew%Init(interpPtr,meshNew%nElem)

    allocate(srcIdx(1:forest%nLeaves),dstIdx(1:forest%nLeaves),genIdx(1:forest%nLeaves))
    nCopy = 0
    nGen = 0
    do li = 1,forest%nLeaves
      lo = 0
      do k = 1,nOld
        if(oldLeaf(k) == forest%leaf(li)) then
          lo = k
          exit
        endif
      enddo
      if(lo > 0) then
        nCopy = nCopy+1
        srcIdx(nCopy) = lo
        dstIdx(nCopy) = li
      else
        nGen = nGen+1
        genIdx(nGen) = li
      endif
    enddo

    ! Generate the changed elements, compacted, so the generation loops run over nGen elements
    ! instead of all of them. Their geometry is then scattered into place.
    if(nGen > 0) then
      allocate(genCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nGen))
      do k = 1,nGen
        genCoords(1:3,:,:,:,k) = meshNew%nodeCoords(1:3,:,:,:,genIdx(k))
      enddo

      call genGeom%Init(interpPtr,nGen)
      call genGeom%GenerateFromNodeCoords(genCoords,nGeo,meshNew%quadrature,nGen)
      call genGeom%x%UpdateDevice()
      call genGeom%x%BoundaryInterp()
      call genGeom%x%UpdateHost()
      call genGeom%CalculateMetricTerms()

      do k = 1,nGen
        srcIdx(nCopy+k) = k
        dstIdx(nCopy+k) = genIdx(k)
      enddo
      call geomNew%CopyElements(genGeom,srcIdx(nCopy+1:),dstIdx(nCopy+1:),nGen)
      deallocate(genCoords)
    endif

    ! Carry the unchanged elements across from the previous epoch's geometry.
    if(nCopy > 0) then
      call geomNew%CopyElements(geomOld,srcIdx,dstIdx,nCopy)
    endif

    call geomNew%UploadGeometry()

    ! ---- (a) kept elements must be bit-for-bit reproducible from the previous epoch ----
    nBad = 0
    do k = 1,nCopy
      li = dstIdx(k)
      if(.not. SameElement(geomNew,li,geomRef,li)) then
        nBad = nBad+1
        if(nBad <= 5) then
          print*,"FAIL: kept element",li,"differs from the full-generation reference", &
            " (epoch-1 source",srcIdx(k),")"
        endif
      endif
    enddo

    print*,"kept (copied) elements =",nCopy," of ",meshNew%nElem, &
      " ; mismatched =",nBad
    if(nCopy == 0) then
      print*,"FAIL: the two-epoch setup produced no kept elements, so nothing was tested"
      r = 1
    endif
    if(nBad /= 0) then
      print*,"FAIL: geometry of kept elements is not reproducible; reuse would be incorrect"
      r = 1
    endif

    ! ---- (b) positive Jacobians ----
    minJ = minval(geomNew%J%interior)
    print*,"min(J) =",minJ
    if(minJ <= 0.0_prec) then
      print*,"FAIL: non-positive Jacobian on the assembled geometry"
      r = 1
    endif

    ! ---- (c) discrete metric identity, per element ----
    ! sum_i d(J dsdx(i,j))/dxi_i = 0 for each j. On a straight-sided mesh the metric terms are
    ! element-wise constant so this is exact up to differentiation round-off.
    nBad = 0
    do iel = 1,geomNew%nElem
      do j = 1,3
        if(MetricResidual(geomNew,iel,j) > metricTol) then
          nBad = nBad+1
          if(nBad <= 5) then
            print*,"FAIL: metric identity violated, element",iel,"direction",j, &
              MetricResidual(geomNew,iel,j)
          endif
        endif
      enddo
    enddo
    if(nBad /= 0) then
      print*,"FAIL: discrete metric identity violated on",nBad,"element/direction pairs"
      r = 1
    endif

    if(r == 0) then
      print*,"PASS: kept-element geometry is bit-for-bit reproducible; J > 0; metric identity holds"
    endif

    deallocate(srcIdx,dstIdx,genIdx)
    deallocate(oldLeaf)
    if(nGen > 0) call genGeom%Free()
    call geomNew%Free()
    call geomRef%Free()
    call geomOld%Free()
    call meshNew%Free()
    call meshOld%Free()
    call forest%Free()
    call interp%Free()
    call baseMesh%Free()

  endfunction geometry_3d_reuse

  logical function SameElement(a,ia,b,ib) result(same)
    !! Exact (bit-for-bit) comparison of every geometry quantity of element ia of a against
    !! element ib of b. dxds is included even though nothing outside the geometry module
    !! consumes it, so that a divergence there is caught rather than tolerated.
    use SELF_Constants
    use SELF_Geometry_3D
    implicit none
    type(SEMHex),intent(in) :: a
    integer,intent(in) :: ia
    type(SEMHex),intent(in) :: b
    integer,intent(in) :: ib

    same = all(a%x%interior(:,:,:,ia,:,:) == b%x%interior(:,:,:,ib,:,:)) .and. &
           all(a%x%boundary(:,:,:,ia,:,:) == b%x%boundary(:,:,:,ib,:,:)) .and. &
           all(a%dxds%interior(:,:,:,ia,:,:,:) == b%dxds%interior(:,:,:,ib,:,:,:)) .and. &
           all(a%dxds%boundary(:,:,:,ia,:,:,:) == b%dxds%boundary(:,:,:,ib,:,:,:)) .and. &
           all(a%dsdx%interior(:,:,:,ia,:,:,:) == b%dsdx%interior(:,:,:,ib,:,:,:)) .and. &
           all(a%dsdx%boundary(:,:,:,ia,:,:,:) == b%dsdx%boundary(:,:,:,ib,:,:,:)) .and. &
           all(a%nHat%boundary(:,:,:,ia,:,:) == b%nHat%boundary(:,:,:,ib,:,:)) .and. &
           all(a%nScale%boundary(:,:,:,ia,:) == b%nScale%boundary(:,:,:,ib,:)) .and. &
           all(a%J%interior(:,:,:,ia,:) == b%J%interior(:,:,:,ib,:)) .and. &
           all(a%J%boundary(:,:,:,ia,:) == b%J%boundary(:,:,:,ib,:))

  endfunction SameElement

  real(prec) function MetricResidual(g,iel,jdir) result(res)
    !! max_ijk | sum_i d(J dsdx(i,jdir))/dxi_i | on element iel, using the interpolant's
    !! derivative matrix. A geometry assembled from inconsistent pieces breaks this even when
    !! global volume and conservation checks still balance.
    use SELF_Constants
    use SELF_Geometry_3D
    implicit none
    type(SEMHex),intent(in) :: g
    integer,intent(in) :: iel
    integer,intent(in) :: jdir
    ! Local
    integer :: i,j,k,ii,Np
    real(prec) :: dxi,deta,dzeta,scale

    Np = g%J%interp%N+1
    res = 0.0_prec
    scale = 0.0_prec

    do k = 1,Np
      do j = 1,Np
        do i = 1,Np
          dxi = 0.0_prec
          deta = 0.0_prec
          dzeta = 0.0_prec
          do ii = 1,Np
            dxi = dxi+g%J%interp%dMatrix(ii,i)* &
                  g%J%interior(ii,j,k,iel,1)*g%dsdx%interior(ii,j,k,iel,1,1,jdir)
            deta = deta+g%J%interp%dMatrix(ii,j)* &
                   g%J%interior(i,ii,k,iel,1)*g%dsdx%interior(i,ii,k,iel,1,2,jdir)
            dzeta = dzeta+g%J%interp%dMatrix(ii,k)* &
                    g%J%interior(i,j,ii,iel,1)*g%dsdx%interior(i,j,ii,iel,1,3,jdir)
          enddo
          res = max(res,abs(dxi+deta+dzeta))
          scale = max(scale,abs(g%J%interior(i,j,k,iel,1)))
        enddo
      enddo
    enddo

    ! Report relative to the element's Jacobian magnitude so the tolerance is scale free.
    if(scale > 0.0_prec) res = res/scale

  endfunction MetricResidual

endprogram test
