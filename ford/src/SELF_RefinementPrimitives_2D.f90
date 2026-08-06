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

module SELF_RefinementPrimitives_2D
!! Element-local primitives for h-refinement of 2-D quadrilateral spectral element meshes
!! (AMR Stage 2). These routines are deliberately free of the Mesh2D container, MPI, and GPU
!! dependencies so that the numerically and topologically delicate pieces - isoparametric child
!! geometry and refined-mesh connectivity - are pure, deterministic, and unit-testable in
!! isolation. The mesh-level driver (SELF_MeshRefinement_2D) assembles a Mesh2D_t from them.
!!
!! ## Child ordering and quadrant convention
!!
!! A parent quadrilateral is split into four children indexed 1..4, each covering one quadrant
!! of the parent reference square [-1,1]^2 and attached to the parent corner of the same index
!! (SELF's CGNS corner order):
!!
!!     eta
!!      ^   4 (NW)  |  3 (NE)         corner 1 = SW = (-1,-1)   corner 3 = NE = (+1,+1)
!!      |  ---------+---------        corner 2 = SE = (+1,-1)   corner 4 = NW = (-1,+1)
!!      |   1 (SW)  |  2 (SE)
!!      +-------------------> xi      child c covers xi in [ax-1, ax], eta in [ay-1, ay]
!!
!! with (ax,ay) = (0,0),(1,0),(1,1),(0,1) for c = 1,2,3,4 (0 = lower/left half, 1 = upper/right).
!!
!! Local side numbering follows SELF: 1=South, 2=East, 3=North, 4=West. A child's local side k
!! that lies on the parent boundary lies on parent side k (subdivision is orientation
!! preserving), which is what makes the connectivity below inherit the parent's neighbor/flip
!! data directly.

  use SELF_Constants
  use SELF_Lagrange

  implicit none

  ! Child-of-side table: childOfSide(t,s) is the child covering sub-position t (1,2) along the
  ! positive direction of parent side s. The child's local side that lies on parent side s is s.
  !   side 1 (South, +xi) : children 1(SW), 2(SE)
  !   side 2 (East,  +eta): children 2(SE), 3(NE)
  !   side 3 (North, +xi) : children 4(NW), 3(NE)
  !   side 4 (West,  +eta): children 1(SW), 4(NW)
  integer,parameter :: childOfSide(1:2,1:4) = reshape([1,2,2,3,4,3,1,4],[2,4])

contains

  subroutine SubdivideNodeCoords(geomInterp,nGeo,parentCoords,childCoords)
    !! Isoparametric subdivision of one element's geometry node coordinates into its four
    !! children. The parent geometry is the degree-nGeo Lagrange interpolant through
    !! parentCoords (sampled at geomInterp's control points, i.e. the mesh geometry nodes); each
    !! child node coordinate is that interpolant evaluated at the corresponding point of the
    !! parent reference square. Exact for any polynomial geometry of degree <= nGeo and for any
    !! control-node type, so straight-sided and curved (isoparametric) elements are both handled
    !! without approximation.
    !!
    !!   childCoords(d,i,j,c) = sum_{ii,jj} H(ii,i,ax) H(jj,j,ay) parentCoords(d,ii,jj)
    !!
    !! where H(ii,i,a) = l_ii( halfmap(s_i,a) ), halfmap(s,a) = (s + (2a-1))/2 maps the child
    !! node reference coordinate s_i into the parent half [a-1,a], and l_ii is the parent's
    !! Lagrange basis. d = 1,2 are the physical coordinate components.
    implicit none
    type(Lagrange),intent(in) :: geomInterp
    integer,intent(in) :: nGeo
    real(prec),intent(in) :: parentCoords(1:2,1:nGeo+1,1:nGeo+1)
    real(prec),intent(out) :: childCoords(1:2,1:nGeo+1,1:nGeo+1,1:4)
    ! Local
    integer :: i,ii,j,jj,c,a,ax,ay
    integer,parameter :: axc(1:4) = [0,1,1,0]
    integer,parameter :: ayc(1:4) = [0,0,1,1]
    real(prec) :: s
    real(prec) :: H(1:nGeo+1,1:nGeo+1,0:1)
    real(prec) :: rowsum(1:2),acc(1:2)

    ! Precompute the half-interval interpolation weights H(ii,i,a) for the lower/left (a=0) and
    ! upper/right (a=1) halves. H(:,i,a) are the parent basis values at the mapped position of
    ! child node i.
    do a = 0,1
      do i = 1,nGeo+1
        s = geomInterp%controlPoints(i)
        H(1:nGeo+1,i,a) = geomInterp%CalculateLagrangePolynomials(0.5_prec*(s+real(2*a-1,prec)))
      enddo
    enddo

    do c = 1,4
      ax = axc(c)
      ay = ayc(c)
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          acc(1:2) = 0.0_prec
          do jj = 1,nGeo+1
            rowsum(1:2) = 0.0_prec
            do ii = 1,nGeo+1
              rowsum(1:2) = rowsum(1:2)+H(ii,i,ax)*parentCoords(1:2,ii,jj)
            enddo
            acc(1:2) = acc(1:2)+H(jj,j,ay)*rowsum(1:2)
          enddo
          childCoords(1:2,i,j,c) = acc(1:2)
        enddo
      enddo
    enddo

  endsubroutine SubdivideNodeCoords

  subroutine RefineConnectivity(nElem,baseSideInfo,baseCorner,nodeOffset,nUniqueSidesBase, &
                                refSideInfo,refCorner,nUniqueSidesRef)
    !! Build the connectivity of a uniformly (2:1 in each direction) refined 2-D mesh from its
    !! base connectivity. Every base element p is split into four children (global id
    !! 4*(p-1)+c). The scheme is fully deterministic integer bookkeeping - no coordinate
    !! hashing, no flip recomputation, no MPI:
    !!
    !!   - Sibling faces interior to a parent are same-orientation (flip 0).
    !!   - A child face on a parent boundary inherits the parent face's neighbor element (the
    !!     appropriate child of the parent's neighbor), neighbor side, and flip. The sub-position
    !!     pairing across the face uses the parent flip: sub-position t maps to t on the neighbor
    !!     when flip = 0 and to 3-t when flip = 1.
    !!
    !! Because refinement is orientation preserving, this reproduces exactly the flips that a
    !! corner-node matching pass (RecalculateFlip) would compute, but without needing globally
    !! consistent node ids, which structured meshes do not guarantee.
    !!
    !! Refined corner node ids (returned in refCorner, SW/SE/NE/NW per child) reuse the base
    !! corner ids for retained corners, place one shared midpoint per base side
    !! (nodeOffset + base global side id) and one center per base element
    !! (nodeOffset + nUniqueSidesBase + p). These feed globalNodeIDs for mesh I/O; the solver
    !! reads flips from refSideInfo directly.
    implicit none
    integer,intent(in) :: nElem
    integer,intent(in) :: baseSideInfo(1:5,1:4,1:nElem)
    integer,intent(in) :: baseCorner(1:4,1:nElem)
    integer,intent(in) :: nodeOffset
    integer,intent(in) :: nUniqueSidesBase
    integer,intent(out) :: refSideInfo(1:5,1:4,1:4*nElem)
    integer,intent(out) :: refCorner(1:4,1:4*nElem)
    integer,intent(out) :: nUniqueSidesRef
    ! Local
    ! Interior (sibling) face pairs within a parent: (childA,sideA) <-> (childB,sideB), flip 0.
    integer,parameter :: nInner = 4
    integer,parameter :: innerA(1:2,1:nInner) = reshape([1,2,1,3,2,3,4,2],[2,nInner])
    integer,parameter :: innerB(1:2,1:nInner) = reshape([2,4,4,1,3,1,3,4],[2,nInner])
    integer :: p,s,t,tq,k,e,nbr,nbrSide
    integer :: q,s2,f,childP,childQ,eBase
    integer :: m(1:4),z,cc(1:4)

    refSideInfo = 0
    refCorner = 0

    do p = 1,nElem
      eBase = 4*(p-1)

      ! ---- Corner node ids for the four children ----
      cc(1:4) = baseCorner(1:4,p)
      do s = 1,4
        m(s) = nodeOffset+abs(baseSideInfo(2,s,p)) ! shared midpoint on base side s
      enddo
      z = nodeOffset+nUniqueSidesBase+p ! element center
      refCorner(1:4,eBase+1) = [cc(1),m(1),z,m(4)] ! child 1 (SW)
      refCorner(1:4,eBase+2) = [m(1),cc(2),m(2),z] ! child 2 (SE)
      refCorner(1:4,eBase+3) = [z,m(2),cc(3),m(3)] ! child 3 (NE)
      refCorner(1:4,eBase+4) = [m(4),z,m(3),cc(4)] ! child 4 (NW)

      ! ---- Interior sibling faces (flip 0) ----
      do k = 1,nInner
        childP = innerA(1,k); s = innerA(2,k)
        refSideInfo(3,s,eBase+childP) = eBase+innerB(1,k) ! neighbor element
        refSideInfo(4,s,eBase+childP) = 10*innerB(2,k) ! 10*neighbor side + flip(0)
        refSideInfo(5,s,eBase+childP) = 0 ! interior: no BC
        childP = innerB(1,k); s = innerB(2,k)
        refSideInfo(3,s,eBase+childP) = eBase+innerA(1,k)
        refSideInfo(4,s,eBase+childP) = 10*innerA(2,k)
        refSideInfo(5,s,eBase+childP) = 0
      enddo

      ! ---- Exterior faces on each parent side ----
      do s = 1,4
        q = baseSideInfo(3,s,p) ! neighbor base element (0 = physical boundary)
        s2 = baseSideInfo(4,s,p)/10 ! neighbor base side
        f = baseSideInfo(4,s,p)-10*s2 ! base flip (0 or 1)
        do t = 1,2
          childP = childOfSide(t,s) ! child on this side; its local side is s
          if(q == 0) then
            ! Physical boundary: inherit the parent side's BC, no neighbor.
            refSideInfo(3,s,eBase+childP) = 0
            refSideInfo(4,s,eBase+childP) = 0
            refSideInfo(5,s,eBase+childP) = baseSideInfo(5,s,p)
          else
            if(f == 0) then
              tq = t
            else
              tq = 3-t
            endif
            childQ = childOfSide(tq,s2)
            refSideInfo(3,s,eBase+childP) = 4*(q-1)+childQ
            refSideInfo(4,s,eBase+childP) = 10*s2+f
            refSideInfo(5,s,eBase+childP) = 0
          endif
        enddo
      enddo
    enddo

    ! ---- Assign a global side id to each unique refined side (count each shared side once) ----
    nUniqueSidesRef = 0
    do e = 1,4*nElem
      do s = 1,4
        if(refSideInfo(2,s,e) /= 0) cycle ! already numbered from its partner
        nbr = refSideInfo(3,s,e)
        nbrSide = refSideInfo(4,s,e)/10
        nUniqueSidesRef = nUniqueSidesRef+1
        refSideInfo(2,s,e) = nUniqueSidesRef
        if(nbr /= 0) refSideInfo(2,nbrSide,nbr) = nUniqueSidesRef
      enddo
    enddo

    ! Side type field is unused in SELF's 2-D path; leave refSideInfo(1,:,:) = 0.

  endsubroutine RefineConnectivity

endmodule SELF_RefinementPrimitives_2D
