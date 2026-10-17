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

module SELF_RefinementPrimitives_3D
!! Element-local primitives for h-refinement of 3-D hexahedral spectral element meshes.
!! These routines are deliberately free of the Mesh3D container, MPI, and GPU dependencies
!! so that the numerically and topologically delicate pieces - isoparametric child geometry
!! and refined-mesh connectivity - are pure, deterministic, and unit-testable in isolation.
!! The mesh-level driver (SELF_MeshRefinement_3D) assembles a Mesh3D_t from them. Direct
!! 3-D analogue of SELF_RefinementPrimitives_2D.
!!
!! ## Child ordering and octant convention
!!
!! A parent hexahedron is split into eight children indexed 1..8, each covering one octant
!! of the parent reference cube [-1,1]^3 and attached to the parent corner of the same
!! index (SELF's CGNS corner order: counter-clockwise looking in -xi3, bottom then top):
!!
!!   child c covers xi1 in [ax-1,ax], xi2 in [ay-1,ay], xi3 in [az-1,az], with
!!   (ax,ay,az) = (0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)
!!   for c = 1..8 (0 = lower half, 1 = upper half per direction).
!!
!! Local face numbering follows SELF: 1=Bottom, 2=South, 3=East, 4=North, 5=West, 6=Top.
!! A child's local face k that lies on the parent boundary lies on parent face k
!! (subdivision is orientation preserving), which is what makes the connectivity below
!! inherit the parent's neighbor/flip data directly.
!!
!! ## Face quadrant convention
!!
!! A parent face is split into four sub-faces ("quadrants") indexed in the face's trace
!! coordinates (SELF's BoundaryInterp convention: faces 1,6 -> (xi1,xi2); faces 2,4 ->
!! (xi1,xi3); faces 3,5 -> (xi2,xi3)) as q = kx + 2*(ky-1), where kx, ky in {1,2} are the
!! half-interval indices along the face's first and second coordinate.

  use SELF_Constants
  use SELF_Lagrange

  implicit none

  ! Half-interval index of each child along the three reference directions.
  integer,parameter :: octAxc(1:8) = [0,1,1,0,0,1,1,0]
  integer,parameter :: octAyc(1:8) = [0,0,1,1,0,0,1,1]
  integer,parameter :: octAzc(1:8) = [0,0,0,0,1,1,1,1]

  ! Child-of-face table: childOfFace(q,s) is the child whose face covers quadrant q of
  ! parent face s. The child's local face that lies on parent face s is s.
  integer,parameter :: childOfFace(1:4,1:6) = reshape([ &
                                                      1,2,4,3, & ! face 1 (Bottom) : (kx,ky) = (xi1,xi2)
                                                      1,2,5,6, & ! face 2 (South)  : (xi1,xi3)
                                                      2,3,6,7, & ! face 3 (East)   : (xi2,xi3)
                                                      4,3,8,7, & ! face 4 (North)  : (xi1,xi3)
                                                      1,4,5,8, & ! face 5 (West)   : (xi2,xi3)
                                                      5,6,8,7],[4,6]) ! face 6 (Top) : (xi1,xi2)

  ! Quadrant permutation under the eight SELF face flips: faceQuadPerm(q,f) is the
  ! donor-face quadrant that coincides with receiver-face quadrant q when the flip f maps
  ! receiver-face indices to donor-face indices (the sideInfo(4) convention; see
  ! MortarFaceMap in SELF_Mesh_3D_t). Obtained by applying the flip index maps to the
  ! half-interval indices.
  integer,parameter :: faceQuadPerm(1:4,0:7) = reshape([ &
                                                       1,2,3,4, & ! flip 0 : i2 = i1,   j2 = j1
                                                       2,1,4,3, & ! flip 1 : i2 = N-i1, j2 = j1
                                                       4,3,2,1, & ! flip 2 : i2 = N-i1, j2 = N-j1
                                                       3,4,1,2, & ! flip 3 : i2 = i1,   j2 = N-j1
                                                       1,3,2,4, & ! flip 4 : i2 = j1,   j2 = i1
                                                       2,4,1,3, & ! flip 5 : i2 = N-j1, j2 = i1
                                                       4,2,3,1, & ! flip 6 : i2 = N-j1, j2 = N-i1
                                                       3,1,4,2],[4,8]) ! flip 7 : i2 = j1, j2 = N-i1

contains

  subroutine SubdivideNodeCoords(geomInterp,nGeo,parentCoords,childCoords)
    !! Isoparametric subdivision of one element's geometry node coordinates into its eight
    !! children. The parent geometry is the degree-nGeo Lagrange interpolant through
    !! parentCoords (sampled at geomInterp's control points, i.e. the mesh geometry nodes);
    !! each child node coordinate is that interpolant evaluated at the corresponding point
    !! of the parent reference cube. Exact for any polynomial geometry of degree <= nGeo
    !! and for any control-node type, so straight-sided and curved (isoparametric) elements
    !! are both handled without approximation.
    !!
    !!   childCoords(d,i,j,k,c) =
    !!     sum_{ii,jj,kk} H(ii,i,ax) H(jj,j,ay) H(kk,k,az) parentCoords(d,ii,jj,kk)
    !!
    !! where H(ii,i,a) = l_ii( halfmap(s_i,a) ), halfmap(s,a) = (s + (2a-1))/2 maps the
    !! child node reference coordinate s_i into the parent half [a-1,a], and l_ii is the
    !! parent's Lagrange basis. d = 1,2,3 are the physical coordinate components.
    implicit none
    type(Lagrange),intent(in) :: geomInterp
    integer,intent(in) :: nGeo
    real(prec),intent(in) :: parentCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1)
    real(prec),intent(out) :: childCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:8)
    ! Local
    integer :: i,ii,j,jj,k,kk,c,a,ax,ay,az
    real(prec) :: s
    real(prec) :: H(1:nGeo+1,1:nGeo+1,0:1)
    real(prec) :: rowsum(1:3),colsum(1:3),acc(1:3)

    ! Precompute the half-interval interpolation weights H(ii,i,a) for the lower (a=0) and
    ! upper (a=1) halves. H(:,i,a) are the parent basis values at the mapped position of
    ! child node i.
    do a = 0,1
      do i = 1,nGeo+1
        s = geomInterp%controlPoints(i)
        H(1:nGeo+1,i,a) = geomInterp%CalculateLagrangePolynomials(0.5_prec*(s+real(2*a-1,prec)))
      enddo
    enddo

    do c = 1,8
      ax = octAxc(c)
      ay = octAyc(c)
      az = octAzc(c)
      do k = 1,nGeo+1
        do j = 1,nGeo+1
          do i = 1,nGeo+1
            acc(1:3) = 0.0_prec
            do kk = 1,nGeo+1
              colsum(1:3) = 0.0_prec
              do jj = 1,nGeo+1
                rowsum(1:3) = 0.0_prec
                do ii = 1,nGeo+1
                  rowsum(1:3) = rowsum(1:3)+H(ii,i,ax)*parentCoords(1:3,ii,jj,kk)
                enddo
                colsum(1:3) = colsum(1:3)+H(jj,j,ay)*rowsum(1:3)
              enddo
              acc(1:3) = acc(1:3)+H(kk,k,az)*colsum(1:3)
            enddo
            childCoords(1:3,i,j,k,c) = acc(1:3)
          enddo
        enddo
      enddo
    enddo

  endsubroutine SubdivideNodeCoords

  subroutine RefineConnectivity(nElem,baseSideInfo,refSideInfo,nUniqueSidesRef)
    !! Build the connectivity of a uniformly (2:1 in each direction) refined 3-D mesh from
    !! its base connectivity. Every base element p is split into eight children (global id
    !! 8*(p-1)+c). The scheme is fully deterministic integer bookkeeping - no coordinate
    !! hashing, no flip recomputation, no MPI:
    !!
    !!   - Sibling faces interior to a parent are same-orientation (flip 0).
    !!   - A child face on a parent boundary inherits the parent face's neighbor element
    !!     (the appropriate child of the parent's neighbor), neighbor face, and flip. The
    !!     quadrant pairing across the face applies the parent flip's quadrant permutation
    !!     (faceQuadPerm).
    !!
    !! Because refinement is orientation preserving, this reproduces exactly the flips that
    !! a corner-node matching pass (RecalculateFlip) would compute, but without needing
    !! globally consistent node ids, which structured meshes do not guarantee. Refined
    !! globalNodeIDs are not produced: the solver reads flips from refSideInfo directly and
    !! no 3-D consumer requires them (the adaptive-mesh emitter sets them to zero as well).
    implicit none
    integer,intent(in) :: nElem
    integer,intent(in) :: baseSideInfo(1:5,1:6,1:nElem)
    integer,intent(out) :: refSideInfo(1:5,1:6,1:8*nElem)
    integer,intent(out) :: nUniqueSidesRef
    ! Local
    ! Interior (sibling) face pairs within a parent: (childA,faceA) <-> (childB,faceB),
    ! flip 0. Four pairs per direction : East/West, North/South, Top/Bottom.
    integer,parameter :: nInner = 12
    integer,parameter :: innerA(1:2,1:nInner) = reshape([ &
                                                        1,3,4,3,5,3,8,3, & ! East faces
                                                        1,4,2,4,5,4,6,4, & ! North faces
                                                        1,6,2,6,3,6,4,6],[2,nInner]) ! Top faces
    integer,parameter :: innerB(1:2,1:nInner) = reshape([ &
                                                        2,5,3,5,6,5,7,5, & ! West faces
                                                        4,2,3,2,8,2,7,2, & ! South faces
                                                        5,1,6,1,7,1,8,1],[2,nInner]) ! Bottom faces
    integer :: p,s,t,tq,k,e,nbr,nbrSide
    integer :: q,s2,f,childP,childQ,eBase

    refSideInfo = 0

    do p = 1,nElem
      eBase = 8*(p-1)

      ! ---- Interior sibling faces (flip 0) ----
      do k = 1,nInner
        childP = innerA(1,k)
        s = innerA(2,k)
        refSideInfo(3,s,eBase+childP) = eBase+innerB(1,k) ! neighbor element
        refSideInfo(4,s,eBase+childP) = 10*innerB(2,k) ! 10*neighbor face + flip(0)
        refSideInfo(5,s,eBase+childP) = 0 ! interior: no BC
        childP = innerB(1,k)
        s = innerB(2,k)
        refSideInfo(3,s,eBase+childP) = eBase+innerA(1,k)
        refSideInfo(4,s,eBase+childP) = 10*innerA(2,k)
        refSideInfo(5,s,eBase+childP) = 0
      enddo

      ! ---- Exterior faces on each parent face ----
      do s = 1,6
        q = baseSideInfo(3,s,p) ! neighbor base element (0 = physical boundary)
        s2 = baseSideInfo(4,s,p)/10 ! neighbor base face
        f = baseSideInfo(4,s,p)-10*s2 ! base flip (0..7)
        do t = 1,4
          childP = childOfFace(t,s) ! child on this face; its local face is s
          if(q == 0) then
            ! Physical boundary: inherit the parent face's BC, no neighbor.
            refSideInfo(3,s,eBase+childP) = 0
            refSideInfo(4,s,eBase+childP) = 0
            refSideInfo(5,s,eBase+childP) = baseSideInfo(5,s,p)
          else
            tq = faceQuadPerm(t,f)
            childQ = childOfFace(tq,s2)
            refSideInfo(3,s,eBase+childP) = 8*(q-1)+childQ
            refSideInfo(4,s,eBase+childP) = 10*s2+f
            refSideInfo(5,s,eBase+childP) = 0
          endif
        enddo
      enddo
    enddo

    ! ---- Assign a global side id to each unique refined face (count each shared face once) ----
    nUniqueSidesRef = 0
    do e = 1,8*nElem
      do s = 1,6
        if(refSideInfo(2,s,e) /= 0) cycle ! already numbered from its partner
        nbr = refSideInfo(3,s,e)
        nbrSide = refSideInfo(4,s,e)/10
        nUniqueSidesRef = nUniqueSidesRef+1
        refSideInfo(2,s,e) = nUniqueSidesRef
        if(nbr /= 0) refSideInfo(2,nbrSide,nbr) = nUniqueSidesRef
      enddo
    enddo

    ! Side type field is only used to mark mortar faces; uniform refinement of a
    ! conforming mesh is conforming, so leave refSideInfo(1,:,:) = 0.

  endsubroutine RefineConnectivity

endmodule SELF_RefinementPrimitives_3D
