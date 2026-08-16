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

module SELF_SolutionTransfer_3D
!! Solution transfer between an element and its eight h-refinement children (AMR Stage 3).
!!
!! When a 3-D hexahedral element is refined into eight children (SELF child ordering = CGNS
!! corner order, child c covering the octant with half-interval indices
!! (octAxc(c),octAyc(c),octAzc(c)) of the parent reference cube; see
!! SELF_RefinementPrimitives_3D), the prognostic solution must move with the mesh:
!!
!!   - **Prolongation** (parent -> 8 children): sample the parent's degree-N nodal polynomial at
!!     the children's nodes. Each child occupies a reference half-interval in each direction, so
!!     this is exactly a triple tensor product of the 2:1 mortar restriction operator
!!     Lagrange%mortarR (built and identity-checked for the nonconforming interface machinery).
!!     It is an exact interpolation: the child fields reproduce the parent polynomial with no
!!     loss.
!!
!!   - **Restriction** (8 children -> parent): the L2 projection of the (generally discontinuous
!!     across children) fine solution back onto the parent's degree-N polynomial space, a triple
!!     tensor product of the mortar projection operator Lagrange%mortarP (the M-inner-product
!!     adjoint of mortarR). Each 1-D mortarP carries the 1/2 per-direction sub-interval
!!     Jacobian, so the triple product carries the 1/8 sub-cell Jacobian appropriate for
!!     projecting *solution* fields: a conservative L2 projection.
!!
!! These inherit the discrete identities established for the 1-D mortar operators
!! (test/mortarprojection_identity.f90):
!!
!!   * Consistency / reversibility:  RestrictFromChildren(ProlongToChildren(u)) = u exactly, from
!!     sum_k P_k R_k = I applied in each tensor direction. Refining then immediately coarsening
!!     does not perturb the solution.
!!   * Conservation:  sum_ijk w_i w_j w_k u_parent(i,j,k) =
!!     (1/8) sum_c sum_ijk w_i w_j w_k u_child_c(i,j,k), i.e. the reference-cell integral of the
!!     restricted parent equals the sum of the children's reference-cell integrals (each child
!!     is an eighth of the parent). Weighted by the geometry Jacobian this is conservation of
!!     the cell-integrated prognostic quantity.
!!
!! The routines are element-local and portable (host `do concurrent` over children/variables);
!! the AMR driver maps forest parent/child relationships onto the element index ranges it passes
!! in. Transfer runs between time steps, so it is not a per-step hot path; on GPU backends the
!! transferred field is re-uploaded as part of the Stage-6 device re-allocation.

  use SELF_Constants
  use SELF_Lagrange

  implicit none

  ! Child octant -> (x-half, y-half, z-half) with 0 = lower ([-1,0]), 1 = upper ([0,1]).
  ! Half index h maps to mortar sub-interval k = h+1 (sub-interval 1 = [-1,0], 2 = [0,1]).
  ! Identical to octAxc/octAyc/octAzc in SELF_RefinementPrimitives_3D.
  integer,parameter :: transferAxc(1:8) = [0,1,1,0,0,1,1,0]
  integer,parameter :: transferAyc(1:8) = [0,0,1,1,0,0,1,1]
  integer,parameter :: transferAzc(1:8) = [0,0,0,0,1,1,1,1]

contains

  subroutine ProlongToChildren(interp,nVar,uParent,uChildren)
    !! Prolong (interpolate) a parent element's nodal solution onto its eight children.
    !! uChildren(:,:,:,:,c) is the field on child c (SELF child ordering = CGNS corner order).
    !! Exact for the parent's degree-N polynomial representation.
    implicit none
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: nVar
    real(prec),intent(in) :: uParent(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nVar)
    real(prec),intent(out) :: uChildren(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nVar,1:8)
    ! Local
    integer :: c,v,Np

    Np = interp%N+1

    do concurrent(c=1:8,v=1:nVar)
      block
        integer :: i,j,k,ii,jj,kk,kx,ky,kz
        real(prec) :: acc
        real(prec) :: tmp1(1:Np,1:Np,1:Np) ! tmp1(childX, parentY, parentZ) after the x pass
        real(prec) :: tmp2(1:Np,1:Np,1:Np) ! tmp2(childX, childY, parentZ) after the y pass

        kx = transferAxc(c)+1
        ky = transferAyc(c)+1
        kz = transferAzc(c)+1

        ! x-direction: tmp1(i,jj,kk) = sum_ii mortarR(ii,i,kx) * uParent(ii,jj,kk)
        do k = 1,Np
          do j = 1,Np
            do i = 1,Np
              acc = 0.0_prec
              do ii = 1,Np
                acc = acc+interp%mortarR(ii,i,kx)*uParent(ii,j,k,v)
              enddo
              tmp1(i,j,k) = acc
            enddo
          enddo
        enddo
        ! y-direction: tmp2(i,j,kk) = sum_jj mortarR(jj,j,ky) * tmp1(i,jj,kk)
        do k = 1,Np
          do j = 1,Np
            do i = 1,Np
              acc = 0.0_prec
              do jj = 1,Np
                acc = acc+interp%mortarR(jj,j,ky)*tmp1(i,jj,k)
              enddo
              tmp2(i,j,k) = acc
            enddo
          enddo
        enddo
        ! z-direction: uChild(i,j,k) = sum_kk mortarR(kk,k,kz) * tmp2(i,j,kk)
        do k = 1,Np
          do j = 1,Np
            do i = 1,Np
              acc = 0.0_prec
              do kk = 1,Np
                acc = acc+interp%mortarR(kk,k,kz)*tmp2(i,j,kk)
              enddo
              uChildren(i,j,k,v,c) = acc
            enddo
          enddo
        enddo
      endblock
    enddo

  endsubroutine ProlongToChildren

  subroutine RestrictFromChildren(interp,nVar,uChildren,uParent)
    !! Restrict (L2-project) the solution on eight children back onto their parent element.
    !! Conservative, and the exact left inverse of ProlongToChildren.
    implicit none
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: nVar
    real(prec),intent(in) :: uChildren(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nVar,1:8)
    real(prec),intent(out) :: uParent(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nVar)
    ! Local
    integer :: v,Np

    Np = interp%N+1

    do concurrent(v=1:nVar)
      block
        integer :: c,i,j,k,ii,jj,kk,kx,ky,kz
        real(prec) :: acc
        real(prec) :: tmp1(1:Np,1:Np,1:Np) ! tmp1(parentX, childY, childZ) after the x projection
        real(prec) :: tmp2(1:Np,1:Np,1:Np) ! tmp2(parentX, parentY, childZ) after the y projection
        real(prec) :: up(1:Np,1:Np,1:Np)

        up(1:Np,1:Np,1:Np) = 0.0_prec
        do c = 1,8
          kx = transferAxc(c)+1
          ky = transferAyc(c)+1
          kz = transferAzc(c)+1
          ! x-direction: tmp1(i,jj,kk) = sum_ii mortarP(ii,i,kx) * uChild(ii,jj,kk,c)
          do kk = 1,Np
            do jj = 1,Np
              do i = 1,Np
                acc = 0.0_prec
                do ii = 1,Np
                  acc = acc+interp%mortarP(ii,i,kx)*uChildren(ii,jj,kk,v,c)
                enddo
                tmp1(i,jj,kk) = acc
              enddo
            enddo
          enddo
          ! y-direction: tmp2(i,j,kk) = sum_jj mortarP(jj,j,ky) * tmp1(i,jj,kk)
          do kk = 1,Np
            do j = 1,Np
              do i = 1,Np
                acc = 0.0_prec
                do jj = 1,Np
                  acc = acc+interp%mortarP(jj,j,ky)*tmp1(i,jj,kk)
                enddo
                tmp2(i,j,kk) = acc
              enddo
            enddo
          enddo
          ! z-direction and accumulate over the eight children
          do k = 1,Np
            do j = 1,Np
              do i = 1,Np
                acc = 0.0_prec
                do kk = 1,Np
                  acc = acc+interp%mortarP(kk,k,kz)*tmp2(i,j,kk)
                enddo
                up(i,j,k) = up(i,j,k)+acc
              enddo
            enddo
          enddo
        enddo
        uParent(1:Np,1:Np,1:Np,v) = up(1:Np,1:Np,1:Np)
      endblock
    enddo

  endsubroutine RestrictFromChildren

endmodule SELF_SolutionTransfer_3D
