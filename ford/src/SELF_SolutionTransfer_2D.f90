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

module SELF_SolutionTransfer_2D
!! Solution transfer between an element and its four h-refinement children (AMR Stage 3).
!!
!! When a 2-D quadrilateral element is refined into four children (SELF child ordering
!! 1=SW, 2=SE, 3=NE, 4=NW, each covering one quadrant of the parent reference square), the
!! prognostic solution must move with the mesh:
!!
!!   - **Prolongation** (parent -> 4 children): sample the parent's degree-N nodal polynomial at
!!     the children's nodes. Each child occupies a reference half-interval in each direction, so
!!     this is exactly a tensor product of the 2:1 mortar restriction operator Lagrange%mortarR
!!     (built and identity-checked for the nonconforming interface machinery). It is an exact
!!     interpolation: the child fields reproduce the parent polynomial with no loss.
!!
!!   - **Restriction** (4 children -> parent): the L2 projection of the (generally discontinuous
!!     across children) fine solution back onto the parent's degree-N polynomial space, a tensor
!!     product of the mortar projection operator Lagrange%mortarP (the M-inner-product adjoint of
!!     mortarR). mortarP already carries the 1/2 per-direction sub-edge Jacobian appropriate for
!!     projecting *solution* traces.
!!
!! These inherit the discrete identities established for the 1-D mortar operators
!! (test/mortarprojection_identity.f90):
!!
!!   * Consistency / reversibility:  RestrictFromChildren(ProlongToChildren(u)) = u exactly, from
!!     sum_k P_k R_k = I applied in each tensor direction. Refining then immediately coarsening
!!     does not perturb the solution.
!!   * Conservation:  sum_ij w_i w_j u_parent(i,j) = (1/4) sum_c sum_ij w_i w_j u_child_c(i,j),
!!     i.e. the reference-cell integral of the restricted parent equals the sum of the children's
!!     reference-cell integrals (each child is a quarter of the parent). Weighted by the geometry
!!     Jacobian this is conservation of the cell-integrated prognostic quantity.
!!
!! The routines are element-local and portable (host `do concurrent` over children/variables);
!! the AMR driver maps forest parent/child relationships onto the element index ranges it passes
!! in. Transfer runs between time steps, so it is not a per-step hot path; on GPU backends the
!! transferred field is re-uploaded as part of the Stage-6 device re-allocation.

  use SELF_Constants
  use SELF_Lagrange

  implicit none

  ! Child quadrant -> (x-half, y-half) with 0 = lower/left ([-1,0]), 1 = upper/right ([0,1]).
  ! Half index h maps to mortar sub-edge k = h+1 (sub-edge 1 = [-1,0], sub-edge 2 = [0,1]).
  integer,parameter :: transferAxc(1:4) = [0,1,1,0]
  integer,parameter :: transferAyc(1:4) = [0,0,1,1]

contains

  subroutine ProlongToChildren(interp,nVar,uParent,uChildren)
    !! Prolong (interpolate) a parent element's nodal solution onto its four children.
    !! uChildren(:,:,:,c) is the field on child c (ordering 1=SW,2=SE,3=NE,4=NW). Exact for the
    !! parent's degree-N polynomial representation.
    implicit none
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: nVar
    real(prec),intent(in) :: uParent(1:interp%N+1,1:interp%N+1,1:nVar)
    real(prec),intent(out) :: uChildren(1:interp%N+1,1:interp%N+1,1:nVar,1:4)
    ! Local
    integer :: c,v,Np

    Np = interp%N+1

    do concurrent(c=1:4,v=1:nVar)
      block
        integer :: i,j,ii,jj,kx,ky
        real(prec) :: acc
        real(prec) :: tmp(1:Np,1:Np) ! tmp(childX, parentY) after the x-direction pass

        kx = transferAxc(c)+1
        ky = transferAyc(c)+1

        ! x-direction: tmp(i,jj) = sum_ii mortarR(ii,i,kx) * uParent(ii,jj)
        do j = 1,Np
          do i = 1,Np
            acc = 0.0_prec
            do ii = 1,Np
              acc = acc+interp%mortarR(ii,i,kx)*uParent(ii,j,v)
            enddo
            tmp(i,j) = acc
          enddo
        enddo
        ! y-direction: uChild(i,j) = sum_jj mortarR(jj,j,ky) * tmp(i,jj)
        do j = 1,Np
          do i = 1,Np
            acc = 0.0_prec
            do jj = 1,Np
              acc = acc+interp%mortarR(jj,j,ky)*tmp(i,jj)
            enddo
            uChildren(i,j,v,c) = acc
          enddo
        enddo
      endblock
    enddo

  endsubroutine ProlongToChildren

  subroutine RestrictFromChildren(interp,nVar,uChildren,uParent)
    !! Restrict (L2-project) the solution on four children back onto their parent element.
    !! Conservative, and the exact left inverse of ProlongToChildren.
    implicit none
    type(Lagrange),intent(in) :: interp
    integer,intent(in) :: nVar
    real(prec),intent(in) :: uChildren(1:interp%N+1,1:interp%N+1,1:nVar,1:4)
    real(prec),intent(out) :: uParent(1:interp%N+1,1:interp%N+1,1:nVar)
    ! Local
    integer :: v,Np

    Np = interp%N+1

    do concurrent(v=1:nVar)
      block
        integer :: c,i,j,ii,jj,kx,ky
        real(prec) :: acc
        real(prec) :: tmp(1:Np,1:Np) ! tmp(parentX, childY) after the x-direction projection
        real(prec) :: up(1:Np,1:Np)

        up(1:Np,1:Np) = 0.0_prec
        do c = 1,4
          kx = transferAxc(c)+1
          ky = transferAyc(c)+1
          ! x-direction: tmp(i,jj) = sum_ii mortarP(ii,i,kx) * uChild(ii,jj,c)
          do jj = 1,Np
            do i = 1,Np
              acc = 0.0_prec
              do ii = 1,Np
                acc = acc+interp%mortarP(ii,i,kx)*uChildren(ii,jj,v,c)
              enddo
              tmp(i,jj) = acc
            enddo
          enddo
          ! y-direction and accumulate over the four children
          do j = 1,Np
            do i = 1,Np
              acc = 0.0_prec
              do jj = 1,Np
                acc = acc+interp%mortarP(jj,j,ky)*tmp(i,jj)
              enddo
              up(i,j) = up(i,j)+acc
            enddo
          enddo
        enddo
        uParent(1:Np,1:Np,v) = up(1:Np,1:Np)
      endblock
    enddo

  endsubroutine RestrictFromChildren

endmodule SELF_SolutionTransfer_2D
