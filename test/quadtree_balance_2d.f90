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

  implicit none
  integer :: exit_code

  exit_code = quadtree_balance_2d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function quadtree_balance_2d() result(r)
    !! Validates the AMR Stage 2b/4a quad-forest on a real structured base mesh:
    !!
    !!   1. Uniform refinement of every leaf stays conforming (MaxLevelJump = 0) and quadruples
    !!      the leaf count.
    !!   2. Adaptive refinement that drives a two-level difference across element faces is reduced
    !!      to a single level by Balance2to1 (the coarse face-neighbours are refined, rippling to
    !!      a fixed point), i.e. MaxLevelJump goes 2 -> 1.
    !!
    !! Exercises Init/AdaptFromFlags/FaceNeighbor/Balance2to1/MaxLevelJump against the real Mesh2D
    !! connectivity produced by StructuredMesh.
    use SELF_Constants
    use SELF_Mesh_2D
    use SELF_QuadTreeMesh_2D

    implicit none

    type(Mesh2D),target :: mesh
    type(QuadTreeMesh2D) :: forest
    integer :: bcids(1:4)
    integer :: nBefore
    integer,allocatable :: flag(:)

    r = 0
    bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    ! 2 x 2 element structured base mesh (4 corner elements, each sharing two interior faces).
    call mesh%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,bcids)

    ! ---- 1. Uniform refinement is conforming ----
    call forest%Init(mesh)
    if(forest%nLeaves /= mesh%nElem) then
      print*,"FAIL: initial leaf count",forest%nLeaves
      r = 1
    endif
    allocate(flag(1:forest%nLeaves))
    flag = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    if(forest%nLeaves /= 4*mesh%nElem) then
      print*,"FAIL: uniform-refine leaf count",forest%nLeaves
      r = 1
    endif
    if(forest%MaxLevelJump() /= 0) then
      print*,"FAIL: uniform refinement not conforming, jump =",forest%MaxLevelJump()
      r = 1
    endif

    ! Calling RefineNode on an already-refined node is a no-op (with a warning).
    call forest%RefineNode(1)
    if(forest%nLeaves /= 4*mesh%nElem) then
      print*,"FAIL: RefineNode on a non-leaf changed the leaf count"
      r = 1
    endif

    ! Coarsening: flag every leaf COARSEN; each family of four leaf siblings merges back, so the
    ! forest returns to the base element count.
    allocate(flag(1:forest%nLeaves))
    flag = QUADTREE_COARSEN
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    if(forest%nLeaves /= mesh%nElem .or. forest%MaxLevel() /= 0) then
      print*,"FAIL: coarsening did not return to the base mesh, nLeaves =",forest%nLeaves
      r = 1
    endif
    call forest%Free()

    ! ---- 2. Adaptive refinement + 2:1 balance ----
    call forest%Init(mesh)
    ! Refine the first root element ...
    allocate(flag(1:forest%nLeaves))
    flag = QUADTREE_KEEP
    flag(1) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    ! ... then refine all four of its children, making that element's region two levels finer
    ! than its still-level-0 face-neighbours.
    allocate(flag(1:forest%nLeaves))
    flag = QUADTREE_KEEP
    flag(1:4) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)

    if(forest%MaxLevelJump() /= 2) then
      print*,"FAIL: expected a two-level jump before balancing, got",forest%MaxLevelJump()
      r = 1
    endif

    nBefore = forest%nLeaves
    call forest%Balance2to1()
    if(forest%MaxLevelJump() /= 1) then
      print*,"FAIL: 2:1 balance did not reduce the jump to 1, got",forest%MaxLevelJump()
      r = 1
    endif
    if(forest%nLeaves <= nBefore) then
      print*,"FAIL: balancing should have refined coarse neighbours"
      r = 1
    endif
    print*,"leaves before/after balance :",nBefore,forest%nLeaves, &
      "  max level :",forest%MaxLevel()

    ! Final consistency: no leaf may see a leaf neighbour more than one level away.
    if(forest%MaxLevelJump() > 1) r = 1

    if(r == 0) print*,"QUADTREE BALANCE CHECKS PASSED"

    call forest%Free()
    call mesh%Free()

  endfunction quadtree_balance_2d

endprogram test
