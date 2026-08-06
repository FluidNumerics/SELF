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

! Guard test (expected to abort): ApplyTransferPlanRange must reject an element range outside
! 1..nNew. Registered with WILL_FAIL in CMake, so `stop 1` is the passing outcome.
program test
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_QuadTreeMesh_2D
  use SELF_TransferPlan_2D
  implicit none
  type(Mesh2D),target :: mesh
  type(QuadTreeMesh2D) :: forest
  type(TransferPlan2D) :: plan
  type(Lagrange),target :: interp
  integer :: bcids(1:4)
  integer :: flag(1:4)
  integer,allocatable :: oldLeaf(:)
  real(prec),allocatable :: uOld(:,:,:,:),uNew(:,:,:,:)

  bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]
  call mesh%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,bcids)
  call interp%Init(N=2,controlNodeType=GAUSS,M=2,targetNodeType=UNIFORM)
  call forest%Init(mesh)

  allocate(oldLeaf(1:forest%nLeaves))
  oldLeaf = forest%leaf(1:forest%nLeaves)
  flag = QUADTREE_REFINE
  call forest%AdaptFromFlags(flag)
  call BuildTransferPlan(forest,4,oldLeaf,plan)

  allocate(uOld(1:3,1:3,1:plan%nOld,1:1),uNew(1:3,1:3,1:plan%nNew,1:1))
  uOld = 0.0_prec
  call ApplyTransferPlanRange(plan,interp,1,uOld,1,plan%nNew+1,uNew) ! must stop 1

  print*,"ERROR: ApplyTransferPlanRange out-of-range guard did not trigger"
endprogram test
