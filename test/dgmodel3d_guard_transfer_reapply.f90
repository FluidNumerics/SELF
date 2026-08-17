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

! Guard test (expected to abort): ApplyTransferPlan must reject a SECOND apply that is not
! preceded by its own StageSolutionForTransfer.
!
! The companion test dgmodel3d_guard_transfer_unstaged covers a fresh model, where any guard
! catches the mistake. This one covers the case that only a consuming guard catches: one full
! stage -> regrid -> apply cycle has already succeeded, so on a GPU build the device staging
! buffer is still allocated (it is a persistent capacity buffer, deliberately not released
! between epochs) and the previous epoch's field is still sitting in it. A guard that merely
! tested whether that buffer exists would pass the unpaired call through and silently transfer
! stale data. Registered with WILL_FAIL in CMake, so the stop 1 is the passing outcome.
program test
  use SELF_Constants
  use SELF_Lagrange
  use self_lineareuler3d
  use self_mesh_3d
  use SELF_Geometry_3D
  use SELF_OctreeMesh_3D
  use SELF_AdaptiveMesh_3D
  use SELF_TransferPlan_3D
  implicit none
  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: baseMesh,oldMesh,newMesh
  type(SEMHex),target :: oldGeom,newGeom
  type(OctreeMesh3D) :: forest
  type(TransferPlan3D) :: plan
  integer :: bcids(1:6)
  integer :: nOld
  integer,allocatable :: flag(:),oldLeaf(:)

  bcids(1:6) = SELF_BC_PRESCRIBED
  call baseMesh%StructuredMesh(2,2,2,1,1,1,0.5_prec,0.5_prec,0.5_prec,bcids)
  call interp%Init(N=3,controlNodeType=GAUSS,M=3,targetNodeType=UNIFORM)

  call forest%Init(baseMesh)
  call EmitMesh(forest,baseMesh,oldMesh)
  call oldGeom%Init(interp,oldMesh%nElem)
  call oldGeom%GenerateFromMesh(oldMesh)
  call modelobj%Init(oldMesh,oldGeom)

  nOld = forest%nLeaves
  allocate(oldLeaf(1:nOld))
  oldLeaf(1:nOld) = forest%leaf(1:nOld)

  ! One legitimate epoch: refine a root, transfer onto the emitted mesh.
  allocate(flag(1:nOld))
  flag = OCTREE_KEEP
  flag(1) = OCTREE_REFINE
  call forest%AdaptFromFlags(flag)
  call forest%Balance2to1()
  call BuildTransferPlan(forest,nOld,oldLeaf,plan)
  call EmitMesh(forest,baseMesh,newMesh)
  call newGeom%Init(interp,newMesh%nElem)
  call newGeom%GenerateFromMesh(newMesh)

  call modelobj%StageSolutionForTransfer()
  call modelobj%Regrid(newMesh,newGeom)
  call modelobj%ApplyTransferPlan(plan,interp,1,plan%nNew)

  ! The staged field has now been consumed. Applying the same plan again without re-staging
  ! must stop 1 rather than silently reusing it.
  call modelobj%ApplyTransferPlan(plan,interp,1,plan%nNew)

  print*,"ERROR: ApplyTransferPlan re-apply guard did not trigger"
endprogram test
