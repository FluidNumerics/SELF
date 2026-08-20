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

! Guard test (expected to abort): ApplyTransferPlan must reject a staged local field and a migrated
! window presented together.
!
! The two are alternative sources for one apply - StageSolutionForTransfer preserves the whole local
! field for a single-rank transfer, MigrateOldWindow assembles the window a multi-rank transfer reads
! - and the caller picks exactly one. Having both set means the caller has done two contradictory
! things, so it is rejected rather than silently resolved by branch order, which would make whichever
! branch happened to come first the winner and leave the other source stale and unnoticed.
!
! The GPU override enforces the same condition; this pins the portable half, and the pair is the
! reason the guard exists at all, since the controller never reaches this state.
program test
  use SELF_Constants
  use SELF_Lagrange
  use self_lineareuler3d
  use self_mesh_3d
  use SELF_Geometry_3D
  use SELF_TransferPlan_3D
  use iso_fortran_env,only:int64
  implicit none
  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  type(TransferPlan3D) :: plan ! never built; the source-conflict guard fires first
  integer :: bcids(1:6)
  integer :: winFirst(1:1),winLast(1:1)
  integer(int64) :: nRecv,nSent,nRemote

  bcids(1:6) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED, &
                SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]
  call mesh%StructuredMesh(2,2,2,1,1,1,0.5_prec,0.5_prec,0.5_prec,bcids)
  call interp%Init(N=3,controlNodeType=GAUSS,M=3,targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)
  call modelobj%Init(mesh,geometry)

  ! A window covering the whole local field: at one rank the migration posts no messages and its
  ! whole job is the local copy, which is all that is needed to set the window marker.
  winFirst(1) = 1
  winLast(1) = mesh%nElem
  nRecv = 0
  nSent = 0
  nRemote = 0
  call modelobj%MigrateOldWindow(winFirst,winLast,1,mesh%nElem,nRecv,nSent,nRemote)

  ! ... and now the other source as well: must stop 1
  call modelobj%StageSolutionForTransfer()
  call modelobj%ApplyTransferPlan(plan,interp,1,mesh%nElem)

  print*,"ERROR: ApplyTransferPlan window-and-stage conflict guard did not trigger"
endprogram test
