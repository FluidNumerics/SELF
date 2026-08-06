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

! Guard test (expected to abort): Regrid must reject a model that was never initialized.

! Guard test (expected to abort): ApplyTransferPlan must reject an unstaged solution.
!
! The AMR transfer is a two-step protocol - StageSolutionForTransfer preserves the solution, then
! Regrid is free to release the storage it lived in, then ApplyTransferPlan consumes the staged
! copy. Calling the second step without the first would read from an unallocated buffer, so it is
! rejected rather than left to fail as a memory error somewhere less obvious.
program test
  use SELF_Constants
  use SELF_Lagrange
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Geometry_2D
  use SELF_TransferPlan_2D
  implicit none
  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  type(TransferPlan2D) :: plan ! never built; the staging guard fires first
  integer :: bcids(1:4)

  bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]
  call mesh%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,bcids)
  call interp%Init(N=3,controlNodeType=GAUSS,M=3,targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)
  call modelobj%Init(mesh,geometry)

  ! No StageSolutionForTransfer: must stop 1
  call modelobj%ApplyTransferPlan(plan,interp,1,mesh%nElem)

  print*,"ERROR: ApplyTransferPlan unstaged-solution guard did not trigger"
endprogram test
