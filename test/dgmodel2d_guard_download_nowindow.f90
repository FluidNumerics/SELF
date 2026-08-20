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

! Guard test (expected to abort): DownloadOldWindow must reject a request that no migrated window
! backs.
!
! DownloadOldWindow exists for the SELF_AMR_MIGRATE_VERIFY diagnostic, which compares the migrated
! window against an allgathered reference bit for bit. That comparison is only meaningful if the
! buffer it reads was actually filled by the migration it claims to describe, so asking for a window
! that was never migrated - or one whose bounds do not match what was - is rejected rather than
! answered with whatever the buffer last held. A stale window would make the diagnostic report
! agreement it never established, which is worse than no diagnostic at all.
program test
  use SELF_Constants
  use SELF_Lagrange
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Geometry_2D
  implicit none
  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)
  real(prec),allocatable :: uWin(:,:,:,:)

  bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]
  call mesh%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,bcids)
  call interp%Init(N=3,controlNodeType=GAUSS,M=3,targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)
  call modelobj%Init(mesh,geometry)

  ! Correctly shaped for the window being asked for, so the shape check passes and the marker check
  ! is what fires: no MigrateOldWindow has been called, so nothing backs this request.
  allocate(uWin(1:interp%N+1,1:interp%N+1,1:mesh%nElem,1:modelobj%nvar))
  call modelobj%DownloadOldWindow(1,mesh%nElem,uWin)

  print*,"ERROR: DownloadOldWindow unmigrated-window guard did not trigger"
endprogram test
