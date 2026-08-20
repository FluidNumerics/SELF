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

! Guard test (expected to abort): DownloadOldWindow must reject a buffer that does not match the
! window it is asked to receive.
!
! The download is a single flat copy sized from the window, not an element-by-element loop, so a
! buffer with the wrong element extent would be filled past its end or short of it rather than
! producing a shape error the compiler could catch. The GPU override makes this concrete: there it
! is one hipMemcpy of perElem*nWinElem*nvar reals into whatever address it is handed.
program test
  use SELF_Constants
  use SELF_Lagrange
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Geometry_2D
  use iso_fortran_env,only:int64
  implicit none
  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)
  integer :: winFirst(1:1),winLast(1:1)
  integer(int64) :: nRecv,nSent,nRemote
  real(prec),allocatable :: uWin(:,:,:,:)

  bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]
  call mesh%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,bcids)
  call interp%Init(N=3,controlNodeType=GAUSS,M=3,targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)
  call modelobj%Init(mesh,geometry)

  ! A real migration first, so the marker matches and the SHAPE check is what fires.
  winFirst(1) = 1
  winLast(1) = mesh%nElem
  nRecv = 0
  nSent = 0
  nRemote = 0
  call modelobj%MigrateOldWindow(winFirst,winLast,1,mesh%nElem,nRecv,nSent,nRemote)

  ! One element short of the window: must stop 1
  allocate(uWin(1:interp%N+1,1:interp%N+1,1:mesh%nElem-1,1:modelobj%nvar))
  call modelobj%DownloadOldWindow(1,mesh%nElem,uWin)

  print*,"ERROR: DownloadOldWindow buffer-mismatch guard did not trigger"
endprogram test
