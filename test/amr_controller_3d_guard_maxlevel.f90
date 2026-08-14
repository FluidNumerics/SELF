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

! Guard test (expected to abort): the AMR controller must reject a negative level cap.
! Registered with WILL_FAIL in CMake, so the stop 1 is the passing outcome.
program test
  use SELF_Constants
  use self_lineareuler3d
  use self_mesh_3d
  use SELF_Geometry_3D
  use SELF_AMRController_3D
  implicit none
  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  type(AMRController3D) :: controller
  integer :: bcids(1:6)

  bcids(1:6) = SELF_BC_PRESCRIBED
  call mesh%StructuredMesh(2,2,2,1,1,1,0.5_prec,0.5_prec,0.5_prec,bcids)
  call interp%Init(N=2,controlNodeType=GAUSS,M=2,targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)
  call modelobj%Init(mesh,geometry)

  call controller%Init(modelobj,refineThreshold=-3.0_prec,coarsenThreshold=-8.0_prec, &
                       ivar=4,maxLevel=-1,nHalo=1) ! must stop 1

  print*,"ERROR: AMRController3D maxLevel guard did not trigger"
endprogram test
