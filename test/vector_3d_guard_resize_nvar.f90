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

! Guard test (expected to abort; registered with WILL_FAIL): Resize must reject a change of nVar.
!
! Resize exists to rebind a live object to a new ELEMENT count while reusing its storage. nVar is a
! different matter: it sizes the metadata and equation-parser arrays, which Resize deliberately
! preserves precisely because they do not depend on nElem. A caller asking for a different nVar is
! therefore not performing a regrid and must go through Free + Init, so the request is rejected
! rather than silently producing an object whose metadata disagrees with its data.
program test
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Vector_3D
  implicit none
  type(Lagrange),target :: interp
  type(Vector3D) :: f

  call interp%Init(N=3,controlNodeType=GAUSS,M=3,targetNodeType=UNIFORM)
  call f%Init(interp,2,8) ! nVar = 2

  call f%Resize(interp,3,8) ! must stop 1: nVar changed from 2 to 3

  print*,"ERROR: Vector3D Resize nVar guard did not trigger"
endprogram test
