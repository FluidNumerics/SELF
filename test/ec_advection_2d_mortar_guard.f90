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

program ec_advection_2d_mortar_guard
!! Verifies that the entropy-conserving split-form models refuse nonconforming
!! (mortar) meshes : the plain L2 mortar projection would break their provable
!! entropy estimate, so ECDGModel2D initialization must abort on a mesh with
!! nMortars > 0. This test is registered with WILL_FAIL, so the expected nonzero
!! exit code from the guard is the success criterion.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_Mesh
  use SELF_ECAdvection2D

  implicit none
  type(ECAdvection2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  bcids(1:4) = [SELF_BC_RADIATION, & ! south
                SELF_BC_RADIATION, & ! east
                SELF_BC_RADIATION, & ! north
                SELF_BC_RADIATION] ! west
  call mesh%SimpleMortarMesh(0.1_prec,bcids)

  call interp%Init(N=3, &
                   controlNodeType=GAUSS_LOBATTO, &
                   M=7, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! The EC model guard must error stop here
  call modelobj%Init(mesh,geometry)

  ! Unreachable when the guard works; exit cleanly (test failure under WILL_FAIL)
  print*,"Error: EC model accepted a nonconforming (mortar) mesh."

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram ec_advection_2d_mortar_guard
