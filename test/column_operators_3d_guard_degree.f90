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
!! ColumnOperators3D%Init must reject an interpolant whose degree differs from
!! the one the geometry was generated with.
!! This is a WILL_FAIL test: it passes when the program aborts.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh
  use SELF_Mesh_2D
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_Scalar_2D
  use SELF_Scalar_3D
  use SELF_ColumnOperators_3D

  implicit none

  integer,parameter :: nz = 2
  type(Lagrange),target :: interp
  type(Lagrange),target :: interp4
  type(Mesh2D),target :: surf
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  type(SEMHex),target :: badGeometry
  type(ColumnOperators3D) :: cols
  type(Scalar3D) :: q
  type(Scalar3D) :: r
  type(Scalar2D) :: qbar
  real(prec) :: h(1:2,1:2,1:4)
  real(prec) :: sigma(0:nz)
  real(prec) :: zLevels(1:2,1:2,0:nz,1:4)

  ! A 2x2 surface extruded in two flat layers: 4 columns, 8 elements
  call surf%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec, &
                           [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                            SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW])
  h = 1.0_prec
  sigma = [-1.0_prec,-0.5_prec,0.0_prec]
  call ExtrudedSigmaLevels(h,sigma,zLevels)
  call mesh%ExtrudedMesh(surf,zLevels,SELF_BC_NONORMALFLOW,SELF_BC_PRESCRIBED)
  call interp%Init(N=3,controlNodeType=GAUSS_LOBATTO,M=3,targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call interp4%Init(N=4,controlNodeType=GAUSS_LOBATTO,M=4,targetNodeType=UNIFORM)
  call cols%Init(interp4,mesh,geometry)

  ! Unreachable when the guard works; exiting cleanly here would make
  ! the WILL_FAIL test fail.
  print*,"FAIL: the guard in ColumnOperators3D did not stop the program"

endprogram test
