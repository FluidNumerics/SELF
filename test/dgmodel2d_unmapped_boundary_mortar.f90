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

program DGModel2D_Unmapped_Boundary_Mortar
!! Guard for the unmapped-boundary scan added for issue #180 on a nonconforming mesh.
!!
!! Mortar sides carry sideInfo(3) = 0 exactly like a physical boundary, and sideInfo(5) = 0,
!! so a scan that keyed only on sideInfo(3) would report every mortar side as an unmapped
!! boundary. The scan excludes them via sideInfo(1), which holds the mortar index.
!!
!! Assertions:
!!  (a) the mesh really does carry mortars, so the guard is exercised;
!!  (b) with every domain edge tagged to a registered bcid, nothing is reported unmapped.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d

  implicit none

  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  real(prec),parameter :: dx = 0.1_prec

  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  ! Radiation is registered by AdditionalInit_LinearEuler2D_t
  bcids(1:4) = [SELF_BC_RADIATION, & ! south
                SELF_BC_RADIATION, & ! east
                SELF_BC_RADIATION, & ! north
                SELF_BC_RADIATION] ! west

  ! Three-element mesh with a single 2:1 mortar interface at x = 2*dx
  call mesh%SimpleMortarMesh(dx,bcids)

  ! (a) the guard is only meaningful on a mesh that has mortar sides
  if(mesh%nMortars <= 0) then
    print*,"Error: SimpleMortarMesh produced no mortars; the guard is not exercised."
    stop 1
  endif

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)

  ! (b) mortar sides must not be mistaken for unmapped physical boundaries
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: ",modelobj%nUnmappedBoundaries, &
      " boundary edges reported unmapped on a fully tagged mortar mesh; bcid ", &
      modelobj%unmappedBoundaryID
    stop 1
  endif

  print*,"dgmodel2d_unmapped_boundary_mortar : all assertions passed"

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram DGModel2D_Unmapped_Boundary_Mortar
