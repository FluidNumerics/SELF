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

module SELF_Mesh

  use SELF_Constants
  use SELF_DomainDecomposition
  use iso_c_binding

  implicit none

  type,abstract :: SEMMesh
    integer :: nGeo
    integer :: nElem
    integer :: nGlobalElem
    integer :: nNodes
    integer :: nSides
    integer :: nCornerNodes
    integer :: nUniqueNodes
    integer :: nUniqueSides
    integer :: nBCs
    integer :: quadrature
    type(DomainDecomposition) :: decomp

  contains
    procedure(SELF_FreeMesh),deferred :: Free

  endtype SEMMesh

  interface
    subroutine SELF_FreeMesh(this)
      import SEMMesh
      implicit none
      class(SEMMesh),intent(inout) :: this
    endsubroutine SELF_FreeMesh
  endinterface

  ! Element Types - From Table 4.1 of https://www.hopr-project.org/externals/Meshformat.pdf
  integer,parameter :: selfLineLinear = 1
  integer,parameter :: selfLineNonlinear = 2
  integer,parameter :: selfTriangleLinear = 3
  integer,parameter :: selfQuadLinear = 4
  integer,parameter :: selfQuadBilinear = 14
  integer,parameter :: selfTriangleNonlinear = 23
  integer,parameter :: selfQuadNonlinear = 24
  integer,parameter :: selfTetrahedronLinear = 104
  integer,parameter :: selfPyramidLinear = 105
  integer,parameter :: selfPrismLinear = 106
  integer,parameter :: selfHexahedronLinear = 108
  integer,parameter :: selfPyramidBilinear = 115
  integer,parameter :: selfPrismBilinear = 116
  integer,parameter :: selfHexahedronBilinear = 118
  integer,parameter :: selfTetrahedronNonlinear = 204
  integer,parameter :: selfPyramidNonlinear = 205
  integer,parameter :: selfPrismNonlinear = 206
  integer,parameter :: selfHexahedronNonlinear = 208
  integer,parameter :: self_BCDefault = 1
  integer,parameter :: self_nBCsDefault = 5

  !==============================================!
  ! --------------- File Types------------------ !
  !==============================================!
  integer,parameter :: SELF_MESH_ISM_V2_2D = 1
  integer,parameter :: SELF_MESH_ISM_V2_3D = 2
  integer,parameter :: SELF_MESH_HOPR_2D = 3
  integer,parameter :: SELF_MESH_HOPR_3D = 4

! //////////////////////////////////////////////// !
!   Boundary Condition parameters
!

  ! Conditions on the solution
  integer,parameter :: SELF_BC_PRESCRIBED = 100
  integer,parameter :: SELF_BC_RADIATION = 101
  integer,parameter :: SELF_BC_NONORMALFLOW = 102
  integer,parameter :: SELF_BC_INFLOW = 103
  integer,parameter :: SELF_BC_OUTFLOW = 104
  integer,parameter :: SELF_BC_NOSLIP = 105

  ! Conditions on the solution gradients
  integer,parameter :: SELF_BC_PRESCRIBED_STRESS = 200
  integer,parameter :: SELF_BC_NOSTRESS = 201 ! Defaults to setting the external gradient field as the negative of the internal gradient field so the average is zero.
  integer,parameter :: SELF_BC_STRESS = 202 ! Defaults to prolonging the internal gradient field to the external gradient field

endmodule SELF_Mesh
