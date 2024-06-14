!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
module SELF_Mesh

  use SELF_Constants

  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type :: SEMMesh
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
  endtype SEMMesh

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

endmodule SELF_Mesh
