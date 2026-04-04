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
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_MappedTwoPointVector_2D_t
  !! Geometry-aware two-point vector and split-form divergence for 2-D.
  !!
  !! The MappedDivergence routine follows the construction described in
  !!
  !!   Winters, Kopriva, Gassner, Hindenlang,
  !!   "Construction of Modern Robust Nodal Discontinuous Galerkin Spectral
  !!    Element Methods for the Compressible Navier-Stokes Equations",
  !!   Lecture Notes in Computational Science and Engineering, 2021.
  !!
  !! The key ingredient is that the contravariant two-point flux in the r-th
  !! computational direction is formed by projecting the physical two-point
  !! flux onto *averaged* contravariant basis vectors:
  !!
  !!   F~^r_{(i,n),j} = sum_d  (Ja^r_d(i,j) + Ja^r_d(n,j))/2 * f^d_{(i,n),j}
  !!
  !! where dsdx%interior(i,j,iEl,1,d,r) = J*a^r_d stores the scaled
  !! contravariant basis vectors (metric terms times Jacobian).
  !! The physical-space divergence at (i,j) is then
  !!
  !!   (1/J_{i,j}) * 2 * sum_n [ D_{n,i} F~^1_{(i,n),j} + D_{n,j} F~^2_{i,(j,n)} ]
  !!
  !! The interior array stores the physical-space two-point flux
  !! interior(n,i,j,iEl,iVar,d) = f^d_{(i or j, n)} where d is the physical
  !! direction.  Metric averaging is performed inside MappedDivergence.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Geometry_2D
  use SELF_TwoPointVector_2D
  use iso_c_binding

  implicit none

  type,extends(TwoPointVector2D),public :: MappedTwoPointVector2D_t

    logical :: geometry_associated = .false.
    type(SEMQuad),pointer :: geometry => null()

  contains

    procedure,public :: AssociateGeometry => AssociateGeometry_MappedTwoPointVector2D_t
    procedure,public :: DissociateGeometry => DissociateGeometry_MappedTwoPointVector2D_t

    generic,public :: MappedDivergence => MappedDivergence_MappedTwoPointVector2D_t
    procedure,private :: MappedDivergence_MappedTwoPointVector2D_t

  endtype MappedTwoPointVector2D_t

contains

  subroutine AssociateGeometry_MappedTwoPointVector2D_t(this,geometry)
    implicit none
    class(MappedTwoPointVector2D_t),intent(inout) :: this
    type(SEMQuad),target,intent(in) :: geometry

    if(.not. associated(this%geometry)) then
      this%geometry => geometry
      this%geometry_associated = .true.
    endif

  endsubroutine AssociateGeometry_MappedTwoPointVector2D_t

  subroutine DissociateGeometry_MappedTwoPointVector2D_t(this)
    implicit none
    class(MappedTwoPointVector2D_t),intent(inout) :: this

    if(associated(this%geometry)) then
      this%geometry => null()
      this%geometry_associated = .false.
    endif

  endsubroutine DissociateGeometry_MappedTwoPointVector2D_t

  subroutine MappedDivergence_MappedTwoPointVector2D_t(this,df)
    !! Computes the physical-space divergence of a 2-D split-form vector field
    !! on a curvilinear mesh.
    !!
    !! Convention (following Trixi.jl for curved meshes):
    !! interior(n,i,j,iEl,iVar,r) holds the pre-projected SCALAR contravariant
    !! two-point flux for the r-th computational direction:
    !!
    !!   interior(n,i,j,iEl,iVar,1) = avg(Ja^1) . F_EC(s(i,j), s(n,j))
    !!   interior(n,i,j,iEl,iVar,2) = avg(Ja^2) . F_EC(s(i,j), s(i,n))
    !!
    !! The metric averaging and direction-correct pairing are the caller's
    !! responsibility (e.g. TwoPointFluxMethod in ECDGModel).
    !! MappedDivergence applies the reference-element split-form sum and
    !! divides by J.
    implicit none
    class(MappedTwoPointVector2D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%nElem,1:this%nVar)
    ! Local
    integer :: i,j,iEl,iVar

    call this%Divergence(df)

    do concurrent(i=1:this%N+1,j=1:this%N+1,iEl=1:this%nElem,iVar=1:this%nVar)
      df(i,j,iEl,iVar) = df(i,j,iEl,iVar)/this%geometry%J%interior(i,j,iEl,1)
    enddo

  endsubroutine MappedDivergence_MappedTwoPointVector2D_t

endmodule SELF_MappedTwoPointVector_2D_t
