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

module SELF_MappedTwoPointVector_3D_t
  !! Geometry-aware two-point vector and split-form divergence for 3-D.
  !!
  !! The MappedDivergence routine follows the construction described in
  !!
  !!   Winters, Kopriva, Gassner, Hindenlang,
  !!   "Construction of Modern Robust Nodal Discontinuous Galerkin Spectral
  !!    Element Methods for the Compressible Navier-Stokes Equations",
  !!   Lecture Notes in Computational Science and Engineering, 2021.
  !!
  !! The contravariant two-point flux in the r-th computational direction is
  !!
  !!   F~^r_{(i,n),j,k} = sum_d  (Ja^r_d(i,j,k) + Ja^r_d(n,j,k))/2 * f^d_{(i,n),j,k}
  !!
  !! and similarly for the xi^2 and xi^3 sums (averaging (i,j,k)-(i,n,k) and
  !! (i,j,k)-(i,j,n) respectively).  The physical divergence is
  !!
  !!   (1/J_{i,j,k}) * 2 * sum_n [ D_{n,i} F~^1 + D_{n,j} F~^2 + D_{n,k} F~^3 ]
  !!
  !! interior(n,i,j,k,iEl,iVar,d) stores the physical-space two-point flux
  !! f^d in the d-th physical direction.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Geometry_3D
  use SELF_TwoPointVector_3D
  use iso_c_binding

  implicit none

  type,extends(TwoPointVector3D),public :: MappedTwoPointVector3D_t

    logical :: geometry_associated = .false.
    type(SEMHex),pointer :: geometry => null()

  contains

    procedure,public :: AssociateGeometry => AssociateGeometry_MappedTwoPointVector3D_t
    procedure,public :: DissociateGeometry => DissociateGeometry_MappedTwoPointVector3D_t

    generic,public :: MappedDivergence => MappedDivergence_MappedTwoPointVector3D_t
    procedure,private :: MappedDivergence_MappedTwoPointVector3D_t

  endtype MappedTwoPointVector3D_t

contains

  subroutine AssociateGeometry_MappedTwoPointVector3D_t(this,geometry)
    implicit none
    class(MappedTwoPointVector3D_t),intent(inout) :: this
    type(SEMHex),target,intent(in) :: geometry

    if(.not. associated(this%geometry)) then
      this%geometry => geometry
      this%geometry_associated = .true.
    endif

  endsubroutine AssociateGeometry_MappedTwoPointVector3D_t

  subroutine DissociateGeometry_MappedTwoPointVector3D_t(this)
    implicit none
    class(MappedTwoPointVector3D_t),intent(inout) :: this

    if(associated(this%geometry)) then
      this%geometry => null()
      this%geometry_associated = .false.
    endif

  endsubroutine DissociateGeometry_MappedTwoPointVector3D_t

  subroutine MappedDivergence_MappedTwoPointVector3D_t(this,df)
    !! Computes the physical-space divergence of a 3-D split-form vector field
    !! on a curvilinear mesh following Winters, Kopriva, Gassner, Hindenlang.
    !!
    !! interior(n,i,j,k,iEl,iVar,d) must contain the physical-space two-point
    !! flux f^d between the node pair in each coordinate direction.
    !! Metric terms are averaged between the two nodes of each pair to satisfy
    !! the discrete metric identities required for entropy conservation.
    implicit none
    class(MappedTwoPointVector3D_t),intent(in) :: this
    real(prec),intent(out) :: &
      df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nElem,1:this%nVar)
    ! Local
    integer :: i,j,k,nn,d,iEl,iVar
    real(prec) :: dfLoc,Fc1,Fc2,Fc3

    do concurrent(i=1:this%N+1,j=1:this%N+1,k=1:this%N+1, &
                  iEl=1:this%nElem,iVar=1:this%nVar)

      dfLoc = 0.0_prec
      do nn = 1,this%N+1

        ! Contravariant two-point flux in the xi^1 direction.
        ! Metric terms are averaged between nodes (i,j,k) and (nn,j,k).
        Fc1 = 0.0_prec
        do d = 1,3
          Fc1 = Fc1+0.5_prec*( &
                this%geometry%dsdx%interior(i,j,k,iEl,1,d,1)+ &
                this%geometry%dsdx%interior(nn,j,k,iEl,1,d,1))* &
                this%interior(nn,i,j,k,iEl,iVar,d)
        enddo

        ! Contravariant two-point flux in the xi^2 direction.
        ! Metric terms are averaged between nodes (i,j,k) and (i,nn,k).
        Fc2 = 0.0_prec
        do d = 1,3
          Fc2 = Fc2+0.5_prec*( &
                this%geometry%dsdx%interior(i,j,k,iEl,1,d,2)+ &
                this%geometry%dsdx%interior(i,nn,k,iEl,1,d,2))* &
                this%interior(nn,i,j,k,iEl,iVar,d)
        enddo

        ! Contravariant two-point flux in the xi^3 direction.
        ! Metric terms are averaged between nodes (i,j,k) and (i,j,nn).
        Fc3 = 0.0_prec
        do d = 1,3
          Fc3 = Fc3+0.5_prec*( &
                this%geometry%dsdx%interior(i,j,k,iEl,1,d,3)+ &
                this%geometry%dsdx%interior(i,j,nn,iEl,1,d,3))* &
                this%interior(nn,i,j,k,iEl,iVar,d)
        enddo

        dfLoc = dfLoc+ &
                this%interp%dMatrix(nn,i)*Fc1+ &
                this%interp%dMatrix(nn,j)*Fc2+ &
                this%interp%dMatrix(nn,k)*Fc3

      enddo

      df(i,j,k,iEl,iVar) = 2.0_prec*dfLoc/this%geometry%J%interior(i,j,k,iEl,1)

    enddo

  endsubroutine MappedDivergence_MappedTwoPointVector3D_t

endmodule SELF_MappedTwoPointVector_3D_t
