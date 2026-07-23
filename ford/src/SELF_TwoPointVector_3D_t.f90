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

module SELF_TwoPointVector_3D_t

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_Data

  use iso_c_binding

  implicit none

  type,extends(SELF_DataObj),public :: TwoPointVector3D_t
    !! A two-point vector field for use in split-form DGSEM in 3-D.
    !!
    !! Memory layout: interior(n, i, j, k, nEl, nVar, 1:3)
    !!   dim 1 (n)   : two-point index; looping over n is equivalent to looping
    !!                 over n in the split-form divergence sum (0:N in 0-based)
    !!   dim 2 (i)   : first computational coordinate direction (xi^1)
    !!   dim 3 (j)   : second computational coordinate direction (xi^2)
    !!   dim 4 (k)   : third computational coordinate direction (xi^3)
    !!   dim 5       : element index
    !!   dim 6       : variable index
    !!   dim 7 (idir): vector direction (1, 2, or 3)

    real(prec),pointer,contiguous,dimension(:,:,:,:,:,:,:) :: interior

  contains

    procedure,public :: Init => Init_TwoPointVector3D_t
    procedure,public :: Free => Free_TwoPointVector3D_t

    procedure,public :: UpdateHost => UpdateHost_TwoPointVector3D_t
    procedure,public :: UpdateDevice => UpdateDevice_TwoPointVector3D_t

    generic,public :: Divergence => Divergence_TwoPointVector3D_t
    procedure,private :: Divergence_TwoPointVector3D_t

  endtype TwoPointVector3D_t

contains

  subroutine Init_TwoPointVector3D_t(this,interp,nVar,nElem)
    !! Allocate the interior array for a 3-D two-point vector field.
    !! The interior array has rank 7 with layout (n,i,j,k,nEl,nVar,idir).
    !!
    !! Requires Gauss-Lobatto quadrature nodes (controlNodeType=GAUSS_LOBATTO).
    implicit none
    class(TwoPointVector3D_t),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer :: i

    if(interp%controlNodeType /= GAUSS_LOBATTO) then
      print*,__FILE__//" : TwoPointVector3D requires Gauss-Lobatto quadrature nodes."
      stop 1
    endif

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:interp%N+1, &
                           1:nElem,1:nVar,1:3))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:3*nVar))

    ! Initialize equation parser to prevent segmentation faults with amdflang
    ! when the parser functions are not allocated (see SELF_Vector_2D_t.f90)
    do i = 1,3*nVar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    this%interior = 0.0_prec

  endsubroutine Init_TwoPointVector3D_t

  subroutine Free_TwoPointVector3D_t(this)
    implicit none
    class(TwoPointVector3D_t),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%meta)
    deallocate(this%eqn)

  endsubroutine Free_TwoPointVector3D_t

  subroutine UpdateHost_TwoPointVector3D_t(this)
    implicit none
    class(TwoPointVector3D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateHost_TwoPointVector3D_t

  subroutine UpdateDevice_TwoPointVector3D_t(this)
    implicit none
    class(TwoPointVector3D_t),intent(inout) :: this
    if(.false.) this%N = this%N ! CPU stub; suppress unused-dummy-argument warning
  endsubroutine UpdateDevice_TwoPointVector3D_t

  subroutine Divergence_TwoPointVector3D_t(this,df)
    !! Computes the split-form (two-point) divergence of a 3-D vector field
    !! in the reference element (computational coordinates).
    !!
    !! The split-form divergence at node (i,j,k) is
    !!
    !!   (nabla.F)_{i,j,k} = 2 sum_n [ D_{n,i} F^1(n,i,j,k)
    !!                                + D_{n,j} F^2(n,i,j,k)
    !!                                + D_{n,k} F^3(n,i,j,k) ]
    !!
    !! where D is the standard derivative matrix (dMatrix) and
    !! F^idir(n,i,j,k,...) = interior(n,i,j,k,iEl,iVar,idir) stores the
    !! two-point flux between nodes i (j or k) and n in the idir-th direction.
    !!
    !! The interior array is assumed to hold contravariant (J-scaled) two-point
    !! fluxes.  To obtain the physical divergence, divide the result by the
    !! element Jacobian J(i,j,k,iEl).
    implicit none
    class(TwoPointVector3D_t),intent(in) :: this
    real(prec),intent(out) :: df(1:this%N+1,1:this%N+1,1:this%N+1,1:this%nElem,1:this%nVar)
    ! Local
    integer :: i,j,k,nn,iEl,iVar
    real(prec) :: dfLoc

    do concurrent(i=1:this%N+1,j=1:this%N+1,k=1:this%N+1, &
                  iEl=1:this%nElem,iVar=1:this%nVar)

      dfLoc = 0.0_prec
      do nn = 1,this%N+1
        dfLoc = dfLoc+ &
                this%interp%dSplitMatrix(nn,i)*this%interior(nn,i,j,k,iEl,iVar,1)+ &
                this%interp%dSplitMatrix(nn,j)*this%interior(nn,i,j,k,iEl,iVar,2)+ &
                this%interp%dSplitMatrix(nn,k)*this%interior(nn,i,j,k,iEl,iVar,3)
      enddo
      df(i,j,k,iEl,iVar) = 2.0_prec*dfLoc

    enddo

  endsubroutine Divergence_TwoPointVector3D_t

endmodule SELF_TwoPointVector_3D_t
