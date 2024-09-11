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

module SELF_MappedScalar_1D

  use SELF_MappedScalar_1D_t
  use SELF_GPU
  use iso_c_binding

  implicit none

  type,extends(MappedScalar1D_t),public :: MappedScalar1D

  contains

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar1D

    procedure,public :: SideExchange => SideExchange_MappedScalar1D
    generic,public :: MappedDerivative => MappedDerivative_MappedScalar1D
    procedure,private :: MappedDerivative_MappedScalar1D

    generic,public :: MappedDGDerivative => MappedDGDerivative_MappedScalar1D
    procedure,private :: MappedDGDerivative_MappedScalar1D

  endtype MappedScalar1D

  interface
    subroutine JacobianWeight_1D_gpu(scalar,dxds,N,nVar,nEl) &
      bind(c,name="JacobianWeight_1D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,dxds
      integer(c_int),value :: N,nVar,nEl
    endsubroutine JacobianWeight_1D_gpu
  endinterface

  interface
    subroutine DGDerivative_BoundaryContribution_1D_gpu(bMatrix,qWeights,bf,df,N,nVar,nEl) &
      bind(c,name="DGDerivative_BoundaryContribution_1D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix,qWeights,bf,df
      integer(c_int),value :: N,nVar,nEl
    endsubroutine DGDerivative_BoundaryContribution_1D_gpu
  endinterface

contains

  subroutine SetInteriorFromEquation_MappedScalar1D(this,time)
    !!  Sets the this % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    real(prec),intent(in) :: time
    ! Local
    integer :: iVar

    do ivar = 1,this%nvar
      this%interior(:,:,ivar) = this%eqn(ivar)%evaluate(this%geometry%x%interior)
    enddo
    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine SetInteriorFromEquation_MappedScalar1D

  subroutine SideExchange_MappedScalar1D(this,mesh)
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    type(Mesh1D),intent(inout) :: mesh
    ! Local
    integer :: e1,e2,s1,s2
    integer :: ivar

    call gpuCheck(hipMemcpy(c_loc(this%boundary),this%boundary_gpu,sizeof(this%boundary),hipMemcpyDeviceToHost))

    do ivar = 1,this%nvar
      do e1 = 1,mesh%nElem

        if(e1 == 1) then

          s1 = 2
          e2 = e1+1
          s2 = 1
          this%extBoundary(s1,e1,ivar) = this%boundary(s2,e2,ivar)

        elseif(e1 == mesh%nElem) then

          s1 = 1
          e2 = e1-1
          s2 = 2
          this%extBoundary(s1,e1,ivar) = this%boundary(s2,e2,ivar)

        else

          s1 = 1
          e2 = e1-1
          s2 = 2
          this%extBoundary(s1,e1,ivar) = this%boundary(s2,e2,ivar)

          s1 = 2
          e2 = e1+1
          s2 = 1
          this%extBoundary(s1,e1,ivar) = this%boundary(s2,e2,ivar)

        endif

      enddo
    enddo
    call gpuCheck(hipMemcpy(this%extboundary_gpu,c_loc(this%extboundary),sizeof(this%extboundary),hipMemcpyHostToDevice))

  endsubroutine SideExchange_MappedScalar1D

  subroutine MappedDerivative_MappedScalar1D(this,dF)
    implicit none
    class(MappedScalar1D),intent(in) :: this
    type(c_ptr),intent(inout) :: df
    ! Local
    integer :: iEl,iVar,i,ii
    real(prec) :: dfloc

    call this%Derivative(df)
    call JacobianWeight_1D_gpu(df,this%geometry%dxds%interior_gpu,this%N,this%nVar,this%nelem)

  endsubroutine MappedDerivative_MappedScalar1D

  subroutine MappedDGDerivative_MappedScalar1D(this,dF)
    implicit none
    class(MappedScalar1D),intent(in) :: this
    type(c_ptr),intent(inout) :: df
    ! Local
    integer :: iEl,iVar,i,ii
    real(prec) :: dfloc

    call self_blas_matrixop_1d(this%interp%dgMatrix_gpu, &
                               this%interior_gpu, &
                               df,this%N+1,this%N+1, &
                               this%nvar*this%nelem,this%blas_handle)

    call DGDerivative_BoundaryContribution_1D_gpu(this%interp%bMatrix_gpu, &
                                                  this%interp%qWeights_gpu, &
                                                  this%boundarynormal_gpu,df, &
                                                  this%N,this%nVar,this%nelem)

    call JacobianWeight_1D_gpu(df,this%geometry%dxds%interior_gpu,this%N,this%nVar,this%nelem)

  endsubroutine MappedDGDerivative_MappedScalar1D

endmodule SELF_MappedScalar_1D
