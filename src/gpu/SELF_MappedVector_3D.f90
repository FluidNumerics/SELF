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

module SELF_MappedVector_3D

  use SELF_MappedVector_3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use SELF_GPUBLAS
  use iso_c_binding

  implicit none

  type,extends(MappedVector3D_t),public :: MappedVector3D

  contains

    procedure,public :: SideExchange => SideExchange_MappedVector3D

    generic,public :: MappedDivergence => MappedDivergence_MappedVector3D
    procedure,private :: MappedDivergence_MappedVector3D

    generic,public :: MappedDGDivergence => MappedDGDivergence_MappedVector3D
    procedure,private :: MappedDGDivergence_MappedVector3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector3D

  endtype MappedVector3D

  interface
    subroutine ContravariantProjection_3D_gpu(f,dsdx,N,nvar,nel) &
      bind(c,name="ContravariantProjection_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,dsdx
      integer(c_int),value :: N,nvar,nel
    endsubroutine ContravariantProjection_3D_gpu
  endinterface

contains
  subroutine SetInteriorFromEquation_MappedVector3D(this,geometry,time)
    !!  Sets the this % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,k,iEl,iVar
    real(prec) :: x
    real(prec) :: y
    real(prec) :: z

    do iVar = 1,this%nVar
      do iEl = 1,this%nElem
        do k = 1,this%interp%N+1
          do j = 1,this%interp%N+1
            do i = 1,this%interp%N+1

              ! Get the mesh positions
              x = geometry%x%interior(i,j,k,iEl,1,1)
              y = geometry%x%interior(i,j,k,iEl,1,2)
              z = geometry%x%interior(i,j,k,iEl,1,3)

              this%interior(i,j,k,iEl,iVar,1) = &
                this%eqn(1+3*(iVar-1))%Evaluate((/x,y,z,time/))

              this%interior(i,j,k,iEl,iVar,2) = &
                this%eqn(2+3*(iVar-1))%Evaluate((/x,y,z,time/))

              this%interior(i,j,k,iEl,iVar,3) = &
                this%eqn(3+3*(iVar-1))%Evaluate((/x,y,z,time/))

            enddo
          enddo
        enddo
      enddo
    enddo

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine SetInteriorFromEquation_MappedVector3D

  subroutine SideExchange_MappedVector3D(this,mesh)
    implicit none
    class(MappedVector3D),intent(inout) :: this
    type(Mesh3D),intent(inout) :: mesh
    ! Local
    integer :: rankId,offset

    rankId = mesh%decomp%rankId
    offset = mesh%decomp%offsetElem(rankId+1)

    !call this%MPIExchangeAsync(mesh%decomp,mesh,resetCount=.true.)

    call SideExchange_3D_gpu(this%extboundary_gpu, &
                             this%boundary_gpu,mesh%sideinfo_gpu,mesh%decomp%elemToRank_gpu, &
                             mesh%decomp%rankid,offset,this%interp%N,3*this%nvar,this%nelem)

    !call decomp%FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    !call this%ApplyFlip(decomp,mesh)

  endsubroutine SideExchange_MappedVector3D

  subroutine MappedDivergence_MappedVector3D(this,df)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector3D),intent(in) :: this
    type(c_ptr),intent(inout) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:,:)
    type(c_ptr) :: fc

    ! Contravariant projection
    call ContravariantProjection_3D_gpu(this%interior_gpu, &
                                        this%geometry%dsdx%interior_gpu,this%interp%N,this%nvar,this%nelem)

    call c_f_pointer(this%interior_gpu,f_p, &
                     [this%interp%N+1,this%interp%N+1,this%interp%N+1,this%nelem,this%nvar,3])

    fc = c_loc(f_p(1,1,1,1,1,1))
    call self_blas_matrixop_dim1_3d(this%interp%dMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,2))
    call self_blas_matrixop_dim2_3d(this%interp%dMatrix_gpu,fc,df,1.0_c_prec, &
                                    this%interp%N,this%interp%N,this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,3))
    call self_blas_matrixop_dim3_3d(this%interp%dMatrix_gpu,fc,df,1.0_c_prec, &
                                    this%interp%N,this%interp%N,this%nvar,this%nelem,this%blas_handle)
    f_p => null()

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu,this%interp%N,this%nVar,this%nelem)

  endsubroutine MappedDivergence_MappedVector3D

  subroutine MappedDGDivergence_MappedVector3D(this,df)
      !! Computes the divergence of a 3-D vector using the weak form
      !! On input, the  attribute of the vector
      !! is assigned and the  attribute is set to the physical
      !! directions of the vector. This method will project the vector
      !! onto the contravariant basis vectors.
    implicit none
    class(MappedVector3D),intent(in) :: this
    type(c_ptr),intent(inout) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:,:)
    type(c_ptr) :: fc

    ! Contravariant projection
    call ContravariantProjection_3D_gpu(this%interior_gpu, &
                                        this%geometry%dsdx%interior_gpu,this%interp%N,this%nvar,this%nelem)

    call c_f_pointer(this%interior_gpu,f_p, &
                     [this%interp%N+1,this%interp%N+1,this%interp%N+1,this%nelem,this%nvar,3])

    fc = c_loc(f_p(1,1,1,1,1,1))
    call self_blas_matrixop_dim1_3d(this%interp%dgMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,2))
    call self_blas_matrixop_dim2_3d(this%interp%dgMatrix_gpu,fc,df,1.0_c_prec, &
                                    this%interp%N,this%interp%N,this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,3))
    call self_blas_matrixop_dim3_3d(this%interp%dgMatrix_gpu,fc,df,1.0_c_prec, &
                                    this%interp%N,this%interp%N,this%nvar,this%nelem,this%blas_handle)
    f_p => null()

    ! Boundary terms
    call DG_BoundaryContribution_3D_gpu(this%interp%bmatrix_gpu,this%interp%qweights_gpu, &
                                        this%boundarynormal_gpu,df,this%interp%N,this%nvar,this%nelem)

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu,this%interp%N,this%nVar,this%nelem)

  endsubroutine MappedDGDivergence_MappedVector3D

endmodule SELF_MappedVector_3D
