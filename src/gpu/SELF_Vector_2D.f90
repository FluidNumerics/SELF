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

module SELF_Vector_2D

  use SELF_Constants
  use SELF_Vector_2D_t
  use SELF_GPU
  use SELF_GPUBLAS
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(Vector2D_t),public :: Vector2D
    character(3) :: backend = "gpu"
    type(c_ptr) :: blas_handle
    type(c_ptr) :: interior_gpu
    type(c_ptr) :: boundary_gpu
    type(c_ptr) :: extBoundary_gpu
    type(c_ptr) :: avgBoundary_gpu
    type(c_ptr) :: boundaryNormal_gpu
    type(c_ptr) :: interpWork

  contains

    procedure,public :: Init => Init_Vector2D
    procedure,public :: Free => Free_Vector2D

    procedure,public :: UpdateHost => UpdateHost_Vector2D
    procedure,public :: UpdateDevice => UpdateDevice_Vector2D

    procedure,public :: BoundaryInterp => BoundaryInterp_Vector2D
    procedure,public :: AverageSides => AverageSides_Vector2D

    generic,public :: GridInterp => GridInterp_Vector2D
    procedure,private :: GridInterp_Vector2D

    generic,public :: Gradient => Gradient_Vector2D
    procedure,private :: Gradient_Vector2D

    generic,public :: Divergence => Divergence_Vector2D
    procedure,private :: Divergence_Vector2D

  endtype Vector2D

contains

  subroutine Init_Vector2D(this,interp,nVar,nElem)
    implicit none
    class(Vector2D),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! local
    integer :: i
    integer(c_size_t) :: workSize

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:nelem,1:nvar,1:2), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%avgBoundary(1:interp%N+1,1:4,1:nelem,1:nvar,1:2), &
             this%boundaryNormal(1:interp%N+1,1:4,1:nelem,1:nvar))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:2*nVar))

    ! Initialize equation parser
    ! This is done to prevent segmentation faults that arise
    ! when building with amdflang that are traced back to
    ! feqparse_functions.f90 : finalize routine
    ! When the equation parser is not initialized, the
    ! functions are not allocated, which I think are the
    ! source of the segfault - joe@fluidnumerics.com
    do i = 1,2*nvar
      this%eqn(i) = EquationParser('f=0',(/'x','y','z','t'/))
    enddo

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%boundarynormal = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec

    call gpuCheck(hipMalloc(this%interior_gpu,sizeof(this%interior)))
    call gpuCheck(hipMalloc(this%boundary_gpu,sizeof(this%boundary)))
    call gpuCheck(hipMalloc(this%extBoundary_gpu,sizeof(this%extBoundary)))
    call gpuCheck(hipMalloc(this%avgBoundary_gpu,sizeof(this%avgBoundary)))
    call gpuCheck(hipMalloc(this%boundaryNormal_gpu,sizeof(this%boundaryNormal)))
    workSize = (interp%N+1)*(interp%M+1)*nelem*nvar*2*prec
    call gpuCheck(hipMalloc(this%interpWork,workSize))

    call this%UpdateDevice()

    call hipblasCheck(hipblasCreate(this%blas_handle))

  endsubroutine Init_Vector2D

  subroutine Free_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    this%interp => null()
    this%nVar = 0
    this%nElem = 0

    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%boundaryNormal)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)

    deallocate(this%meta)
    deallocate(this%eqn)

    call gpuCheck(hipFree(this%interior_gpu))
    call gpuCheck(hipFree(this%boundary_gpu))
    call gpuCheck(hipFree(this%extBoundary_gpu))
    call gpuCheck(hipFree(this%avgBoundary_gpu))
    call gpuCheck(hipFree(this%boundaryNormal_gpu))
    call gpuCheck(hipFree(this%interpWork))
    call hipblasCheck(hipblasDestroy(this%blas_handle))

  endsubroutine Free_Vector2D

  subroutine UpdateHost_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(c_loc(this%interior),this%interior_gpu,sizeof(this%interior),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%boundary),this%boundary_gpu,sizeof(this%boundary),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%extboundary),this%extboundary_gpu,sizeof(this%extboundary),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%avgboundary),this%avgboundary_gpu,sizeof(this%avgboundary),hipMemcpyDeviceToHost))
    call gpuCheck(hipMemcpy(c_loc(this%boundaryNormal),this%boundaryNormal_gpu,sizeof(this%boundaryNormal),hipMemcpyDeviceToHost))

  endsubroutine UpdateHost_Vector2D

  subroutine UpdateDevice_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%boundary_gpu,c_loc(this%boundary),sizeof(this%boundary),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%extboundary_gpu,c_loc(this%extboundary),sizeof(this%extboundary),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%avgboundary_gpu,c_loc(this%avgboundary),sizeof(this%avgboundary),hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%boundaryNormal_gpu,c_loc(this%boundaryNormal),sizeof(this%boundaryNormal),hipMemcpyHostToDevice))

  endsubroutine UpdateDevice_Vector2D

  subroutine GridInterp_Vector2D(this,f)
    implicit none
    class(Vector2D),intent(inout) :: this
    type(c_ptr),intent(inout) :: f

    call self_blas_matrixop_dim1_2d(this%interp%iMatrix_gpu,this%interior_gpu, &
                                    this%interpWork,this%N,this%M,2*this%nvar,this%nelem, &
                                    this%blas_handle)

    call self_blas_matrixop_dim2_2d(this%interp%iMatrix_gpu,this%interpWork,f, &
                                    0.0_c_prec,this%N,this%M,2*this%nvar,this%nelem, &
                                    this%blas_handle)

  endsubroutine GridInterp_Vector2D

  subroutine AverageSides_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    call Average_gpu(this%avgBoundary_gpu,this%boundary_gpu,this%extBoundary_gpu,size(this%boundary))

  endsubroutine AverageSides_Vector2D

  subroutine BoundaryInterp_Vector2D(this)
    implicit none
    class(Vector2D),intent(inout) :: this

    call BoundaryInterp_2D_gpu(this%interp%bMatrix_gpu,this%interior_gpu,this%boundary_gpu, &
                               this%interp%N,2*this%nvar,this%nelem)

  endsubroutine BoundaryInterp_Vector2D

  subroutine Gradient_Vector2D(this,df)
    implicit none
    class(Vector2D),intent(in) :: this
    type(c_ptr),intent(inout) :: df
    !Local
    real(prec),pointer :: df_p(:,:,:,:,:,:)
    real(prec),pointer :: dfloc(:,:,:,:)
    type(c_ptr) :: dfc

    call c_f_pointer(df,df_p,[this%interp%N+1,this%interp%N+1,this%nelem,this%nvar,2,2])

    dfloc(1:,1:,1:,1:) => df_p(1:,1:,1:,1:,1,1)
    dfc = c_loc(dfloc)
    call self_blas_matrixop_dim1_2d(this%interp%dMatrix_gpu,this%interior_gpu,dfc, &
                                    this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    dfloc(1:,1:,1:,1:) => df_p(1:,1:,1:,1:,1,2)
    dfc = c_loc(dfloc)
    call self_blas_matrixop_dim2_2d(this%interp%dMatrix_gpu,this%interior_gpu,dfc,0.0_c_prec, &
                                    this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    dfloc => null()
    df_p => null()

  endsubroutine Gradient_Vector2D

  subroutine Divergence_Vector2D(this,df)
    implicit none
    class(Vector2D),intent(in) :: this
    type(c_ptr),intent(inout) :: df
    !Local
    !real(prec),pointer :: f_p(:,:,:,:,:)
    !type(c_ptr) :: fc

    call Divergence_2D_gpu(this%interior_gpu,df,this%interp%dMatrix_gpu, &
                           this%interp%N,this%nvar,this%nelem)
    ! call c_f_pointer(this%interior_gpu,f_p,[this%interp%N+1,this%interp%N+1,this%nelem,this%nvar,2])

    ! fc = c_loc(f_p(1,1,1,1,1))
    ! call self_blas_matrixop_dim1_2d(this%interp%dMatrix_gpu,fc,df, &
    !                                 this%interp%N,this%interp%N,this%nvar,this%nelem,this%blas_handle)

    ! fc = c_loc(f_p(1,1,1,1,2))
    ! call self_blas_matrixop_dim2_2d(this%interp%dMatrix_gpu,fc,df, &
    !                                 1.0_c_prec,this%interp%N,this%interp%N,this%nvar,this%nelem,this%blas_handle)

    ! f_p => null()

  endsubroutine Divergence_Vector2D

endmodule SELF_Vector_2D
