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
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(Vector2D_t),public :: Vector2D
    character(3) :: backend = "gpu"
    type(c_ptr) :: interior_gpu
    type(c_ptr) :: boundary_gpu
    type(c_ptr) :: extBoundary_gpu
    type(c_ptr) :: avgBoundary_gpu
    type(c_ptr) :: boundaryNormal_gpu

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

    call this%UpdateDevice()

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

    call GridInterp_2D_gpu(this%interp%iMatrix_gpu,this%interior_gpu, &
                           f,this%N,this%M,2*this%nvar,this%nelem)

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

    ! The vector gradient is the (row,col) tensor df(i,j,e,v,row,col) = d(f_row)/dxi^col.
    ! Treating the two vector components as 2*nvar scalar fields, this is the scalar
    ! gradient of each, with the direction index (col) written to separate slots.
    call ScalarGradient_2D_gpu(this%interp%dMatrix_gpu,this%interior_gpu,df, &
                               this%interp%N,2*this%nvar,this%nelem)

  endsubroutine Gradient_Vector2D

  subroutine Divergence_Vector2D(this,df)
    implicit none
    class(Vector2D),intent(in) :: this
    type(c_ptr),intent(inout) :: df

    call Divergence_2D_gpu(this%interior_gpu,df,this%interp%dMatrix_gpu, &
                           this%interp%N,this%nvar,this%nelem)

  endsubroutine Divergence_Vector2D

endmodule SELF_Vector_2D
