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

module SELF_MappedScalar_3D

  use SELF_MappedScalar_3D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(MappedScalar3D_t),public :: MappedScalar3D

    type(c_ptr) :: jas_gpu ! jacobian weighted scalar for gradient calculation

  contains
    procedure,public :: Init => Init_MappedScalar3D
    procedure,public :: Free => Free_MappedScalar3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar3D

    procedure,public :: SideExchange => SideExchange_MappedScalar3D

    generic,public :: MappedGradient => MappedGradient_MappedScalar3D
    procedure,private :: MappedGradient_MappedScalar3D

    generic,public :: MappedDGGradient => MappedDGGradient_MappedScalar3D
    procedure,private :: MappedDGGradient_MappedScalar3D

  endtype MappedScalar3D

  interface
    subroutine ContravariantWeight_3D_gpu(f,dsdx,jaf,N,nvar,nel) &
      bind(c,name="ContravariantWeight_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,dsdx,jaf
      integer(c_int),value :: N,nvar,nel
    endsubroutine ContravariantWeight_3D_gpu
  endinterface

  interface
    subroutine NormalWeight_3D_gpu(fb,nhat,nscale,fbn,N,nvar,nel) &
      bind(c,name="NormalWeight_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: fb,nhat,nscale,fbn
      integer(c_int),value :: N,nvar,nel
    endsubroutine NormalWeight_3D_gpu
  endinterface

contains

  subroutine Init_MappedScalar3D(this,interp,nVar,nElem)
    implicit none
    class(MappedScalar3D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem
    ! Local
    integer(c_size_t) :: workSize

    this%interp => interp
    this%nVar = nVar
    this%nElem = nElem
    this%N = interp%N
    this%M = interp%M

    allocate(this%interior(1:interp%N+1,1:interp%N+1,1:interp%N+1,1:nelem,1:nvar), &
             this%boundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%extBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%avgBoundary(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:nvar), &
             this%boundarynormal(1:interp%N+1,1:interp%N+1,1:6,1:nelem,1:3*nvar))

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec
    this%boundarynormal = 0.0_prec

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

    call gpuCheck(hipMalloc(this%interior_gpu,sizeof(this%interior)))
    call gpuCheck(hipMalloc(this%boundary_gpu,sizeof(this%boundary)))
    call gpuCheck(hipMalloc(this%extBoundary_gpu,sizeof(this%extBoundary)))
    call gpuCheck(hipMalloc(this%avgBoundary_gpu,sizeof(this%avgBoundary)))
    call gpuCheck(hipMalloc(this%boundarynormal_gpu,sizeof(this%boundarynormal)))
    workSize = (interp%N+1)*(interp%N+1)*(interp%M+1)*nelem*nvar*prec
    call gpuCheck(hipMalloc(this%interpWork1,workSize))
    workSize = (interp%N+1)*(interp%M+1)*(interp%M+1)*nelem*nvar*prec
    call gpuCheck(hipMalloc(this%interpWork2,workSize))
    workSize = (interp%N+1)*(interp%N+1)*(interp%N+1)*nelem*nvar*9*prec
    call gpuCheck(hipMalloc(this%jas_gpu,workSize))

    call this%UpdateDevice()

    call hipblasCheck(hipblasCreate(this%blas_handle))

  endsubroutine Init_MappedScalar3D

  subroutine Free_MappedScalar3D(this)
    implicit none
    class(MappedScalar3D),intent(inout) :: this

    this%nVar = 0
    this%nElem = 0
    this%interp => null()
    deallocate(this%interior)
    deallocate(this%boundary)
    deallocate(this%extBoundary)
    deallocate(this%avgBoundary)
    deallocate(this%boundarynormal)
    deallocate(this%meta)
    deallocate(this%eqn)

    call gpuCheck(hipFree(this%interior_gpu))
    call gpuCheck(hipFree(this%boundary_gpu))
    call gpuCheck(hipFree(this%extBoundary_gpu))
    call gpuCheck(hipFree(this%avgBoundary_gpu))
    call gpuCheck(hipFree(this%boundarynormal_gpu))
    call gpuCheck(hipFree(this%interpWork1))
    call gpuCheck(hipFree(this%interpWork2))
    call gpuCheck(hipFree(this%jas_gpu))
    call hipblasCheck(hipblasDestroy(this%blas_handle))

  endsubroutine Free_MappedScalar3D

  subroutine SetInteriorFromEquation_MappedScalar3D(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
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

              this%interior(i,j,k,iEl,iVar) = &
                this%eqn(iVar)%Evaluate((/x,y,z,time/))

            enddo
          enddo
        enddo
      enddo
    enddo

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine SetInteriorFromEquation_MappedScalar3D

  subroutine SideExchange_MappedScalar3D(this,mesh,decomp)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: offset

    offset = decomp%offsetElem(decomp%rankId+1)

    !call this%MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    call SideExchange_3D_gpu(this%extboundary_gpu, &
                             this%boundary_gpu,mesh%sideinfo_gpu,decomp%elemToRank_gpu, &
                             decomp%rankid,offset,this%interp%N,this%nvar,this%nelem)

    !call decomp%FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    !call this%ApplyFlip(decomp,mesh)

  endsubroutine SideExchange_MappedScalar3D

  subroutine MappedGradient_MappedScalar3D(this,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(c_ptr),intent(out) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:,:)
    type(c_ptr) :: fc

    call ContravariantWeight_3D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! From Vector divergence
    call c_f_pointer(this%jas_gpu,f_p, &
                     [this%interp%N+1,this%interp%N+1,this%interp%N+1,this%nelem,3*this%nvar,3])

    fc = c_loc(f_p(1,1,1,1,1,1))
    call self_blas_matrixop_dim1_3d(this%interp%dMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,2))
    call self_blas_matrixop_dim2_3d(this%interp%dMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,3))
    call self_blas_matrixop_dim3_3d(this%interp%dMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    f_p => null()

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu,this%N,3*this%nVar,this%nelem)

  endsubroutine MappedGradient_MappedScalar3D

  subroutine MappedDGGradient_MappedScalar3D(this,df)
    !! Calculates the gradient of a function using the weak form of the gradient
    !! and the average boundary state.
    !! This method will compute the average boundary state from the
    !! and  attributes of
    implicit none
    class(MappedScalar3D),intent(in) :: this
    type(c_ptr),intent(inout) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:,:)
    type(c_ptr) :: fc

    call ContravariantWeight_3D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! From Vector divergence
    call c_f_pointer(this%jas_gpu,f_p, &
                     [this%interp%N+1,this%interp%N+1,this%interp%N+1,this%nelem,3*this%nvar,3])

    fc = c_loc(f_p(1,1,1,1,1,1))
    call self_blas_matrixop_dim1_3d(this%interp%dgMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,2))
    call self_blas_matrixop_dim2_3d(this%interp%dgMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,1,3))
    call self_blas_matrixop_dim3_3d(this%interp%dgMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,3*this%nvar,this%nelem,this%blas_handle)

    f_p => null()

    ! Do the boundary terms
    call NormalWeight_3D_gpu(this%avgBoundary_gpu, &
                             this%geometry%nhat%boundary_gpu,this%geometry%nscale%boundary_gpu, &
                             this%boundarynormal_gpu, &
                             this%interp%N,this%nvar,this%nelem)

    call DG_BoundaryContribution_3D_gpu(this%interp%bmatrix_gpu,this%interp%qweights_gpu, &
                                        this%boundarynormal_gpu,df,this%interp%N,3*this%nvar,this%nelem)

    call JacobianWeight_3D_gpu(df,this%geometry%J%interior_gpu,this%N,3*this%nVar,this%nelem)

  endsubroutine MappedDGGradient_MappedScalar3D

endmodule SELF_MappedScalar_3D