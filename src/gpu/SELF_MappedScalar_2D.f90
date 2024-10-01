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

module SELF_MappedScalar_2D

  use SELF_MappedScalar_2D_t
  use SELF_GPU
  use SELF_GPUInterfaces
  use iso_c_binding

  implicit none

  type,extends(MappedScalar2D_t),public :: MappedScalar2D

    type(c_ptr) :: jas_gpu ! jacobian weighted scalar for gradient calculation

  contains
    procedure,public :: Init => Init_MappedScalar2D
    procedure,public :: Free => Free_MappedScalar2D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar2D

    procedure,public :: SideExchange => SideExchange_MappedScalar2D
    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar2D

    generic,public :: MappedGradient => MappedGradient_MappedScalar2D
    procedure,private :: MappedGradient_MappedScalar2D

    generic,public :: MappedDGGradient => MappedDGGradient_MappedScalar2D
    procedure,private :: MappedDGGradient_MappedScalar2D

  endtype MappedScalar2D

  interface
    subroutine ContravariantWeight_2D_gpu(f,dsdx,jaf,N,nvar,nel) &
      bind(c,name="ContravariantWeight_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,dsdx,jaf
      integer(c_int),value :: N,nvar,nel
    endsubroutine ContravariantWeight_2D_gpu
  endinterface

  interface
    subroutine NormalWeight_2D_gpu(fb,nhat,nscale,fbn,N,nvar,nel) &
      bind(c,name="NormalWeight_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: fb,nhat,nscale,fbn
      integer(c_int),value :: N,nvar,nel
    endsubroutine NormalWeight_2D_gpu
  endinterface

contains

  subroutine Init_MappedScalar2D(this,interp,nVar,nElem)
    implicit none
    class(MappedScalar2D),intent(out) :: this
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

    allocate(this%interior(1:interp%N+1,interp%N+1,nelem,nvar), &
             this%boundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%extBoundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%avgBoundary(1:interp%N+1,1:4,1:nelem,1:nvar), &
             this%boundarynormal(1:interp%N+1,1:4,1:nelem,1:2*nvar))

    allocate(this%meta(1:nVar))
    allocate(this%eqn(1:nVar))

    this%interior = 0.0_prec
    this%boundary = 0.0_prec
    this%extBoundary = 0.0_prec
    this%avgBoundary = 0.0_prec
    this%boundarynormal = 0.0_prec

    call gpuCheck(hipMalloc(this%interior_gpu,sizeof(this%interior)))
    call gpuCheck(hipMalloc(this%boundary_gpu,sizeof(this%boundary)))
    call gpuCheck(hipMalloc(this%extBoundary_gpu,sizeof(this%extBoundary)))
    call gpuCheck(hipMalloc(this%avgBoundary_gpu,sizeof(this%avgBoundary)))
    call gpuCheck(hipMalloc(this%boundarynormal_gpu,sizeof(this%boundarynormal)))
    workSize = (interp%N+1)*(interp%M+1)*nelem*nvar*prec
    call gpuCheck(hipMalloc(this%interpWork,workSize))
    workSize = (interp%N+1)*(interp%N+1)*nelem*nvar*4*prec
    call gpuCheck(hipMalloc(this%jas_gpu,workSize))

    call hipblasCheck(hipblasCreate(this%blas_handle))

  endsubroutine Init_MappedScalar2D

  subroutine Free_MappedScalar2D(this)
    implicit none
    class(MappedScalar2D),intent(inout) :: this

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
    call gpuCheck(hipFree(this%interpWork))
    call gpuCheck(hipFree(this%jas_gpu))
    call hipblasCheck(hipblasDestroy(this%blas_handle))

  endsubroutine Free_MappedScalar2D

  subroutine SetInteriorFromEquation_MappedScalar2D(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: x
    real(prec) :: y

    do iVar = 1,this%nVar
      do iEl = 1,this%nElem
        do j = 1,this%interp%N+1
          do i = 1,this%interp%N+1

            ! Get the mesh positions
            x = geometry%x%interior(i,j,iEl,1,1)
            y = geometry%x%interior(i,j,iEl,1,2)

            this%interior(i,j,iEl,iVar) = &
              this%eqn(iVar)%Evaluate((/x,y,0.0_prec,time/))

          enddo
        enddo
      enddo
    enddo

    call gpuCheck(hipMemcpy(this%interior_gpu,c_loc(this%interior),sizeof(this%interior),hipMemcpyHostToDevice))

  endsubroutine SetInteriorFromEquation_MappedScalar2D

  subroutine MPIExchangeAsync_MappedScalar2D(this,mesh)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: e1,s1,e2,s2,ivar
    integer :: globalSideId,r2,tag
    integer :: iError
    integer :: msgCount
    real(prec),pointer :: boundary(:,:,:,:)
    real(prec),pointer :: extboundary(:,:,:,:)

    msgCount = 0
    call c_f_pointer(this%boundary_gpu,boundary,[this%interp%N+1,4,this%nelem,this%nvar])
    call c_f_pointer(this%extboundary_gpu,extboundary,[this%interp%N+1,4,this%nelem,this%nvar])

    do ivar = 1,this%nvar
      do e1 = 1,this%nElem
        do s1 = 1,4

          e2 = mesh%sideInfo(3,s1,e1) ! Neighbor Element
          if(e2 > 0) then
            r2 = mesh%decomp%elemToRank(e2) ! Neighbor Rank

            if(r2 /= mesh%decomp%rankId) then

              s2 = mesh%sideInfo(4,s1,e1)/10
              globalSideId = abs(mesh%sideInfo(2,s1,e1))
              ! create unique tag for each side and each variable
              tag = globalsideid+mesh%nUniqueSides*(ivar-1)

              msgCount = msgCount+1
              call MPI_IRECV(extBoundary(:,s1,e1,ivar), &
                             (this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             r2,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)

              msgCount = msgCount+1
              call MPI_ISEND(boundary(:,s1,e1,ivar), &
                             (this%interp%N+1), &
                             mesh%decomp%mpiPrec, &
                             r2,tag, &
                             mesh%decomp%mpiComm, &
                             mesh%decomp%requests(msgCount),iError)
            endif
          endif

        enddo
      enddo
    enddo

    mesh%decomp%msgCount = msgCount

  endsubroutine MPIExchangeAsync_MappedScalar2D

  subroutine SideExchange_MappedScalar2D(this,mesh)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(Mesh2D),intent(inout) :: mesh
    ! Local
    integer :: offset

    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)
    if(mesh%decomp%mpiEnabled) then
      call this%MPIExchangeAsync(mesh)
    endif

    ! Do the side exchange internal to this mpi process
    call SideExchange_2D_gpu(this%extboundary_gpu, &
                             this%boundary_gpu,mesh%sideinfo_gpu,mesh%decomp%elemToRank_gpu, &
                             mesh%decomp%rankid,offset,this%interp%N,this%nvar,this%nelem)

    if(mesh%decomp%mpiEnabled) then
      call mesh%decomp%FinalizeMPIExchangeAsync()
      ! Apply side flips for data exchanged with MPI
      call ApplyFlip_2D_gpu(this%extboundary_gpu,mesh%sideInfo_gpu, &
                            mesh%decomp%elemToRank_gpu,mesh%decomp%rankId, &
                            offset,this%interp%N,this%nVar,this%nElem)
    endif

  endsubroutine SideExchange_MappedScalar2D

  subroutine MappedGradient_MappedScalar2D(this,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(c_ptr),intent(out) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:)
    type(c_ptr) :: fc

    call ContravariantWeight_2D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! From Vector divergence
    call c_f_pointer(this%jas_gpu,f_p,[this%interp%N+1,this%interp%N+1,this%nelem,2*this%nvar,2])

    fc = c_loc(f_p(1,1,1,1,1))
    call self_blas_matrixop_dim1_2d(this%interp%dMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,2))
    call self_blas_matrixop_dim2_2d(this%interp%dMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    f_p => null()

    call JacobianWeight_2D_gpu(df,this%geometry%J%interior_gpu,this%N,2*this%nVar,this%nelem)

  endsubroutine MappedGradient_MappedScalar2D

  subroutine MappedDGGradient_MappedScalar2D(this,df)
    !!
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(c_ptr),intent(inout) :: df
    ! Local
    real(prec),pointer :: f_p(:,:,:,:,:)
    type(c_ptr) :: fc

    call ContravariantWeight_2D_gpu(this%interior_gpu, &
                                    this%geometry%dsdx%interior_gpu,this%jas_gpu, &
                                    this%interp%N,this%nvar,this%nelem)

    ! From Vector divergence
    call c_f_pointer(this%jas_gpu,f_p,[this%interp%N+1,this%interp%N+1,this%nelem,2*this%nvar,2])

    fc = c_loc(f_p(1,1,1,1,1))
    call self_blas_matrixop_dim1_2d(this%interp%dgMatrix_gpu,fc,df, &
                                    this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    fc = c_loc(f_p(1,1,1,1,2))
    call self_blas_matrixop_dim2_2d(this%interp%dgMatrix_gpu,fc,df, &
                                    1.0_c_prec,this%interp%N,this%interp%N,2*this%nvar,this%nelem,this%blas_handle)

    f_p => null()

    ! Do the boundary terms
    call NormalWeight_2D_gpu(this%avgBoundary_gpu, &
                             this%geometry%nhat%boundary_gpu,this%geometry%nscale%boundary_gpu, &
                             this%boundarynormal_gpu, &
                             this%interp%N,this%nvar,this%nelem)

    call DG_BoundaryContribution_2D_gpu(this%interp%bmatrix_gpu,this%interp%qweights_gpu, &
                                        this%boundarynormal_gpu,df,this%interp%N,2*this%nvar,this%nelem)

    call JacobianWeight_2D_gpu(df,this%geometry%J%interior_gpu,this%N,2*this%nVar,this%nelem)

  endsubroutine MappedDGGradient_MappedScalar2D

endmodule SELF_MappedScalar_2D
