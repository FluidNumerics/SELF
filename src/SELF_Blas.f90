! SELF_Blas.f90
!
! Copyright 2024 Fluid Numerics LLC
! Authors : Joseph Schoonover (joe@fluidnumerics.com)
!           Garrett Byrd      (garrett@fluidnumerics.com)
!
!   This module provides wrappers for hipBLAS built specifically for SELF. Only
!   subroutines necessary for SELF are included: (gemvstridedbatched, gemm).
!   These wrappers are intended to reduce the complexity of implementation as
!   the respective subroutines are called throughout the library. I.e., instead
!   of needing to specify all arguments needed for BLAS subroutines (e.g., 
!   `lda`, `incx`, `handle`, `status`, etc.), the user only needs to pass in
!   the relevant scalars, vectors*, matrices*, and batchcount. (*by reference)
!   
!   Interfaces are provided so that one does not need to distinguish between
!   single and double precision in implementation. Note: the user cannot mix
!   and match precision, but the provided interfaces still reduce complexity in
!   implementation.
!
!   Also provided are `SUBROUTINE_op` subroutines, which allow the user to
!   specify if a particular matrix is to be passed with `HIPBLAS_OP_N` or 
!   `HIPBLAS_OP_T`. 
!

module self_blas
    use hipfort_hipblas
    use hipfort
    use hipfort_check
    use iso_c_binding
    use iso_fortran_env

    implicit none

    interface gemvstridedbatched
        module procedure sgemvstridedbatched
        module procedure sgemvstridedbatched_op
        module procedure dgemvstridedbatched
        module procedure dgemvstridedbatched_op
    end interface gemvstridedbatched
    
    interface gemm
        module procedure sgemm
        module procedure sgemm_op
        module procedure dgemm
        module procedure dgemm_op
    end interface gemm

    contains

    !!!!!!!!!!!
    ! LEVEL 2 !
    !!!!!!!!!!!

    ! gemvstridedbatched
    subroutine sgemvstridedbatched(alpha, A, strideA, x, stridex, beta, y, stridey, batchCount)
        implicit none

        real(kind=real32), intent(in) :: alpha
        real(kind=real32), pointer, intent(in) :: A(:,:)
        real(kind=real32), pointer, intent(in) :: x(:)
        real(kind=real32), intent(in) :: beta
        real(kind=real32), pointer, intent(inout) :: y(:)
        integer(c_int64_t), intent(in) :: strideA, stridex, stridey
        integer(c_int), intent(in) :: batchCount

        type(c_ptr) :: handle
        integer :: operation = HIPBLAS_OP_T
        integer :: m, n
        integer :: lda
        integer :: incx = 1
        integer :: incy = 1
        type(c_ptr) :: A_gpu, x_gpu, y_gpu
        integer :: status

        call sgemvstridedbatched_op(alpha, A, strideA, x, stridex, beta, y, stridey, batchCount, operation)

    end subroutine sgemvstridedbatched

    subroutine sgemvstridedbatched_op(alpha, A, strideA, x, stridex, beta, y, stridey, batchCount, operation)
        implicit none

        real(kind=real32), intent(in) :: alpha
        real(kind=real32), pointer, intent(in) :: A(:,:)
        real(kind=real32), pointer, intent(in) :: x(:)
        real(kind=real32), intent(in) :: beta
        real(kind=real32), pointer, intent(inout) :: y(:)
        integer(c_int64_t), intent(in) :: strideA, stridex, stridey
        integer(c_int), intent(in) :: batchCount
        integer(kind(HIPBLAS_OP_N)), intent(in) :: operation

        type(c_ptr) :: handle
        integer :: m, n
        integer :: lda
        integer :: incx = 1
        integer :: incy = 1
        type(c_ptr) :: A_gpu, x_gpu, y_gpu
        integer :: status

        m = stridey
        n = stridex

        call hipcheck(hipmalloc(A_gpu, sizeof(A)))
        call hipcheck(hipmalloc(x_gpu, sizeof(x)))
        call hipcheck(hipmalloc(y_gpu, sizeof(y)))
        call hipblasCheck(hipblasCreate(handle))

        call hipcheck(hipmemcpy(A_gpu, c_loc(A), sizeof(A), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(x_gpu, c_loc(x), sizeof(x), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(y_gpu, c_loc(y), sizeof(y), hipmemcpyhosttodevice))

        if( operation == HIPBLAS_OP_N) then
            status = hipblasSgemvStridedBatched(handle, operation, m, n, alpha, A_gpu, m, strideA, x_gpu, incx, stridex, beta, y_gpu, incy, stridey, batchcount)
        else if( operation == HIPBLAS_OP_T) then 
            status = hipblasSgemvStridedBatched(handle, operation, n, m, alpha, A_gpu, n, strideA, x_gpu, incx, stridex, beta, y_gpu, incy, stridey, batchcount)
        ! else
        !     throw an error
        end if

        call hipcheck(hipmemcpy(c_loc(y), y_gpu, sizeof(y), hipmemcpydevicetohost))

        call hipcheck(hipfree(A_gpu))
        call hipcheck(hipfree(x_gpu))
        call hipcheck(hipfree(y_gpu))
        call hipblasCheck(hipblasDestroy(handle))

    end subroutine sgemvstridedbatched_op

    subroutine dgemvstridedbatched(alpha, A, strideA, x, stridex, beta, y, stridey, batchCount)
        implicit none

        real(kind=real64), intent(in) :: alpha
        real(kind=real64), pointer, intent(in) :: A(:,:)
        real(kind=real64), pointer, intent(in) :: x(:)
        real(kind=real64), intent(in) :: beta
        real(kind=real64), pointer, intent(inout) :: y(:)
        integer(c_int64_t), intent(in) :: strideA, stridex, stridey
        integer(c_int), intent(in) :: batchCount

        type(c_ptr) :: handle
        integer :: operation = HIPBLAS_OP_T
        integer :: m, n
        integer :: lda
        integer :: incx = 1
        integer :: incy = 1
        type(c_ptr) :: A_gpu, x_gpu, y_gpu
        integer :: status

        call dgemvstridedbatched_op(alpha, A, strideA, x, stridex, beta, y, stridey, batchCount, operation)

    end subroutine dgemvstridedbatched

    subroutine dgemvstridedbatched_op(alpha, A, strideA, x, stridex, beta, y, stridey, batchCount, operation)
        implicit none

        real(kind=real64), intent(in) :: alpha
        real(kind=real64), pointer, intent(in) :: A(:,:)
        real(kind=real64), pointer, intent(in) :: x(:)
        real(kind=real64), intent(in) :: beta
        real(kind=real64), pointer, intent(inout) :: y(:)
        integer(c_int64_t), intent(in) :: strideA, stridex, stridey
        integer(c_int), intent(in) :: batchCount
        integer, intent(in) :: operation

        type(c_ptr) :: handle
        integer :: m, n
        integer :: lda
        integer :: incx = 1
        integer :: incy = 1
        type(c_ptr) :: A_gpu, x_gpu, y_gpu
        integer :: status

        m = stridey
        n = stridex

        call hipcheck(hipmalloc(A_gpu, sizeof(A)))
        call hipcheck(hipmalloc(x_gpu, sizeof(x)))
        call hipcheck(hipmalloc(y_gpu, sizeof(y)))
        call hipblasCheck(hipblasCreate(handle))

        call hipcheck(hipmemcpy(A_gpu, c_loc(A), sizeof(A), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(x_gpu, c_loc(x), sizeof(x), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(y_gpu, c_loc(y), sizeof(y), hipmemcpyhosttodevice))

        if( operation == HIPBLAS_OP_N) then
            status = hipblasDgemvStridedBatched(handle, operation, m, n, alpha, A_gpu, m, strideA, x_gpu, incx, stridex, beta, y_gpu, incy, stridey, batchcount)
        else if( operation == HIPBLAS_OP_T) then 
            status = hipblasDgemvStridedBatched(handle, operation, n, m, alpha, A_gpu, n, strideA, x_gpu, incx, stridex, beta, y_gpu, incy, stridey, batchcount)
        ! else
        !     throw an error
        end if

        call hipcheck(hipmemcpy(c_loc(y), y_gpu, sizeof(y), hipmemcpydevicetohost))

        call hipcheck(hipfree(A_gpu))
        call hipcheck(hipfree(x_gpu))
        call hipcheck(hipfree(y_gpu))
        call hipblasCheck(hipblasDestroy(handle))

    end subroutine dgemvstridedbatched_op

    !!!!!!!!!!!
    ! LEVEL 3 !
    !!!!!!!!!!!

    ! gemm
    subroutine sgemm(alpha, A, B, beta, C)
        implicit none

        real(kind=real32), intent(in) :: alpha
        real(kind=real32), pointer, intent(in) :: A(:,:)
        real(kind=real32), pointer, intent(in) :: B(:,:)
        real(kind=real32), intent(in) :: beta
        real(kind=real32), pointer, intent(inout) :: C(:,:)

        type(c_ptr) :: handle
        integer :: operationA = HIPBLAS_OP_T
        integer :: m, n, k
        integer :: lda, ldb, ldc
        type(c_ptr) :: A_gpu, B_gpu, C_gpu
        integer :: status

        call sgemm_op(alpha, A, B, beta, C, operationA)

    end subroutine sgemm

    subroutine sgemm_op(alpha, A, B, beta, C, operationA)
        implicit none

        real(kind=real32), intent(in) :: alpha
        real(kind=real32), pointer, intent(in) :: A(:,:)
        real(kind=real32), pointer, intent(in) :: B(:,:)
        real(kind=real32), intent(in) :: beta
        real(kind=real32), pointer, intent(inout) :: C(:,:)
        integer, intent(in) :: operationA

        type(c_ptr) :: handle
        integer :: operationB = HIPBLAS_OP_N
        integer :: m, n, k
        integer :: lda, ldb, ldc
        type(c_ptr) :: A_gpu, B_gpu, C_gpu
        integer :: status

        ! Assign m,n,k
        ! A is mxm, so 
        m = size(A, 1)
        n = size(B, 2)
        k = m
        ! Assign lda, ldb, ldc
        lda = m
        ldb = m
        ldc = m

        call hipcheck(hipmalloc(A_gpu, sizeof(A)))
        call hipcheck(hipmalloc(B_gpu, sizeof(B)))
        call hipcheck(hipmalloc(C_gpu, sizeof(C)))        
        call hipblasCheck(hipblasCreate(handle))

        call hipcheck(hipmemcpy(A_gpu, c_loc(A), sizeof(A), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(B_gpu, c_loc(B), sizeof(B), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(C_gpu, c_loc(C), sizeof(C), hipmemcpyhosttodevice))

        status = hipblassgemm(handle, operationA, operationB, m, n, k, alpha, A_gpu, lda, B_gpu, ldb, beta, C_gpu, ldc)

        call hipcheck(hipmemcpy(c_loc(C), C_gpu, sizeof(C), hipmemcpydevicetohost))

        call hipcheck(hipfree(A_gpu))
        call hipcheck(hipfree(B_gpu))
        call hipcheck(hipfree(C_gpu))
        call hipblasCheck(hipblasDestroy(handle))

    end subroutine sgemm_op

    subroutine dgemm(alpha, A, B, beta, C)
        implicit none

        real(kind=real64), intent(in) :: alpha
        real(kind=real64), pointer, intent(in) :: A(:,:)
        real(kind=real64), pointer, intent(in) :: B(:,:)
        real(kind=real64), intent(in) :: beta
        real(kind=real64), pointer, intent(inout) :: C(:,:)

        type(c_ptr) :: handle
        integer :: operationA = HIPBLAS_OP_T
        integer :: m, n, k
        integer :: lda, ldb, ldc
        integer :: batchCount
        type(c_ptr) :: A_gpu, B_gpu, C_gpu
        integer :: status

        ! Assign m,n,k
        ! A is mxm, so
        m = size(A, 1)
        n = size(B, 2)
        k = m
        ! Assign lda, ldb, ldc
        lda = k
        ldb = n
        ldc = m

        call dgemm_op(alpha, A, B, beta, C, operationA)

    end subroutine dgemm

    subroutine dgemm_op(alpha, A, B, beta, C, operationA)
        implicit none

        real(kind=real64), intent(in) :: alpha
        real(kind=real64), pointer, intent(in) :: A(:,:)
        real(kind=real64), pointer, intent(in) :: B(:,:)
        real(kind=real64), intent(in) :: beta
        real(kind=real64), pointer, intent(inout) :: C(:,:)
        integer, intent(in) :: operationA

        type(c_ptr) :: handle
        integer :: operationB = HIPBLAS_OP_N
        integer :: m, n, k
        integer :: lda, ldb, ldc
        type(c_ptr) :: A_gpu, B_gpu, C_gpu
        integer :: status

        ! Assign m,n,k
        ! A is mxm, so 
        m = size(A, 1)
        n = size(B, 2)
        k = m
        ! Assign lda, ldb, ldc
        lda = m
        ldb = m
        ldc = m

        call hipcheck(hipmalloc(A_gpu, sizeof(A)))
        call hipcheck(hipmalloc(B_gpu, sizeof(B)))
        call hipcheck(hipmalloc(C_gpu, sizeof(C)))        
        call hipblasCheck(hipblasCreate(handle))

        call hipcheck(hipmemcpy(A_gpu, c_loc(A), sizeof(A), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(B_gpu, c_loc(B), sizeof(B), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(C_gpu, c_loc(C), sizeof(C), hipmemcpyhosttodevice))

        status = hipblasdgemm(handle, operationA, operationB, m, n, k, alpha, A_gpu, lda, B_gpu, ldb, beta, C_gpu, ldc)

        call hipcheck(hipmemcpy(c_loc(C), C_gpu, sizeof(C), hipmemcpydevicetohost))

        call hipcheck(hipfree(A_gpu))
        call hipcheck(hipfree(B_gpu))
        call hipcheck(hipfree(C_gpu))
        call hipblasCheck(hipblasDestroy(handle))

    end subroutine dgemm_op

end module self_blas