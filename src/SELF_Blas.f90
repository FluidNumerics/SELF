! SELF_Blas.f90
!
! Copyright 2024 Fluid Numerics LLC
! Authors : Joseph Schoonover (joe@fluidnumerics.com)
!           Garrett Byrd      (garrett@fluidnumerics.com)
!

module self_blas
    use hipfort_hipblas
    use self_hip
    use self_hip_enums
    use self_constants
    use iso_c_binding
    use iso_fortran_env

    implicit none

    ! instantiate handle

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
        integer :: m
        integer :: n
        integer :: lda
        integer :: incx = 1
        integer :: incy = 1
        type(c_ptr) :: A_gpu
        type(c_ptr) :: x_gpu
        type(c_ptr) :: y_gpu
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
        integer :: m
        integer :: n
        integer :: lda
        integer :: incx = 1
        integer :: incy = 1
        type(c_ptr) :: A_gpu
        type(c_ptr) :: x_gpu
        type(c_ptr) :: y_gpu
        integer :: status

        ! check if all dimensions are fine
        m = stridey
        n = stridex

        call hipcheck(hipmalloc(A_gpu, sizeof(A)))
        call hipcheck(hipmalloc(x_gpu, sizeof(x)))
        call hipcheck(hipmalloc(y_gpu, sizeof(y)))
        call hipblasCheck(hipblasCreate(handle))

        call hipcheck(hipmemcpy(A_gpu, c_loc(A), sizeof(A), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(x_gpu, c_loc(x), sizeof(x), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(y_gpu, c_loc(y), sizeof(y), hipmemcpyhosttodevice))

        ! HIPBLAS_OP_N = no transpose
        ! HIPBLAS_OP_T = transpose
        ! come back to this
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
        integer :: m
        integer :: n
        integer :: lda
        integer :: incx = 1
        integer :: incy = 1
        type(c_ptr) :: A_gpu
        type(c_ptr) :: x_gpu
        type(c_ptr) :: y_gpu
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
        integer :: m
        integer :: n
        integer :: lda
        integer :: incx = 1
        integer :: incy = 1
        type(c_ptr) :: A_gpu
        type(c_ptr) :: x_gpu
        type(c_ptr) :: y_gpu
        integer :: status

        ! check if all dimensions are fine
        m = stridey
        n = stridex

        call hipcheck(hipmalloc(A_gpu, sizeof(A)))
        call hipcheck(hipmalloc(x_gpu, sizeof(x)))
        call hipcheck(hipmalloc(y_gpu, sizeof(y)))
        call hipblasCheck(hipblasCreate(handle))

        call hipcheck(hipmemcpy(A_gpu, c_loc(A), sizeof(A), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(x_gpu, c_loc(x), sizeof(x), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(y_gpu, c_loc(y), sizeof(y), hipmemcpyhosttodevice))

        ! HIPBLAS_OP_N = no transpose
        ! HIPBLAS_OP_T = transpose
        ! come back to this
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
    subroutine sgemm(scalar1, A, matrix2, scalar2, matrix3)
        ! prefer scalarX, matrixX, etc, or alpha, beta, A, B, C?
        implicit none

        real(kind=real32), intent(in) :: scalar1
        real(kind=real32), pointer, intent(in) :: A(:,:)
        real(kind=real32), pointer, intent(in) :: matrix2(:,:)
        real(kind=real32), intent(in) :: scalar2
        real(kind=real32), pointer, intent(inout) :: matrix3(:,:)

        type(c_ptr) :: handle
        integer :: operation1 = HIPBLAS_OP_T
        integer :: operation2 = HIPBLAS_OP_T
        ! prefer integer :: m, n, k ?
        integer :: m
        integer :: n
        integer :: k
        integer :: lda
        integer :: ldb
        integer :: ldc
        type(c_ptr) :: A_gpu
        type(c_ptr) :: B_gpu
        type(c_ptr) :: C_gpu
        integer :: status

        ! check if dimensions are fine

        m = size(matrix3, 1)
        n = size(matrix3, 2)
        k = size(A, 1) ! this changes for non-transposed matrix A
        lda = k
        ldb = n
        ldc = m

        call hipcheck(hipmalloc(A_gpu, sizeof(A)))
        call hipcheck(hipmalloc(B_gpu, sizeof(matrix2)))
        call hipcheck(hipmalloc(C_gpu, sizeof(matrix3)))        
        call hipblasCheck(hipblasCreate(handle))

        call hipcheck(hipmemcpy(A_gpu, c_loc(A), sizeof(A), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(B_gpu, c_loc(matrix2), sizeof(matrix2), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(C_gpu, c_loc(matrix3), sizeof(matrix3), hipmemcpyhosttodevice))

        status = hipblassgemm(handle, operation1, operation2, m, n, k, scalar1, A_gpu, lda, B_gpu, ldb, scalar2, C_gpu, ldc)

        call hipcheck(hipmemcpy(c_loc(matrix3), C_gpu, sizeof(matrix3), hipmemcpydevicetohost))

        call hipcheck(hipfree(A_gpu))
        call hipcheck(hipfree(B_gpu))
        call hipcheck(hipfree(C_gpu))
        call hipblasCheck(hipblasDestroy(handle))

    end subroutine sgemm
    subroutine sgemm_op(alpha, A, B, beta, C, operation_A)
        implicit none

        real(kind=real32), intent(in) :: alpha
        real(kind=real32), pointer, intent(in) :: A(:,:)
        real(kind=real32), pointer, intent(in) :: B(:,:)
        real(kind=real32), intent(in) :: beta
        real(kind=real32), pointer, intent(inout) :: C(:,:)
        integer, intent(in) :: operation_A

        type(c_ptr) :: handle
        integer :: operation_B = HIPBLAS_OP_N
        integer :: m, n, k
        integer :: lda, ldb, ldc
        type(c_ptr) :: A_gpu, B_gpu, C_gpu
        integer :: status

        ! check if dimensions are fine

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

        status = hipblassgemm(handle, operation_A, operation_B, m, n, k, alpha, A_gpu, lda, B_gpu, ldb, beta, C_gpu, ldc)

        call hipcheck(hipmemcpy(c_loc(C), C_gpu, sizeof(C), hipmemcpydevicetohost))

        call hipcheck(hipfree(A_gpu))
        call hipcheck(hipfree(B_gpu))
        call hipcheck(hipfree(C_gpu))
        call hipblasCheck(hipblasDestroy(handle))

    end subroutine sgemm_op


    subroutine dgemm(scalar1, A, matrix2, scalar2, matrix3)
        implicit none

        real(kind=real64), intent(in) :: scalar1
        real(kind=real64), pointer, intent(in) :: A(:,:)
        real(kind=real64), pointer, intent(in) :: matrix2(:,:)
        real(kind=real64), intent(in) :: scalar2
        real(kind=real64), pointer, intent(inout) :: matrix3(:,:)

        type(c_ptr) :: handle
        integer :: operation1 = HIPBLAS_OP_T
        integer :: operation2 = HIPBLAS_OP_T
        ! prefer integer :: m, n, k ?
        integer :: m
        integer :: n
        integer :: k
        integer :: lda
        integer :: ldb
        integer :: ldc
        integer :: batchCount
        type(c_ptr) :: A_gpu
        type(c_ptr) :: B_gpu
        type(c_ptr) :: C_gpu
        integer :: status

        ! check if dimensions are fine

        m = size(matrix3, 1)
        n = size(matrix3, 2)
        k = size(A, 1) ! this changes for non-transposed matrix A
        lda = k
        ldb = n
        ldc = m

        call hipcheck(hipmalloc(A_gpu, sizeof(A)))
        call hipcheck(hipmalloc(B_gpu, sizeof(matrix2)))
        call hipcheck(hipmalloc(C_gpu, sizeof(matrix3)))        
        call hipblasCheck(hipblasCreate(handle))

        call hipcheck(hipmemcpy(A_gpu, c_loc(A), sizeof(A), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(B_gpu, c_loc(matrix2), sizeof(matrix2), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(C_gpu, c_loc(matrix3), sizeof(matrix3), hipmemcpyhosttodevice))

        status = hipblasdgemm(handle, operation1, operation2, m, n, k, scalar1, A_gpu, lda, B_gpu, ldb, scalar2, C_gpu, ldc)

        call hipcheck(hipmemcpy(c_loc(matrix3), C_gpu, sizeof(matrix3), hipmemcpydevicetohost))

        call hipcheck(hipfree(A_gpu))
        call hipcheck(hipfree(B_gpu))
        call hipcheck(hipfree(C_gpu))
        call hipblasCheck(hipblasDestroy(handle))

    end subroutine dgemm
    subroutine dgemm_op(alpha, A, B, beta, C, operation_A)
        implicit none

        real(kind=real64), intent(in) :: alpha
        real(kind=real64), pointer, intent(in) :: A(:,:)
        real(kind=real64), pointer, intent(in) :: B(:,:)
        real(kind=real64), intent(in) :: beta
        real(kind=real64), pointer, intent(inout) :: C(:,:)
        integer, intent(in) :: operation_A

        type(c_ptr) :: handle
        integer :: operation_B = HIPBLAS_OP_N
        integer :: m, n, k
        integer :: lda, ldb, ldc
        type(c_ptr) :: A_gpu, B_gpu, C_gpu
        integer :: status

        ! check if dimensions are fine

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

        status = hipblasdgemm(handle, operation_A, operation_B, m, n, k, alpha, A_gpu, lda, B_gpu, ldb, beta, C_gpu, ldc)

        call hipcheck(hipmemcpy(c_loc(C), C_gpu, sizeof(C), hipmemcpydevicetohost))

        call hipcheck(hipfree(A_gpu))
        call hipcheck(hipfree(B_gpu))
        call hipcheck(hipfree(C_gpu))
        call hipblasCheck(hipblasDestroy(handle))

    end subroutine dgemm_op


end module self_blas