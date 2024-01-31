! SELF_Blas.f90
!
! Copyright 2024 Fluid Numerics LLC
! Authors : Joseph Schoonover (joe@fluidnumerics.com)
!           Garrett Byrd      (garrett@fluidnumerics.com)
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
! This module seeks to provider generic implementations for the real32 (single, S) and 
! real64 (double, D) floating point operations in the BLAS library.
!
! Below is an outline of all operations provided by BLAS (specifically hipBLAS). Each operation
! is an interface that maps to specific functions that accept singles and doubles. E.g., `axpy` maps 
! to `Saxpy` and `Daxpy`, and applies the appropriate function. Any complex-specific operations are
! NOT implemented and are marked with an asterisk (*). Eventually, all 
! batch/strided/stridedbatched versions of each function will also be implemented.
!
!   |- Level 1      |- Level 2      |- Level 3
!   | |- amax       | |- gbmv       | |- gemm
!   | |- amin       | |- gemv       | |- herk*
!   | |- asu        | |- ger        | |- herkx*
!   | |- axpy       | |- hbmv*      | |- her2k*
!   | |- copy       | |- hemv*      | |- symm
!   | |- dot        | |- her*       | |- syrk
!   | |- nrm2       | |- her2*      | |- syr2k
!   | |- rot        | |- hpmv*      | |- geam
!   | |- rotg       | |- hpr*       | |- hemm*
!   | |- rotm       | |- hpr2*      | |- trmm
!   | |- rotmg      | |- sbmv       | |- trsm
!   | |- scal       | |- spmv       | |- trtri
!   | |- swap       | |- spr        | |- dgmm
!                   | |- spr2
!                   | |- symv
!                   | |- syr
!                   | |- syr2
!                   | |- tbmv
!                   | |- tbsv
!                   | |- tpmv
!                   | |- tpsv
!                   | |- trmv
!                   | |- trsv

module self_blas
    use hipfort_hipblas
    use self_hip
    use self_hip_enums
    use self_constants
    use iso_c_binding
    use iso_fortran_env

    implicit none

    interface axpy
        module procedure saxpy
        module procedure daxpy
    end interface axpy

    contains

    !!!!!!!!!!!
    ! LEVEL 1 !
    !!!!!!!!!!!

    ! axpy
    function saxpy(alpha, vector1, vector2) result(z)
        implicit none
        
        real(kind=real32), intent(in) :: alpha
        real(kind=real32), intent(in) :: vector1(:), vector2(:)

        type(c_ptr) :: handle
        integer :: n
        ! real(real32) :: alpha = alpha
        real(kind=real32), pointer, dimension(:) :: x
        integer :: incx = 1
        real(kind=real32), pointer, dimension(:) :: y
        integer :: incy = 1
        type(c_ptr) :: x_gpu
        type(c_ptr) :: y_gpu
        integer :: status
        real(kind=real32), pointer, dimension(:) :: z

        ! Check if vectors are the same length
        if (size(vector1) /= size(vector2)) then
            print*,"Error: saxpy input vectors are not the same length."
            print*,"       size(vector1): ", size(vector1)
            print*,"       size(vector2): ", size(vector2)
            stop 1
        endif

        n = size(vector1)

        allocate(x(n), y(n), z(n))
        call hipcheck(hipmalloc(x_gpu, sizeof(x)))
        call hipcheck(hipmalloc(y_gpu, sizeof(y)))
        call hipblasCheck(hipblasCreate(handle))

        x = vector1
        y = vector2

        call hipcheck(hipmemcpy(x_gpu, c_loc(x), sizeof(x), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(y_gpu, c_loc(y), sizeof(y), hipmemcpyhosttodevice))

        status = hipblasSaxpy_(handle, n, alpha, x_gpu, incx, y_gpu, incy)
        call hipcheck(hipmemcpy(c_loc(y), y_gpu, sizeof(y), hipmemcpydevicetohost))

        z = y

        deallocate(x, y)
        call hipcheck(hipfree(x_gpu))
        call hipcheck(hipfree(y_gpu))
        call hipblasCheck(hipblasDestroy(handle))
    end function saxpy

    function daxpy(alpha, vector1, vector2) result(z)
        implicit none

        real(kind=real64), intent(in) :: alpha
        real(kind=real64), intent(in) :: vector1(:), vector2(:)

        type(c_ptr) :: handle
        integer :: n
        ! real(real64) :: alpha = alpha
        real(kind=real64), pointer, dimension(:) :: x
        integer :: incx = 1
        real(kind=real64), pointer, dimension(:) :: y
        integer :: incy = 1
        type(c_ptr) :: x_gpu
        type(c_ptr) :: y_gpu
        integer :: status
        real(kind=real64), pointer, dimension(:) :: z

        ! Check if vectors are the same length
        if (size(vector1) /= size(vector2)) then
            print*,"Error: saxpy input vectors are not the same length."
            print*,"       size(vector1): ", size(vector1)
            print*,"       size(vector2): ", size(vector2)
            stop 1
        endif

        n = size(vector1)

        allocate(x(n), y(n), z(n))
        call hipcheck(hipmalloc(x_gpu, sizeof(x)))
        call hipcheck(hipmalloc(y_gpu, sizeof(y)))
        call hipblasCheck(hipblasCreate(handle))

        x = vector1
        y = vector2

        call hipcheck(hipmemcpy(x_gpu, c_loc(x), sizeof(x), hipmemcpyhosttodevice))
        call hipcheck(hipmemcpy(y_gpu, c_loc(y), sizeof(y), hipmemcpyhosttodevice))

        status = hipblasDaxpy_(handle, n, alpha, x_gpu, incx, y_gpu, incy)
        
        call hipcheck(hipmemcpy(c_loc(y), y_gpu, sizeof(y), hipmemcpydevicetohost))

        z = y

        deallocate(x, y)
        call hipcheck(hipfree(x_gpu))
        call hipcheck(hipfree(y_gpu))
        call hipblasCheck(hipblasDestroy(handle))
    end function daxpy

    !!!!!!!!!!!
    ! LEVEL 2 !
    !!!!!!!!!!!

    !!!!!!!!!!!
    ! LEVEL 3 !
    !!!!!!!!!!!

end module self_blas