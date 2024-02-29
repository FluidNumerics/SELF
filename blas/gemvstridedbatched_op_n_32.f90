program test_program
    implicit none

    integer :: exit_code
    exit_code = test_function()
    print*,exit_code
    stop exit_code

    contains
    integer function test_function() result(r)
        use self_blas
        implicit none

        ! general
        logical :: vectors_equal = .true.
        integer(c_int) :: i
        ! real32
        real(real32) :: alpha = 1.0_real32
        real(real32) :: beta = 1.0_real32
        real(real32), pointer, dimension(:,:) :: A
        real(real32), pointer, dimension(:) :: x
        real(real32), pointer, dimension(:) :: y
        real(real32), pointer, dimension(:) :: expected
        real(real32),parameter :: tolerance = 10.0_real32**(-7)
        integer(c_int),parameter :: m = 10
        integer(c_int),parameter :: n = 20
        integer(c_int),parameter :: lda = m
        ! strided batched stuff
        integer(c_int64_t) :: stride_A = 0 ! m*n
        integer(c_int64_t) :: stride_x = n
        integer(c_int64_t) :: stride_y = m
        integer(c_int) :: batchCount = 1000
        integer(c_int) :: operation = HIPBLAS_OP_N

        ! gemv test for real32
        allocate(A(m,n), x(m*batchCount), y(n*batchCount), expected(m))

        A = 1.0_real32
        x = 1.0_real32
        y = 0.0_real32
        expected = dble(n)
        
        call sgemvstridedbatched_op(alpha, A, stride_A, x, stride_x, beta, y, stride_y, batchCount, operation)
        
        ! print*,y(1:m)

        do i=1,m
            if(abs(expected(i) - y(i)) > tolerance) then
                vectors_equal = .false.
                exit
            end if
        end do
        
        if(vectors_equal)then
            r = 0
        else
            r = 1
        endif

        deallocate(A, x, y, expected)

    end function test_function

end program test_program