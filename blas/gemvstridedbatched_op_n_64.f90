program test_program
    implicit none

    integer :: exit_code
    integer :: i
    ! passing argument stuff
    integer :: arg_count, int_arg1, int_arg2,int_arg3, ierr
    character(len=100) :: arg1, arg2, arg3
    call getarg(1, arg1)
    call getarg(2, arg2)
    call getarg(3, arg3)

    read(arg1, *, iostat=ierr) int_arg1
    if (ierr /= 0) then
        print *, "Error: Argument must be an integer"
        exit_code = 1
        stop exit_code
    end if
    read(arg2, *, iostat=ierr) int_arg2
    if (ierr /= 0) then
        print *, "Error: Argument must be an integer"
        exit_code = 1
        stop exit_code
    end if
    read(arg3, *, iostat=ierr) int_arg3
    if (ierr /= 0) then
        print *, "Error: Argument must be an integer"
        exit_code = 1
        stop exit_code
    end if

    do i=1,100
        exit_code = test_function(int_arg1,int_arg2,int_arg3)
    end do
    print*,exit_code
    stop exit_code

    contains
    integer function test_function(m,n,batchCount) result(r)
        use self_blas
        implicit none

        ! general
        integer(c_int), intent(in):: m
        integer(c_int), intent(in):: n
        integer(c_int), intent(in) :: batchCount
        logical :: vectors_equal = .true.
        integer(c_int) :: i
        ! real64
        real(real64) :: alpha = 1.0_real64
        real(real64) :: beta = 1.0_real64
        real(real64), pointer, dimension(:,:) :: A
        real(real64), pointer, dimension(:) :: x
        real(real64), pointer, dimension(:) :: y
        real(real64), pointer, dimension(:) :: expected
        ! import prec
        real(real64),parameter :: tolerance = 10.0_real64**(-7)
        integer(c_int):: lda
        ! strided batched stuff
        integer(c_int64_t) :: stride_A = 0 ! m*n
        integer(c_int64_t) :: stride_x
        integer(c_int64_t) :: stride_y
        integer(c_int) :: operation = HIPBLAS_OP_N

        lda = m
        stride_x = n
        stride_y = m

        ! gemv test for real64
        allocate(A(m,n), x(n*batchCount), y(m*batchCount), expected(m))

        A = 1.0_real64
        x = 1.0_real64
        y = 0.0_real64
        expected = dble(n)
        
        call dgemvstridedbatched_op(alpha, A, stride_A, x, stride_x, beta, y, stride_y, batchCount, operation)
        
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
