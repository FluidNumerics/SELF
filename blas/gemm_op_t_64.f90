program test_program
    implicit none

    integer :: exit_code
    integer :: i
    
    ! passing argument stuff
    integer :: arg_count, int_arg1, int_arg2, ierr
    character(len=100) :: arg1, arg2
    call getarg(1, arg1)
    call getarg(2, arg2)

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

    ! main loop
    do i=1,100
        exit_code = test_function(int_arg1,int_arg2)
    end do

    ! print exit code
    print*,exit_code
    stop exit_code

    contains
    integer function test_function(m,n) result(r)
        use self_blas
        implicit none

        ! general
        integer(c_int), intent(in) :: m, n
        logical :: matrices_equal = .true.
        integer(c_int) :: i,j
        ! real64
        real(real64) :: alpha = 1.0_real64
        real(real64) :: beta = 1.0_real64
        real(real64), pointer, dimension(:,:) :: A, B, C
        real(real64), pointer, dimension(:,:) :: expected
        real(real64),parameter :: tolerance = 10.0_real64**(-7)
        integer(c_int) :: lda, ldb, ldc
        integer(c_int) :: operation = HIPBLAS_OP_T

        lda = m
        ldb = m
        ldc = m

        allocate(A(m,m), B(m,n), C(m,n), expected(m,n))

        A = 1.0_real64
        B = 1.0_real64
        C = 0.0_real64
        expected = dble(m)

        call dgemm_op(alpha, A, B, beta, C, operation)

        do i=1,m
            do j=1,n
                if(abs(expected(i,j) - C(i,j)) > tolerance) then
                    matrices_equal = .false.
                    exit
                end if
            end do
        end do
        
        if(matrices_equal)then
            r = 0
        else
            r = 1
        endif

        deallocate(A, B, C, expected)

    end function test_function
end program test_program