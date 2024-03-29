program test_program
    implicit none

    integer :: exit_code
    integer :: i
    
    ! passing argument stuff
    integer :: arg_count, int_arg1, ierr
    character(len=100) :: arg1
    call getarg(1, arg1)

    read(arg1, *, iostat=ierr) int_arg1
    if (ierr /= 0) then
        print *, "Error: Argument must be an integer"
        exit_code = 1
        stop exit_code
    end if

    ! main loop
    do i=1,100
        exit_code = test_function(int_arg1)
    end do

    ! print exit code
    print*,exit_code
    stop exit_code

    contains
    integer function test_function(n) result(r)
        use self_blas
        implicit none

        ! general
        integer(c_int), intent(in) :: n
        logical :: vectors_equal = .true.
        integer(c_int) :: j

        real(real32) :: alpha = 1.0_real32
        real(real32), pointer, dimension(:) :: x, y
        real(real32), pointer, dimension(:) :: expected
        real(real32),parameter :: tolerance = 10.0_real32**(-7)
        
        allocate(x(n), y(n), expected(n))

        x = 1.0_real32
        y = 0.0_real32
        expected = 1.0_real32

        call saxpy(alpha, x, y)

        do j=1,n
            if(abs(expected(j) - x(j)) > tolerance) then
                vectors_equal = .false.
                exit
            end if
        end do

        if(vectors_equal)then
            r = 0
        else
            r = 1
        endif

        deallocate(x, y, expected)

    end function test_function
end program test_program