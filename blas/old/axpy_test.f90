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

        logical :: vectors_equal = .true.
        integer :: i
        ! real32
        real(real32) :: alpha_input_32 = 2.0_real32
        real(real32), pointer, dimension(:) :: x_input_32
        real(real32), pointer, dimension(:) :: y_input_32
        real(real32), pointer, dimension(:) :: z_output_32
        real(real32), pointer, dimension(:) :: expected_output_32
        real(real32),parameter :: tolerance_32 = 10.0_real32**(-3)
        integer,parameter :: n_32 = 3
        ! real64
        real(real64) :: alpha_input_64 = 2.0_real64
        real(real64), pointer, dimension(:) :: x_input_64
        real(real64), pointer, dimension(:) :: y_input_64
        real(real64), pointer, dimension(:) :: z_output_64
        real(real64), pointer, dimension(:) :: expected_output_64
        real(real64),parameter :: tolerance_64 = 10.0_real64**(-7)
        integer,parameter :: n_64 = 3

        ! axpy test for real32
        allocate(x_input_32(n_32), y_input_32(n_32), z_output_32(n_32), expected_output_32(n_32))

        x_input_32 = (/1.0_real32, 2.0_real32, 3.0_real32/)
        y_input_32 = (/4.0_real32, 5.0_real32, 6.0_real32/)
        expected_output_32 = (/6.0_real32, 9.0_real32, 12.0_real32/)

        z_output_32 = axpy(alpha_input_32, x_input_32, y_input_32)

        do i=1,n_32
            if(abs(expected_output_32(i) - z_output_32(i)) > tolerance_32) then
                vectors_equal = .false.
                exit
            end if
        end do
          
        if(vectors_equal)then
            r = 0
        else
            r = 1
        endif
    
        deallocate(x_input_32, y_input_32, z_output_32, expected_output_32)

        ! axpy test for real64
        allocate(x_input_64(n_64), y_input_64(n_64), z_output_64(n_64), expected_output_64(n_64))

        x_input_64 = (/1.0_real64, 2.0_real64, 3.0_real64/)
        y_input_64 = (/4.0_real64, 5.0_real64, 6.0_real64/)
        expected_output_64 = (/6.0_real64, 9.0_real64, 12.0_real64/)

        z_output_64 = axpy(alpha_input_64, x_input_64, y_input_64)

        do i=1,n_64
            if(abs(expected_output_64(i) - z_output_64(i)) > tolerance_64) then
                vectors_equal = .false.
                exit
            end if
        end do
          
        if(vectors_equal)then
            r = 0
        else
            r = 1
        endif
    
        deallocate(x_input_64, y_input_64, z_output_64, expected_output_64)

    end function test_function

end program test_program