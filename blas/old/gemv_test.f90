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
        integer :: i, j
        ! real32
        ! real(real32) :: alpha_input_32 = 2.0_real32
        ! real(real32), pointer, dimension(:) :: x_input_32
        ! real(real32), pointer, dimension(:) :: y_input_32
        ! real(real32), pointer, dimension(:) :: z_output_32
        ! real(real32), pointer, dimension(:) :: expected_output_32
        ! real(real32),parameter :: tolerance_32 = 10.0_real32**(-3)
        ! integer,parameter :: n_32 = 3
        ! real64
        real(real64) :: alpha_input_64 = 2.0_real64
        real(real64) :: beta_input_64 = 3.0_real64
        real(real64), pointer, dimension(:,:) :: A_input_64
        real(real64), pointer, dimension(:) :: x_input_64
        real(real64), pointer, dimension(:) :: y_input_64
        real(real64), pointer, dimension(:) :: z_output_64
        real(real64), pointer, dimension(:) :: expected_output_64
        real(real64),parameter :: tolerance_64 = 10.0_real64**(-7)
        integer,parameter :: m_64 = 2
        integer,parameter :: n_64 = 3

        ! gemv test for real64
        ! allocate(A_input_64(m_64,n_64), x_input_64(n_64), y_input_64(m_64), z_output_64(m_64), expected_output_64(m_64))
        allocate(A_input_64(n_64,m_64), x_input_64(n_64), y_input_64(m_64), z_output_64(m_64), expected_output_64(m_64))

        A_input_64 = (reshape((/1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64/), (/n_64, m_64/)))
        print*,A_input_64
        print*,shape(A_input_64)
        x_input_64 = (/1.0_real64, 2.0_real64, 3.0_real64/)
        y_input_64 = (/4.0_real64, 5.0_real64/)
        expected_output_64 = (/40.0_real64, 79.0_real64/)

        z_output_64 = gemv(alpha_input_64, A_input_64, x_input_64, beta_input_64, y_input_64)

        print*,z_output_64

        do i=1,m_64
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
    
        deallocate(A_input_64, x_input_64, y_input_64, z_output_64, expected_output_64)

    end function test_function

end program test_program