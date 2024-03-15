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
        logical :: matrices_equal = .true.
        integer i, j
        ! real32
        real(real32) :: alpha_input_32 = 2.0_real32
        real(real32) :: beta_input_32 = 3.0_real32
        real(real32), pointer, dimension(:,:) :: A_input_32
        real(real32), pointer, dimension(:,:) :: B_input_32
        real(real32), pointer, dimension(:,:) :: C_input_32
        real(real32), pointer, dimension(:,:) :: expected_output_32
        real(real32),parameter :: tolerance_32 = 10.0_real32**(-3)
        integer,parameter :: m_32 = 2
        integer,parameter :: n_32 = 3
        integer,parameter :: k_32 = 1
        ! real64
        real(real64) :: alpha_input_64 = 2.0_real64
        real(real64) :: beta_input_64 = 3.0_real64
        real(real64), pointer, dimension(:,:) :: A_input_64
        real(real64), pointer, dimension(:,:) :: B_input_64
        real(real64), pointer, dimension(:,:) :: C_input_64
        real(real64), pointer, dimension(:,:) :: expected_output_64
        real(real64),parameter :: tolerance_64 = 10.0_real64**(-7)
        integer,parameter :: m_64 = 2
        integer,parameter :: n_64 = 3
        integer,parameter :: k_64 = 1

        ! gemm test for real32
        allocate(A_input_32(k_32,m_32), B_input_32(n_32,k_32), C_input_32(m_32,n_32), expected_output_32(m_32,n_32))

        A_input_32 = (reshape((/1.0_real32, 2.0_real32/), (/k_32, m_32/)))
        B_input_32 = (reshape((/3.0_real32, 4.0_real32, 5.0_real32/), (/n_32, k_32/)))
        C_input_32 = transpose(reshape((/1.0_real32, 2.0_real32, 3.0_real32, 4.0_real32, 5.0_real32, 6.0_real32/), (/n_32, m_32/)))
        expected_output_32 = transpose(reshape((/9.0_real32, 14.0_real32, 19.0_real32, 24.0_real32, 31.0_real32, 38.0_real32/),(/n_32,m_32/)))

        call gemm(alpha_input_32, A_input_32, B_input_32, beta_input_32, C_input_32)

        print*,C_input_32
        print*,expected_output_32

        do i=1,m_32
            do j=1,n_32
                if(abs(expected_output_32(i,j) - C_input_32(i,j)) > tolerance_32) then
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
    
        deallocate(A_input_32, B_input_32, C_input_32, expected_output_32)

        ! gemm test for real64
        allocate(A_input_64(k_64,m_64), B_input_64(n_64,k_64), C_input_64(m_64,n_64), expected_output_64(m_64,n_64))

        A_input_64 = (reshape((/1.0_real64, 2.0_real64/), (/k_64, m_64/)))
        B_input_64 = (reshape((/3.0_real64, 4.0_real64, 5.0_real64/), (/n_64, k_64/)))
        C_input_64 = transpose(reshape((/1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64/), (/n_64, m_64/)))
        expected_output_64 = transpose(reshape((/9.0_real64, 14.0_real64, 19.0_real64, 24.0_real64, 31.0_real64, 38.0_real64/),(/n_64,m_64/)))

        call gemm(alpha_input_64, A_input_64, B_input_64, beta_input_64, C_input_64)

        print*,C_input_64
        print*,expected_output_64

        do i=1,m_64
            do j=1,n_64
                if(abs(expected_output_64(i,j) - C_input_64(i,j)) > tolerance_64) then
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
    
        deallocate(A_input_64, B_input_64, C_input_64, expected_output_64)

    end function test_function

end program test_program