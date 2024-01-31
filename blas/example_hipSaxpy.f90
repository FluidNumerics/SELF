program blas_program
    implicit none
  
    integer :: exit_code
    exit_code = example_hipSaxpy()
    print*,exit_code
  
    stop exit_code
  
    contains
    integer function example_hipSaxpy() result(r)
      use hipfort_hipblas
      use self_hip
      use self_hip_enums
      use self_constants
      use iso_c_binding
      use iso_fortran_env
      implicit none
  
      ! hipblasDaxpy inputs
      ! use type(c_ptr) instead of hipbasHandle_t
      type(c_ptr) :: handle
      integer,parameter :: n = 3
      real(real32) :: alpha = 2.0_real32
      real, pointer, dimension(:) :: x
      integer :: incx = 1
      real, pointer, dimension(:) :: y
      integer :: incy = 1
  
      type(c_ptr) :: x_gpu
      type(c_ptr) :: y_gpu
  
      ! expected results
      real, pointer, dimension(:) :: expected_y
  
      ! tolerance
      real(prec),parameter :: tolerance = 10.0_real32**(-3)
  
      ! compare vectors declarations
      logical :: vectors_equal = .true.
      integer :: i
  
      integer :: status
  
      ! execute block
      allocate(x(n), y(n), expected_y(n))
      call hipcheck(hipmalloc(x_gpu, sizeof(x)))
      call hipcheck(hipmalloc(y_gpu, sizeof(y)))
  
      call hipblasCheck(hipblasCreate(handle))
  
      x = (/1.0_real32, 2.0_real32, 3.0_real32/)
      y = (/4.0_real32, 5.0_real32, 6.0_real32/)
      expected_y = (/6.0_real32, 9.0_real32, 12.0_real32/)    
  
      call hipcheck(hipmemcpy(x_gpu, c_loc(x), sizeof(x), hipmemcpyhosttodevice))
      call hipcheck(hipmemcpy(y_gpu, c_loc(y), sizeof(y), hipmemcpyhosttodevice))
  
      ! alpha * x + y
      ! 2 * (1,2,3) + (4,5,6)
      ! = (6,9,12)
      status = hipblasSaxpy_(handle, n, alpha, x_gpu, incx, y_gpu, incy)
      
      print*,status
  
      call hipcheck(hipmemcpy(c_loc(y), y_gpu, sizeof(y), hipmemcpydevicetohost))
  
      do i=1,n
        if(abs(expected_y(i) - y(i)) > tolerance) then
          vectors_equal = .false.
          exit
        end if
      end do
      
      if(vectors_equal)then
        r = 0
      else
        r = 1
      endif
  
      deallocate(x, y, expected_y)
      call hipcheck(hipfree(x_gpu))
      call hipcheck(hipfree(y_gpu))
  
      call hipblasCheck(hipblasDestroy(handle))
  
    end function example_hipSaxpy
  
  end program blas_program
  