module SELF_HIP_Support

  use iso_c_binding
  use iso_fortran_env
  use hipfort
  use hipfort_check

  logical, private :: acquired = .false.
  
contains

  subroutine init_gpu(dev_id)
    integer,intent(in) :: dev_id
        !! The integer id of the GPU to use
    ! Number of compute devices
    integer :: ndevices

    ! Initialise resources the best practice way
    if (.not. acquired) then
      ! Initialise HIP
      call hipCheck(hipinit(0))
    end if

    ! Get the number of compute devices
    call hipCheck(hipgetdevicecount(ndevices))

    if ((dev_id .ge. 0) .and. (dev_id .lt. ndevices)) then
      ! Choose a compute device
      call hipCheck(hipsetdevice(dev_id))
    else
      write (error_unit,*) 'Error, dev_id was not inside the range of available devices.'
      stop 1
    end if

    ! Reset the GPU and all resources allocated on it
    call reset_gpu

    ! Set the device id for the GPU
    device_id = dev_id

  end subroutine init_gpu

  subroutine reset_gpu

    use hipfort
    use hipfort_check

    ! Release all resources on the gpu
    if (acquired) then
      call hipCheck(hipdevicereset())
    end if
  end subroutine reset_gpu

end module SELF_HIP_Support
