! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_GPU

  use iso_c_binding
  use SELF_GPU_enums

  implicit none

  interface hipGetDeviceCount
#ifdef HAVE_HIP
    function hipGetDeviceCount_(count) bind(c,name="hipGetDeviceCount")
#elif HAVE_CUDA
      function hipGetDeviceCount_(count) bind(c,name="cudaGetDeviceCount")
#endif
        use iso_c_binding
        use SELF_GPU_enums
        implicit none
        integer(kind(hipSuccess)) :: hipGetDeviceCount_
        integer(c_int) :: count
      endfunction
      endinterface

      interface hipMalloc
#ifdef HAVE_HIP
        function hipMalloc_(ptr,mySize) bind(c,name="hipMalloc")
#elif HAVE_CUDA
          function hipMalloc_(ptr,mySize) bind(c,name="cudaMalloc")
#endif
            use iso_c_binding
            use SELF_GPU_enums
            implicit none
            integer(kind(hipSuccess)) :: hipMalloc_
            type(c_ptr) :: ptr
            integer(c_size_t),value :: mySize
          endfunction
          endinterface hipMalloc

          interface hipFree
#ifdef HAVE_HIP
            function hipFree_(ptr) bind(c,name="hipFree")
#elif HAVE_CUDA
              function hipFree_(ptr) bind(c,name="cudaFree")
#endif
                use iso_c_binding
                use SELF_GPU_enums
                implicit none
                integer(kind(hipSuccess)) :: hipFree_
                type(c_ptr),value :: ptr
              endfunction
              endinterface hipFree

              interface hipMemcpy
#ifdef HAVE_HIP
                function hipMemcpy_(dest,src,sizeBytes,myKind) bind(c,name="hipMemcpy")
#elif HAVE_CUDA
                  function hipMemcpy_(dest,src,sizeBytes,myKind) bind(c,name="cudaMemcpy")
#endif
                    use iso_c_binding
                    use SELF_GPU_enums
                    implicit none
                    integer(kind(hipSuccess)) :: hipMemcpy_
                    type(c_ptr),value :: dest
                    type(c_ptr),value :: src
                    integer(c_size_t),value :: sizeBytes
                    integer(kind(hipMemcpyHostToHost)),value :: myKind
                  endfunction hipMemcpy_
                  endinterface hipMemcpy

                  contains

                  subroutine gpuCheck(gpuError_t)
                    implicit none
                    integer(kind(hipSuccess)) :: gpuError_t

                    if(gpuError_t /= hipSuccess) then
                      write(*,*) "GPU ERROR: Error code = ",gpuError_t
                      call exit(gpuError_t)
                    endif
                  endsubroutine gpuCheck

                  function GPUAvailable() result(avail)
                    implicit none
                    logical :: avail
                    ! Local
                    integer(c_int) :: gpuCount
                    integer(kind(hipSuccess)) :: err

                    err = hipGetDeviceCount(gpuCount)
                    if(gpuCount > 0 .and. err == hipSuccess) then
                      avail = .true.
                    else
                      avail = .false.
                    endif

                  endfunction GPUAvailable

                  ! subroutine HostToDevice(fsource,fdest)
                  !   implicit none
                  !   type(c_ptr), intent(inout) :: fsource
                  !   type(c_ptr), intent(inout) :: fdest

                  !   call gpuCheck(hipMemcpy(fdest,fsource,sizeof(fsource),hipMemcpyHostToDevice))

                  ! endsubroutine HostToDevice

                  ! subroutine DeviceToHost(fsource,fdest)
                  !   implicit none
                  !   type(c_ptr), intent(inout)    :: fsource
                  !   type(c_ptr), intent(inout) :: fdest

                  !   call gpuCheck(hipMemcpy(fdest,fsource,sizeof(fsource),hipMemcpyDeviceToHost))

                  ! endsubroutine DeviceToHost

                  endmodule SELF_GPU
