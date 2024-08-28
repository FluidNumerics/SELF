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

module SELF_GPU_enums

  use iso_c_binding

  implicit none

  enum,bind(c)
    enumerator :: hipSuccess = 0
  endenum

  enum,bind(c)
    enumerator :: hipMemcpyHostToHost = 0
    enumerator :: hipMemcpyHostToDevice = 1
    enumerator :: hipMemcpyDeviceToHost = 2
    enumerator :: hipMemcpyDeviceToDevice = 3
    enumerator :: hipMemcpyDefault = 4
  endenum

  enum,bind(c)
    enumerator :: HIPBLAS_STATUS_SUCCESS = 0
    enumerator :: HIPBLAS_STATUS_NOT_INITIALIZED = 1
    enumerator :: HIPBLAS_STATUS_ALLOC_FAILED = 2
    enumerator :: HIPBLAS_STATUS_INVALID_VALUE = 3
    enumerator :: HIPBLAS_STATUS_MAPPING_ERROR = 4
    enumerator :: HIPBLAS_STATUS_EXECUTION_FAILED = 5
    enumerator :: HIPBLAS_STATUS_INTERNAL_ERROR = 6
    enumerator :: HIPBLAS_STATUS_NOT_SUPPORTED = 7
    enumerator :: HIPBLAS_STATUS_ARCH_MISMATCH = 8
    enumerator :: HIPBLAS_STATUS_HANDLE_IS_NULLPTR = 9
    enumerator :: HIPBLAS_STATUS_INVALID_ENUM = 10
    enumerator :: HIPBLAS_STATUS_UNKNOWN = 11
  endenum

  enum,bind(c)
#ifdef HAVE_CUDA
    enumerator :: HIPBLAS_OP_N = 0
#else
    enumerator :: HIPBLAS_OP_N = 111
#endif
#ifdef HAVE_CUDA
    enumerator :: HIPBLAS_OP_T = 1
#else
    enumerator :: HIPBLAS_OP_T = 112
#endif
  endenum

endmodule SELF_GPU_enums
