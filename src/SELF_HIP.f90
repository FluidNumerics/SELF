MODULE SELF_HIP

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ==============================================================================
! hipfort: FORTRAN INTERFACEs for GPU kernels
! ==============================================================================
! Copyright (c) 2020-2022 Advanced Micro Devices, Inc. All rights reserved.
! [MITx11 License]
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE ISO_C_BINDING

  IMPLICIT NONE

  INTERFACE hipGetDeviceCount
    FUNCTION hipGetDeviceCount_(count) BIND(c,name="hipGetDeviceCount")
      USE ISO_C_BINDING
      USE SELF_HIP_enums
      IMPLICIT NONE
      INTEGER(KIND(hipSuccess)) :: hipGetDeviceCount_
      INTEGER(C_INT) :: count
    END FUNCTION
  END INTERFACE

  INTERFACE hipMalloc
    FUNCTION hipMalloc_(ptr,mySize) BIND(c,name="hipMalloc")
      USE ISO_C_BINDING
      USE SELF_HIP_enums
      IMPLICIT NONE
      INTEGER(KIND(hipSuccess)) :: hipMalloc_
      TYPE(C_PTR) :: ptr
      INTEGER(C_SIZE_T),VALUE :: mySize
    END FUNCTION
  END INTERFACE hipMalloc

  INTERFACE hipFree
    FUNCTION hipFree_(ptr) BIND(c,name="hipFree")
      USE ISO_C_BINDING
      USE SELF_HIP_enums
      IMPLICIT NONE
      INTEGER(KIND(hipSuccess)) :: hipFree_
      TYPE(C_PTR),VALUE :: ptr
    END FUNCTION
  END INTERFACE hipFree

  INTERFACE hipMemcpy
    FUNCTION hipMemcpy_(dest,src,sizeBytes,myKind) BIND(c,name="hipMemcpy")
      USE ISO_C_BINDING
      USE SELF_HIP_enums
      IMPLICIT NONE
      INTEGER(KIND(hipSuccess)) :: hipMemcpy_
      TYPE(C_PTR),VALUE :: dest
      TYPE(C_PTR),VALUE :: src
      INTEGER(C_SIZE_T),VALUE :: sizeBytes
      INTEGER(KIND(hipMemcpyHostToHost)),VALUE :: myKind
    END FUNCTION hipMemcpy_
  END INTERFACE hipMemcpy

  INTERFACE hipSetDevice
    FUNCTION hipSetDevice_(deviceId) BIND(c,name="hipSetDevice")
      USE ISO_C_BINDING
      USE SELF_HIP_enums
      IMPLICIT NONE
      INTEGER(KIND(hipSuccess)) :: hipSetDevice_
      INTEGER(C_INT),VALUE :: deviceId
    END FUNCTION hipSetDevice_
  END INTERFACE hipSetDevice

CONTAINS

  SUBROUTINE hipCheck(hipError_t)
    USE SELF_HIP_enums
    IMPLICIT NONE
    INTEGER(KIND(hipSuccess)) :: hipError_t

    IF (hipError_t /= hipSuccess) THEN
      WRITE (*,*) "HIP ERROR: Error code = ",hipError_t
      CALL EXIT(hipError_t)
    END IF
  END SUBROUTINE hipCheck

  subroutine hipblasCheck(hipblasError_t)
    use hipfort_hipblas_enums

    implicit none

    integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasError_t

    if(hipblasError_t /= HIPBLAS_STATUS_SUCCESS)then
       write(*,*) "HIPBLAS ERROR: Error code = ", hipblasError_t
       call exit(hipblasError_t)
    end if
  end subroutine hipblasCheck

END MODULE SELF_HIP
