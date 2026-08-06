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

! Device-memory conservation across repeated Init/Free of the 2-D data classes.
!
! The adaptive loop performs one Free/Init cycle per epoch on every mesh-sized field the model
! owns (see Regrid_DGModel2D_t), so a single device pointer that Init allocates and Free omits
! is enough to exhaust HBM over a long adapting run. That is not a hypothetical: boundarynormal_gpu
! was allocated in Init_MappedScalar2D and never released, leaking
! (N+1)*4*nElem*2*nvar reals per cycle per field, which aborted a 640-epoch maxLevel-8
! ultrasound run with an out-of-memory error after ~555 adaptations on a 192 GiB MI300X.
!
! Nothing in the suite observed device memory, so the leak was invisible. This test closes that
! gap for every 2-D class the adaptive loop cycles: it records free device memory, performs many
! Init/Free cycles at two alternating element counts, and requires the free byte count to return
! to its starting value. The alternating sizes also exercise the path where a reallocation
! changes size rather than repeating it.
!
! On CPU builds there is no device to interrogate and the class hierarchy has no _gpu members, so
! the body compiles out and the test trivially succeeds.
program test

#ifdef ENABLE_GPU
  use iso_c_binding
  use SELF_Constants
  use SELF_Lagrange
  use SELF_MappedScalar_2D
  use SELF_MappedVector_2D
  use SELF_Scalar_2D
  use SELF_Vector_2D
  use SELF_Tensor_2D
  use SELF_GPU
  use SELF_GPU_enums
  implicit none

  integer,parameter :: N = 4 ! polynomial degree
  integer,parameter :: nvar = 3
  integer,parameter :: nElemA = 64 ! the two alternating element counts
  integer,parameter :: nElemB = 96
  integer,parameter :: nWarm = 3 ! cycles before the baseline is taken
  integer,parameter :: nCycle = 24 ! measured cycles

  type(Lagrange),target :: interp
  integer(c_size_t) :: freeBaseline,freeNow,totalBytes
  integer :: cycle,nel

  call interp%Init(N=N,controlNodeType=GAUSS,M=N,targetNodeType=UNIFORM)

  ! Warm-up cycles first: the very first allocation on a device can move the free-memory
  ! figure through runtime bookkeeping unrelated to these classes, so the baseline is taken
  ! once that has settled.
  do cycle = 1,nWarm
    call CycleAll(interp,nvar,nElemA)
    call CycleAll(interp,nvar,nElemB)
  enddo

  call gpuCheck(hipMemGetInfo(freeBaseline,totalBytes))

  do cycle = 1,nCycle
    if(mod(cycle,2) == 0) then
      nel = nElemA
    else
      nel = nElemB
    endif
    call CycleAll(interp,nvar,nel)
  enddo

  call gpuCheck(hipMemGetInfo(freeNow,totalBytes))

  print*,"device free bytes: baseline =",freeBaseline,", after",nCycle,"cycles =",freeNow

  if(freeNow < freeBaseline) then
    print*,"ERROR: device memory leaked across Init/Free cycles."
    print*,"       leaked bytes =",freeBaseline-freeNow, &
      " over",nCycle,"cycles =",(freeBaseline-freeNow)/int(nCycle,c_size_t),"bytes/cycle"
    stop 1
  endif

  call interp%Free()

  print*,"PASS: device memory is conserved across repeated Init/Free of the 2-D data classes."

contains

  subroutine CycleAll(interp,nvar,nElem)
    !! One Init/Free cycle of every 2-D data class the adaptive loop re-creates per epoch.
    implicit none
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nvar
    integer,intent(in) :: nElem
    ! Local
    type(MappedScalar2D) :: ms
    type(MappedVector2D) :: mv
    type(Scalar2D) :: s
    type(Vector2D) :: v
    type(Tensor2D) :: t

    call ms%Init(interp,nvar,nElem)
    call ms%Free()

    call mv%Init(interp,nvar,nElem)
    call mv%Free()

    call s%Init(interp,nvar,nElem)
    call s%Free()

    call v%Init(interp,nvar,nElem)
    call v%Free()

    call t%Init(interp,nvar,nElem)
    call t%Free()

  endsubroutine CycleAll

#else
  implicit none
  print*,"SKIP: CPU build has no device memory to account for."
#endif

endprogram test
