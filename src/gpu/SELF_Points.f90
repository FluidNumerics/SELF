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

module SELF_Points
!! GPU-resident backend for the Points class.
!!
!! Construction, point location, and basis-cache computation happen on the
!! host (inherited from Points_t). After LocatePoints, all per-point state
!! (elements, coordinates, lS/lT/lU caches) is mirrored to the device via
!! UpdateDevice. EvaluateScalar then runs entirely on the device using the
!! precomputed cached basis values — no Lagrange polynomial evaluation
!! happens in the kernel, only the tensor-product contraction against a
!! scalar field already resident on the device.

  use SELF_Constants
  use SELF_Points_t
  use SELF_Geometry_2D
  use SELF_Geometry_3D
  use SELF_MappedScalar_2D
  use SELF_MappedScalar_3D
  use SELF_GPU
  use iso_c_binding

  implicit none

  type,extends(Points_t),public :: Points
    character(3) :: backend = "gpu"
    integer :: nPointsAlloc = 0
      !! Capacity of elements_gpu / coordinates_gpu (set in Init).
    integer :: nCachedAlloc = 0
      !! Capacity (polynomial degree) of l*_cache_gpu; 0 = unallocated.
    type(c_ptr) :: elements_gpu = c_null_ptr
    type(c_ptr) :: coordinates_gpu = c_null_ptr
    type(c_ptr) :: lS_cache_gpu = c_null_ptr
    type(c_ptr) :: lT_cache_gpu = c_null_ptr
    type(c_ptr) :: lU_cache_gpu = c_null_ptr

  contains
    procedure,public :: Init => Init_Points
    procedure,public :: Free => Free_Points
    procedure,public :: UpdateDevice => UpdateDevice_Points

    ! Override the inherited specifics so LocatePoints automatically syncs
    ! device state after the host search.
    procedure,public :: LocatePoints_2D_Points_t => LocatePoints_2D_Points
    procedure,public :: LocatePoints_3D_Points_t => LocatePoints_3D_Points

    ! Add device-output specifics under the same EvaluateScalar generic. The
    ! inherited (host-output) specifics remain available for callers that
    ! want results in a Fortran array.
    procedure,private :: EvalScalar_2D_dev_Points
    procedure,private :: EvalScalar_3D_dev_Points
    generic,public :: EvaluateScalar => EvalScalar_2D_dev_Points,EvalScalar_3D_dev_Points

  endtype Points

  interface
    subroutine EvalScalarPoints_2D_gpu(values,elements,lS,lT,scalar,N,nPoints,nElem,nVar) &
      bind(c,name="EvalScalarPoints_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: values,elements,lS,lT,scalar
      integer(c_int),value :: N,nPoints,nElem,nVar
    endsubroutine EvalScalarPoints_2D_gpu
  endinterface

  interface
    subroutine EvalScalarPoints_3D_gpu(values,elements,lS,lT,lU,scalar,N,nPoints,nElem,nVar) &
      bind(c,name="EvalScalarPoints_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: values,elements,lS,lT,lU,scalar
      integer(c_int),value :: N,nPoints,nElem,nVar
    endsubroutine EvalScalarPoints_3D_gpu
  endinterface

contains

  subroutine Init_Points(this,nPoints,nDim)
    !! Allocate host + device storage for nPoints in nDim dimensions. The
    !! basis-cache device buffers are allocated lazily in UpdateDevice once
    !! the polynomial degree N is known (it is set by LocatePoints).
    implicit none
    class(Points),intent(out) :: this
    integer,intent(in) :: nPoints
    integer,intent(in) :: nDim
    ! Local
    integer(c_size_t) :: nBytes

    if(nDim /= 2 .and. nDim /= 3) then
      print*,"SELF_Points (gpu)::Init: nDim must be 2 or 3, got ",nDim
      stop 1
    endif

    this%nPoints = nPoints
    this%nDim = nDim
    this%nPointsAlloc = nPoints
    this%nCachedAlloc = 0
    allocate(this%x(1:nPoints,1:nDim))
    allocate(this%elements(1:nPoints))
    allocate(this%coordinates(1:nPoints,1:nDim))
    this%x = 0.0_prec
    this%elements = 0
    this%coordinates = 0.0_prec

    nBytes = int(nPoints,c_size_t)*c_sizeof(0_c_int)
    call gpuCheck(hipMalloc(this%elements_gpu,nBytes))

    nBytes = int(nPoints*nDim,c_size_t)*int(prec,c_size_t)
    call gpuCheck(hipMalloc(this%coordinates_gpu,nBytes))

  endsubroutine Init_Points

  subroutine Free_Points(this)
    implicit none
    class(Points),intent(inout) :: this

    if(allocated(this%x)) deallocate(this%x)
    if(allocated(this%elements)) deallocate(this%elements)
    if(allocated(this%coordinates)) deallocate(this%coordinates)
    if(allocated(this%lS_cache)) deallocate(this%lS_cache)
    if(allocated(this%lT_cache)) deallocate(this%lT_cache)
    if(allocated(this%lU_cache)) deallocate(this%lU_cache)

    if(c_associated(this%elements_gpu)) call gpuCheck(hipFree(this%elements_gpu))
    if(c_associated(this%coordinates_gpu)) call gpuCheck(hipFree(this%coordinates_gpu))
    if(c_associated(this%lS_cache_gpu)) call gpuCheck(hipFree(this%lS_cache_gpu))
    if(c_associated(this%lT_cache_gpu)) call gpuCheck(hipFree(this%lT_cache_gpu))
    if(c_associated(this%lU_cache_gpu)) call gpuCheck(hipFree(this%lU_cache_gpu))
    this%elements_gpu = c_null_ptr
    this%coordinates_gpu = c_null_ptr
    this%lS_cache_gpu = c_null_ptr
    this%lT_cache_gpu = c_null_ptr
    this%lU_cache_gpu = c_null_ptr

    this%nPoints = 0
    this%nDim = 0
    this%nCached = 0
    this%nPointsAlloc = 0
    this%nCachedAlloc = 0

  endsubroutine Free_Points

  subroutine LocatePoints_2D_Points(this,geometry)
    !! Run the host spatial-hash + Newton search (inherited from Points_t),
    !! then mirror per-point state to the device.
    implicit none
    class(Points),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry

    call this%Points_t%LocatePoints_2D_Points_t(geometry)
    call this%UpdateDevice()

  endsubroutine LocatePoints_2D_Points

  subroutine LocatePoints_3D_Points(this,geometry)
    implicit none
    class(Points),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry

    call this%Points_t%LocatePoints_3D_Points_t(geometry)
    call this%UpdateDevice()

  endsubroutine LocatePoints_3D_Points

  subroutine UpdateDevice_Points(this)
    !! Copy elements, coordinates, and the per-point Lagrange basis cache
    !! from host to device. Lazily (re)allocates the cache buffers if their
    !! degree changed since the previous call.
    implicit none
    class(Points),intent(inout) :: this
    ! Local
    integer(c_size_t) :: nBytes

    if(this%nPoints == 0) return

    nBytes = int(this%nPoints,c_size_t)*c_sizeof(0_c_int)
    call gpuCheck(hipMemcpy(this%elements_gpu,c_loc(this%elements),nBytes, &
                            hipMemcpyHostToDevice))

    nBytes = int(this%nPoints*this%nDim,c_size_t)*int(prec,c_size_t)
    call gpuCheck(hipMemcpy(this%coordinates_gpu,c_loc(this%coordinates),nBytes, &
                            hipMemcpyHostToDevice))

    if(this%nCached <= 0) return
    if(.not. allocated(this%lS_cache) .or. .not. allocated(this%lT_cache)) return

    ! (Re)allocate per-point basis caches on the device if needed.
    if(this%nCachedAlloc /= this%nCached) then
      if(c_associated(this%lS_cache_gpu)) call gpuCheck(hipFree(this%lS_cache_gpu))
      if(c_associated(this%lT_cache_gpu)) call gpuCheck(hipFree(this%lT_cache_gpu))
      if(c_associated(this%lU_cache_gpu)) call gpuCheck(hipFree(this%lU_cache_gpu))
      this%lS_cache_gpu = c_null_ptr
      this%lT_cache_gpu = c_null_ptr
      this%lU_cache_gpu = c_null_ptr
      nBytes = int((this%nCached+1)*this%nPoints,c_size_t)*int(prec,c_size_t)
      call gpuCheck(hipMalloc(this%lS_cache_gpu,nBytes))
      call gpuCheck(hipMalloc(this%lT_cache_gpu,nBytes))
      if(this%nDim == 3) call gpuCheck(hipMalloc(this%lU_cache_gpu,nBytes))
      this%nCachedAlloc = this%nCached
    endif

    nBytes = int((this%nCached+1)*this%nPoints,c_size_t)*int(prec,c_size_t)
    call gpuCheck(hipMemcpy(this%lS_cache_gpu,c_loc(this%lS_cache),nBytes, &
                            hipMemcpyHostToDevice))
    call gpuCheck(hipMemcpy(this%lT_cache_gpu,c_loc(this%lT_cache),nBytes, &
                            hipMemcpyHostToDevice))
    if(this%nDim == 3 .and. allocated(this%lU_cache)) then
      call gpuCheck(hipMemcpy(this%lU_cache_gpu,c_loc(this%lU_cache),nBytes, &
                              hipMemcpyHostToDevice))
    endif

  endsubroutine UpdateDevice_Points

  subroutine EvalScalar_2D_dev_Points(this,scalar,values_dev)
    !! Evaluate a 2D MappedScalar at the cached points on the device. The
    !! caller is responsible for hipMalloc-ing values_dev with capacity
    !! nPoints*nVar*prec bytes. Layout is column-major (nPoints, nVar).
    !!
    !! Requires: the cache must have been populated by LocatePoints (or a
    !! subsequent UpdateDevice). Mismatch with scalar%interp%N is fatal.
    implicit none
    class(Points),intent(in) :: this
    class(MappedScalar2D),intent(in) :: scalar
    type(c_ptr),intent(inout) :: values_dev

    if(this%nDim /= 2) then
      print*,"SELF_Points (gpu)::EvaluateScalar (2D): nDim must be 2"
      stop 1
    endif
    if(this%nCached /= scalar%interp%N .or. .not. c_associated(this%lS_cache_gpu)) then
      print*,"SELF_Points (gpu)::EvaluateScalar (2D): basis cache not synchronized; ", &
        "call LocatePoints (or UpdateDevice) first. nCached=", &
        this%nCached," scalar%N=",scalar%interp%N
      stop 1
    endif

    call EvalScalarPoints_2D_gpu(values_dev, &
                                 this%elements_gpu, &
                                 this%lS_cache_gpu,this%lT_cache_gpu, &
                                 scalar%interior_gpu, &
                                 scalar%interp%N,this%nPoints,scalar%nElem,scalar%nVar)

  endsubroutine EvalScalar_2D_dev_Points

  subroutine EvalScalar_3D_dev_Points(this,scalar,values_dev)
    implicit none
    class(Points),intent(in) :: this
    class(MappedScalar3D),intent(in) :: scalar
    type(c_ptr),intent(inout) :: values_dev

    if(this%nDim /= 3) then
      print*,"SELF_Points (gpu)::EvaluateScalar (3D): nDim must be 3"
      stop 1
    endif
    if(this%nCached /= scalar%interp%N .or. .not. c_associated(this%lS_cache_gpu) .or. &
       .not. c_associated(this%lU_cache_gpu)) then
      print*,"SELF_Points (gpu)::EvaluateScalar (3D): basis cache not synchronized; ", &
        "call LocatePoints (or UpdateDevice) first. nCached=", &
        this%nCached," scalar%N=",scalar%interp%N
      stop 1
    endif

    call EvalScalarPoints_3D_gpu(values_dev, &
                                 this%elements_gpu, &
                                 this%lS_cache_gpu,this%lT_cache_gpu,this%lU_cache_gpu, &
                                 scalar%interior_gpu, &
                                 scalar%interp%N,this%nPoints,scalar%nElem,scalar%nVar)

  endsubroutine EvalScalar_3D_dev_Points

endmodule SELF_Points
