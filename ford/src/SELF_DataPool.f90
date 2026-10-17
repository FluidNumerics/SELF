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

module SELF_DataPool
!! High-water-mark storage for the mesh-sized arrays of the data classes (AMR Stage 6b).
!!
!! Why this exists
!! ---------------
!! Adaptive mesh refinement changes the element count almost every epoch, and the adaptive loop
!! responded by freeing and re-initializing every mesh-sized field: seven in the model plus six
!! in the geometry. Profiling one MI300X showed that cycle to be the single largest component of
!! an adaptation - larger than the solution transfer and larger than geometry regeneration - and
!! that the cost is NOT the device allocator (hipMalloc plus hipFree together are only ~8% of an
!! adaptation) but the host-side work that accompanies a fresh allocation: allocating, zeroing
!! ~140 MB of arrays per epoch at modest resolution, and then uploading those zeros to the
!! device only for the solution transfer to overwrite them.
!!
!! How it avoids the stride problem
!! --------------------------------
!! A naive high-water-mark scheme allocates the element dimension larger than the logical
!! element count. That does not work here: the device kernels take nElem as BOTH the launch
!! bound and the array stride (see SC_2D_INDEX in src/gpu/SELF_GPU_Macros.h), the variable index
!! is the slowest-varying dimension so a logical element prefix is not a contiguous prefix of
!! the buffer, HDF5 hyperslabs are taken from shape(), and several routines assert on nElem
!! equality.
!!
!! Instead, each array is backed by a rank-1 pool that is grown monotonically, and the array
!! itself is a pointer remapped onto the leading part of that pool at the EXACT logical shape.
!! nElem therefore keeps meaning the logical element count everywhere, every stride, bound,
!! shape() and equality assertion stays correct with no change, and no padding is ever computed
!! over or transferred. Only the pool is oversized, and nothing indexes the pool directly.
!!
!! The device side needs no shape at all - a device pointer is an address - so a device buffer is
!! reused whenever the bytes required fit the bytes already allocated.
!!
!! The pools are POINTER rather than ALLOCATABLE components deliberately: a derived-type
!! component cannot carry the TARGET attribute, and pointing at a component of an object that is
!! not itself a target is not conforming. An allocated pointer, by contrast, is always a valid
!! target, so the public arrays can be remapped onto it safely regardless of how the owning
!! object was declared.
!!
!! Growth policy
!! -------------
!! poolGrowth is the over-allocation factor applied when a pool must grow. 1.0 reproduces exact
!! sizing, i.e. the pre-6b behaviour, and is the escape hatch if capacity reuse is ever
!! suspected of misbehaving: set it to 1.0 and every Resize reallocates exactly as Init used to.
!! Above 1.0, a monotonically growing element count (the usual case while a wavefront expands)
!! stops reallocating almost immediately.
!!
!! This mirrors the amortized-capacity pattern already used for the quadtree forest in
!! EnsureCapacity_QuadTreeMesh2D (src/SELF_QuadTreeMesh_2D.f90).

  use SELF_Constants

  implicit none

  real(prec),parameter :: poolGrowth = 1.25_prec
  !! Over-allocation factor for pool growth; 1.0 disables capacity reuse (see above).

  public :: EnsurePool,poolGrowth

contains

  subroutine EnsurePool(pool,needed)
    !! Guarantee that pool has at least `needed` elements, growing it by poolGrowth when it must
    !! reallocate. Existing contents are NOT preserved: every caller re-establishes the array
    !! contents after resizing (the solution through the AMR transfer, the geometry through
    !! GenerateFromMesh, everything else by being written before it is read), so copying the old
    !! data would be wasted bandwidth.
    implicit none
    real(prec),pointer,contiguous,intent(inout) :: pool(:)
    integer,intent(in) :: needed
    ! Local
    integer :: capacity

    if(associated(pool)) then
      if(size(pool) >= needed) return
      deallocate(pool)
    endif

    capacity = max(needed,int(real(needed,prec)*poolGrowth))
    allocate(pool(1:capacity))

  endsubroutine EnsurePool

endmodule SELF_DataPool
