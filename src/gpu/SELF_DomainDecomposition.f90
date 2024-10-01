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

module SELF_DomainDecomposition

  use SELF_DomainDecomposition_t
  use mpi
  use iso_c_binding

  implicit none

  type,extends(DomainDecomposition_t) :: DomainDecomposition
    type(c_ptr) :: elemToRank_gpu

  contains

    procedure :: Free => Free_DomainDecomposition

    procedure :: SetElemToRank => SetElemToRank_DomainDecomposition

  endtype DomainDecomposition

contains

  subroutine Free_DomainDecomposition(this)
    implicit none
    class(DomainDecomposition),intent(inout) :: this
    ! Local
    integer :: ierror

    if(associated(this%offSetElem)) then
      deallocate(this%offSetElem)
    endif
    if(associated(this%elemToRank)) then
      deallocate(this%elemToRank)
      call gpuCheck(hipFree(this%elemToRank_gpu))
    endif

    if(allocated(this%requests)) deallocate(this%requests)
    if(allocated(this%stats)) deallocate(this%stats)

    if(this%mpiEnabled) then
      print*,__FILE__," : Rank ",this%rankId+1,"/",this%nRanks," checking out."
      call MPI_FINALIZE(ierror)
    endif

  endsubroutine Free_DomainDecomposition

  subroutine SetElemToRank_DomainDecomposition(this,nElem)
    implicit none
    class(DomainDecomposition),intent(inout) :: this
    integer,intent(in) :: nElem
    ! Local
    integer :: iel

    this%nElem = nElem

    allocate(this%elemToRank(1:nelem))
    call gpuCheck(hipMalloc(this%elemToRank_gpu,sizeof(this%elemToRank)))

    call DomainDecomp(nElem, &
                      this%nRanks, &
                      this%offSetElem)

    do iel = 1,nElem
      call ElemToRank(this%nRanks, &
                      this%offSetElem, &
                      iel, &
                      this%elemToRank(iel))
    enddo
    call gpuCheck(hipMemcpy(this%elemToRank_gpu,c_loc(this%elemToRank),sizeof(this%elemToRank),hipMemcpyHostToDevice))

  endsubroutine SetElemToRank_DomainDecomposition

endmodule SELF_DomainDecomposition
