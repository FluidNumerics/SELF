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

    procedure :: Init => Init_DomainDecomposition
    procedure :: Free => Free_DomainDecomposition

    procedure :: SetElemToRank => SetElemToRank_DomainDecomposition

  endtype DomainDecomposition

  interface
    function check_gpu_aware_support() bind(c,name="check_gpu_aware_support")
      use iso_c_binding
      integer(c_int) :: check_gpu_aware_support
    endfunction check_gpu_aware_support
  endinterface
contains

  subroutine Init_DomainDecomposition(this,enableMPI)
    implicit none
    class(DomainDecomposition),intent(out) :: this
    logical,intent(in) :: enableMPI
    ! Local
    integer       :: ierror
    integer(c_int) :: gpuaware

    this%mpiComm = 0
    this%mpiPrec = prec
    this%rankId = 0
    this%nRanks = 1
    this%nElem = 0
    this%mpiEnabled = enableMPI

    if(enableMPI) then
      this%mpiComm = MPI_COMM_WORLD
      print*,__FILE__," : Initializing MPI"
      call MPI_INIT(ierror)

      if(check_gpu_aware_support() == 0) then
        print*,__FILE__" : Error! GPU Aware support is not detected. Stopping."
        call MPI_FINALIZE(ierror)
        stop 1
      endif

      call MPI_COMM_RANK(this%mpiComm,this%rankId,ierror)
      call MPI_COMM_SIZE(this%mpiComm,this%nRanks,ierror)
      print*,__FILE__," : Rank ",this%rankId+1,"/",this%nRanks," checking in."
    else
      print*,__FILE__," : MPI not initialized. No domain decomposition used."
    endif

    if(prec == real32) then
      this%mpiPrec = MPI_FLOAT
    else
      this%mpiPrec = MPI_DOUBLE
    endif

    allocate(this%offsetElem(1:this%nRanks+1))

  endsubroutine Init_DomainDecomposition

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
