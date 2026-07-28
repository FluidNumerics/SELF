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

program read_hopr_external_comm_2d
!! Reads a HOPr mesh through an externally managed communicator.
!!
!! Read_HOPr takes the same optional `comm` as the structured constructors,
!! but it is the only entry point that forwards the communicator into parallel
!! HDF5 (h5pset_fapl_mpio_f). That makes it the highest-risk consumer of an
!! external handle and it had no coverage: a communicator that MPI accepts can
!! still be rejected -- or silently ignored -- by the HDF5 file access property
!! list.
!!
!! The caller owns the MPI lifecycle here, exactly as mpi4py does.

  use self_data
  use SELF_Mesh_2D
  use mpi

  implicit none
  type(Mesh2D),target :: mesh
  character(LEN=255) :: WORKSPACE
  integer :: dupComm,dupRank,dupSize
  integer :: cmpResult,ierror
  integer :: nElemLocal,nElemTotal
  logical :: mpiIsFinalized

  call mpi_init(ierror)
  call MPI_Comm_dup(MPI_COMM_WORLD,dupComm,ierror)
  call MPI_Comm_rank(dupComm,dupRank,ierror)
  call MPI_Comm_size(dupComm,dupSize,ierror)

  call get_environment_variable("WORKSPACE",WORKSPACE)
  call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5",comm=dupComm)

  call MPI_Comm_compare(mesh%decomp%mpiComm,dupComm,cmpResult,ierror)
  if(cmpResult /= MPI_IDENT .and. cmpResult /= MPI_CONGRUENT) then
    print*,"Error: Read_HOPr did not use the provided communicator"
    stop 1
  endif
  ! A dup is MPI_CONGRUENT with MPI_COMM_WORLD, so the check that actually
  ! detects `comm` being dropped is that the handle is not MPI_COMM_WORLD itself.
  call MPI_Comm_compare(mesh%decomp%mpiComm,MPI_COMM_WORLD,cmpResult,ierror)
  if(cmpResult == MPI_IDENT) then
    print*,"Error: Read_HOPr fell back to MPI_COMM_WORLD instead of the provided communicator"
    stop 1
  endif
  if(mesh%decomp%nRanks /= dupSize .or. mesh%decomp%rankId /= dupRank) then
    print*,"Error: decomposition rank/size disagree with the provided communicator", &
      mesh%decomp%rankId,dupRank,mesh%decomp%nRanks,dupSize
    stop 1
  endif

  ! Every element of the global mesh must be owned by exactly one rank.
  nElemLocal = mesh%nElem
  if(nElemLocal <= 0) then
    print*,"Error: rank",dupRank,"was assigned no elements"
    stop 1
  endif
  call MPI_Allreduce(nElemLocal,nElemTotal,1,MPI_INTEGER,MPI_SUM,dupComm,ierror)
  if(nElemTotal /= mesh%decomp%nElem) then
    print*,"Error: local element counts do not sum to the global element count", &
      nElemTotal,mesh%decomp%nElem
    stop 1
  endif

  call mesh%free() ! must NOT finalize MPI: the communicator is externally managed

  call MPI_Finalized(mpiIsFinalized,ierror)
  if(mpiIsFinalized) then
    print*,"Error: SELF finalized MPI despite an externally provided communicator"
    stop 1
  endif

  call MPI_Comm_free(dupComm,ierror)
  call mpi_finalize(ierror)

endprogram read_hopr_external_comm_2d
