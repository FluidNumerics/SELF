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

program domaindecomposition_two_meshes
!! Two meshes coexisting in one process, with SELF owning the MPI lifecycle.
!!
!! SELF initializes MPI when the first mesh is created and must keep it alive
!! until the last mesh is gone. Freeing the first mesh while the second is
!! still in use must NOT finalize MPI -- otherwise every subsequent operation
!! on the survivor runs against a finalized library.
!!
!! This is the no-comm counterpart to structuredmesh_external_comm_2d: there
!! the caller owns MPI, here SELF does.

  use self_data
  use SELF_Mesh_2D
  use mpi

  implicit none
  type(Mesh2D),target :: meshA,meshB
  integer :: bcids(1:4)
  integer :: ierror
  integer :: survivorRank,survivorSize
  logical :: mpiIsFinalized

  bcids(1:4) = [SELF_BC_NONORMALFLOW, & ! South
                SELF_BC_NONORMALFLOW, & ! East
                SELF_BC_NONORMALFLOW, & ! North
                SELF_BC_NONORMALFLOW] ! West

  ! First mesh: MPI is not yet up, so SELF initializes it and takes ownership.
  call meshA%StructuredMesh(5,5,2,2,0.1_prec,0.1_prec,bcids)

  ! Second mesh: MPI is already up. SELF must reuse it and must NOT re-init.
  call meshB%StructuredMesh(5,5,2,2,0.1_prec,0.1_prec,bcids)

  ! Retire the mesh that owns MPI while the other is still live.
  call meshA%free()

  call MPI_Finalized(mpiIsFinalized,ierror)
  if(mpiIsFinalized) then
    print*,"Error: freeing the first mesh finalized MPI while a second mesh is still live"
    stop 1
  endif

  ! The survivor's communicator must still be usable.
  call MPI_Comm_rank(meshB%decomp%mpiComm,survivorRank,ierror)
  if(ierror /= MPI_SUCCESS) then
    print*,"Error: the surviving mesh's communicator is no longer valid",ierror
    stop 1
  endif
  call MPI_Comm_size(meshB%decomp%mpiComm,survivorSize,ierror)
  if(ierror /= MPI_SUCCESS) then
    print*,"Error: MPI_Comm_size failed on the surviving mesh's communicator",ierror
    stop 1
  endif
  if(survivorRank /= meshB%decomp%rankId .or. survivorSize /= meshB%decomp%nRanks) then
    print*,"Error: surviving mesh decomposition is inconsistent with its communicator", &
      survivorRank,meshB%decomp%rankId,survivorSize,meshB%decomp%nRanks
    stop 1
  endif

  ! Freeing the last mesh releases MPI.
  call meshB%free()

  call MPI_Finalized(mpiIsFinalized,ierror)
  if(.not. mpiIsFinalized) then
    print*,"Error: freeing the last SELF-owned mesh did not finalize MPI"
    stop 1
  endif

endprogram domaindecomposition_two_meshes
