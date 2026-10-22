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

program test
!! Verifies that a conforming side exchange stops before posting MPI messages
!! when the mesh's side-id tag span exceeds the communicator's MPI_TAG_UB.
!!
!! The mesh-side count is changed after construction instead of allocating a
!! mesh large enough to reach a real implementation's tag limit. SELF requires
!! global side ids to lie in 1..nUniqueSides, so this is the same tag-space
!! contract as a genuinely large mesh and exercises the real exchange entry
!! point without a large run.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_MappedScalar_2D_t
  use mpi

  implicit none

  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(MappedScalar2D_t) :: f
  integer :: bcids(1:4)
  integer :: tagUpperBound
  integer :: ierror
  logical :: flag

  bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                SELF_BC_PRESCRIBED, & ! East
                SELF_BC_PRESCRIBED, & ! North
                SELF_BC_PRESCRIBED] ! West

  call mesh%StructuredMesh(3,3,1,1,0.1_prec,0.1_prec,bcids)
  call interp%Init(N=1,controlNodeType=GAUSS,M=2,targetNodeType=UNIFORM)
  call f%Init(interp,2,mesh%nElem)

  ! Two tag blocks at the full MPI_TAG_UB side-id stride require tags past the
  ! limit without multiplying the two large values in the test itself.
  call MPI_Comm_get_attr(mesh%decomp%mpiComm,MPI_TAG_UB,tagUpperBound,flag,ierror)
  if(ierror /= MPI_SUCCESS .or. .not. flag) then
    error stop "SELF side-exchange tag guard test could not query MPI_TAG_UB."
  endif
  mesh%nUniqueSides = tagUpperBound
  call f%SideExchange(mesh)

  error stop "SELF side-exchange tag guard did not stop the exchange."

endprogram test