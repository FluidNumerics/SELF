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

program structuredmesh_subcomm_2d
!! Builds a mesh on a communicator that is a strict SUBSET of MPI_COMM_WORLD.
!!
!! structuredmesh_external_comm_2d passes a dup of MPI_COMM_WORLD, which has
!! the same size and rank ordering as MPI_COMM_WORLD -- so it cannot tell a
!! working implementation apart from one that quietly ignores `comm` and uses
!! MPI_COMM_WORLD instead. This test can: the last world rank is excluded from
!! the mesh communicator and never calls into SELF, so any collective SELF
!! performs on MPI_COMM_WORLD deadlocks, and any size/rank SELF reads from
!! MPI_COMM_WORLD disagrees with the subcommunicator.
!!
!! This is the shape a Python caller produces when it hands SELF a communicator
!! from mpi4py that covers only the ranks assigned to the solver.

  use self_data
  use self_LinearEuler2D
  use mpi

  implicit none
  real(prec),parameter :: rho0 = 1.225_prec
  real(prec),parameter :: rhoprime = 0.01_prec
  real(prec),parameter :: c = 1.0_prec
  real(prec),parameter :: Lr = 0.06_prec
  real(prec),parameter :: x0 = 0.5_prec
  real(prec),parameter :: y0 = 0.5_prec
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 15
  real(prec),parameter :: dt = 2.0_prec*10.0_prec**(-4)
  real(prec),parameter :: endtime = 5.0_prec*10.0_prec**(-3)
  real(prec),parameter :: iointerval = 5.0_prec*10.0_prec**(-3)
  real(prec) :: ef
  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)
  integer :: worldRank,worldSize
  integer :: subComm,subRank,subSize,color
  integer :: cmpResult,ierror

  call mpi_init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD,worldRank,ierror)
  call MPI_Comm_size(MPI_COMM_WORLD,worldSize,ierror)

  ! Exclude the highest world rank from the solver communicator. With the
  ! default two-rank test launch this leaves a single-rank subcommunicator;
  ! with more ranks it leaves a genuine multi-rank decomposition.
  if(worldRank < worldSize-1) then
    color = 0
  else
    color = MPI_UNDEFINED
  endif
  call MPI_Comm_split(MPI_COMM_WORLD,color,worldRank,subComm,ierror)

  if(color /= MPI_UNDEFINED) then

    call MPI_Comm_rank(subComm,subRank,ierror)
    call MPI_Comm_size(subComm,subSize,ierror)

    bcids(1:4) = [SELF_BC_NONORMALFLOW, & ! South
                  SELF_BC_NONORMALFLOW, & ! East
                  SELF_BC_NONORMALFLOW, & ! North
                  SELF_BC_NONORMALFLOW] ! West

    call mesh%StructuredMesh(5,5,2,2,0.1_prec,0.1_prec,bcids,comm=subComm)

    ! The decomposition must describe the subcommunicator, not MPI_COMM_WORLD.
    call MPI_Comm_compare(mesh%decomp%mpiComm,subComm,cmpResult,ierror)
    if(cmpResult /= MPI_IDENT .and. cmpResult /= MPI_CONGRUENT) then
      print*,"Error: mesh decomposition is not using the provided communicator"
      stop 1
    endif
    if(mesh%decomp%nRanks /= subSize .or. mesh%decomp%rankId /= subRank) then
      print*,"Error: decomposition rank/size disagree with the subcommunicator", &
        mesh%decomp%rankId,subRank,mesh%decomp%nRanks,subSize
      stop 1
    endif
    if(mesh%decomp%nRanks == worldSize) then
      print*,"Error: decomposition spans MPI_COMM_WORLD rather than the subcommunicator"
      stop 1
    endif

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call modelobj%Init(mesh,geometry)
    modelobj%prescribed_bcs_enabled = .false.
    modelobj%tecplot_enabled = .false.
    modelobj%rho0 = rho0

    call modelobj%SphericalSoundWave(rhoprime,Lr,x0,y0,c)

    call modelobj%SetTimeIntegrator(integrator)
    call modelobj%ForwardStep(endtime,dt,iointerval)

    ef = modelobj%entropy
    if(ef /= ef) then
      print*,"Error: Final entropy is inf or nan",ef
      stop 1
    endif

    call modelobj%free()
    call mesh%free() ! must NOT finalize MPI: the communicator is externally managed
    call geometry%free()
    call interp%free()

    call MPI_Comm_free(subComm,ierror)

  endif

  ! Ranks outside the solver communicator rejoin here. Reaching this barrier
  ! at all is the proof that SELF ran no MPI_COMM_WORLD collectives.
  call MPI_Barrier(MPI_COMM_WORLD,ierror)

  call mpi_finalize(ierror)

endprogram structuredmesh_subcomm_2d
