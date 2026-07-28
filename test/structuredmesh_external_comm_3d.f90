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

program structuredmesh_external_comm_3d
!! Exercises the externally-managed MPI lifecycle used by Python (mpi4py)
!! callers: the program initializes MPI itself, hands a duplicated
!! communicator to StructuredMesh, runs a LinearEuler3D model, frees the
!! mesh, and verifies SELF did NOT finalize MPI before finalizing itself.

  use self_data
  use self_LinearEuler3D
  use mpi

  implicit none
  real(prec),parameter :: rho0 = 1.225_prec
  real(prec),parameter :: rhoprime = 0.01_prec
  real(prec),parameter :: c = 1.0_prec
  real(prec),parameter :: Lr = 0.06_prec
  real(prec),parameter :: x0 = 0.5_prec
  real(prec),parameter :: y0 = 0.5_prec
  real(prec),parameter :: z0 = 0.5_prec
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk3'
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 15
  real(prec),parameter :: dt = 5.0_prec*10.0_prec**(-4)
  real(prec),parameter :: endtime = 5.0_prec*10.0_prec**(-3)
  real(prec),parameter :: iointerval = 5.0_prec*10.0_prec**(-3)
  real(prec) :: ef
  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  integer :: dupComm,dupRank,dupSize
  integer :: cmpResult,ierror
  logical :: mpiIsFinalized

  ! The caller owns the MPI lifecycle (mpi4py does exactly this on import).
  call mpi_init(ierror)
  call MPI_Comm_dup(MPI_COMM_WORLD,dupComm,ierror)

  bcids(1:6) = [SELF_BC_RADIATION, & ! Bottom
                SELF_BC_RADIATION, & ! South
                SELF_BC_RADIATION, & ! East
                SELF_BC_RADIATION, & ! North
                SELF_BC_RADIATION, & ! West
                SELF_BC_RADIATION] ! Top

  call mesh%StructuredMesh(5,5,5,1,1,1, &
                           0.2_prec,0.2_prec,0.2_prec,bcids,comm=dupComm)

  ! The decomposition must be running on the communicator we handed over. A dup
  ! of MPI_COMM_WORLD has the same size and rank ordering, so rank/size alone
  ! cannot catch `comm` being dropped -- compare the handles instead. Note that
  ! a dup is MPI_CONGRUENT with MPI_COMM_WORLD, so the telling check is that the
  ! decomposition is not *identical* to MPI_COMM_WORLD.
  call MPI_Comm_rank(dupComm,dupRank,ierror)
  call MPI_Comm_size(dupComm,dupSize,ierror)
  call MPI_Comm_compare(mesh%decomp%mpiComm,dupComm,cmpResult,ierror)
  if(cmpResult /= MPI_IDENT .and. cmpResult /= MPI_CONGRUENT) then
    print*,"Error: mesh decomposition is not using the provided communicator"
    stop 1
  endif
  call MPI_Comm_compare(mesh%decomp%mpiComm,MPI_COMM_WORLD,cmpResult,ierror)
  if(cmpResult == MPI_IDENT) then
    print*,"Error: mesh decomposition fell back to MPI_COMM_WORLD instead of the provided communicator"
    stop 1
  endif
  if(mesh%decomp%nRanks /= dupSize .or. mesh%decomp%rankId /= dupRank) then
    print*,"Error: decomposition rank/size disagree with the provided communicator", &
      mesh%decomp%rankId,dupRank,mesh%decomp%nRanks,dupSize
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
  modelobj%c = c

  call modelobj%SphericalSoundWave(rhoprime,Lr,x0,y0,z0)

  call modelobj%SetTimeIntegrator(integrator)
  call modelobj%ForwardStep(endtime,dt,iointerval)

  ef = modelobj%entropy
  if(ef /= ef) then
    print*,"Error: Final entropy is inf or nan",ef
    stop 1
  endif

  call modelobj%free()
  call mesh%free() ! must NOT finalize MPI: the communicator is externally managed

  call MPI_Finalized(mpiIsFinalized,ierror)
  if(mpiIsFinalized) then
    print*,"Error: SELF finalized MPI despite an externally provided communicator"
    stop 1
  endif

  call geometry%free()
  call interp%free()

  call MPI_Comm_free(dupComm,ierror)
  call mpi_finalize(ierror)

endprogram structuredmesh_external_comm_3d
