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

program bench_lineareuler3d_scaling
!! Multi-process / multi-GPU scaling benchmark for the 3-D linear Euler solver.
!!
!! Solves the linear Euler equations (spherical sound wave initial condition,
!! radiation boundary conditions) on a uniform structured mesh of
!! nex x ney x nez elements with degree `cdegree` polynomials, using the
!! low-storage RK3 integrator. No file IO is performed; the program reports
!! wall-time per time step and degrees-of-freedom throughput so that runs at
!! 1, 2, 4, and 8 MPI ranks measure parallel scaling.
!!
!! The mesh is tiled so that each MPI rank owns exactly one cuboid tile
!! (SELF's linear element decomposition maps tiles to contiguous element
!! ranges). The number of tiles (ntx*nty*ntz) must equal the number of ranks,
!! and each global element count must be divisible by its tile count.
!!
!! Usage:
!!   bench_lineareuler3d_scaling [nex ney nez ntx nty ntz nsteps cdegree]
!!
!! Defaults (no arguments) run a small single-rank problem: 8^3 elements,
!! 1x1x1 tiles, 5 steps, degree 7.
!!
!! The time step is set from a CFL-like rule dt = cfl*dx/(c*(N+1)^2) with
!! cfl = 0.5, c = 1, so the run is stable for any resolution; a NaN check on
!! the final entropy guards against silent blowup.

  use iso_fortran_env,only:real64
  use self_data
  use self_LinearEuler3D
  use mpi

  implicit none

  real(prec),parameter :: rho0 = 1.225_prec ! Background density
  real(prec),parameter :: rhoprime = 0.01_prec ! Sound wave density anomaly
  real(prec),parameter :: c = 1.0_prec ! Speed of sound
  real(prec),parameter :: cfl = 0.5_prec ! CFL-like number for time step selection

  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  integer :: nex,ney,nez ! Global element counts in each direction
  integer :: ntx,nty,ntz ! Number of tiles (= ranks) in each direction
  integer :: nsteps ! Number of benchmark time steps
  integer :: nwarmup ! Number of untimed warm-up time steps
  integer :: cdegree ! Control polynomial degree
  integer :: istep
  integer :: rankid,nranks,ierror
  logical :: mpiInitialized
  real(prec) :: dx,dt
  real(real64) :: t1,t2,tstep
  real(real64) :: ndof
  character(len=32) :: arg

  nex = 8
  ney = 8
  nez = 8
  ntx = 1
  nty = 1
  ntz = 1
  nsteps = 5
  nwarmup = 2
  cdegree = 7

  if(command_argument_count() >= 6) then
    call get_command_argument(1,arg); read(arg,*) nex
    call get_command_argument(2,arg); read(arg,*) ney
    call get_command_argument(3,arg); read(arg,*) nez
    call get_command_argument(4,arg); read(arg,*) ntx
    call get_command_argument(5,arg); read(arg,*) nty
    call get_command_argument(6,arg); read(arg,*) ntz
  endif
  if(command_argument_count() >= 7) then
    call get_command_argument(7,arg); read(arg,*) nsteps
  endif
  if(command_argument_count() >= 8) then
    call get_command_argument(8,arg); read(arg,*) cdegree
  endif

  if(mod(nex,ntx) /= 0 .or. mod(ney,nty) /= 0 .or. mod(nez,ntz) /= 0) then
    print*,'bench_lineareuler3d_scaling: element counts must be divisible by tile counts'
    stop 1
  endif

  ! Radiation boundary conditions on all sides
  bcids(1:6) = [SELF_BC_RADIATION,SELF_BC_RADIATION,SELF_BC_RADIATION, &
                SELF_BC_RADIATION,SELF_BC_RADIATION,SELF_BC_RADIATION]

  ! Uniform element width for a unit cube domain in x; dy,dz chosen equal to dx
  ! so that elements are isotropic.
  dx = 1.0_prec/real(nex,prec)

  ! Create a structured mesh with one tile per rank. Mesh generation
  ! initializes the domain decomposition (and MPI).
  call mesh%StructuredMesh(nex/ntx,ney/nty,nez/ntz, &
                           ntx,nty,ntz,dx,dx,dx,bcids)

  rankid = mesh%decomp%rankId
  nranks = mesh%decomp%nRanks

  if(nranks /= ntx*nty*ntz) then
    if(rankid == 0) then
      print*,'bench_lineareuler3d_scaling: number of tiles must equal number of ranks: ', &
        ntx*nty*ntz,nranks
    endif
    stop 1
  endif

  call interp%Init(N=cdegree, &
                   controlNodeType=GAUSS, &
                   M=cdegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%tecplot_enabled = .false.
  modelobj%rho0 = rho0
  modelobj%c = c

  ! Spherical sound wave centered in the unit cube; length scale set relative
  ! to the domain so the initial condition is resolution independent.
  call modelobj%SphericalSoundWave(rhoprime,0.06_prec,0.5_prec,0.5_prec,0.5_prec)

  call modelobj%SetTimeIntegrator('rk3')

  ! Stable explicit time step for degree N spectral elements
  dt = cfl*dx/(c*real((cdegree+1)*(cdegree+1),prec))
  modelobj%dt = dt

  if(rankid == 0) then
    print*,'bench_lineareuler3d_scaling : global elements :',nex,ney,nez
    print*,'bench_lineareuler3d_scaling : tiles           :',ntx,nty,ntz
    print*,'bench_lineareuler3d_scaling : ranks           :',nranks
    print*,'bench_lineareuler3d_scaling : degree          :',cdegree
    print*,'bench_lineareuler3d_scaling : local elements  :',mesh%nElem
    print*,'bench_lineareuler3d_scaling : dt              :',dt
    print*,'bench_lineareuler3d_scaling : nsteps          :',nsteps
  endif

  ! Untimed warm-up steps (kernel compilation, MPI wire-up, allocator warmup)
  do istep = 1,nwarmup
    call modelobj%timeIntegrator(modelobj%t+dt)
  enddo

  ! Synchronize all ranks; CalculateEntropy performs a device synchronization
  ! (device-to-host copy) and an MPI reduction, so all warm-up work is complete
  ! before the timers start.
  call modelobj%CalculateEntropy()
  call MPI_Barrier(mesh%decomp%mpiComm,ierror)

  t1 = MPI_Wtime()
  do istep = 1,nsteps
    call modelobj%timeIntegrator(modelobj%t+dt)
  enddo
  ! CalculateEntropy forces completion of all device work and includes the
  ! final MPI reduction; its cost is amortized over nsteps.
  call modelobj%CalculateEntropy()
  t2 = MPI_Wtime()
  call MPI_Barrier(mesh%decomp%mpiComm,ierror)

  tstep = (t2-t1)/real(nsteps,real64)
  ndof = real(nex,real64)*real(ney,real64)*real(nez,real64)* &
         real((cdegree+1)**3,real64)

  if(rankid == 0) then
    print*,'bench_lineareuler3d_scaling : total time (s)          :',t2-t1
    print*,'bench_lineareuler3d_scaling : time per step (s)       :',tstep
    print*,'bench_lineareuler3d_scaling : DOF                     :',ndof
    print*,'bench_lineareuler3d_scaling : DOF/s (per prognostic)  :',ndof/tstep
  endif

  if(modelobj%entropy /= modelobj%entropy) then
    print*,'Error: Final entropy is inf or nan',modelobj%entropy
    stop 1
  endif

  if(rankid == 0) then
    print*,'bench_lineareuler3d_scaling : final entropy           :',modelobj%entropy
  endif

  call modelobj%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram bench_lineareuler3d_scaling
