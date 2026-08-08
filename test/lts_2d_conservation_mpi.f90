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

program LTS_2D_Conservation_MPI
!! Multi-rank local time stepping (AMR Stage 7) across 2:1 level interfaces.
!!
!! This is the check that the LTS recursion is collective in the right way. Element lists
!! are rank-local, but the mortar table, mortarLevel and hence the per-level mortar lists
!! are replicated, so every rank walks the SAME recursion and reaches every MortarExchange
!! and MortarFluxCollect at the same point - including ranks that own no element at the
!! level currently being advanced, which must still participate so that the sends posted by
!! the other side of a split interface are matched. If that were wrong the run would either
!! deadlock or mismatch messages, and if the per-rank register bookkeeping were wrong the
!! conservation check below would drift.
!!
!! DoubleMortarMesh is used precisely because its two-rank decomposition is adversarial by
!! construction (see its docstring): mortar 1 splits the big side from its small sides across
!! the rank boundary, while mortar 2 puts BOTH small elements on a rank remote from the big
!! element's, so the big-side trace is received independently per sub-edge. One of its small
!! sides also carries flip = 1, so the trace reorientation is exercised at the same time.
!! Three-level nesting of the recursion is covered serially by lts_2d_conservation.
!!
!! Assertions:
!!  (a) the run completes on >= 2 ranks (no deadlock, no message mismatch);
!!  (b) the GLOBAL domain integral of pressure is conserved to round-off - the rank-local
!!      integrals are summed with an allreduce, so this is a statement about the whole
!!      domain including the interfaces that straddle a rank boundary;
!!  (c) entropy stays finite and non-increasing, and the solution is NaN-free.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_Model
  use SELF_LocalTimeStepping_2D
  use mpi

  implicit none

  integer,parameter :: controlDegree = 5
  integer,parameter :: targetDegree = 10
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-4_prec
  real(prec),parameter :: Lr = 0.05_prec
  real(prec),parameter :: dtBase = 5.0e-4_prec
  real(prec),parameter :: endtime = 0.01_prec
  real(prec),parameter :: iointerval = 0.01_prec

  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  type(LTSSchedule2D) :: sched
  integer :: bcids(1:4)
  real(prec) :: m0,mf,e0,ef,scale,defect

  bcids(1:4) = [SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW, &
                SELF_BC_NONORMALFLOW,SELF_BC_NONORMALFLOW]

  ! Six elements, two 2:1 mortars; big elements at level 0, small elements at level 1.
  call mesh%DoubleMortarMesh(dx,bcids)

  if(mesh%maxElemLevel /= 1 .or. mesh%nMortars /= 2) then
    print*,"Error: expected a two-level mesh with two mortars, got level cap ", &
      mesh%maxElemLevel," and ",mesh%nMortars," mortars."
    stop 1
  endif

  call interp%Init(N=controlDegree,controlNodeType=GAUSS, &
                   M=targetDegree,targetNodeType=UNIFORM)
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%tecplot_enabled = .false.
  modelobj%rho0 = rho0

  ! Centred on the first mortar so the interface carries a real signal from step one.
  call modelobj%SphericalSoundWave(amp,Lr,2.0_prec*dx,dx,c0)

  m0 = GlobalPressureIntegral()
  call modelobj%CalculateEntropy()
  e0 = modelobj%entropy

  call ForwardStepLTS(modelobj,sched,tn=endtime,dtBase=dtBase,ioInterval=iointerval)

  mf = GlobalPressureIntegral()
  call modelobj%CalculateEntropy()
  ef = modelobj%entropy

  if(ef /= ef) then
    print*,"Error: entropy is NaN after multi-rank local time stepping."
    stop 1
  endif
  if(ef > e0) then
    print*,"Error: entropy grew under multi-rank local time stepping :",e0,ef
    stop 1
  endif

  scale = amp*c0*c0*(12.0_prec*dx*dx)
  defect = abs(mf-m0)/scale
  if(mesh%decomp%rankId == 0) then
    print*,"LTS MPI conservation : integral before, after, relative defect :",m0,mf,defect
  endif
  if(defect > 1.0e-10_prec) then
    print*,"Error: the global pressure integral drifted under multi-rank local time ", &
      "stepping across the level interfaces."
    stop 1
  endif

  call sched%Free()
  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

contains

  real(prec) function GlobalPressureIntegral() result(total)
    !! Jacobian-and-quadrature-weighted integral of the pressure over the WHOLE domain:
    !! summed locally, then reduced, so a defect hiding on one side of a rank boundary
    !! cannot cancel against the other side's.
    implicit none
    ! Local
    integer :: iel,i,j,iError
    real(prec) :: localSum

    localSum = 0.0_prec
    do iel = 1,geometry%nelem
      do j = 1,modelobj%solution%interp%N+1
        do i = 1,modelobj%solution%interp%N+1
          localSum = localSum+modelobj%solution%interior(i,j,iel,3)* &
                     abs(geometry%J%interior(i,j,iel,1))* &
                     modelobj%solution%interp%qWeights(i)* &
                     modelobj%solution%interp%qWeights(j)
        enddo
      enddo
    enddo

    if(mesh%decomp%mpiEnabled) then
      call mpi_allreduce(localSum,total,1,mesh%decomp%mpiPrec,MPI_SUM, &
                         mesh%decomp%mpiComm,iError)
    else
      total = localSum
    endif

  endfunction GlobalPressureIntegral

endprogram LTS_2D_Conservation_MPI
