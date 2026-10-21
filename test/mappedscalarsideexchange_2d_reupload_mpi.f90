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
!! Regression test for FluidNumerics/SELF#201: UpdateDevice must invalidate the
!! "static traces already exchanged" state of the aggregated GPU halo exchange.
!!
!! SideExchangeStart(mesh,nactive) exchanges only the leading nactive variables
!! once a full-width exchange has carried every variable; the trailing
!! (static, non-stepped) variables are assumed to keep valid neighbor traces in
!! extboundary_gpu. UpdateDevice overwrites extboundary_gpu with the host
!! extBoundary, so after an upload the next exchange must again carry every
!! variable, or the static variables keep the host values on every rank-cut edge.
!!
!! The test reproduces the sequence a model follows when it resets its initial
!! condition after a first tendency evaluation:
!!
!!   1. set six variables from linear functions, interpolate to the boundaries,
!!      and run a full side exchange (this marks the static traces as carried);
!!   2. set six different linear functions on the host, fill the host
!!      extBoundary with a sentinel that no valid trace can equal, and upload;
!!   3. interpolate to the boundaries and exchange with nactive = 3.
!!
!! The functions are continuous, so on every interior edge the neighbor trace
!! (extBoundary) must equal the local trace (boundary) for all six variables, to
!! round-off. Without the fix, variables 4-6 keep the sentinel on rank-cut edges.
!! Same-rank edges are filled by the local exchange kernel for every variable
!! either way, so the test needs at least two ranks to bite; it asserts that some
!! edges are rank-remote when run on more than one rank.
!!
!! On CPU builds there is no prefix exchange; the test runs a plain SideExchange
!! in step 3 and checks the same assertion.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_MappedScalar_2D
  use mpi

  implicit none
  integer :: exit_code

  exit_code = mappedscalarsideexchange_2d_reupload_mpi()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mappedscalarsideexchange_2d_reupload_mpi() result(r)

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 6
    integer,parameter :: nstepped = 3
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-10)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-4)
#endif
    real(prec),parameter :: sentinel = 10.0_prec**(30)
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    integer :: bcids(1:4)
    integer :: iel,iside,ivar,i,e2,r2
    integer :: nRemoteSides,nRemoteSidesGlobal,ierror
    real(prec) :: maxError,maxErrorGlobal

    bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                  SELF_BC_PRESCRIBED, & ! East
                  SELF_BC_PRESCRIBED, & ! North
                  SELF_BC_PRESCRIBED] ! West

    ! Split over two ranks, a 3x3 mesh has rank-remote edges on every rank.
    ! dx and dy differ so that a trace taken from the wrong neighbor cannot
    ! coincide with the right one.
    call mesh%StructuredMesh(3,3,1,1,0.1_prec,0.13_prec,bcids)

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nElem)
    call f%AssociateGeometry(geometry)

    ! Step 1: a full-width exchange, which carries the static variables once.
    call f%SetEquation(1,'f = x + 2.7*y')
    call f%SetEquation(2,'f = 3.1*x - 1.3*y')
    call f%SetEquation(3,'f = 0.9*x + 4.2*y')
    call f%SetEquation(4,'f = -2.2*x + 0.4*y')
    call f%SetEquation(5,'f = 5.6*x + 1.1*y')
    call f%SetEquation(6,'f = 1.7*x - 3.8*y')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    call f%BoundaryInterp()
    call f%SideExchange(mesh)

    ! Step 2: new data for every variable, uploaded from the host together with a
    ! host extBoundary that holds no valid trace.
    call f%SetEquation(1,'f = 1.5 - 0.6*x + 3.3*y')
    call f%SetEquation(2,'f = -0.8 + 2.4*x + 0.5*y')
    call f%SetEquation(3,'f = 2.1 + 1.2*x - 2.9*y')
    call f%SetEquation(4,'f = -1.9 + 4.4*x + 2.2*y')
    call f%SetEquation(5,'f = 0.7 - 3.6*x + 1.8*y')
    call f%SetEquation(6,'f = 3.2 + 0.3*x - 4.1*y')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    f%extBoundary = sentinel
    call f%UpdateDevice()

    ! Step 3: exchange only the leading nstepped variables, as a model's tendency
    ! evaluation does. The first exchange after an upload must still carry all of
    ! them.
    call f%BoundaryInterp()
#ifdef ENABLE_GPU
    call f%SideExchangeStart(mesh,nstepped)
    call f%SideExchangeFinish(mesh)
#else
    call f%SideExchange(mesh)
#endif
    call f%UpdateHost()

    ! Compare the two traces of the interface. Domain boundaries are skipped:
    ! extBoundary is untouched there (it is the boundary condition's to set).
    maxError = 0.0_prec
    nRemoteSides = 0
    do iel = 1,mesh%nElem
      do iside = 1,4

        e2 = mesh%sideInfo(3,iside,iel) ! Neighboring element (global id)
        if(e2 == 0) cycle

        r2 = mesh%decomp%elemToRank(e2)
        if(r2 /= mesh%decomp%rankId) nRemoteSides = nRemoteSides+1

        do ivar = 1,nvar
          do i = 1,controlDegree+1
            maxError = max(maxError, &
                           abs(f%extBoundary(i,iside,iel,ivar)-f%boundary(i,iside,iel,ivar)))
          enddo
        enddo

      enddo
    enddo

    ! Reduce before tearing the mesh down: freeing the mesh releases MPI.
    if(mesh%decomp%mpiEnabled) then
      call MPI_Allreduce(maxError,maxErrorGlobal,1,mesh%decomp%mpiPrec, &
                         MPI_MAX,mesh%decomp%mpiComm,ierror)
      call MPI_Allreduce(nRemoteSides,nRemoteSidesGlobal,1,MPI_INTEGER, &
                         MPI_SUM,mesh%decomp%mpiComm,ierror)
    else
      maxErrorGlobal = maxError
      nRemoteSidesGlobal = nRemoteSides
    endif

    print*,"rank ",mesh%decomp%rankId," : n remote edges = ",nRemoteSides, &
      " max |extBoundary-boundary| = ",maxError

    r = 0
    if(maxErrorGlobal > tolerance) then
      print*,"FAIL: max |extBoundary-boundary| ",maxErrorGlobal," exceeds tolerance ",tolerance
      r = 1
    endif
    ! With more than one rank the mesh above always straddles a rank boundary. If
    ! no edge is remote the exchange under test never ran, and the assertion
    ! above passed vacuously.
    if(mesh%decomp%nRanks > 1 .and. nRemoteSidesGlobal == 0) then
      print*,"FAIL: no edges were exchanged over MPI on ",mesh%decomp%nRanks," ranks"
      r = 1
    endif

    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

  endfunction mappedscalarsideexchange_2d_reupload_mpi
endprogram test
