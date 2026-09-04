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
!! Exercises the scalar side exchange -- and nothing else -- across ranks in 2-D.
!!
!! A MappedScalar2D is set from continuous analytical functions of (x,y). The
!! trace of a continuous function is single valued on an interface, so after
!! SideExchange the neighbor trace (extBoundary) must equal the local trace
!! (boundary) at every interior edge, to round-off. The functions are linear, so
!! the nodal interpolant is exact in both elements and the comparison is a
!! round-off level assertion rather than an accuracy one.
!!
!! What this covers that the DG tests do not: the exchange is checked directly,
!! rather than through a derivative that averages the two traces (and would
!! therefore be insensitive to the neighbor trace being wrong in a way that is
!! symmetric). Several variables are carried, each with a different function, so
!! that a message pairing the wrong variables shows up as a trace mismatch. That
!! is what covers the per-variable offset of the MPI message tag,
!! tag = |globalSideId| + nUniqueSides*(ivar-1): a stride that fails to separate
!! the variables lands one variable's trace in another's buffer, and the two
!! functions disagree.
!!
!! The mesh and variable count are also chosen so that the exchange posts more
!! messages than there are unique edges in the mesh (see the comment on the mesh
!! below): the decomposition sizes its non-blocking request scratch from the mesh
!! alone, and an exchange that posts 2*nvar messages per remote edge overruns it.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_MappedScalar_2D
  use mpi

  implicit none
  integer :: exit_code

  exit_code = scalarsideexchange_2d_mpi()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function scalarsideexchange_2d_mpi() result(r)

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 6
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-10)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-4)
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    integer :: bcids(1:4)
    integer :: iel,iside,ivar,i,e2,r2
    integer :: nRemoteSides,nRemoteSidesGlobal,ierror
    integer :: msgCountMax
    real(prec) :: maxError,maxErrorGlobal

    bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                  SELF_BC_PRESCRIBED, & ! East
                  SELF_BC_PRESCRIBED, & ! North
                  SELF_BC_PRESCRIBED] ! West

    ! A 3x3 mesh has 24 unique edges, which is what the decomposition sizes its
    ! request scratch from. Split over two ranks it has 4 rank-remote edges per
    ! rank, so the exchange posts 2*4*nvar = 48 messages -- twice what the scratch
    ! was sized for. dx and dy differ so that a trace taken from the wrong
    ! neighbor cannot coincide with the right one.
    call mesh%StructuredMesh(3,3,1,1,0.1_prec,0.13_prec,bcids)

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nElem)
    call f%AssociateGeometry(geometry)

    ! Continuous (indeed linear) analytical functions. The coefficients are
    ! chosen so that no variable is symmetric under a swap of the coordinate
    ! directions, and so that no two variables are multiples of each other.
    call f%SetEquation(1,'f = x + 2.7*y')
    call f%SetEquation(2,'f = 3.1*x - 1.3*y')
    call f%SetEquation(3,'f = 0.9*x + 4.2*y')
    call f%SetEquation(4,'f = -2.2*x + 0.4*y')
    call f%SetEquation(5,'f = 5.6*x + 1.1*y')
    call f%SetEquation(6,'f = 1.7*x - 3.8*y')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    call f%BoundaryInterp()
    call f%SideExchange(mesh)
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
#ifndef ENABLE_GPU
    ! Keep this a regression test for the request-scratch sizing: the point of the
    ! mesh and nvar above is that the exchange posts more messages than the mesh
    ! alone would have reserved room for.
    !
    ! CPU builds only. The GPU exchange packs all of a neighbor's edges into a
    ! single message and keeps its own persistent request array, so it neither
    ! fills decomp%msgCount nor has a per-side scratch to overrun; the trace
    ! assertion above is what covers the exchange there.
    !
    ! Reduce over the ranks rather than testing each one: which rank carries the
    ! overrun depends on how the elements were dealt out, and SELF_MPIEXEC_NUMPROCS
    ! is configurable, so a given rank may legitimately post fewer messages.
    if(mesh%decomp%mpiEnabled) then
      call MPI_Allreduce(mesh%decomp%msgCount,msgCountMax,1,MPI_INTEGER, &
                         MPI_MAX,mesh%decomp%mpiComm,ierror)
      if(msgCountMax <= mesh%nUniqueSides) then
        print*,"FAIL: the busiest rank posted ",msgCountMax," messages for ", &
          mesh%nUniqueSides," unique edges; the test no longer exercises the request scratch"
        r = 1
      endif
    endif
#endif

    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

  endfunction scalarsideexchange_2d_mpi
endprogram test
