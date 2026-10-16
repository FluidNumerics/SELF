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
!! Verifies that the halo exchange graph contains exactly the rank-remote
!! neighbors present in the mesh connectivity.
  use SELF_Constants
  use SELF_Mesh_3D
  use mpi

  implicit none
  integer :: exit_code

  exit_code = domaindecomposition_graph_mpi()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function domaindecomposition_graph_mpi() result(r)
    implicit none
    type(Mesh3D) :: mesh
    integer :: bcids(1:6)
    integer :: e1,s1,e2,neighborRank
    integer :: i,ierror,nExpected
    integer :: graphSize,indegree,outdegree
    logical :: weighted
    integer :: failures,globalFailures
    logical,allocatable :: expected(:)
    integer,allocatable :: sources(:),destinations(:)
    integer,allocatable :: sourceWeights(:),destinationWeights(:)

    bcids(1:6) = SELF_BC_PRESCRIBED
    call mesh%StructuredMesh(2,2,2,1,1,1,0.1_prec,0.13_prec,0.17_prec,bcids)
    call mesh%decomp%EnsureExchangeGraph(mesh%sideInfo,mesh%nElem,6, &
                                         mesh%mortarInfo,mesh%nMortars,4)

    allocate(expected(0:mesh%decomp%nRanks-1))
    expected = .false.
    do e1 = 1,mesh%nElem
      do s1 = 1,6
        e2 = mesh%sideInfo(3,s1,e1)
        if(e2 > 0) then
          neighborRank = mesh%decomp%elemToRank(e2)
          if(neighborRank /= mesh%decomp%rankId) expected(neighborRank) = .true.
        endif
      enddo
    enddo

    nExpected = count(expected)
    failures = 0
    if(mesh%decomp%nExchangeNeighbors /= nExpected) then
      print*,"FAIL: graph neighbor count ",mesh%decomp%nExchangeNeighbors, &
        " does not match remote-side count ",nExpected
      failures = failures+1
    endif
    if(allocated(mesh%decomp%exchangeNeighbors)) then
      do i = 1,mesh%decomp%nExchangeNeighbors
        if(.not. expected(mesh%decomp%exchangeNeighbors(i))) then
          print*,"FAIL: graph contains unexpected rank ",mesh%decomp%exchangeNeighbors(i)
          failures = failures+1
        endif
        if(i > 1) then
          if(mesh%decomp%exchangeNeighbors(i-1) >= mesh%decomp%exchangeNeighbors(i)) then
            print*,"FAIL: graph neighbor ranks are not strictly increasing"
            failures = failures+1
          endif
        endif
      enddo
    else
      print*,"FAIL: graph neighbor table was not allocated"
      failures = failures+1
    endif

    call MPI_Comm_size(mesh%decomp%exchangeComm,graphSize,ierror)
    if(ierror /= MPI_SUCCESS) then
      print*,"FAIL: MPI_Comm_size returned ",ierror
      failures = failures+1
    endif
    if(graphSize /= mesh%decomp%nRanks) then
      print*,"FAIL: graph communicator size ",graphSize," differs from world size ",mesh%decomp%nRanks
      failures = failures+1
    endif
    call MPI_Dist_graph_neighbors_count(mesh%decomp%exchangeComm,indegree, &
                                        outdegree,weighted,ierror)
    if(ierror /= MPI_SUCCESS) then
      print*,"FAIL: MPI_Dist_graph_neighbors_count returned ",ierror
      failures = failures+1
    endif
    if(indegree /= nExpected .or. outdegree /= nExpected .or. weighted) then
      print*,"FAIL: MPI graph degree mismatch: in=",indegree," out=",outdegree, &
        " expected=",nExpected," weighted=",weighted
      failures = failures+1
    endif

    allocate(sources(1:max(1,indegree)),destinationWeights(1:max(1,outdegree)))
    allocate(destinations(1:max(1,outdegree)),sourceWeights(1:max(1,indegree)))
    sourceWeights = 1
    destinationWeights = 1
    call MPI_Dist_graph_neighbors(mesh%decomp%exchangeComm,indegree,sources, &
                                  sourceWeights,outdegree,destinations, &
                                  destinationWeights,ierror)
    if(ierror /= MPI_SUCCESS) then
      print*,"FAIL: MPI_Dist_graph_neighbors returned ",ierror
      failures = failures+1
    endif
    if(indegree > 0) then
      do i = 1,indegree
        if(sources(i) /= mesh%decomp%exchangeNeighbors(i)) then
          print*,"FAIL: MPI source neighbor ",sources(i)," differs from SELF rank ", &
            mesh%decomp%exchangeNeighbors(i)
          failures = failures+1
        endif
      enddo
    endif
    if(outdegree > 0) then
      do i = 1,outdegree
        if(destinations(i) /= mesh%decomp%exchangeNeighbors(i)) then
          print*,"FAIL: MPI destination neighbor ",destinations(i)," differs from SELF rank ", &
            mesh%decomp%exchangeNeighbors(i)
          failures = failures+1
        endif
      enddo
    endif

    call MPI_Allreduce(failures,globalFailures,1,MPI_INTEGER,MPI_SUM, &
                       mesh%decomp%mpiComm,ierror)
    if(ierror /= MPI_SUCCESS) then
      print*,"FAIL: final MPI_Allreduce returned ",ierror
      globalFailures = globalFailures+1
    endif
    r = merge(0,1,globalFailures == 0)
    if(r /= 0) then
      print*,"FAIL: domain-decomposition graph validation failed on ",globalFailures," checks"
    else if(mesh%decomp%rankId == 0) then
      print*,"PASS: domain-decomposition graph validation neighbors=",nExpected, &
        " ranks=",mesh%decomp%nRanks
    endif
    if(mesh%decomp%rankId == 0) print*,"RESULT: domaindecomposition_graph_mpi r=",r
    flush(6)

    deallocate(expected,sources,destinations,sourceWeights,destinationWeights)
    call mesh%Free()

  endfunction domaindecomposition_graph_mpi
endprogram test
