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

program read_hohqmesh_mpi
!! The distributed reader must preserve the serial mesh geometry and exchange
!! continuous linear traces across rank boundaries, including oriented faces.
!! The caller owns MPI; separate two-rank groups also exercise a subcommunicator.
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_MappedScalar_3D
  use mpi
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite

  implicit none
  integer :: comm,groupComm,rank,ierror

  call MPI_Init(ierror)
  call MPI_Comm_dup(MPI_COMM_WORLD,comm,ierror)
  call MPI_Comm_rank(comm,rank,ierror)
  call CheckMesh("MultiMaterial3D/ReaderFixture_ManyMaterials.mesh",comm)
  call CheckMesh("MultiMaterial3D/InsulatedWire3D.mesh",comm)
  call MPI_Comm_split(comm,rank/2,rank,groupComm,ierror)
  call CheckMesh("MultiMaterial3D/ReaderFixture_RotatedStack.mesh",groupComm)
  call MPI_Comm_free(groupComm,ierror)
  call MPI_Comm_free(comm,ierror)
  call MPI_Finalize(ierror)

contains

  subroutine CheckMesh(fixture,comm)
    !! Compares each owned element with the serial reader, then checks that
    !! independent linear fields have single-valued traces after SideExchange.
    implicit none
    character(*),intent(in) :: fixture
    integer,intent(in) :: comm
    ! Local
    integer,parameter :: degree = 7,nvar = 3
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-10_prec
#else
    real(prec),parameter :: tolerance = 2.0e-4_prec
#endif
    type(Mesh3D),target :: mesh,reference
    type(Lagrange),target :: interp
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: field
    character(LEN=1024) :: workspace
    integer :: rank,nRanks,ierror,comparison,offset,nOwned,nTotal,e,g,s,e2
    integer :: nRemote,nRemoteGlobal
    real(prec) :: error,errorGlobal

    call MPI_Comm_rank(comm,rank,ierror)
    call MPI_Comm_size(comm,nRanks,ierror)
    call get_environment_variable("WORKSPACE",workspace)
    call reference%Read_HOHQMesh(trim(workspace)//"/share/mesh/"//fixture,comm=MPI_COMM_SELF)
    if(reference%nElem < nRanks) then
      call reference%Free()
      return
    endif
    call mesh%Read_HOHQMesh(trim(workspace)//"/share/mesh/"//fixture,comm=comm)
    call MPI_Comm_compare(mesh%decomp%mpiComm,comm,comparison,ierror)
    if(comparison /= MPI_IDENT) call Fail("reader did not retain the caller's communicator",comm)

    offset = mesh%decomp%offsetElem(rank+1)
    nOwned = mesh%decomp%offsetElem(rank+2)-offset
    if(mesh%nElem /= nOwned) call Fail("mesh retains elements owned by other ranks",comm)
    call MPI_Allreduce(mesh%nElem,nTotal,1,MPI_INTEGER,MPI_SUM,comm,ierror)
    if(nTotal /= reference%nElem) call Fail("elements are not owned exactly once",comm)
    if(mesh%nGlobalElem /= reference%nElem) call Fail("global element count changed",comm)
    if(mesh%nUniqueSides /= reference%nUniqueSides) call Fail("global side numbering changed",comm)
    if(mesh%nMaterials /= reference%nMaterials .or. mesh%nBCs /= reference%nBCs) &
      call Fail("material or boundary table size changed",comm)
    if(any(mesh%materialNames /= reference%materialNames)) call Fail("material names changed",comm)
    if(any(mesh%BCNames /= reference%BCNames)) call Fail("boundary names changed",comm)
    if(any(mesh%BCType /= reference%BCType)) call Fail("boundary types changed",comm)

    error = 0.0_prec
    do e = 1,mesh%nElem
      g = offset+e
      if(mesh%elemMaterial(e) /= reference%elemMaterial(g)) call Fail("material assignment changed",comm)
      if(any(mesh%elemInfo(:,e) /= reference%elemInfo(:,g))) call Fail("element metadata changed",comm)
      if(any(mesh%sideInfo(:,:,e) /= reference%sideInfo(:,:,g))) call Fail("side connectivity changed",comm)
      if(any(mesh%globalNodeIDs(:,:,:,e) /= reference%globalNodeIDs(:,:,:,g))) &
        call Fail("corner node identifiers changed",comm)
      error = max(error,maxval(abs(mesh%nodeCoords(:,:,:,:,e)-reference%nodeCoords(:,:,:,:,g))))
    enddo
    call MPI_Allreduce(error,errorGlobal,1,mesh%decomp%mpiPrec,MPI_MAX,comm,ierror)
    if(.not. ieee_is_finite(errorGlobal) .or. errorGlobal > tolerance) &
      call Fail("distributed geometry differs from serial geometry",comm)
    call reference%Free()

    call interp%Init(N=degree,controlNodeType=GAUSS,M=degree,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)
    call field%Init(interp,nvar,mesh%nElem)
    call field%AssociateGeometry(geometry)
    call field%SetEquation(1,'f = x + 2.7*y + 5.3*z')
    call field%SetEquation(2,'f = 3.1*x - 1.3*y + 0.7*z')
    call field%SetEquation(3,'f = 0.9*x + 4.2*y - 2.5*z')
    call field%SetInteriorFromEquation(geometry,0.0_prec)
    call field%BoundaryInterp()
    call field%SideExchange(mesh)
    call field%UpdateHost()

    nRemote = 0
    error = 0.0_prec
    do e = 1,mesh%nElem
      do s = 1,6
        e2 = mesh%sideInfo(3,s,e)
        if(e2 == 0) cycle
        if(mesh%decomp%elemToRank(e2) /= rank) nRemote = nRemote+1
        if(any(.not. ieee_is_finite(field%extBoundary(:,:,s,e,:)))) &
          call Fail("neighbor trace is not finite",comm)
        error = max(error,maxval(abs(field%extBoundary(:,:,s,e,:)-field%boundary(:,:,s,e,:))))
      enddo
    enddo
    call MPI_Allreduce(error,errorGlobal,1,mesh%decomp%mpiPrec,MPI_MAX,comm,ierror)
    call MPI_Allreduce(nRemote,nRemoteGlobal,1,MPI_INTEGER,MPI_SUM,comm,ierror)
    if(nRanks > 1 .and. nRemoteGlobal == 0) call Fail("no rank-remote face was exercised",comm)
    if(.not. ieee_is_finite(errorGlobal) .or. errorGlobal > tolerance) call Fail("linear traces disagree",comm)
    if(rank == 0) print*,"PASS ",fixture," ranks=",nRanks," elements=",nTotal, &
      " remote_faces=",nRemoteGlobal," trace_error=",errorGlobal

    call field%DissociateGeometry()
    call field%Free()
    call geometry%Free()
    call interp%Free()
    call mesh%Free()
  endsubroutine CheckMesh

  subroutine Fail(message,comm)
    implicit none
    character(*),intent(in) :: message
    integer,intent(in) :: comm
    ! Local
    integer :: ierror

    print*,"FAIL: ",message
    call MPI_Abort(comm,1,ierror)
  endsubroutine Fail

endprogram read_hohqmesh_mpi
