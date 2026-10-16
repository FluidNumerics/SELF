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
!! Exercises SideExchange on extruded meshes (Mesh3D%ExtrudedMesh) built from the
!! curved Circle surface mesh, read both from HOHQMesh (Circle.mesh) and from HOPr
!! (Circle_mesh.h5), with terrain-following levels over a sloping bathymetry.
!!
!! A MappedScalar3D carries several linear functions of (x,y,z). On an extruded element
!! z is the interpolant of the levels, a polynomial of degree nGeo in (xi,eta) and
!! linear in zeta, and x and y are degree nGeo in (xi,eta), so with N >= nGeo each
!! function is interpolated exactly and its trace on a face is single valued. After
!! SideExchange the neighbor trace (extBoundary) must therefore equal the local trace
!! (boundary) node by node, to round-off, on every interior face.
!!
!! The two traces are only compared node by node correctly if the face flip recorded
!! for a lateral face is the right one; a wrong flip pairs mirror-image nodes and the
!! traces differ at O(1). The Circle surface has edges whose two sides run in
!! opposite directions (2-D flip 1), and so does a structured surface on which every
!! other element is turned through 180 degrees; the test asserts that such flipped
!! lateral faces are present and that the exchange is exact across them. On the
!! turned surface the flipped faces also straddle the two ranks of a parallel run,
!! which exercises the flip applied after an MPI exchange (ApplyFlip). Several variables
!! with unrelated coefficients are carried so that a message pairing the wrong
!! variables or faces shows up as a mismatch.
!!
!! This is the MPI variant. It uses the HOPr surface only, because the HOHQMesh reader
!! keeps every element on every rank instead of decomposing the mesh, and it asserts
!! that rank-remote faces, flipped ones among them, are exchanged.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_MappedScalar_3D
  use mpi

  implicit none
  integer :: exit_code

  exit_code = mappedscalarsideexchange_3d_extruded()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  logical function IsRotated(eGlobal) result(rotated)
    !! The surfaces built by RotateAlternateElements turn every element with an odd
    !! global id through 180 degrees.
    implicit none
    integer,intent(in) :: eGlobal

    rotated = modulo(eGlobal,2) == 1

  endfunction IsRotated

  subroutine RotateAlternateElements(surf)
    !! Turns every element of a conforming structured surface with an odd global id
    !! through 180 degrees in its reference coordinates, (xi,eta) -> (-xi,-eta). The
    !! geometry is unchanged, but the element's south and north sides, and its east
    !! and west sides, trade places, and each of its edges now runs in the opposite
    !! direction. An edge between a turned and an unturned element therefore gets
    !! flip 1, which a structured mesh never has on its own. Every rank derives the
    !! same orientation for every global element id, so the sides of an edge that
    !! straddles two ranks agree.
    implicit none
    type(Mesh2D),intent(inout) :: surf
    ! Local
    integer,parameter :: rot(1:4) = [3,4,1,2] ! side s of the original is side rot(s) after the turn
    integer :: e,eGlobal,s,sNew,e2,s2,s2New,flip,i,j,n
    integer :: si(1:5,1:4)
    real(prec),allocatable :: xy(:,:,:)

    n = surf%nGeo+1
    allocate(xy(1:2,1:n,1:n))
    do e = 1,surf%nElem
      eGlobal = surf%decomp%offsetElem(surf%decomp%rankId+1)+e
      si = surf%sideInfo(1:5,1:4,e)
      do s = 1,4
        sNew = s
        if(IsRotated(eGlobal)) sNew = rot(s)
        surf%sideInfo(1:5,sNew,e) = si(1:5,s)
        e2 = si(3,s)
        if(e2 > 0) then
          s2 = si(4,s)/10
          s2New = s2
          if(IsRotated(e2)) s2New = rot(s2)
          flip = 0
          if(IsRotated(eGlobal) .neqv. IsRotated(e2)) flip = 1
          surf%sideInfo(4,sNew,e) = 10*s2New+flip
        endif
      enddo
      if(IsRotated(eGlobal)) then
        xy = surf%nodeCoords(1:2,1:n,1:n,e)
        do j = 1,n
          do i = 1,n
            surf%nodeCoords(1:2,i,j,e) = xy(1:2,n+1-i,n+1-j)
          enddo
        enddo
      endif
    enddo
    deallocate(xy)
    call surf%UpdateDevice()

  endsubroutine RotateAlternateElements

  integer function CheckExchange(surf,label,requireRemoteFlip) result(r)
    implicit none
    type(Mesh2D),intent(in) :: surf
    character(*),intent(in) :: label
    logical,intent(in) :: requireRemoteFlip
    ! Local
    integer,parameter :: N = 7
    integer,parameter :: nvar = 4
    integer,parameter :: nz = 3
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-10_prec
#else
    real(prec),parameter :: tolerance = 1.0e-4_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: f
    real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
    integer :: nGeo,e,i,j,iel,iside,ivar,e2,s2,flip,r2,ierror
    integer :: counts(1:4),countsGlobal(1:4) ! interior, remote, flipped, remote flipped faces
    real(prec) :: x,y,maxError,maxErrorGlobal

    r = 0
    nGeo = surf%nGeo
    allocate(h(1:nGeo+1,1:nGeo+1,1:surf%nElem))
    allocate(zLevels(1:nGeo+1,1:nGeo+1,0:nz,1:surf%nElem))
    do e = 1,surf%nElem
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          x = surf%nodeCoords(1,i,j,e)
          y = surf%nodeCoords(2,i,j,e)
          h(i,j,e) = 1.0_prec+0.3_prec*x+0.2_prec*y*y
        enddo
      enddo
    enddo
    call ExtrudedSigmaLevels(h,[-1.0_prec,-0.55_prec,-0.2_prec,0.0_prec],zLevels)
    call mesh%ExtrudedMesh(surf,zLevels,SELF_BC_NONORMALFLOW,SELF_BC_PRESCRIBED)

    call interp%Init(N=N,controlNodeType=GAUSS,M=N+3,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nElem)
    call f%AssociateGeometry(geometry)
    call f%SetEquation(1,'f = x + 2.7*y + 5.3*z')
    call f%SetEquation(2,'f = 3.1*x - 1.3*y + 0.7*z')
    call f%SetEquation(3,'f = 0.9*x + 4.2*y - 2.5*z')
    call f%SetEquation(4,'f = -2.2*x + 0.4*y + 1.9*z')
    call f%SetInteriorFromEquation(geometry,0.0_prec)
    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%UpdateHost()

    maxError = 0.0_prec
    counts = 0
    do iel = 1,mesh%nElem
      do iside = 1,6
        e2 = mesh%sideInfo(3,iside,iel)
        if(e2 == 0) cycle
        s2 = mesh%sideInfo(4,iside,iel)/10
        flip = mesh%sideInfo(4,iside,iel)-10*s2
        r2 = mesh%decomp%elemToRank(e2)
        counts(1) = counts(1)+1
        if(r2 /= mesh%decomp%rankId) counts(2) = counts(2)+1
        if(flip /= 0) then
          counts(3) = counts(3)+1
          if(r2 /= mesh%decomp%rankId) counts(4) = counts(4)+1
        endif
        do ivar = 1,nvar
          do j = 1,N+1
            do i = 1,N+1
              maxError = max(maxError,abs(f%extBoundary(i,j,iside,iel,ivar)-f%boundary(i,j,iside,iel,ivar)))
            enddo
          enddo
        enddo
      enddo
    enddo

    call MPI_Allreduce(maxError,maxErrorGlobal,1,mesh%decomp%mpiPrec,MPI_MAX,mesh%decomp%mpiComm,ierror)
    call MPI_Allreduce(counts,countsGlobal,4,MPI_INTEGER,MPI_SUM,mesh%decomp%mpiComm,ierror)

    if(mesh%decomp%rankId == 0) then
      print*,label,": ",mesh%nGlobalElem," elements on ",mesh%decomp%nRanks," rank(s)"
      print*,"  interior face slots ",countsGlobal(1),", rank-remote ",countsGlobal(2)
      print*,"  flipped face slots  ",countsGlobal(3),", rank-remote ",countsGlobal(4)
      print*,"  max |extBoundary - boundary| = ",maxErrorGlobal," (tolerance ",tolerance,")"
    endif

    if(maxErrorGlobal > tolerance) then
      if(mesh%decomp%rankId == 0) print*,"FAIL (",label,"): neighbor trace does not match"
      r = 1
    endif
    if(countsGlobal(3) == 0) then
      if(mesh%decomp%rankId == 0) print*,"FAIL (",label,"): no flipped lateral face to exercise"
      r = 1
    endif
    if(mesh%decomp%nRanks > 1 .and. countsGlobal(2) == 0) then
      if(mesh%decomp%rankId == 0) print*,"FAIL (",label,"): no rank-remote face to exercise"
      r = 1
    endif
    if(requireRemoteFlip .and. mesh%decomp%nRanks > 1 .and. countsGlobal(4) == 0) then
      if(mesh%decomp%rankId == 0) print*,"FAIL (",label,"): no rank-remote flipped face to exercise"
      r = 1
    endif

    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call interp%Free()
    call mesh%Free()
    deallocate(h,zLevels)

  endfunction CheckExchange

  integer function mappedscalarsideexchange_3d_extruded() result(r)
    implicit none
    ! Local
    type(Mesh2D) :: turned,circleHOPr
    character(LEN=255) :: WORKSPACE

    call get_environment_variable("WORKSPACE",WORKSPACE)
    ! 8 x 4 elements in two tiles side by side, so that two ranks meet on a vertical
    ! grid line, where the ids of the two elements of each edge differ by an odd number
    call turned%StructuredMesh(4,4,2,1,0.25_prec,0.2_prec, &
                               [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED])
    call RotateAlternateElements(turned)
    call circleHOPr%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Circle/Circle_mesh.h5")

    r = 0
    r = max(r,CheckExchange(turned,"StructuredMesh, odd elements turned",.true.))
    r = max(r,CheckExchange(circleHOPr,"HOPr Circle_mesh.h5",.false.))

    call circleHOPr%Free()
    call turned%Free()

  endfunction mappedscalarsideexchange_3d_extruded

endprogram test
