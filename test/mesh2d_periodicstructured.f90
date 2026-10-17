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
!! Verifies the connectivity of the per-axis periodic structured 2-D mesh
!! (Mesh2D%PeriodicStructuredMesh) for all four combinations of the periodic flags.
!!
!! Every element is rank-local in serial, so the pairing of the two sides of each
!! edge can be checked directly. For each mesh the test asserts that
!!
!!   1. nUniqueSides is nLx*nY + nLy*nX, with nLx = nX in a periodic x direction and
!!      nX+1 otherwise (likewise nLy), and the global side ids are dense in
!!      1..nUniqueSides: an interior (or periodic) id is used by exactly two side slots
!!      and a boundary id by exactly one,
!!   2. the neighbor of an interior side names this element and side back, carries the
!!      same global id, and both flips are 0,
!!   3. the two sides of an edge coincide geometrically, up to a shift by the domain
!!      length Lx or Ly across a periodic wrap, so the wiring joins the right elements,
!!   4. boundary sides occur only on the walls of a non-periodic direction and carry
!!      that wall's boundary id,
!!   5. with neither direction periodic, the neighbors and boundary ids are identical to
!!      those of Mesh2D%StructuredMesh.

  use SELF_Constants
  use SELF_Mesh_2D

  implicit none
  integer :: exit_code

  exit_code = mesh2d_periodicstructured()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function CheckPeriodicMesh(periodicX,periodicY) result(r)
    implicit none
    logical,intent(in) :: periodicX
    logical,intent(in) :: periodicY
    ! Local
    integer,parameter :: nxPerTile = 3,nyPerTile = 2,nTileX = 2,nTileY = 2
    integer,parameter :: nX = nxPerTile*nTileX,nY = nyPerTile*nTileY
    real(prec),parameter :: dx = 0.25_prec,dy = 0.4_prec
    real(prec),parameter :: Lx = dx*real(nX,prec),Ly = dy*real(nY,prec)
    real(prec),parameter :: tol = 10.0_prec*epsilon(1.0_prec)
    integer,parameter :: bcids(1:4) = [11,12,13,14]
    type(Mesh2D) :: mesh,reference
    integer :: e1,s1,e2,s2,g,nLx,nLy,nBoundary,expectedBoundary
    integer,allocatable :: idcount(:)
    real(prec) :: a1(1:2),b1(1:2),a2(1:2),b2(1:2),shift(1:2)
    logical :: ok
    character(len=64) :: label

    r = 1
    write(label,'(A,L1,A,L1)') "periodicX=",periodicX,", periodicY=",periodicY

    call mesh%PeriodicStructuredMesh(nxPerTile,nyPerTile,nTileX,nTileY,dx,dy, &
                                     periodicX,periodicY,bcids)

    nLx = nX+1
    if(periodicX) nLx = nX
    nLy = nY+1
    if(periodicY) nLy = nY

    if(mesh%nElem /= nX*nY) then
      print*,"FAIL (",trim(label),"): nElem ",mesh%nElem," /= ",nX*nY
      return
    endif
    if(mesh%nUniqueSides /= nLx*nY+nLy*nX) then
      print*,"FAIL (",trim(label),"): nUniqueSides ",mesh%nUniqueSides," /= ",nLx*nY+nLy*nX
      return
    endif

    allocate(idcount(1:mesh%nUniqueSides))
    idcount = 0
    nBoundary = 0

    do e1 = 1,mesh%nElem
      do s1 = 1,4
        g = mesh%sideInfo(2,s1,e1)
        if(g < 1 .or. g > mesh%nUniqueSides) then
          print*,"FAIL (",trim(label),"): global side id ",g," out of range on element ",e1," side ",s1
          return
        endif
        idcount(g) = idcount(g)+1

        call EdgeEndpoints(mesh,e1,s1,a1,b1)
        e2 = mesh%sideInfo(3,s1,e1)

        if(e2 == 0) then
          nBoundary = nBoundary+1
          ! Only a wall of a non-periodic direction may be a boundary
          if(((s1 == 2 .or. s1 == 4) .and. periodicX) .or. &
             ((s1 == 1 .or. s1 == 3) .and. periodicY)) then
            print*,"FAIL (",trim(label),"): boundary side in a periodic direction, element ",e1," side ",s1
            return
          endif
          if(mesh%sideInfo(5,s1,e1) /= bcids(s1)) then
            print*,"FAIL (",trim(label),"): boundary id ",mesh%sideInfo(5,s1,e1)," /= ",bcids(s1)
            return
          endif
          ok = .false.
          if(s1 == 1) ok = abs(a1(2)) < tol .and. abs(b1(2)) < tol
          if(s1 == 2) ok = abs(a1(1)-Lx) < tol .and. abs(b1(1)-Lx) < tol
          if(s1 == 3) ok = abs(a1(2)-Ly) < tol .and. abs(b1(2)-Ly) < tol
          if(s1 == 4) ok = abs(a1(1)) < tol .and. abs(b1(1)) < tol
          if(.not. ok) then
            print*,"FAIL (",trim(label),"): boundary side is not on its domain wall, element ",e1," side ",s1
            return
          endif
          cycle
        endif

        s2 = mesh%sideInfo(4,s1,e1)/10
        if(mesh%sideInfo(4,s1,e1)-10*s2 /= 0) then
          print*,"FAIL (",trim(label),"): non-zero flip on element ",e1," side ",s1
          return
        endif
        if(mesh%sideInfo(5,s1,e1) /= 0) then
          print*,"FAIL (",trim(label),"): interior side carries a boundary id, element ",e1," side ",s1
          return
        endif
        if(mesh%sideInfo(3,s2,e2) /= e1 .or. mesh%sideInfo(4,s2,e2) /= 10*s1 .or. &
           mesh%sideInfo(2,s2,e2) /= g) then
          print*,"FAIL (",trim(label),"): neighbor does not name element ",e1," side ",s1," back"
          return
        endif

        ! The two sides must coincide, after undoing a periodic wrap
        call EdgeEndpoints(mesh,e2,s2,a2,b2)
        shift = 0.0_prec
        if(s1 == 2 .and. a2(1) < a1(1)-0.5_prec*dx) shift(1) = Lx
        if(s1 == 4 .and. a2(1) > a1(1)+0.5_prec*dx) shift(1) = -Lx
        if(s1 == 3 .and. a2(2) < a1(2)-0.5_prec*dy) shift(2) = Ly
        if(s1 == 1 .and. a2(2) > a1(2)+0.5_prec*dy) shift(2) = -Ly
        if(maxval(abs(a2+shift-a1)) > tol .or. maxval(abs(b2+shift-b1)) > tol) then
          print*,"FAIL (",trim(label),"): sides do not coincide, element ",e1," side ",s1
          return
        endif
      enddo
    enddo

    expectedBoundary = 0
    if(.not. periodicX) expectedBoundary = expectedBoundary+2*nY
    if(.not. periodicY) expectedBoundary = expectedBoundary+2*nX
    if(nBoundary /= expectedBoundary) then
      print*,"FAIL (",trim(label),"): ",nBoundary," boundary sides, expected ",expectedBoundary
      return
    endif
    if(any(idcount == 0) .or. any(idcount > 2)) then
      print*,"FAIL (",trim(label),"): global side ids are not each used once or twice"
      return
    endif
    if(count(idcount == 1) /= expectedBoundary) then
      print*,"FAIL (",trim(label),"): ",count(idcount == 1)," ids used once, expected ",expectedBoundary
      return
    endif

    if(.not.(periodicX .or. periodicY)) then
      call reference%StructuredMesh(nxPerTile,nyPerTile,nTileX,nTileY,dx,dy,bcids)
      if(any(reference%sideInfo(3:5,:,:) /= mesh%sideInfo(3:5,:,:))) then
        print*,"FAIL (",trim(label),"): connectivity differs from StructuredMesh"
        return
      endif
      if(reference%nUniqueSides /= mesh%nUniqueSides) then
        print*,"FAIL (",trim(label),"): nUniqueSides differs from StructuredMesh"
        return
      endif
      if(any(reference%nodeCoords /= mesh%nodeCoords)) then
        print*,"FAIL (",trim(label),"): node coordinates differ from StructuredMesh"
        return
      endif
      call reference%Free()
    endif

    deallocate(idcount)
    call mesh%Free()
    print*,"PASS (",trim(label),")"
    r = 0

  endfunction CheckPeriodicMesh

  subroutine EdgeEndpoints(mesh,e,s,a,b)
    !! Corner coordinates of side s of element e, in the side's own direction
    !! (CGNSSideMap order), for a bilinear (nGeo = 1) element.
    implicit none
    type(Mesh2D),intent(in) :: mesh
    integer,intent(in) :: e
    integer,intent(in) :: s
    real(prec),intent(out) :: a(1:2)
    real(prec),intent(out) :: b(1:2)
    ! Local
    integer :: c1,c2

    c1 = mesh%CGNSSideMap(1,s)
    c2 = mesh%CGNSSideMap(2,s)
    a = mesh%nodeCoords(1:2,mesh%CGNSCornerMap(1,c1),mesh%CGNSCornerMap(2,c1),e)
    b = mesh%nodeCoords(1:2,mesh%CGNSCornerMap(1,c2),mesh%CGNSCornerMap(2,c2),e)

  endsubroutine EdgeEndpoints

  integer function mesh2d_periodicstructured() result(r)
    implicit none
    ! Local
    type(Mesh2D) :: keepalive

    ! SELF finalizes MPI when the last live decomposition is freed, and MPI cannot be
    ! initialized twice. This mesh keeps one decomposition alive while the meshes
    ! under test are built and freed in turn.
    call keepalive%StructuredMesh(2,1,1,1,1.0_prec,1.0_prec,[1,1,1,1])

    r = 0
    r = max(r,CheckPeriodicMesh(.true.,.true.))
    r = max(r,CheckPeriodicMesh(.true.,.false.))
    r = max(r,CheckPeriodicMesh(.false.,.true.))
    r = max(r,CheckPeriodicMesh(.false.,.false.))

    call keepalive%Free()

  endfunction mesh2d_periodicstructured

endprogram test
