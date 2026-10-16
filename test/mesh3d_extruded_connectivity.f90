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
!! Verifies the connectivity, numbering, decomposition and node coordinates of meshes
!! built with Mesh3D%ExtrudedMesh from five surface meshes: a structured mesh, the same
!! kind of mesh with every other element turned through 180 degrees (so that its edges
!! carry 2-D flip 1), a channel that is periodic in x, the curved HOHQMesh Circle mesh,
!! and the HOPr Block2D mesh (whose global edge ids carry a master/slave sign).
!!
!! In serial every element is rank-local, so both sides of each face can be read
!! directly. For each extruded mesh the test asserts that
!!
!!   1. element counts, nLayers, nGlobalElem and nUniqueSides follow the documented
!!      formulas, and the decomposition is the surface one expanded by nz,
!!   2. x and y of every node equal the surface node coordinates, and the bottom and
!!      top node layers of element k equal the levels z_{k-1} and z_k bitwise,
!!   3. every face has a global id within 1..nUniqueSides, and the ids are unique: an
!!      interior face id is used by exactly two face slots, a boundary id by one,
!!   4. the neighbor of an interior face names this element and face back with the
!!      same id, in the same layer for a lateral face and in the adjacent layer for a
!!      horizontal one, and both sides agree on the flip,
!!   5. a lateral boundary face carries the surface edge's boundary id, and the bottom
!!      and top faces of a column carry bottomBC and topBC,
!!   6. on the non-periodic meshes, the face node coordinates of the two sides of every
!!      interior face coincide under the recorded flip (MortarFaceMap, the index map
!!      SideExchange applies), which pins the lateral flips down geometrically.

  use SELF_Constants
  use SELF_Mesh_2D
  use SELF_Mesh_3D

  implicit none
  integer :: exit_code

  exit_code = mesh3d_extruded_connectivity()
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

  subroutine FaceNode(mesh,e,s,i,j,p)
    !! Node coordinates at geometry-node index (i,j) of face s of element e, in the face
    !! trace coordinates of BoundaryInterp: faces 1,6 -> (xi1,xi2); faces 2,4 ->
    !! (xi1,xi3); faces 3,5 -> (xi2,xi3).
    implicit none
    type(Mesh3D),intent(in) :: mesh
    integer,intent(in) :: e,s,i,j
    real(prec),intent(out) :: p(1:3)
    ! Local
    integer :: n

    n = mesh%nGeo+1
    select case(s)
    case(1)
      p = mesh%nodeCoords(1:3,i,j,1,e)
    case(2)
      p = mesh%nodeCoords(1:3,i,1,j,e)
    case(3)
      p = mesh%nodeCoords(1:3,n,i,j,e)
    case(4)
      p = mesh%nodeCoords(1:3,i,n,j,e)
    case(5)
      p = mesh%nodeCoords(1:3,1,i,j,e)
    case default
      p = mesh%nodeCoords(1:3,i,j,n,e)
    endselect

  endsubroutine FaceNode

  integer function CheckExtruded(surf,label,checkGeometry) result(r)
    implicit none
    type(Mesh2D),intent(in) :: surf
    character(*),intent(in) :: label
    logical,intent(in) :: checkGeometry
    ! Local
    integer,parameter :: nz = 3
    integer,parameter :: bottomBC = 21,topBC = 22
    real(prec),parameter :: sigma(0:nz) = [-1.0_prec,-0.6_prec,-0.25_prec,0.0_prec]
    real(prec),parameter :: tol = 1.0e-10_prec
    type(Mesh3D) :: mesh
    real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
    integer,allocatable :: idcount(:)
    integer :: nGeo,e2,e3,k,s1,s2,e1g,e2g,e3n,i,j,i2,j2,g,flip,flip2,r2D
    integer :: nGlobal2D,nFlipped,nMismatch
    real(prec) :: p1(1:3),p2(1:3),x,y

    r = 1
    nGeo = surf%nGeo

    ! Smooth, positive bathymetry sampled at the surface geometry nodes
    allocate(h(1:nGeo+1,1:nGeo+1,1:surf%nElem))
    allocate(zLevels(1:nGeo+1,1:nGeo+1,0:nz,1:surf%nElem))
    do e2 = 1,surf%nElem
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          x = surf%nodeCoords(1,i,j,e2)
          y = surf%nodeCoords(2,i,j,e2)
          h(i,j,e2) = 1.0_prec+0.3_prec*x+0.2_prec*y*y
        enddo
      enddo
    enddo
    call ExtrudedSigmaLevels(h,sigma,zLevels)

    call mesh%ExtrudedMesh(surf,zLevels,bottomBC,topBC)

    nGlobal2D = surf%decomp%offsetElem(surf%decomp%nRanks+1)

    ! 1. Counts and decomposition
    if(mesh%nLayers /= nz .or. mesh%nElem /= nz*surf%nElem .or. mesh%nGlobalElem /= nz*nGlobal2D) then
      print*,"FAIL (",label,"): element counts ",mesh%nLayers,mesh%nElem,mesh%nGlobalElem
      return
    endif
    if(mesh%nUniqueSides /= nz*surf%nUniqueSides+nGlobal2D*(nz+1)) then
      print*,"FAIL (",label,"): nUniqueSides ",mesh%nUniqueSides
      return
    endif
    if(mesh%nGeo /= nGeo .or. mesh%quadrature /= surf%quadrature) then
      print*,"FAIL (",label,"): nGeo or quadrature not inherited"
      return
    endif
    if(any(mesh%decomp%offsetElem /= nz*surf%decomp%offsetElem)) then
      print*,"FAIL (",label,"): offsetElem is not nz times the surface offsets"
      return
    endif
    do e3 = 1,mesh%nGlobalElem
      r2D = surf%decomp%elemToRank((e3-1)/nz+1)
      if(mesh%decomp%elemToRank(e3) /= r2D) then
        print*,"FAIL (",label,"): elemToRank is not the surface one expanded, element ",e3
        return
      endif
    enddo

    ! 2. Node coordinates
    do e2 = 1,surf%nElem
      do k = 1,nz
        e3 = (e2-1)*nz+k
        do j = 1,nGeo+1
          do i = 1,nGeo+1
            if(any(mesh%nodeCoords(1,i,j,:,e3) /= surf%nodeCoords(1,i,j,e2)) .or. &
               any(mesh%nodeCoords(2,i,j,:,e3) /= surf%nodeCoords(2,i,j,e2))) then
              print*,"FAIL (",label,"): horizontal coordinates depend on zeta, element ",e3
              return
            endif
            if(mesh%nodeCoords(3,i,j,1,e3) /= zLevels(i,j,k-1,e2) .or. &
               mesh%nodeCoords(3,i,j,nGeo+1,e3) /= zLevels(i,j,k,e2)) then
              print*,"FAIL (",label,"): bottom/top node layer is not the level surface, element ",e3
              return
            endif
          enddo
        enddo
      enddo
    enddo

    ! 3-6. Faces
    allocate(idcount(1:mesh%nUniqueSides))
    idcount = 0
    nFlipped = 0
    nMismatch = 0
    do e3 = 1,mesh%nElem
      e1g = e3 ! serial: global id = local id
      e2 = (e3-1)/nz+1
      k = e3-(e2-1)*nz
      do s1 = 1,6
        g = abs(mesh%sideInfo(2,s1,e3))
        if(g < 1 .or. g > mesh%nUniqueSides) then
          print*,"FAIL (",label,"): face id ",mesh%sideInfo(2,s1,e3)," out of range"
          return
        endif
        idcount(g) = idcount(g)+1
        if(mesh%sideInfo(1,s1,e3) /= 0) then
          print*,"FAIL (",label,"): non-zero side type"
          return
        endif

        e2g = mesh%sideInfo(3,s1,e3)
        if(e2g == 0) then
          if(s1 == 1) then
            if(k /= 1 .or. mesh%sideInfo(5,s1,e3) /= bottomBC) then
              print*,"FAIL (",label,"): bottom boundary face misplaced or mislabelled, element ",e3
              return
            endif
          elseif(s1 == 6) then
            if(k /= nz .or. mesh%sideInfo(5,s1,e3) /= topBC) then
              print*,"FAIL (",label,"): top boundary face misplaced or mislabelled, element ",e3
              return
            endif
          else
            if(surf%sideInfo(3,s1-1,e2) /= 0 .or. &
               mesh%sideInfo(5,s1,e3) /= surf%sideInfo(5,s1-1,e2)) then
              print*,"FAIL (",label,"): lateral boundary face does not match its surface edge, element ",e3
              return
            endif
          endif
          cycle
        endif

        s2 = mesh%sideInfo(4,s1,e3)/10
        flip = mesh%sideInfo(4,s1,e3)-10*s2
        e3n = e2g ! serial: local id = global id
        if(mesh%sideInfo(3,s2,e3n) /= e1g .or. mesh%sideInfo(4,s2,e3n)/10 /= s1 .or. &
           abs(mesh%sideInfo(2,s2,e3n)) /= g) then
          print*,"FAIL (",label,"): neighbor does not name element ",e3," face ",s1," back"
          return
        endif
        flip2 = mesh%sideInfo(4,s2,e3n)-10*s1
        if(flip2 /= flip .or. (flip /= 0 .and. flip /= 1)) then
          print*,"FAIL (",label,"): flips ",flip,flip2," are not a consistent lateral or vertical pair"
          return
        endif
        if(s1 == 1 .or. s1 == 6) then
          if(flip /= 0 .or. (s1 == 1 .and. e3n /= e3-1) .or. (s1 == 6 .and. e3n /= e3+1)) then
            print*,"FAIL (",label,"): vertical neighbor is not the adjacent layer, element ",e3
            return
          endif
        else
          if(e3n-((e3n-1)/nz)*nz /= k) then
            print*,"FAIL (",label,"): lateral neighbor is in another layer, element ",e3
            return
          endif
          if(mesh%sideInfo(5,s1,e3) /= 0) then
            print*,"FAIL (",label,"): interior lateral face carries a boundary id"
            return
          endif
          if(flip == 1) nFlipped = nFlipped+1
        endif

        if(checkGeometry) then
          do j = 1,nGeo+1
            do i = 1,nGeo+1
              call MortarFaceMap(i,j,nGeo,flip,i2,j2)
              call FaceNode(mesh,e3,s1,i,j,p1)
              call FaceNode(mesh,e3n,s2,i2,j2,p2)
              if(maxval(abs(p1-p2)) > tol) nMismatch = nMismatch+1
            enddo
          enddo
        endif
      enddo
    enddo

    if(nMismatch > 0) then
      print*,"FAIL (",label,"): ",nMismatch," face nodes do not coincide under the recorded flip"
      return
    endif
    ! An interior face id must be used by exactly its two face slots, a boundary face
    ! id by exactly one
    do e3 = 1,mesh%nElem
      do s1 = 1,6
        g = abs(mesh%sideInfo(2,s1,e3))
        if(mesh%sideInfo(3,s1,e3) == 0) then
          if(idcount(g) /= 1) then
            print*,"FAIL (",label,"): boundary face id ",g," is used ",idcount(g)," times"
            return
          endif
        elseif(idcount(g) /= 2) then
          print*,"FAIL (",label,"): interior face id ",g," is used ",idcount(g)," times"
          return
        endif
      enddo
    enddo

    print*,"PASS (",label,"): ",mesh%nElem," elements, ",nFlipped," flipped lateral face slots"
    r = 0

    call mesh%Free()
    deallocate(h,zLevels,idcount)

  endfunction CheckExtruded

  integer function mesh3d_extruded_connectivity() result(r)
    implicit none
    ! Local
    type(Mesh2D) :: structured,turned,channel,circle,block
    character(LEN=255) :: WORKSPACE
    integer :: bcids(1:4)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    bcids = [11,12,13,14]

    call structured%StructuredMesh(3,2,2,1,0.25_prec,0.4_prec,bcids)
    call turned%StructuredMesh(4,4,2,1,0.25_prec,0.2_prec,bcids)
    call RotateAlternateElements(turned)
    call channel%PeriodicStructuredMesh(3,2,2,1,0.25_prec,0.4_prec,.true.,.false.,bcids)
    call circle%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Circle/Circle.mesh")
    call block%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5")

    r = 0
    r = max(r,CheckExtruded(structured,"StructuredMesh",.true.))
    r = max(r,CheckExtruded(turned,"StructuredMesh, odd elements turned 180 degrees",.true.))
    r = max(r,CheckExtruded(channel,"PeriodicStructuredMesh, periodic in x",.false.))
    r = max(r,CheckExtruded(circle,"HOHQMesh Circle.mesh",.true.))
    r = max(r,CheckExtruded(block,"HOPr Block2D_mesh.h5",.true.))

    call structured%Free()
    call turned%Free()
    call channel%Free()
    call circle%Free()
    call block%Free()

  endfunction mesh3d_extruded_connectivity

endprogram test
