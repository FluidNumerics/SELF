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
!! Verifies that the global face ids -- sideInfo(2,:,:) -- are valid for the 3-D
!! meshes SELF can build or read.
!!
!! The global face id is the only thing that pairs the two sides of an interface
!! across ranks: MappedScalar3D/MappedVector3D build their MPI message tag from
!! it (tag = |globalSideId| + nUniqueSides*(ivar-1)). An unset or mismatched id
!! therefore does not produce a wrong answer locally -- it mispairs, or hangs, a
!! parallel run. This test asserts the invariants in serial, where every element
!! is rank-local and the pairing can be checked directly.

  use SELF_Constants
  use SELF_Mesh_3D

  implicit none
  integer :: exit_code

  exit_code = mesh3d_globalsideids()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function CheckGlobalSideIds(mesh,label,contiguousIds) result(r)
    !! Asserts, for a mesh whose elements are all rank-local:
    !!
    !!   1. every side slot carries a non-zero global face id,
    !!   2. an interior face id is used by exactly two side slots and a boundary
    !!      face id by exactly one,
    !!   3. the neighbor's side slot names this element back and carries the same
    !!      face id (up to sign; HOPr signs the id to mark master/slave),
    !!   4. when `contiguousIds` is set, the ids are exactly 1..nUniqueSides.
    implicit none
    type(Mesh3D),intent(in) :: mesh
    character(*),intent(in) :: label
    logical,intent(in) :: contiguousIds
    ! Local
    integer :: e1,e2,s1,s2,i,maxid
    integer,allocatable :: idcount(:)

    r = 1

    do e1 = 1,mesh%nElem
      do s1 = 1,6
        if(mesh%sideInfo(2,s1,e1) == 0) then
          print*,"FAIL ("//label//"): global face id is not set at element ",e1," side ",s1
          return
        endif
      enddo
    enddo

    maxid = maxval(abs(mesh%sideInfo(2,1:6,1:mesh%nElem)))
    if(contiguousIds .and. maxid /= mesh%nUniqueSides) then
      print*,"FAIL ("//label//"): largest global face id ",maxid, &
        " does not match nUniqueSides ",mesh%nUniqueSides
      return
    endif

    allocate(idcount(1:maxid))
    idcount(1:maxid) = 0
    do e1 = 1,mesh%nElem
      do s1 = 1,6
        i = abs(mesh%sideInfo(2,s1,e1))
        idcount(i) = idcount(i)+1
      enddo
    enddo

    do e1 = 1,mesh%nElem
      do s1 = 1,6

        i = abs(mesh%sideInfo(2,s1,e1))
        e2 = mesh%sideInfo(3,s1,e1)

        if(e2 == 0) then ! boundary face

          if(idcount(i) /= 1) then
            print*,"FAIL ("//label//"): boundary face id ",i," is used ",idcount(i), &
              " times (expected 1) at element ",e1," side ",s1
            return
          endif

        else ! interior face

          if(idcount(i) /= 2) then
            print*,"FAIL ("//label//"): interior face id ",i," is used ",idcount(i), &
              " times (expected 2) at element ",e1," side ",s1
            return
          endif
          if(e2 < 1 .or. e2 > mesh%nElem) then
            print*,"FAIL ("//label//"): neighbor element ",e2," is out of range at element ", &
              e1," side ",s1
            return
          endif
          s2 = mesh%sideInfo(4,s1,e1)/10
          if(s2 < 1 .or. s2 > 6) then
            print*,"FAIL ("//label//"): neighbor side ",s2," is out of range at element ", &
              e1," side ",s1
            return
          endif
          if(mesh%sideInfo(3,s2,e2) /= e1) then
            print*,"FAIL ("//label//"): neighbor (",e2,",",s2,") points back at element ", &
              mesh%sideInfo(3,s2,e2)," instead of ",e1
            return
          endif
          if(abs(mesh%sideInfo(2,s2,e2)) /= i) then
            print*,"FAIL ("//label//"): face id ",i," at element ",e1," side ",s1, &
              " does not match the neighbor's face id ",abs(mesh%sideInfo(2,s2,e2))
            return
          endif

        endif

      enddo
    enddo

    if(contiguousIds) then
      do i = 1,maxid
        if(idcount(i) == 0) then
          print*,"FAIL ("//label//"): global face id ",i," is never used"
          return
        endif
      enddo
    endif

    print*,"PASS ("//label//"): nElem = ",mesh%nElem," nUniqueSides = ",mesh%nUniqueSides

    r = 0

  endfunction CheckGlobalSideIds

  integer function mesh3d_globalsideids() result(r)
    !! The meshes are held simultaneously and freed at the end: SELF owns the MPI
    !! lifecycle here, and freeing the last live mesh finalizes MPI, which would
    !! leave any mesh created afterwards without one.
    implicit none
    type(Mesh3D),target :: meshTiled,meshSingleTile,meshBlock,meshHopr
    character(LEN=255) :: WORKSPACE
    integer :: bcids(1:6)
    integer :: ri

    call get_environment_variable("WORKSPACE",WORKSPACE)

    bcids(1:6) = [SELF_BC_PRESCRIBED, & ! Bottom
                  SELF_BC_PRESCRIBED, & ! South
                  SELF_BC_PRESCRIBED, & ! East
                  SELF_BC_PRESCRIBED, & ! North
                  SELF_BC_PRESCRIBED, & ! West
                  SELF_BC_PRESCRIBED] ! Top

    r = 0

    ! Structured mesh. More than one tile in each direction is required: the tile
    ! boundaries are where a side takes its face id from a neighbor in another
    ! tile, and those are the sides that become MPI interfaces under a
    ! decomposition.
    call meshTiled%StructuredMesh(2,2,2,2,2,2,0.1_prec,0.13_prec,0.17_prec,bcids)
    ri = CheckGlobalSideIds(meshTiled,"StructuredMesh 2x2x2 per tile, 2x2x2 tiles",.true.)
    r = max(r,ri)

    ! Single tile: every face id is assigned in place rather than copied from a
    ! neighboring tile.
    call meshSingleTile%StructuredMesh(3,2,2,1,1,1,0.1_prec,0.13_prec,0.17_prec,bcids)
    ri = CheckGlobalSideIds(meshSingleTile,"StructuredMesh 3x2x2, single tile",.true.)
    r = max(r,ri)

    ! HOHQMesh (ISM) meshes assign their face ids while matching sides.
    call meshBlock%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D.mesh")
    ri = CheckGlobalSideIds(meshBlock,"HOHQMesh Block3D.mesh",.true.)
    r = max(r,ri)

    ! HOPr meshes carry their face ids in the file.
    call meshHopr%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D_mesh.h5")
    ri = CheckGlobalSideIds(meshHopr,"HOPr Block3D_mesh.h5",.true.)
    r = max(r,ri)

    call meshTiled%Free()
    call meshSingleTile%Free()
    call meshBlock%Free()
    call meshHopr%Free()

  endfunction mesh3d_globalsideids

endprogram test
