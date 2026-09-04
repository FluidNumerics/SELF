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
!! Verifies that the global edge ids -- sideInfo(2,:,:) -- are valid for the 2-D
!! meshes SELF can build or read.
!!
!! The global edge id is the only thing that pairs the two sides of an interface
!! across ranks: MappedScalar2D/MappedVector2D build their MPI message tag from
!! it (tag = |globalSideId| + nUniqueSides*(ivar-1)). An unset or mismatched id
!! therefore does not produce a wrong answer locally -- it mispairs, or hangs, a
!! parallel run. This test asserts the invariants in serial, where every element
!! is rank-local and the pairing can be checked directly.

  use SELF_Constants
  use SELF_Mesh_2D

  implicit none
  integer :: exit_code

  exit_code = mesh2d_globalsideids()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function CheckGlobalSideIds(mesh,label,contiguousIds) result(r)
    !! Asserts, for a mesh whose elements are all rank-local:
    !!
    !!   1. every side slot carries a non-zero global edge id,
    !!   2. an interior edge id is used by exactly two side slots and a boundary
    !!      edge id by exactly one,
    !!   3. the neighbor's side slot names this element back and carries the same
    !!      edge id (up to sign; HOPr signs the id to mark master/slave),
    !!   4. when `contiguousIds` is set, the ids are exactly 1..nUniqueSides.
    !!
    !! `contiguousIds` is false only for HOPr 2-D meshes: those are read from an
    !! extruded 3-D file and keep the 3-D file's edge ids, which run past the
    !! 2-D nUniqueSides that Read_HOPr derives.
    implicit none
    type(Mesh2D),intent(in) :: mesh
    character(*),intent(in) :: label
    logical,intent(in) :: contiguousIds
    ! Local
    integer :: e1,e2,s1,s2,i,maxid
    integer,allocatable :: idcount(:)

    r = 1

    do e1 = 1,mesh%nElem
      do s1 = 1,4
        if(mesh%sideInfo(2,s1,e1) == 0) then
          print*,"FAIL ("//label//"): global edge id is not set at element ",e1," side ",s1
          return
        endif
      enddo
    enddo

    maxid = maxval(abs(mesh%sideInfo(2,1:4,1:mesh%nElem)))
    if(contiguousIds .and. maxid /= mesh%nUniqueSides) then
      print*,"FAIL ("//label//"): largest global edge id ",maxid, &
        " does not match nUniqueSides ",mesh%nUniqueSides
      return
    endif

    allocate(idcount(1:maxid))
    idcount(1:maxid) = 0
    do e1 = 1,mesh%nElem
      do s1 = 1,4
        i = abs(mesh%sideInfo(2,s1,e1))
        idcount(i) = idcount(i)+1
      enddo
    enddo

    do e1 = 1,mesh%nElem
      do s1 = 1,4

        i = abs(mesh%sideInfo(2,s1,e1))
        e2 = mesh%sideInfo(3,s1,e1)

        if(e2 == 0) then ! boundary edge

          if(idcount(i) /= 1) then
            print*,"FAIL ("//label//"): boundary edge id ",i," is used ",idcount(i), &
              " times (expected 1) at element ",e1," side ",s1
            return
          endif

        else ! interior edge

          if(idcount(i) /= 2) then
            print*,"FAIL ("//label//"): interior edge id ",i," is used ",idcount(i), &
              " times (expected 2) at element ",e1," side ",s1
            return
          endif
          if(e2 < 1 .or. e2 > mesh%nElem) then
            print*,"FAIL ("//label//"): neighbor element ",e2," is out of range at element ", &
              e1," side ",s1
            return
          endif
          s2 = mesh%sideInfo(4,s1,e1)/10
          if(s2 < 1 .or. s2 > 4) then
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
            print*,"FAIL ("//label//"): edge id ",i," at element ",e1," side ",s1, &
              " does not match the neighbor's edge id ",abs(mesh%sideInfo(2,s2,e2))
            return
          endif

        endif

      enddo
    enddo

    if(contiguousIds) then
      do i = 1,maxid
        if(idcount(i) == 0) then
          print*,"FAIL ("//label//"): global edge id ",i," is never used"
          return
        endif
      enddo
    endif

    print*,"PASS ("//label//"): nElem = ",mesh%nElem," nUniqueSides = ",mesh%nUniqueSides

    r = 0

  endfunction CheckGlobalSideIds

  integer function mesh2d_globalsideids() result(r)
    !! The meshes are held simultaneously and freed at the end: SELF owns the MPI
    !! lifecycle here, and freeing the last live mesh finalizes MPI, which would
    !! leave any mesh created afterwards without one.
    implicit none
    type(Mesh2D),target :: meshTiled,meshSingleTile,meshBlock,meshCircle,meshHopr
    character(LEN=255) :: WORKSPACE
    integer :: bcids(1:4)
    integer :: ri

    call get_environment_variable("WORKSPACE",WORKSPACE)

    bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                  SELF_BC_PRESCRIBED, & ! East
                  SELF_BC_PRESCRIBED, & ! North
                  SELF_BC_PRESCRIBED] ! West

    r = 0

    ! Structured mesh. More than one tile in each direction is required: the tile
    ! boundaries are where a side takes its edge id from a neighbor in another
    ! tile, and those are the sides that become MPI interfaces under a
    ! decomposition.
    call meshTiled%StructuredMesh(3,3,2,2,0.1_prec,0.13_prec,bcids)
    ri = CheckGlobalSideIds(meshTiled,"StructuredMesh 3x3 per tile, 2x2 tiles",.true.)
    r = max(r,ri)

    ! Single tile: every edge id is assigned in place rather than copied from a
    ! neighboring tile.
    call meshSingleTile%StructuredMesh(4,3,1,1,0.1_prec,0.13_prec,bcids)
    ri = CheckGlobalSideIds(meshSingleTile,"StructuredMesh 4x3, single tile",.true.)
    r = max(r,ri)

    ! HOHQMesh (ISM) meshes assign their edge ids while matching sides.
    call meshBlock%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D.mesh")
    ri = CheckGlobalSideIds(meshBlock,"HOHQMesh Block2D.mesh",.true.)
    r = max(r,ri)

    ! Circle.mesh is unstructured and curved (nGeo > 1), so its sides carry
    ! non-trivial neighbor orientations.
    call meshCircle%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Circle/Circle.mesh")
    ri = CheckGlobalSideIds(meshCircle,"HOHQMesh Circle.mesh",.true.)
    r = max(r,ri)

    ! HOPr meshes carry their edge ids in the file; see the note on
    ! `contiguousIds` above for why the range is not asserted here.
    call meshHopr%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5")
    ri = CheckGlobalSideIds(meshHopr,"HOPr Block2D_mesh.h5",.false.)
    r = max(r,ri)

    call meshTiled%Free()
    call meshSingleTile%Free()
    call meshBlock%Free()
    call meshCircle%Free()
    call meshHopr%Free()

  endfunction mesh2d_globalsideids

endprogram test
