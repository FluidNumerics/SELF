! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2026 Fluid Numerics LLC
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
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

program test

  implicit none
  integer :: exit_code

  exit_code = mesh3d_read_manymaterials()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mesh3d_read_manymaterials() result(r)
    !! Reads share/mesh/MultiMaterial3D/ReaderFixture_ManyMaterials.mesh:
    !! a 3x3x1 block (ISM-MM, polyOrder=1) where every element carries a
    !! unique material name (9 materials) and every physical boundary
    !! face carries a unique boundary name (30 names). This overflows
    !! the reader's initial material-name (8) and boundary-name (16)
    !! tables, validating the incremental table growth.

    use SELF_Constants
    use SELF_Mesh_3D

    implicit none

    type(Mesh3D),target :: mesh
    character(LEN=255) :: WORKSPACE
    character(LEN=16) :: expected
    integer :: e,s
    integer :: nInterior,nBoundary

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOHQMesh(trim(WORKSPACE)// &
                            "/share/mesh/MultiMaterial3D/ReaderFixture_ManyMaterials.mesh")

    r = 1
    if(mesh%nElem /= 9) then
      print*,"FAIL: nElem expected 9, got ",mesh%nElem
      return
    endif
    if(mesh%nGeo /= 1) then
      print*,"FAIL: nGeo expected 1, got ",mesh%nGeo
      return
    endif

    ! Nine unique materials, inserted in element order
    if(mesh%nMaterials /= 9) then
      print*,"FAIL: nMaterials expected 9, got ",mesh%nMaterials
      return
    endif
    do e = 1,mesh%nElem
      write(expected,'(A,I0)') "mat",e
      if(mesh%elemMaterial(e) /= e) then
        print*,"FAIL: element ",e," expected material id ",e, &
          ", got ",mesh%elemMaterial(e)
        return
      endif
      if(trim(mesh%materialNames(e)) /= trim(expected)) then
        print*,"FAIL: material ",e," expected name '",trim(expected), &
          "', got '",trim(mesh%materialNames(e)),"'"
        return
      endif
    enddo

    ! Thirty unique boundary names: 9 bottom + 9 top + 12 perimeter
    if(mesh%nBCs /= 30) then
      print*,"FAIL: nBCs expected 30, got ",mesh%nBCs
      return
    endif

    nInterior = 0
    nBoundary = 0
    do e = 1,mesh%nElem
      do s = 1,6
        if(mesh%sideInfo(3,s,e) == 0) then
          nBoundary = nBoundary+1
          if(mesh%sideInfo(5,s,e) <= 0 .or. mesh%sideInfo(5,s,e) > mesh%nBCs) then
            print*,"FAIL: boundary face has bad bc id at e=",e," side=",s
            return
          endif
        else
          nInterior = nInterior+1
        endif
      enddo
    enddo
    if(nBoundary /= 30 .or. nInterior /= 24) then
      print*,"FAIL: expected 30 boundary/24 interior side slots, got ", &
        nBoundary,nInterior
      return
    endif
    if(mesh%nUniqueSides /= 42) then
      print*,"FAIL: nUniqueSides expected 42, got ",mesh%nUniqueSides
      return
    endif

    print*,"PASS: nElem=",mesh%nElem," nMaterials=",mesh%nMaterials, &
      " nBCs=",mesh%nBCs," nBoundary=",nBoundary," nInterior=",nInterior

    call mesh%Free()

    r = 0
  endfunction mesh3d_read_manymaterials
endprogram test
