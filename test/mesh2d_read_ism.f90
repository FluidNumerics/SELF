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

  exit_code = mesh2d_read_ism()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mesh2d_read_ism() result(r)
    !! Reads share/mesh/Block2D/Block2D.mesh (plain ISM format,
    !! polyOrder=1, 25 quad elements, single "default" material) and
    !! sanity-checks the resulting Mesh2D state.

    use SELF_Constants
    use SELF_Mesh_2D

    implicit none

    type(Mesh2D),target :: mesh
    character(LEN=255) :: WORKSPACE
    integer :: e,iSide,nInterior,nBoundary
    integer :: bcid

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D.mesh")

    r = 1
    if(mesh%nElem /= 25) then
      print*,"FAIL: nElem expected 25, got ",mesh%nElem
      return
    endif
    if(mesh%nGeo /= 1) then
      print*,"FAIL: nGeo expected 1, got ",mesh%nGeo
      return
    endif
    if(mesh%nMaterials /= 1) then
      print*,"FAIL: nMaterials expected 1 (single default), got ",mesh%nMaterials
      return
    endif
    if(trim(mesh%materialNames(1)) /= "default") then
      print*,"FAIL: default material name, got '",trim(mesh%materialNames(1)),"'"
      return
    endif

    ! Every element should be tagged with material id 1
    do e = 1,mesh%nElem
      if(mesh%elemMaterial(e) /= 1) then
        print*,"FAIL: element ",e," has unexpected material id ",mesh%elemMaterial(e)
        return
      endif
    enddo

    ! Count interior vs boundary sides. For a 5x5 structured block
    ! mesh we expect 20 boundary sides (4*5) and 80 interior sides
    ! (5*4 + 4*5).
    nInterior = 0
    nBoundary = 0
    do e = 1,mesh%nElem
      do iSide = 1,4
        if(mesh%sideInfo(3,iSide,e) == 0) then
          nBoundary = nBoundary+1
          bcid = mesh%sideInfo(5,iSide,e)
          if(bcid <= 0 .or. bcid > mesh%nBCs) then
            print*,"FAIL: boundary side has bad bc id ",bcid," (nBCs=",mesh%nBCs,")"
            return
          endif
        else
          nInterior = nInterior+1
        endif
      enddo
    enddo
    if(nBoundary /= 20) then
      print*,"FAIL: nBoundary expected 20, got ",nBoundary
      return
    endif
    if(nInterior /= 80) then
      print*,"FAIL: nInterior expected 80, got ",nInterior
      return
    endif

    if(mesh%nBCs < 4) then
      print*,"FAIL: nBCs expected at least 4 (Bottom,Top,Left,Right), got ",mesh%nBCs
      return
    endif

    call mesh%Free()

    r = 0
  endfunction mesh2d_read_ism
endprogram test
