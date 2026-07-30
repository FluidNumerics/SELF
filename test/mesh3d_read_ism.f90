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

  exit_code = mesh3d_read_ism()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mesh3d_read_ism() result(r)
    !! Reads share/mesh/Block3D/Block3D.mesh (HOHQMesh 3-D plain ISM
    !! format, polyOrder=1, 2x2x2 unit cube, single "default"
    !! material) and sanity-checks the resulting Mesh3D state.

    use SELF_Constants
    use SELF_Mesh_3D

    implicit none

    type(Mesh3D),target :: mesh
    character(LEN=255) :: WORKSPACE
    integer :: e,s,flip
    integer :: nInterior,nBoundary
    integer :: bcid

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D.mesh")

    r = 1
    if(mesh%nElem /= 8) then
      print*,"FAIL: nElem expected 8, got ",mesh%nElem
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

    ! The unit cube carries all six boundary names
    if(mesh%nBCs /= 6) then
      print*,"FAIL: nBCs expected 6, got ",mesh%nBCs
      return
    endif

    ! Count interior vs boundary faces. For a 2x2x2 block every
    ! element has 3 boundary faces (24 total) and 3 interior faces
    ! (24 slots pairing into 12 unique faces).
    nInterior = 0
    nBoundary = 0
    do e = 1,mesh%nElem
      do s = 1,6
        if(mesh%sideInfo(3,s,e) == 0) then
          nBoundary = nBoundary+1
          bcid = mesh%sideInfo(5,s,e)
          if(bcid <= 0 .or. bcid > mesh%nBCs) then
            print*,"FAIL: boundary face has bad bc id ",bcid," (nBCs=",mesh%nBCs,")"
            return
          endif
        else
          nInterior = nInterior+1
          if(mesh%sideInfo(3,s,e) < 1 .or. mesh%sideInfo(3,s,e) > mesh%nElem) then
            print*,"FAIL: interior face has out-of-range neighbor at e=",e," side=",s
            return
          endif
          ! An axis-aligned structured block has identically oriented
          ! neighbors, so every interior face must have flip 0.
          flip = mesh%sideInfo(4,s,e)-10*(mesh%sideInfo(4,s,e)/10)
          if(flip /= 0) then
            print*,"FAIL: expected flip 0 on structured block, got ",flip, &
              " at e=",e," side=",s
            return
          endif
        endif
      enddo
    enddo
    if(nBoundary /= 24 .or. nInterior /= 24) then
      print*,"FAIL: expected 24 boundary/24 interior side slots, got ", &
        nBoundary,nInterior
      return
    endif
    if(mesh%nUniqueSides /= 36) then
      print*,"FAIL: nUniqueSides expected 36, got ",mesh%nUniqueSides
      return
    endif

    print*,"PASS: nElem=",mesh%nElem," nGeo=",mesh%nGeo, &
      " nMaterials=",mesh%nMaterials," nBoundary=",nBoundary," nInterior=",nInterior

    call mesh%Free()

    r = 0
  endfunction mesh3d_read_ism
endprogram test
