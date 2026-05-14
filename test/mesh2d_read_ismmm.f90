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

  exit_code = mesh2d_read_ismmm()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mesh2d_read_ismmm() result(r)
    !! Reads share/mesh/MultiMaterial2D/BoneAndMarrow.mesh
    !! (ISM-MM format, polyOrder=4, curved outer boundary,
    !!  multiple materials including "Muscle" and "Bone").

    use SELF_Constants
    use SELF_Mesh_2D

    implicit none

    type(Mesh2D),target :: mesh
    character(LEN=255) :: WORKSPACE
    integer :: e,iSide,nInterior,nBoundary,m
    logical :: sawMuscle,sawBone

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/MultiMaterial2D/BoneAndMarrow.mesh")

    r = 1
    if(mesh%nElem /= 973) then
      print*,"FAIL: nElem expected 973, got ",mesh%nElem
      return
    endif
    if(mesh%nGeo /= 4) then
      print*,"FAIL: nGeo expected 4, got ",mesh%nGeo
      return
    endif
    if(mesh%quadrature /= CHEBYSHEV_GAUSS_LOBATTO) then
      print*,"FAIL: quadrature expected CHEBYSHEV_GAUSS_LOBATTO (4), got ",mesh%quadrature
      return
    endif

    if(mesh%nMaterials < 2) then
      print*,"FAIL: nMaterials expected >= 2, got ",mesh%nMaterials
      return
    endif

    sawMuscle = .false.
    sawBone = .false.
    do m = 1,mesh%nMaterials
      if(trim(mesh%materialNames(m)) == "Muscle") sawMuscle = .true.
      if(trim(mesh%materialNames(m)) == "Bone") sawBone = .true.
    enddo
    if(.not. sawMuscle) then
      print*,"FAIL: 'Muscle' material missing from table"
      return
    endif
    if(.not. sawBone) then
      print*,"FAIL: 'Bone' material missing from table"
      return
    endif

    ! Every elemMaterial must point at a valid table index
    do e = 1,mesh%nElem
      if(mesh%elemMaterial(e) < 1 .or. mesh%elemMaterial(e) > mesh%nMaterials) then
        print*,"FAIL: element ",e," has out-of-range material id ",mesh%elemMaterial(e)
        return
      endif
    enddo

    ! Every interior side must point to a real neighbor; every boundary
    ! side must carry a positive bc id.
    nInterior = 0
    nBoundary = 0
    do e = 1,mesh%nElem
      do iSide = 1,4
        if(mesh%sideInfo(3,iSide,e) == 0) then
          nBoundary = nBoundary+1
          if(mesh%sideInfo(5,iSide,e) <= 0) then
            print*,"FAIL: boundary side has zero/negative bc id at e=",e," side=",iSide
            return
          endif
        else
          nInterior = nInterior+1
          if(mesh%sideInfo(3,iSide,e) < 1 .or. mesh%sideInfo(3,iSide,e) > mesh%nElem) then
            print*,"FAIL: interior side has out-of-range neighbor at e=",e," side=",iSide
            return
          endif
        endif
      enddo
    enddo
    if(nInterior+nBoundary /= 4*mesh%nElem) then
      print*,"FAIL: side accounting mismatch"
      return
    endif

    print*,"PASS: nElem=",mesh%nElem," nGeo=",mesh%nGeo, &
      " nMaterials=",mesh%nMaterials," nBoundary=",nBoundary," nInterior=",nInterior

    call mesh%Free()

    r = 0
  endfunction mesh2d_read_ismmm
endprogram test
