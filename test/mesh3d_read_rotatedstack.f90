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

  exit_code = mesh3d_read_rotatedstack()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mesh3d_read_rotatedstack() result(r)
    !! Reads share/mesh/MultiMaterial3D/ReaderFixture_RotatedStack.mesh:
    !! a 1x1x3 stack of unit cubes (ISM, polyOrder=2) where the middle
    !! element's corner numbering is rotated 90 degrees about z and the
    !! top element's by 180 degrees. The shared horizontal faces
    !! therefore pair with transposed flips (here flip 7), which an
    !! extruded HOHQMesh mesh can never produce. Selected faces carry
    !! (bilinear) face-point data in flag patterns HOHQMesh's sweep
    !! never writes, exercising the flagged-face-preferred edge
    !! extraction in the transfinite interpolation. All boundary names
    !! are "---", so the reader must fall back to a single "unused"
    !! boundary-condition slot.

    use SELF_Constants
    use SELF_Mesh_3D

    implicit none

#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Mesh3D),target :: mesh
    character(LEN=255) :: WORKSPACE
    integer :: e,e2,s,s2,flip
    integer :: i,j,i2,j2,ng1
    integer :: nInterior,nBoundary
    real(prec) :: x1(1:3),x2(1:3)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOHQMesh(trim(WORKSPACE)// &
                            "/share/mesh/MultiMaterial3D/ReaderFixture_RotatedStack.mesh")

    r = 1
    if(mesh%nElem /= 3) then
      print*,"FAIL: nElem expected 3, got ",mesh%nElem
      return
    endif
    if(mesh%nGeo /= 2) then
      print*,"FAIL: nGeo expected 2, got ",mesh%nGeo
      return
    endif
    if(mesh%nMaterials /= 1 .or. trim(mesh%materialNames(1)) /= "default") then
      print*,"FAIL: plain ISM stack must carry the single default material"
      return
    endif

    ! With no named boundaries the reader must keep one "unused" slot
    if(mesh%nBCs /= 1 .or. trim(mesh%BCNames(1)) /= "unused") then
      print*,"FAIL: expected a single 'unused' BC slot, got nBCs=",mesh%nBCs
      return
    endif

    ! 14 boundary faces (4 sidewalls per element + bottom + top) and 4
    ! interior side slots (2 shared faces)
    nInterior = 0
    nBoundary = 0
    do e = 1,mesh%nElem
      do s = 1,6
        if(mesh%sideInfo(3,s,e) == 0) then
          nBoundary = nBoundary+1
        else
          nInterior = nInterior+1
        endif
      enddo
    enddo
    if(nBoundary /= 14 .or. nInterior /= 4) then
      print*,"FAIL: expected 14 boundary/4 interior side slots, got ", &
        nBoundary,nInterior
      return
    endif
    if(mesh%nUniqueSides /= 16) then
      print*,"FAIL: nUniqueSides expected 16, got ",mesh%nUniqueSides
      return
    endif

    ! The 90-degree relative rotations must be detected as transposed
    ! flips: looking up from each element's top side, the neighbor's
    ! bottom face is traversed with flip 7 (i2 = j1, j2 = N+2-i1).
    do e = 1,2
      if(mesh%sideInfo(3,6,e) /= e+1) then
        print*,"FAIL: element ",e," top neighbor expected ",e+1, &
          ", got ",mesh%sideInfo(3,6,e)
        return
      endif
      s2 = mesh%sideInfo(4,6,e)/10
      flip = mesh%sideInfo(4,6,e)-10*s2
      if(s2 /= 1 .or. flip /= 7) then
        print*,"FAIL: element ",e," top side expected neighbor side 1 flip 7, got ", &
          s2,flip
        return
      endif
    enddo

    ! Watertight interior faces: the stored flips must map each face
    ! point onto the physically coincident neighbor point.
    ng1 = mesh%nGeo+1
    do e = 1,mesh%nElem
      do s = 1,6
        e2 = mesh%sideInfo(3,s,e)
        if(e2 == 0) cycle
        s2 = mesh%sideInfo(4,s,e)/10
        flip = mesh%sideInfo(4,s,e)-10*s2
        do j = 1,ng1
          do i = 1,ng1
            select case(flip)
            case(0); i2 = i; j2 = j
            case(1); i2 = ng1+1-i; j2 = j
            case(2); i2 = ng1+1-i; j2 = ng1+1-j
            case(3); i2 = i; j2 = ng1+1-j
            case(4); i2 = j; j2 = i
            case(5); i2 = ng1+1-j; j2 = i
            case(6); i2 = ng1+1-j; j2 = ng1+1-i
            case(7); i2 = j; j2 = ng1+1-i
            case default
              print*,"FAIL: invalid flip ",flip," at e=",e," side=",s
              return
            endselect
            x1 = face_point(mesh,e,s,i,j)
            x2 = face_point(mesh,e2,s2,i2,j2)
            if(maxval(abs(x1-x2)) > tolerance) then
              print*,"FAIL: face mismatch at e=",e," side=",s," (i,j)=",i,j, &
                " |dx|=",maxval(abs(x1-x2))
              return
            endif
          enddo
        enddo
      enddo
    enddo

    print*,"PASS: nElem=",mesh%nElem," nBCs=",mesh%nBCs, &
      " nBoundary=",nBoundary," nInterior=",nInterior

    call mesh%Free()

    r = 0
  endfunction mesh3d_read_rotatedstack

  function face_point(mesh,e,s,i,j) result(x)
    !! Returns the mesh node coordinate on local side s of element e at
    !! face indices (i,j), following the boundary-interpolation index
    !! convention (side 1: (i,j,1), side 2: (i,1,j), side 3: (N+1,i,j),
    !! side 4: (i,N+1,j), side 5: (1,i,j), side 6: (i,j,N+1)).
    use SELF_Constants
    use SELF_Mesh_3D
    implicit none
    type(Mesh3D),intent(in) :: mesh
    integer,intent(in) :: e,s,i,j
    real(prec) :: x(1:3)
    integer :: ng1

    ng1 = mesh%nGeo+1
    select case(s)
    case(1); x = mesh%nodeCoords(1:3,i,j,1,e)
    case(2); x = mesh%nodeCoords(1:3,i,1,j,e)
    case(3); x = mesh%nodeCoords(1:3,ng1,i,j,e)
    case(4); x = mesh%nodeCoords(1:3,i,ng1,j,e)
    case(5); x = mesh%nodeCoords(1:3,1,i,j,e)
    case(6); x = mesh%nodeCoords(1:3,i,j,ng1,e)
    case default; x = 0.0_prec
    endselect
  endfunction face_point

endprogram test
