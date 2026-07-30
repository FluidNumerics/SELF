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

  exit_code = mesh3d_read_ismmm()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mesh3d_read_ismmm() result(r)
    !! Reads share/mesh/MultiMaterial3D/InsulatedWire3D.mesh
    !! (HOHQMesh 3-D ISM-MM format, polyOrder=4, curved outer
    !! cylinder faces, two materials "Copper" and "Insulator").
    !!
    !! Validation criteria:
    !!  * mesh counts, quadrature, and the parsed material table
    !!  * every interior face pairs with a real neighbor, boundary
    !!    faces carry a valid boundary-condition id, and the unique
    !!    side count is exact for this mesh topology
    !!  * the flip stored for every interior face maps each face
    !!    point of the primary side onto the physically coincident
    !!    point of the secondary side (watertight element faces)
    !!  * the metric terms generated from the mesh have strictly
    !!    positive Jacobians (valid right-handed curved elements)

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D

    implicit none

    integer,parameter :: controlDegree = 7
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    character(LEN=255) :: WORKSPACE
    integer :: e,e2,s,s2,flip,m
    integer :: i,j,k,i2,j2,ng1
    integer :: nInterior,nBoundary
    integer :: nCopper,nInsulator
    logical :: sawCopper,sawInsulator
    real(prec) :: x1(1:3),x2(1:3)
    real(prec) :: minJ

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/MultiMaterial3D/InsulatedWire3D.mesh")

    r = 1
    if(mesh%nElem /= 24) then
      print*,"FAIL: nElem expected 24, got ",mesh%nElem
      return
    endif
    if(mesh%nGeo /= 4) then
      print*,"FAIL: nGeo expected 4, got ",mesh%nGeo
      return
    endif
    if(mesh%quadrature /= CHEBYSHEV_GAUSS_LOBATTO) then
      print*,"FAIL: quadrature expected CHEBYSHEV_GAUSS_LOBATTO, got ",mesh%quadrature
      return
    endif

    if(mesh%nMaterials /= 2) then
      print*,"FAIL: nMaterials expected 2, got ",mesh%nMaterials
      return
    endif
    sawCopper = .false.
    sawInsulator = .false.
    do m = 1,mesh%nMaterials
      if(trim(mesh%materialNames(m)) == "Copper") sawCopper = .true.
      if(trim(mesh%materialNames(m)) == "Insulator") sawInsulator = .true.
    enddo
    if(.not. sawCopper) then
      print*,"FAIL: 'Copper' material missing from table"
      return
    endif
    if(.not. sawInsulator) then
      print*,"FAIL: 'Insulator' material missing from table"
      return
    endif

    ! Material assignment: 4 core elements and 8 ring elements per
    ! layer, two layers.
    nCopper = 0
    nInsulator = 0
    do e = 1,mesh%nElem
      if(mesh%elemMaterial(e) < 1 .or. mesh%elemMaterial(e) > mesh%nMaterials) then
        print*,"FAIL: element ",e," has out-of-range material id ",mesh%elemMaterial(e)
        return
      endif
      if(trim(mesh%materialNames(mesh%elemMaterial(e))) == "Copper") then
        nCopper = nCopper+1
      else
        nInsulator = nInsulator+1
      endif
    enddo
    if(nCopper /= 8 .or. nInsulator /= 16) then
      print*,"FAIL: material counts expected 8 Copper / 16 Insulator, got ", &
        nCopper,nInsulator
      return
    endif

    ! Side accounting. The insulated wire has 16 outer cylinder faces,
    ! 12 bottom faces, and 12 top faces on the physical boundary; the
    ! remaining 104 face slots pair into 52 interior faces.
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
          if(mesh%sideInfo(3,s,e) < 1 .or. mesh%sideInfo(3,s,e) > mesh%nElem) then
            print*,"FAIL: interior face has out-of-range neighbor at e=",e," side=",s
            return
          endif
          if(mesh%sideInfo(5,s,e) /= 0) then
            print*,"FAIL: interior face carries a bc id at e=",e," side=",s
            return
          endif
        endif
      enddo
    enddo
    if(nBoundary /= 40 .or. nInterior /= 104) then
      print*,"FAIL: expected 40 boundary/104 interior side slots, got ", &
        nBoundary,nInterior
      return
    endif
    if(mesh%nUniqueSides /= 92) then
      print*,"FAIL: nUniqueSides expected 92, got ",mesh%nUniqueSides
      return
    endif

    ! Watertight interior faces: for every interior face, the stored
    ! flip must map each face point of this element onto the
    ! physically coincident face point of the neighbor.
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

    ! Metric terms: curved multi-material elements must produce
    ! strictly positive Jacobians.
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=controlDegree+2, &
                     targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    minJ = minval(geometry%J%interior)
    if(.not.(minJ > 0.0_prec)) then
      print*,"FAIL: nonpositive Jacobian, min J = ",minJ
      call geometry%Free()
      call interp%Free()
      return
    endif

    print*,"PASS: nElem=",mesh%nElem," nGeo=",mesh%nGeo, &
      " nMaterials=",mesh%nMaterials," nBoundary=",nBoundary, &
      " nInterior=",nInterior," min J=",minJ

    call geometry%Free()
    call interp%Free()
    call mesh%Free()

    r = 0
  endfunction mesh3d_read_ismmm

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
