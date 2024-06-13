!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
module SELF_Geometry

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use SELF_Data_1D
  use SELF_Data_2D
  use SELF_Data_3D
  use SELF_SupportRoutines
  use SELF_Mesh

  implicit none

#include "SELF_Macros.h"
  type,public :: SEMGeometry
    integer :: nElem
  endtype SEMGeometry

  type,extends(SEMGeometry),public :: Geometry1D
    type(Scalar1D) :: x ! Physical Positions
    type(Scalar1D) :: dxds ! Conversion from computational to physical space

  contains

    procedure,public :: Init => Init_Geometry1D
    procedure,public :: Free => Free_Geometry1D
    procedure,public :: GenerateFromMesh => GenerateFromMesh_Geometry1D
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_Geometry1D

    procedure :: write => Write_Geometry1D

  endtype Geometry1D

  type,extends(SEMGeometry),public :: SEMQuad
    type(Vector2D) :: x ! Physical positions
    type(Tensor2D) :: dxds ! Covariant basis vectors
    type(Tensor2D) :: dsdx ! Contavariant basis vectors
    type(Vector2D) :: nHat ! Normal Vectors pointing across coordinate lines
    type(Scalar2D) :: nScale ! Boundary scale
    type(Scalar2D) :: J ! Jacobian of the transformation

  contains

    procedure,public :: Init => Init_SEMQuad
    procedure,public :: Free => Free_SEMQuad
    procedure,public :: GenerateFromMesh => GenerateFromMesh_SEMQuad
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_SEMQuad
    procedure,private :: CalculateContravariantBasis => CalculateContravariantBasis_SEMQuad

    !PROCEDURE,PUBLIC :: CovariantArcMin => CovariantArcMin_SEMQuad
    !PROCEDURE :: Write => Write_SEMQuad

  endtype SEMQuad

  type,extends(SEMGeometry),public :: SEMHex
    type(Vector3D) :: x ! Physical positions
    type(Tensor3D) :: dxds ! Covariant basis vectors
    type(Tensor3D) :: dsdx ! Contavariant basis vectors
    type(Vector3D) :: nHat ! Normal Vectors pointing across coordinate lines
    type(Scalar3D) :: nScale ! Boundary scale
    type(Scalar3D) :: J ! Jacobian of the transformation

  contains

    procedure,public :: Init => Init_SEMHex
    procedure,public :: Free => Free_SEMHex
    procedure,public :: GenerateFromMesh => GenerateFromMesh_SEMHex
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_SEMHex
    procedure,private :: CalculateContravariantBasis => CalculateContravariantBasis_SEMHex
    procedure,private :: CheckSides => CheckSides_SEMHex

    !PROCEDURE :: Write => Write_SEMHex

  endtype SEMHex

contains

  subroutine Init_Geometry1D(myGeom,interp,nElem)
    implicit none
    class(Geometry1D),intent(out) :: myGeom
    type(Lagrange),pointer,intent(in) :: interp
    integer,intent(in) :: nElem

    myGeom%nElem = nElem

    call myGeom%x%Init(interp=interp, &
                       nVar=1, &
                       nElem=nElem)

    call myGeom%dxds%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

  endsubroutine Init_Geometry1D

  subroutine Free_Geometry1D(myGeom)
    implicit none
    class(Geometry1D),intent(inout) :: myGeom

    call myGeom%x%Free()
    call myGeom%dxds%Free()

  endsubroutine Free_Geometry1D

  subroutine GenerateFromMesh_Geometry1D(myGeom,mesh)
    ! Generates the geometry for a 1-D mesh ( set of line segments )
    ! Assumes that mesh is using Gauss-Lobatto quadrature and the degree is given by mesh % nGeo
    implicit none
    class(Geometry1D),intent(inout) :: myGeom
    type(Mesh1D),intent(in) :: mesh
    ! Local
    integer :: iel,i,nid
    type(Lagrange),target :: meshToModel
    type(Scalar1D) :: xMesh

    call meshToModel%Init(mesh%nGeo,mesh%quadrature, &
                          myGeom%x%interp%N, &
                          myGeom%x%interp%controlNodeType)

    call xMesh%Init(meshToModel, &
                    1,mesh%nElem)

    ! Set the element internal mesh locations
    nid = 1
    do iel = 1,mesh%nElem
      do i = 1,mesh%nGeo+1
        xMesh%interior(i,iel,1) = mesh%nodeCoords(nid)
        nid = nid+1
      enddo
    enddo

    ! Interpolate from the mesh hopr_nodeCoords to the geometry (Possibly not gauss_lobatto quadrature)
    call xMesh%GridInterp(myGeom%x)

    call myGeom%x%BoundaryInterp()

    call myGeom%CalculateMetricTerms()

    call xMesh%Free()

    call meshToModel%Free()

  endsubroutine GenerateFromMesh_Geometry1D

  subroutine CalculateMetricTerms_Geometry1D(myGeom)
    implicit none
    class(Geometry1D),intent(inout) :: myGeom

    call myGeom%x%Derivative(myGeom%dxds)
    call myGeom%dxds%BoundaryInterp()

  endsubroutine CalculateMetricTerms_Geometry1D

  subroutine Write_Geometry1D(myGeom,fileName)
    implicit none
    class(Geometry1D),intent(in) :: myGeom
    character(*),optional,intent(in) :: fileName
    ! Local
    integer(HID_T) :: fileId
    ! Local
    character(LEN=self_FileNameLength) :: pickupFile

    if(present(filename)) then
      pickupFile = filename
    else
      pickupFile = 'mesh.h5'
    endif

    call Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

    call CreateGroup_HDF5(fileId,'/quadrature')

    call WriteArray_HDF5(fileId,'/quadrature/xi', &
                         myGeom%x%interp%controlPoints)

    call WriteArray_HDF5(fileId,'/quadrature/weights', &
                         myGeom%x%interp%qWeights)

    call WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
                         myGeom%x%interp%dgMatrix)

    call WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
                         myGeom%x%interp%dMatrix)

    call CreateGroup_HDF5(fileId,'/mesh')

    call CreateGroup_HDF5(fileId,'/mesh/interior')

    call CreateGroup_HDF5(fileId,'/mesh/boundary')

    call WriteArray_HDF5(fileId,'/mesh/interior/x',myGeom%x%interior)

    call WriteArray_HDF5(fileId,'/mesh/interior/dxds',myGeom%dxds%interior)

    call WriteArray_HDF5(fileId,'/mesh/boundary/x',myGeom%x%boundary)

    call WriteArray_HDF5(fileId,'/mesh/boundary/dxds',myGeom%dxds%boundary)

    call Close_HDF5(fileId)

  endsubroutine Write_Geometry1D

  subroutine Init_SEMQuad(myGeom,interp,nElem)
    implicit none
    class(SEMQuad),intent(out) :: myGeom
    type(Lagrange),pointer,intent(in) :: interp
    integer,intent(in) :: nElem

    myGeom%nElem = nElem

    call myGeom%x%Init(interp=interp, &
                       nVar=1, &
                       nElem=nElem)

    call myGeom%dxds%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%dsdx%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%nHat%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%nScale%Init(interp=interp, &
                            nVar=1, &
                            nElem=nElem)

    call myGeom%J%Init(interp=interp, &
                       nVar=1, &
                       nElem=nElem)

  endsubroutine Init_SEMQuad

  subroutine Free_SEMQuad(myGeom)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom

    call myGeom%x%Free()
    call myGeom%dxds%Free()
    call myGeom%dsdx%Free()
    call myGeom%nHat%Free()
    call myGeom%nScale%Free()
    call myGeom%J%Free()

  endsubroutine Free_SEMQuad

  subroutine GenerateFromMesh_SEMQuad(myGeom,mesh)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: iel
    integer :: i,j,nid
    type(Lagrange),target :: meshToModel
    type(Vector2D) :: xMesh

    call meshToModel%Init(mesh%nGeo, &
                          mesh%quadrature, &
                          myGeom%x%interp%N, &
                          myGeom%x%interp%controlNodeType)

    call xMesh%Init(meshToModel,1,mesh%nElem)

    ! Set the element internal mesh locations
    do iel = 1,mesh%nElem
      do j = 1,mesh%nGeo+1
        do i = 1,mesh%nGeo+1
          xMesh%interior(i,j,iel,1,1:2) = mesh%nodeCoords(1:2,i,j,iel)
        enddo
      enddo
    enddo

    call xMesh%GridInterp(myGeom%x)
    call myGeom%x%BoundaryInterp()
    call myGeom%CalculateMetricTerms()

    call xMesh%Free()
    call meshToModel%Free()

  endsubroutine GenerateFromMesh_SEMQuad

  subroutine CalculateContravariantBasis_SEMQuad(myGeom)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom
    ! Local
    integer :: iEl,i,j,k
    real(prec) :: fac
    real(prec) :: mag

    ! Now calculate the contravariant basis vectors
    ! In this convention, dsdx(j,i) is contravariant vector i, component j
    ! To project onto contravariant vector i, dot vector along the first dimension
    do iEl = 1,myGeom%nElem
      do j = 1,myGeom%dxds%interp%N+1
        do i = 1,myGeom%dxds%interp%N+1

          myGeom%dsdx%interior(i,j,iel,1,1,1) = myGeom%dxds%interior(i,j,iel,1,2,2)
          myGeom%dsdx%interior(i,j,iel,1,2,1) = -myGeom%dxds%interior(i,j,iel,1,1,2)
          myGeom%dsdx%interior(i,j,iel,1,1,2) = -myGeom%dxds%interior(i,j,iel,1,2,1)
          myGeom%dsdx%interior(i,j,iel,1,2,2) = myGeom%dxds%interior(i,j,iel,1,1,1)

        enddo
      enddo
    enddo

    ! Interpolate the contravariant tensor to the boundaries
    call myGeom%dsdx%BoundaryInterp()

    ! Now, modify the sign of dsdx so that
    ! myGeom % dsdx % boundary is equal to the outward pointing normal vector
    do iEl = 1,myGeom%nElem
      do k = 1,4
        do i = 1,myGeom%J%interp%N+1
          if(k == selfSide2D_East .or. k == selfSide2D_North) then
            fac = sign(1.0_prec,myGeom%J%boundary(i,k,iEl,1))
          else
            fac = -sign(1.0_prec,myGeom%J%boundary(i,k,iEl,1))
          endif

          if(k == 1) then ! South

            mag = sqrt(myGeom%dsdx%boundary(i,k,iEl,1,1,2)**2+ &
                       myGeom%dsdx%boundary(i,k,iEl,1,2,2)**2)

            myGeom%nScale%boundary(i,k,iEl,1) = mag

            myGeom%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*myGeom%dsdx%boundary(i,k,iEl,1,1:2,2)/mag

          elseif(k == 2) then ! East

            mag = sqrt(myGeom%dsdx%boundary(i,k,iEl,1,1,1)**2+ &
                       myGeom%dsdx%boundary(i,k,iEl,1,2,1)**2)

            myGeom%nScale%boundary(i,k,iEl,1) = mag

            myGeom%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*myGeom%dsdx%boundary(i,k,iEl,1,1:2,1)/mag

          elseif(k == 3) then ! North

            mag = sqrt(myGeom%dsdx%boundary(i,k,iEl,1,1,2)**2+ &
                       myGeom%dsdx%boundary(i,k,iEl,1,2,2)**2)

            myGeom%nScale%boundary(i,k,iEl,1) = mag

            myGeom%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*myGeom%dsdx%boundary(i,k,iEl,1,1:2,2)/mag

          elseif(k == 4) then ! West

            mag = sqrt(myGeom%dsdx%boundary(i,k,iEl,1,1,1)**2+ &
                       myGeom%dsdx%boundary(i,k,iEl,1,2,1)**2)

            myGeom%nScale%boundary(i,k,iEl,1) = mag

            myGeom%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*myGeom%dsdx%boundary(i,k,iEl,1,1:2,1)/mag

          endif

          ! Set the directionality for dsdx on the boundaries
          myGeom%dsdx%boundary(i,k,iEl,1,1:2,1:2) = &
            myGeom%dsdx%boundary(i,k,iEl,1,1:2,1:2)*fac

        enddo
      enddo
    enddo

  endsubroutine CalculateContravariantBasis_SEMQuad

  subroutine CalculateMetricTerms_SEMQuad(myGeom)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom

    call myGeom%x%Gradient(myGeom%dxds)
    call myGeom%dxds%BoundaryInterp()
    call myGeom%dxds%Determinant(myGeom%J)
    call myGeom%J%BoundaryInterp()

    call myGeom%CalculateContravariantBasis()

  endsubroutine CalculateMetricTerms_SEMQuad

  ! SUBROUTINE Write_SEMQuad(myGeom,fileName)
  !   IMPLICIT NONE
  !   CLASS(SEMQuad),INTENT(in) :: myGeom
  !   CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
  !   ! Local
  !   INTEGER(HID_T) :: fileId
  !   ! Local
  !   CHARACTER(LEN=self_FileNameLength) :: pickupFile

  !   IF( PRESENT(filename) )THEN
  !     pickupFile = filename
  !   ELSE
  !     pickupFile = 'mesh.h5'
  !   ENDIF

  !   CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

  !   CALL CreateGroup_HDF5(fileId,'/quadrature')

  !   CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
  !                        myGeom % x % interp % controlPoints)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
  !                        myGeom % x % interp % qWeights)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
  !                        myGeom % x % interp % dgMatrix)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
  !                        myGeom % x % interp % dMatrix)

  !   CALL CreateGroup_HDF5(fileId,'/mesh')

  !   CALL CreateGroup_HDF5(fileId,'/mesh/interior')

  !   CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/x',myGeom % x % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/dxds',myGeom % dxds % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/dsdx',myGeom % dsdx % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/J',myGeom % J % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/x',myGeom % x % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/dxds',myGeom % dxds % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/dsdx',myGeom % dsdx % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/nHat',myGeom % nHat % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/nScale',myGeom % nScale % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/J',myGeom % J % boundary)

  !   CALL Close_HDF5(fileId)

  ! END SUBROUTINE Write_SEMQuad

  subroutine Init_SEMHex(myGeom,interp,nElem)
    implicit none
    class(SEMHex),intent(out) :: myGeom
    type(Lagrange),pointer,intent(in) :: interp
    integer,intent(in) :: nElem

    myGeom%nElem = nElem

    call myGeom%x%Init(interp=interp, &
                       nVar=1, &
                       nElem=nElem)

    call myGeom%dxds%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%dsdx%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%nHat%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%nScale%Init(interp=interp, &
                            nVar=1, &
                            nElem=nElem)

    call myGeom%J%Init(interp=interp, &
                       nVar=1, &
                       nElem=nElem)

  endsubroutine Init_SEMHex

  subroutine Free_SEMHex(myGeom)
    implicit none
    class(SEMHex),intent(inout) :: myGeom

    call myGeom%x%Free()
    call myGeom%dxds%Free()
    call myGeom%dsdx%Free()
    call myGeom%nHat%Free()
    call myGeom%nScale%Free()
    call myGeom%J%Free()

  endsubroutine Free_SEMHex

  subroutine GenerateFromMesh_SEMHex(myGeom,mesh)
    implicit none
    class(SEMHex),intent(inout) :: myGeom
    type(Mesh3D),intent(in) :: mesh
    ! Local
    integer :: iel
    integer :: i,j,k,nid
    type(Lagrange),target :: meshToModel
    type(Vector3D) :: xMesh

    call meshToModel%Init(mesh%nGeo,mesh%quadrature, &
                          myGeom%x%interp%N, &
                          myGeom%x%interp%controlNodeType)

    call xMesh%Init(meshToModel, &
                    1,mesh%nElem)

    ! Set the element internal mesh locations
    do iel = 1,mesh%nElem
      do k = 1,mesh%nGeo+1
        do j = 1,mesh%nGeo+1
          do i = 1,mesh%nGeo+1
            xMesh%interior(i,j,k,iel,1,1:3) = mesh%nodeCoords(1:3,i,j,k,iel)
          enddo
        enddo
      enddo
    enddo

    call xMesh%GridInterp(myGeom%x)
    call myGeom%x%BoundaryInterp()
    call myGeom%CalculateMetricTerms()

    call xMesh%Free()
    call meshToModel%Free()

  endsubroutine GenerateFromMesh_SEMHex

  subroutine CheckSides_SEMHex(myGeom,mesh)
    implicit none
    class(SEMHex),intent(in) :: myGeom
    type(Mesh3D),intent(in) :: mesh
    !
    integer :: e1,s1
    integer :: e2,s2
    integer :: i1,j1
    integer :: i2,j2
    integer :: flip,bcid
    real(prec) :: rms

    do e1 = 1,mesh%nElem
      do s1 = 1,6

        e2 = mesh%sideInfo(3,s1,e1)
        s2 = mesh%sideInfo(4,s1,e1)/10
        flip = mesh%sideInfo(4,s1,e1)-s2*10
        bcid = mesh%sideInfo(5,s1,e1)

        if(bcid == 0) then ! Interior

          rms = 0.0_prec

          if(flip == 0) then

            do j1 = 1,myGeom%x%interp%N+1
              do i1 = 1,myGeom%x%interp%N+1
                rms = rms+ &
                      sqrt((myGeom%x%boundary(i1,j1,s1,e1,1,1)- &
                            myGeom%x%boundary(i1,j1,s2,e2,1,1))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,2)- &
                            myGeom%x%boundary(i1,j1,s2,e2,1,2))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,3)- &
                            myGeom%x%boundary(i1,j1,s2,e2,1,3))**2)
              enddo
            enddo

          elseif(flip == 1) then

            do j1 = 1,myGeom%x%interp%N+1
              do i1 = 1,myGeom%x%interp%N+1
                i2 = j1
                j2 = myGeom%x%interp%N-i1
                rms = rms+ &
                      sqrt((myGeom%x%boundary(i1,j1,s1,e1,1,1)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,1))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,2)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,2))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,3)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,3))**2)
              enddo
            enddo

          elseif(flip == 2) then

            do j1 = 1,myGeom%x%interp%N+1
              do i1 = 1,myGeom%x%interp%N+1
                i2 = myGeom%x%interp%N-i1
                j2 = myGeom%x%interp%N-j1
                rms = rms+ &
                      sqrt((myGeom%x%boundary(i1,j1,s1,e1,1,1)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,1))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,2)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,2))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,3)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,3))**2)
              enddo
            enddo

          elseif(flip == 3) then

            do j1 = 1,myGeom%x%interp%N+1
              do i1 = 1,myGeom%x%interp%N+1
                i2 = myGeom%x%interp%N-j1
                j2 = i1
                rms = rms+ &
                      sqrt((myGeom%x%boundary(i1,j1,s1,e1,1,1)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,1))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,2)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,2))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,3)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,3))**2)
              enddo
            enddo

          elseif(flip == 4) then

            do j1 = 1,myGeom%x%interp%N+1
              do i1 = 1,myGeom%x%interp%N+1
                i2 = j1
                j2 = i1
                rms = rms+ &
                      sqrt((myGeom%x%boundary(i1,j1,s1,e1,1,1)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,1))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,2)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,2))**2+ &
                           (myGeom%x%boundary(i1,j1,s1,e1,1,3)- &
                            myGeom%x%boundary(i2,j2,s2,e2,1,3))**2)
              enddo
            enddo
          endif

        endif

      enddo
    enddo

  endsubroutine CheckSides_SEMHex

  subroutine CalculateContravariantBasis_SEMHex(myGeom)
    implicit none
    class(SEMHex),intent(inout) :: myGeom
    ! Local
    integer :: iEl,i,j,k
    real(prec) :: fac
    real(prec) :: mag

    ! Now calculate the contravariant basis vectors
    ! In this convention, dsdx(j,i) is contravariant vector i, component j
    ! To project onto contravariant vector i, dot vector along the first dimension
    ! TO DO : Curl Invariant Form
    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1

            ! Ja1
            myGeom%dsdx%interior(i,j,k,iel,1,1,1) = &
              myGeom%dxds%interior(i,j,k,iel,1,2,2)* &
              myGeom%dxds%interior(i,j,k,iel,1,3,3)- &
              myGeom%dxds%interior(i,j,k,iel,1,3,2)* &
              myGeom%dxds%interior(i,j,k,iel,1,2,3)

            myGeom%dsdx%interior(i,j,k,iel,1,2,1) = &
              myGeom%dxds%interior(i,j,k,iel,1,1,3)* &
              myGeom%dxds%interior(i,j,k,iel,1,3,2)- &
              myGeom%dxds%interior(i,j,k,iel,1,3,3)* &
              myGeom%dxds%interior(i,j,k,iel,1,1,2)

            myGeom%dsdx%interior(i,j,k,iel,1,3,1) = &
              myGeom%dxds%interior(i,j,k,iel,1,1,2)* &
              myGeom%dxds%interior(i,j,k,iel,1,2,3)- &
              myGeom%dxds%interior(i,j,k,iel,1,2,2)* &
              myGeom%dxds%interior(i,j,k,iel,1,1,3)

            ! Ja2
            myGeom%dsdx%interior(i,j,k,iel,1,1,2) = &
              myGeom%dxds%interior(i,j,k,iel,1,2,3)* &
              myGeom%dxds%interior(i,j,k,iel,1,3,1)- &
              myGeom%dxds%interior(i,j,k,iel,1,3,3)* &
              myGeom%dxds%interior(i,j,k,iel,1,2,1)

            myGeom%dsdx%interior(i,j,k,iel,1,2,2) = &
              myGeom%dxds%interior(i,j,k,iel,1,1,1)* &
              myGeom%dxds%interior(i,j,k,iel,1,3,3)- &
              myGeom%dxds%interior(i,j,k,iel,1,3,1)* &
              myGeom%dxds%interior(i,j,k,iel,1,1,3)

            myGeom%dsdx%interior(i,j,k,iel,1,3,2) = &
              myGeom%dxds%interior(i,j,k,iel,1,1,3)* &
              myGeom%dxds%interior(i,j,k,iel,1,2,1)- &
              myGeom%dxds%interior(i,j,k,iel,1,2,3)* &
              myGeom%dxds%interior(i,j,k,iel,1,1,1)

            ! Ja3
            myGeom%dsdx%interior(i,j,k,iel,1,1,3) = &
              myGeom%dxds%interior(i,j,k,iel,1,2,1)* &
              myGeom%dxds%interior(i,j,k,iel,1,3,2)- &
              myGeom%dxds%interior(i,j,k,iel,1,3,1)* &
              myGeom%dxds%interior(i,j,k,iel,1,2,2)

            myGeom%dsdx%interior(i,j,k,iel,1,2,3) = &
              myGeom%dxds%interior(i,j,k,iel,1,1,2)* &
              myGeom%dxds%interior(i,j,k,iel,1,3,1)- &
              myGeom%dxds%interior(i,j,k,iel,1,3,2)* &
              myGeom%dxds%interior(i,j,k,iel,1,1,1)

            myGeom%dsdx%interior(i,j,k,iel,1,3,3) = &
              myGeom%dxds%interior(i,j,k,iel,1,1,1)* &
              myGeom%dxds%interior(i,j,k,iel,1,2,2)- &
              myGeom%dxds%interior(i,j,k,iel,1,2,1)* &
              myGeom%dxds%interior(i,j,k,iel,1,1,2)

          enddo
        enddo
      enddo
    enddo

    ! Interpolate the contravariant tensor to the boundaries
    call myGeom%dsdx%BoundaryInterp()

    ! Now, calculate nHat (outward pointing normal)
    do iEl = 1,myGeom%nElem
      do k = 1,6
        do j = 1,myGeom%J%interp%N+1
          do i = 1,myGeom%J%interp%N+1
            if(k == selfSide3D_Top .or. k == selfSide3D_East .or. k == selfSide3D_North) then
              fac = sign(1.0_prec,myGeom%J%boundary(i,j,k,iEl,1))
            else
              fac = -sign(1.0_prec,myGeom%J%boundary(i,j,k,iEl,1))
            endif

            if(k == 1) then ! Bottom

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,3)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)/mag

            elseif(k == 2) then ! South

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,2)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)/mag

            elseif(k == 3) then ! East

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,1)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)/mag

            elseif(k == 4) then ! North

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,2)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)/mag

            elseif(k == 5) then ! West

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,1)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)/mag

            elseif(k == 6) then ! Top

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,3)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)/mag

            endif

            ! Set the directionality for dsdx on the boundaries
            ! This is primarily used for DG gradient calculations,
            ! which do not use nHat for the boundary terms.
            myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1:3) = &
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1:3)*fac

          enddo
        enddo
      enddo
    enddo

  endsubroutine CalculateContravariantBasis_SEMHex

  subroutine CalculateMetricTerms_SEMHex(myGeom)
    implicit none
    class(SEMHex),intent(inout) :: myGeom

    call myGeom%x%Gradient(myGeom%dxds)
    call myGeom%dxds%BoundaryInterp()
    call myGeom%dxds%Determinant(myGeom%J)
    call myGeom%J%BoundaryInterp()

    call myGeom%CalculateContravariantBasis()

  endsubroutine CalculateMetricTerms_SEMHex

  ! SUBROUTINE Write_SEMHex(myGeom,fileName)
  !   IMPLICIT NONE
  !   CLASS(SEMHex),INTENT(in) :: myGeom
  !   CHARACTER(*),OPTIONAL,INTENT(in) :: fileName
  !   ! Local
  !   INTEGER(HID_T) :: fileId
  !   ! Local
  !   CHARACTER(LEN=self_FileNameLength) :: pickupFile

  !   IF( PRESENT(filename) )THEN
  !     pickupFile = filename
  !   ELSE
  !     pickupFile = 'mesh.h5'
  !   ENDIF

  !   CALL Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

  !   CALL CreateGroup_HDF5(fileId,'/quadrature')

  !   CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
  !                        myGeom % x % interp % controlPoints)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
  !                        myGeom % x % interp % qWeights)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
  !                        myGeom % x % interp % dgMatrix)

  !   CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
  !                        myGeom % x % interp % dMatrix)

  !   CALL CreateGroup_HDF5(fileId,'/mesh')

  !   CALL CreateGroup_HDF5(fileId,'/mesh/interior')

  !   CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/x',myGeom % x % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/dxds',myGeom % dxds % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/dsdx',myGeom % dsdx % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/interior/J',myGeom % J % interior)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/x',myGeom % x % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/dxds',myGeom % dxds % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/dsdx',myGeom % dsdx % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/nHat',myGeom % nHat % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/nScale',myGeom % nScale % boundary)

  !   CALL WriteArray_HDF5(fileId,'/mesh/boundary/J',myGeom % J % boundary)

  !   CALL Close_HDF5(fileId)

  ! END SUBROUTINE Write_SEMHex

endmodule SELF_Geometry
