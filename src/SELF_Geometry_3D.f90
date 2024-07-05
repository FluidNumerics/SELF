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

module SELF_Geometry_3D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use SELF_Scalar_3D
  use SELF_Vector_3D
  use SELF_Tensor_3D
  use SELF_SupportRoutines
  use SELF_Mesh_3D

  implicit none

#include "SELF_Macros.h"

  type,public :: SEMHex
    type(Vector3D) :: x ! Physical positions
    type(Tensor3D) :: dxds ! Covariant basis vectors
    type(Tensor3D) :: dsdx ! Contavariant basis vectors
    type(Vector3D) :: nHat ! Normal Vectors pointing across coordinate lines
    type(Scalar3D) :: nScale ! Boundary scale
    type(Scalar3D) :: J ! Jacobian of the transformation
    integer :: nElem

  contains

    procedure,public :: Init => Init_SEMHex
    procedure,public :: Free => Free_SEMHex
    procedure,public :: GenerateFromMesh => GenerateFromMesh_SEMHex
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_SEMHex
    procedure,private :: CalculateContravariantBasis => CalculateContravariantBasis_SEMHex
    procedure,private :: CheckSides => CheckSides_SEMHex
    !PROCEDURE,public :: Write => Write_SEMHex

  endtype SEMHex

contains

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
    integer :: i,j,k
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

    myGeom%x%interior = xMesh%GridInterp()
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
    type(Vector3D) :: xlgradxm,xmgradxl
    type(Vector3D) :: curl_xlgradxm,curl_xmgradxl

    ! Here we use the curl invariant form from Kopriva (2006)
    ! to calculate the contravariant basis vectors
    call xlgradxm%Init(myGeom%x%interp,1,myGeom%x%nElem)
    call xmgradxl%Init(myGeom%x%interp,1,myGeom%x%nElem)

    call curl_xlgradxm%Init(myGeom%x%interp,1,myGeom%x%nElem)
    call curl_xmgradxl%Init(myGeom%x%interp,1,myGeom%x%nElem)

    ! Ja^{1:3}_1 (n=1, m=2, l=3) First component of the contravariant basis vectors
    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,3)*myGeom%dxds%interior(i,j,k,iel,1,2,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=3,m=2
            xmgradxl%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,2)*myGeom%dxds%interior(i,j,k,iel,1,3,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=3,m=2
          enddo
        enddo
      enddo
    enddo

    curl_xlgradxm%interior = xlgradxm%Curl()
    curl_xmgradxl%interior = xmgradxl%Curl()

    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=1
            myGeom%dsdx%interior(i,j,k,iel,1,1,1:3) = 0.5_prec*( &
                                                      curl_xmgradxl%interior(i,j,k,iel,1,1:3)- &
                                                      curl_xlgradxm%interior(i,j,k,iel,1,1:3))
          enddo
        enddo
      enddo
    enddo

    ! Ja^{1:3}_2 (n=2, m=3, l=1) Second component of the contravariant basis vectors
    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,1)*myGeom%dxds%interior(i,j,k,iel,1,3,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=1,m=3
            xmgradxl%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,3)*myGeom%dxds%interior(i,j,k,iel,1,1,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=1,m=3
          enddo
        enddo
      enddo
    enddo

    curl_xlgradxm%interior = xlgradxm%Curl()
    curl_xmgradxl%interior = xmgradxl%Curl()

    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=2
            myGeom%dsdx%interior(i,j,k,iel,1,2,1:3) = 0.5_prec*( &
                                                      curl_xmgradxl%interior(i,j,k,iel,1,1:3)- &
                                                      curl_xlgradxm%interior(i,j,k,iel,1,1:3))
          enddo
        enddo
      enddo
    enddo

    ! Ja^{1:3}_3 (n=3, m=1, l=2) Third component of the contravariant basis vectors
    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,2)*myGeom%dxds%interior(i,j,k,iel,1,1,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=2,m=1
            xmgradxl%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,1)*myGeom%dxds%interior(i,j,k,iel,1,2,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=2,m=1
          enddo
        enddo
      enddo
    enddo

    curl_xlgradxm%interior = xlgradxm%Curl()
    curl_xmgradxl%interior = xmgradxl%Curl()

    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=3
            myGeom%dsdx%interior(i,j,k,iel,1,3,1:3) = 0.5_prec*( &
                                                      curl_xmgradxl%interior(i,j,k,iel,1,1:3)- &
                                                      curl_xlgradxm%interior(i,j,k,iel,1,1:3))
          enddo
        enddo
      enddo
    enddo

    call xlgradxm%Free()
    call xmgradxl%Free()
    call curl_xlgradxm%Free()
    call curl_xmgradxl%Free()

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
              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)*fac

            elseif(k == 2) then ! South

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,2)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)*fac

            elseif(k == 3) then ! East

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,1)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)/mag
              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)*fac

            elseif(k == 4) then ! North

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,2)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)*fac

            elseif(k == 5) then ! West

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,1)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)*fac

            elseif(k == 6) then ! Top

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,3)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)*fac

            endif

          enddo
        enddo
      enddo
    enddo

  endsubroutine CalculateContravariantBasis_SEMHex

  subroutine CalculateMetricTerms_SEMHex(myGeom)
    implicit none
    class(SEMHex),intent(inout) :: myGeom

    myGeom%dxds%interior = myGeom%x%Gradient()
    call myGeom%dxds%BoundaryInterp()
    call myGeom%dxds%Determinant(myGeom%J%interior)
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

endmodule SELF_Geometry_3D
