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

module SELF_Geometry_2D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use SELF_Scalar_2D
  use SELF_Vector_2D
  use SELF_Tensor_2D
  use SELF_SupportRoutines
  use SELF_Mesh_2D

  implicit none

#include "SELF_Macros.h"

  type,public :: SEMQuad
    type(Vector2D) :: x ! Physical positions
    type(Tensor2D) :: dxds ! Covariant basis vectors
    type(Tensor2D) :: dsdx ! Contavariant basis vectors
    type(Vector2D) :: nHat ! Normal Vectors pointing across coordinate lines
    type(Scalar2D) :: nScale ! Boundary scale
    type(Scalar2D) :: J ! Jacobian of the transformation
    integer :: nElem
  contains

    procedure,public :: Init => Init_SEMQuad
    procedure,public :: Free => Free_SEMQuad
    procedure,public :: GenerateFromMesh => GenerateFromMesh_SEMQuad
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_SEMQuad
    procedure,private :: CalculateContravariantBasis => CalculateContravariantBasis_SEMQuad

  endtype SEMQuad

contains

  subroutine Init_SEMQuad(myGeom,interp,nElem)
    implicit none
    class(SEMQuad),intent(out) :: myGeom
    type(Lagrange),pointer,intent(in) :: interp
    integer,intent(in) :: nElem

    myGeom%nElem = nElem

    call myGeom%x%Init(interp=interp, &
                       nVar=1, &
                       nElem=nElem)

    call myGeom%x%meta(1)%SetName("x")

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
    integer :: i,j
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

    call xMesh%GridInterp(myGeom%x%interior)
    call myGeom%x%UpdateDevice()
    call myGeom%x%BoundaryInterp() ! Boundary interp will run on GPU if enabled, hence why we close in update host/device
    call myGeom%x%UpdateHost()
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
    call myGeom%dsdx%BoundaryInterp() ! Tensor boundary interp is not offloaded

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

    call myGeom%dsdx%UpdateDevice()
    call myGeom%nHat%UpdateDevice()
    call myGeom%nScale%UpdateDevice()

  endsubroutine CalculateContravariantBasis_SEMQuad

  subroutine CalculateMetricTerms_SEMQuad(myGeom)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom

    call myGeom%x%Gradient(myGeom%dxds%interior)
    call myGeom%dxds%BoundaryInterp() ! Tensor boundary interp is not offloaded to GPU
    call myGeom%dxds%UpdateDevice()

    call myGeom%dxds%Determinant(myGeom%J%interior)

    call myGeom%J%UpdateDevice()
    call myGeom%J%BoundaryInterp()
    call myGeom%J%UpdateHost()

    call myGeom%CalculateContravariantBasis()

  endsubroutine CalculateMetricTerms_SEMQuad

endmodule SELF_Geometry_2D
