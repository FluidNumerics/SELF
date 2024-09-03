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

module SELF_Geometry_1D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use SELF_Scalar_1D
  use SELF_SupportRoutines
  use SELF_Mesh_1D

  implicit none

#include "SELF_Macros.h"

  type,public :: Geometry1D
    type(Scalar1D) :: x ! Physical Positions
    type(Scalar1D) :: dxds ! Conversion from computational to physical space
    integer :: nElem

  contains

    procedure,public :: Init => Init_Geometry1D
    procedure,public :: Free => Free_Geometry1D
    procedure,public :: GenerateFromMesh => GenerateFromMesh_Geometry1D
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_Geometry1D

    procedure :: write => Write_Geometry1D

  endtype Geometry1D

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
    call xMesh%GridInterp(myGeom%x%interior)
    call myGeom%x%UpdateDevice()
    call myGeom%x%BoundaryInterp()

    call myGeom%CalculateMetricTerms()

    call xMesh%Free()

    call meshToModel%Free()

  endsubroutine GenerateFromMesh_Geometry1D

  subroutine CalculateMetricTerms_Geometry1D(myGeom)
    implicit none
    class(Geometry1D),intent(inout) :: myGeom

    call myGeom%x%Derivative(myGeom%dxds%interior)
    call myGeom%dxds%UpdateDevice()
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

endmodule SELF_Geometry_1D
