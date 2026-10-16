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

program test
!! DG divergence of a linear vector field on extruded meshes (Mesh3D%ExtrudedMesh)
!! with terrain-following levels over a sloping bathymetry.
!!
!! The field is f = (1.3x + 0.4y - 0.7z, -0.6x + 0.8y + 0.5z, 0.2x - 0.9y + 1.1z), whose
!! divergence is 3.2. The interface flux is the average of the two traces, which is
!! exact for a continuous field, and on the domain boundary the exterior state is set
!! to the interior trace. Two surface meshes are extruded:
!!
!!   - a structured (affine, nGeo = 1) mesh. The contravariant fluxes J a^i . f are
!!     polynomials of degree at most 3 per direction there, so with N = 7 the weak
!!     form is integrated exactly and the divergence must be exact to round-off
!!     (`exactTolerance`).
!!   - the curved Circle mesh, from HOHQMesh (nGeo = 7) and from HOPr (nGeo = 6). The
!!     fluxes J a^i . f are no longer resolved at N = 7, so the divergence carries a
!!     truncation error, which must be below `curvedTolerance`.
!!
!! This is the MPI variant: it extrudes the structured and HOPr surfaces only, because
!! the HOHQMesh reader keeps every element on every rank instead of decomposing the
!! mesh. The lateral faces between columns on different ranks are exchanged over MPI.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_MappedScalar_3D
  use SELF_MappedVector_3D
  use mpi

  implicit none
  integer :: exit_code

  exit_code = mappedvectordgdivergence_3d_extruded()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  real(prec) function DivergenceError(surf,N) result(errGlobal)
    !! Maximum nodal error, over all ranks, of the DG divergence of f on surf
    !! extruded through three sigma layers.
    implicit none
    type(Mesh2D),intent(in) :: surf
    integer,intent(in) :: N
    ! Local
    integer,parameter :: nz = 3
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedVector3D) :: f
    type(MappedScalar3D) :: df
    real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
    integer :: nGeo,e,i,j,k,iel,e2,ierror
    real(prec) :: x,y,nhat(1:3),nmag,fx,fy,fz,err

    nGeo = surf%nGeo
    allocate(h(1:nGeo+1,1:nGeo+1,1:surf%nElem))
    allocate(zLevels(1:nGeo+1,1:nGeo+1,0:nz,1:surf%nElem))
    do e = 1,surf%nElem
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          x = surf%nodeCoords(1,i,j,e)
          y = surf%nodeCoords(2,i,j,e)
          h(i,j,e) = 1.0_prec+0.3_prec*x+0.2_prec*y*y
        enddo
      enddo
    enddo
    call ExtrudedSigmaLevels(h,[-1.0_prec,-0.55_prec,-0.2_prec,0.0_prec],zLevels)
    call mesh%ExtrudedMesh(surf,zLevels,SELF_BC_NONORMALFLOW,SELF_BC_PRESCRIBED)

    call interp%Init(N=N,controlNodeType=GAUSS,M=N+3,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,1,mesh%nElem)
    call df%Init(interp,1,mesh%nElem)
    call f%AssociateGeometry(geometry)
    call f%SetEquation(1,1,'f = 1.3*x + 0.4*y - 0.7*z')
    call f%SetEquation(2,1,'f = -0.6*x + 0.8*y + 0.5*z')
    ! feq-parse reads a - b + c as a - (b + c), so a subtraction is written last
    call f%SetEquation(3,1,'f = 0.2*x + 1.1*z - 0.9*y')
    call f%SetInteriorFromEquation(geometry,0.0_prec)
    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%UpdateHost()

    do iel = 1,f%nElem
      do k = 1,6
        e2 = mesh%sideInfo(3,k,iel)
        if(e2 == 0) then
          f%extBoundary(:,:,k,iel,1,1:3) = f%boundary(:,:,k,iel,1,1:3)
        endif
        do j = 1,N+1
          do i = 1,N+1
            nhat(1:3) = geometry%nHat%boundary(i,j,k,iel,1,1:3)
            nmag = geometry%nScale%boundary(i,j,k,iel,1)
            fx = 0.5_prec*(f%boundary(i,j,k,iel,1,1)+f%extBoundary(i,j,k,iel,1,1))
            fy = 0.5_prec*(f%boundary(i,j,k,iel,1,2)+f%extBoundary(i,j,k,iel,1,2))
            fz = 0.5_prec*(f%boundary(i,j,k,iel,1,3)+f%extBoundary(i,j,k,iel,1,3))
            f%boundaryNormal(i,j,k,iel,1) = (fx*nhat(1)+fy*nhat(2)+fz*nhat(3))*nmag
          enddo
        enddo
      enddo
    enddo
    call f%UpdateDevice()

#ifdef ENABLE_GPU
    call f%MappedDGDivergence(df%interior_gpu)
#else
    call f%MappedDGDivergence(df%interior)
#endif
    call df%UpdateHost()

    err = maxval(abs(df%interior-3.2_prec))
    call MPI_Allreduce(err,errGlobal,1,mesh%decomp%mpiPrec,MPI_MAX,mesh%decomp%mpiComm,ierror)

    call f%DissociateGeometry()
    call f%Free()
    call df%Free()
    call geometry%Free()
    call interp%Free()
    call mesh%Free()
    deallocate(h,zLevels)

  endfunction DivergenceError

  integer function mappedvectordgdivergence_3d_extruded() result(r)
    implicit none
    ! Local
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: exactTolerance = 1.0e-10_prec
    real(prec),parameter :: curvedTolerance = 1.0e-4_prec
#else
    real(prec),parameter :: exactTolerance = 1.0e-3_prec
    real(prec),parameter :: curvedTolerance = 1.0e-2_prec
#endif
    type(Mesh2D) :: structured,circleHOPr
    character(LEN=255) :: WORKSPACE
    real(prec) :: errStructured,errCircleHOPr
    integer :: rankId,ierror

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call structured%StructuredMesh(3,2,2,2,0.3_prec,0.25_prec, &
                                   [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED])
    call circleHOPr%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Circle/Circle_mesh.h5")

    errStructured = DivergenceError(structured,7)
    errCircleHOPr = DivergenceError(circleHOPr,7)

    call MPI_Comm_rank(structured%decomp%mpiComm,rankId,ierror)
    if(rankId == 0) then
      print*,"structured surface, sigma levels: max |div f - 3.2| = ",errStructured," (tolerance ",exactTolerance,")"
      print*,"HOPr Circle, sigma levels       : max |div f - 3.2| = ",errCircleHOPr," (tolerance ",curvedTolerance,")"
    endif

    r = 0
    if(errStructured > exactTolerance) r = 1
    if(errCircleHOPr > curvedTolerance) r = 1

    call structured%Free()
    call circleHOPr%Free()

  endfunction mappedvectordgdivergence_3d_extruded

endprogram test
