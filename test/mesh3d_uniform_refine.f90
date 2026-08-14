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

  implicit none
  integer :: exit_code

  exit_code = mesh3d_uniform_refine()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mesh3d_uniform_refine() result(r)
    !! Validates uniform 3-D mesh refinement (SELF_MeshRefinement_3D):
    !!
    !!   1. The refined mesh has 8x the elements and remains conforming (no mortars).
    !!   2. Its geometry has strictly positive Jacobians and integrates to the same total
    !!      volume as the base mesh.
    !!   3. For a globally linear field, SideExchange on the refined mesh reproduces the
    !!      field's trace in extBoundary on every interior face to roundoff, which
    !!      exercises the inherited neighbor/flip connectivity (including the sibling
    !!      face pairs and the quadrant pairing across parent faces) directly.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D
    use SELF_MeshRefinement_3D

    implicit none

    integer,parameter :: controlDegree = 4
    integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-11_prec
#else
    real(prec),parameter :: tolerance = 1.0e-4_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: baseMesh,mesh
    type(SEMHex),target :: baseGeom,geometry
    type(MappedScalar3D) :: f
    integer :: bcids(1:6)
    integer :: i,j,k,iel,is
    real(prec) :: x,y,z,maxErr,baseVolume,volume

    r = 0
    bcids(1:6) = SELF_BC_PRESCRIBED

    call baseMesh%StructuredMesh(2,2,2,1,1,1,0.5_prec,0.5_prec,0.5_prec,bcids)

    call UniformRefineMesh(baseMesh,mesh)

    ! 1. Element count and conformity.
    if(mesh%nElem /= 8*baseMesh%nElem) then
      print*,"FAIL: refined element count",mesh%nElem
      r = 1
    endif
    if(mesh%nMortars /= 0) then
      print*,"FAIL: uniform refinement produced mortars",mesh%nMortars
      r = 1
    endif

    ! 2. Geometry: positive Jacobians, volume conservation.
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree, &
                     targetNodeType=UNIFORM)
    call baseGeom%Init(interp,baseMesh%nElem)
    call baseGeom%GenerateFromMesh(baseMesh)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    if(minval(geometry%J%interior) <= 0.0_prec) then
      print*,"FAIL: non-positive Jacobian on refined mesh"
      r = 1
    endif
    baseVolume = integrateJ(baseGeom,interp)
    volume = integrateJ(geometry,interp)
    print*,"base volume, refined volume :",baseVolume,volume
    if(abs(baseVolume-volume) > tolerance*abs(baseVolume)) then
      print*,"FAIL: refined mesh does not preserve volume"
      r = 1
    endif

    ! 3. Conforming exchange of a globally linear field on the refined connectivity.
    call f%Init(interp,nvar,mesh%nElem)
    call f%AssociateGeometry(geometry)
    do iel = 1,mesh%nElem
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            x = geometry%x%interior(i,j,k,iel,1,1)
            y = geometry%x%interior(i,j,k,iel,1,2)
            z = geometry%x%interior(i,j,k,iel,1,3)
            f%interior(i,j,k,iel,1) = 1.0_prec+0.8_prec*x+0.5_prec*y+0.3_prec*z
          enddo
        enddo
      enddo
    enddo
    call f%UpdateDevice()
    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%UpdateHost()

    maxErr = 0.0_prec
    do iel = 1,mesh%nElem
      do is = 1,6
        if(mesh%sideInfo(3,is,iel) > 0) then
          do j = 1,controlDegree+1
            do i = 1,controlDegree+1
              maxErr = max(maxErr,abs(f%extBoundary(i,j,is,iel,1)- &
                                      f%boundary(i,j,is,iel,1)))
            enddo
          enddo
        endif
      enddo
    enddo
    print*,"refined-mesh exchange absmax error :",maxErr
    if(maxErr > tolerance) then
      print*,"FAIL: exchange did not reproduce the linear field on interior faces"
      r = 1
    endif

    if(r == 0) print*,"MESH3D UNIFORM REFINE CHECKS PASSED"

    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call baseGeom%Free()
    call interp%Free()
    call mesh%Free()
    call baseMesh%Free()

  endfunction mesh3d_uniform_refine

  function integrateJ(geom,interp) result(volume)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Geometry_3D
    implicit none
    type(SEMHex),intent(in) :: geom
    type(Lagrange),intent(in) :: interp
    real(prec) :: volume
    integer :: e,i,j,k

    volume = 0.0_prec
    do e = 1,geom%nElem
      do k = 1,interp%N+1
        do j = 1,interp%N+1
          do i = 1,interp%N+1
            volume = volume+geom%J%interior(i,j,k,e,1)* &
                     interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)
          enddo
        enddo
      enddo
    enddo

  endfunction integrateJ

endprogram test
