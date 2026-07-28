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

  exit_code = adaptive_mortar_2d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function adaptive_mortar_2d() result(r)
    !! End-to-end validation of the AMR pipeline through Stage 4b: adapt a forest, balance it, emit
    !! a nonconforming Mesh2D_t, and exercise the emitted connectivity with the existing solver
    !! machinery.
    !!
    !! A structured 2x2 base mesh is refined adaptively (one corner element to two levels) and
    !! 2:1-balanced, then EmitMesh produces the mesh (2:1 mortar interfaces at the level jumps).
    !! On this mesh:
    !!
    !!   1. The geometry built from the emitted node coordinates has strictly positive Jacobians
    !!      and integrates to the same total area as the base mesh (the refined mesh is the same
    !!      domain).
    !!   2. For a globally linear field, MortarExchange + SideExchange reproduce the field's trace
    !!      in extBoundary on every interior side - both the conforming sides and the 2:1 mortar
    !!      sides - to roundoff. This is the same criterion as the hand-built mortar-mesh tests and
    !!      exercises the emitted sideInfo and mortarInfo tables directly.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_QuadTreeMesh_2D
    use SELF_AdaptiveMesh_2D

    implicit none

    integer,parameter :: controlDegree = 4
    integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-10_prec
#else
    real(prec),parameter :: tolerance = 1.0e-4_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: baseMesh,mesh
    type(SEMQuad),target :: baseGeom,geometry
    type(QuadTreeMesh2D) :: forest
    type(MappedScalar2D) :: f
    integer :: bcids(1:4)
    integer :: i,j,iel,is
    real(prec) :: x,y,maxErr,baseArea,area
    integer,allocatable :: flag(:)

    r = 0
    bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    call baseMesh%StructuredMesh(2,2,1,1,0.5_prec,0.5_prec,bcids)

    ! Adapt: refine one corner element two levels, then 2:1-balance.
    call forest%Init(baseMesh)
    allocate(flag(1:forest%nLeaves)); flag = QUADTREE_KEEP; flag(1) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag); deallocate(flag)
    allocate(flag(1:forest%nLeaves)); flag = QUADTREE_KEEP; flag(1:4) = QUADTREE_REFINE
    call forest%AdaptFromFlags(flag); deallocate(flag)
    call forest%Balance2to1()

    call EmitMesh(forest,baseMesh,mesh)
    print*,"emitted nElem, nMortars :",mesh%nElem,mesh%nMortars
    if(mesh%nElem /= forest%nLeaves .or. mesh%nMortars <= 0) then
      print*,"FAIL: emitted mesh element/mortar count"
      r = 1
    endif

    ! Geometry on the base and emitted meshes.
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree,targetNodeType=UNIFORM)
    call baseGeom%Init(interp,baseMesh%nElem)
    call baseGeom%GenerateFromMesh(baseMesh)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    ! 1. Positive Jacobian and area conservation.
    if(minval(geometry%J%interior) <= 0.0_prec) then
      print*,"FAIL: non-positive Jacobian on emitted mesh"
      r = 1
    endif
    baseArea = integrateJ(baseGeom,interp)
    area = integrateJ(geometry,interp)
    print*,"base area, emitted area :",baseArea,area
    if(abs(baseArea-area) > tolerance*abs(baseArea)) then
      print*,"FAIL: emitted mesh does not preserve area"
      r = 1
    endif

    ! 2. Mortar/conforming exchange of a globally linear field.
    call f%Init(interp,nvar,mesh%nElem)
    call f%AssociateGeometry(geometry)
    do iel = 1,mesh%nElem
      do j = 1,controlDegree+1
        do i = 1,controlDegree+1
          x = geometry%x%interior(i,j,iel,1,1)
          y = geometry%x%interior(i,j,iel,1,2)
          f%interior(i,j,iel,1) = 1.0_prec+0.8_prec*x+0.5_prec*y
        enddo
      enddo
    enddo
    call f%UpdateDevice()
    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%MortarExchange(mesh)
    call f%UpdateHost()

    ! On every interior side (conforming: sideInfo(3)>0, or mortar: sideInfo(1)>0) the external
    ! trace of a linear field must equal this element's own trace.
    maxErr = 0.0_prec
    do iel = 1,mesh%nElem
      do is = 1,4
        if(mesh%sideInfo(3,is,iel) > 0 .or. mesh%sideInfo(1,is,iel) > 0) then
          do i = 1,controlDegree+1
            maxErr = max(maxErr,abs(f%extBoundary(i,is,iel,1)-f%boundary(i,is,iel,1)))
          enddo
        endif
      enddo
    enddo
    print*,"adaptive mortar/conforming exchange absmax error :",maxErr
    if(maxErr > tolerance) then
      print*,"FAIL: exchange did not reproduce the linear field on interior sides"
      r = 1
    endif

    if(r == 0) print*,"ADAPTIVE MORTAR CHECKS PASSED"

    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call baseGeom%Free()
    call interp%Free()
    call forest%Free()
    call mesh%Free()
    call baseMesh%Free()

  endfunction adaptive_mortar_2d

  function integrateJ(geom,interp) result(area)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Geometry_2D
    implicit none
    type(SEMQuad),intent(in) :: geom
    type(Lagrange),intent(in) :: interp
    real(prec) :: area
    integer :: e,i,j

    area = 0.0_prec
    do e = 1,geom%nElem
      do j = 1,interp%N+1
        do i = 1,interp%N+1
          area = area+geom%J%interior(i,j,e,1)*interp%qWeights(i)*interp%qWeights(j)
        enddo
      enddo
    enddo

  endfunction integrateJ

endprogram test
