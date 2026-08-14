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

  exit_code = adaptive_mortar_3d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function adaptive_mortar_3d() result(r)
    !! End-to-end validation of the 3-D adaptive pipeline through mesh emission: adapt a
    !! forest, balance it, emit a nonconforming Mesh3D_t, and exercise the emitted
    !! connectivity with the existing solver machinery.
    !!
    !! A structured 2x2x2 base mesh is refined adaptively (one corner element to two
    !! levels) and 2:1-balanced, then EmitMesh produces the mesh (2:1 mortar interfaces
    !! at the level jumps). On this mesh:
    !!
    !!   1. The geometry built from the emitted node coordinates has strictly positive
    !!      Jacobians and integrates to the same total volume as the base mesh (the
    !!      refined mesh is the same domain).
    !!   2. For a globally linear field, MortarExchange + SideExchange reproduce the
    !!      field's trace in extBoundary on every interior face - both the conforming
    !!      faces and the 2:1 mortar faces - to roundoff. This is the same criterion as
    !!      the hand-built mortar-mesh tests and exercises the emitted sideInfo and
    !!      mortarInfo tables directly.
    !!   3. A uniformly refined forest emits a mesh with no mortars (the conforming
    !!      emission path).
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D
    use SELF_OctreeMesh_3D
    use SELF_AdaptiveMesh_3D

    implicit none

    integer,parameter :: controlDegree = 4
    integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-10_prec
#else
    real(prec),parameter :: tolerance = 1.0e-4_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: baseMesh,mesh
    type(SEMHex),target :: baseGeom,geometry
    type(OctreeMesh3D) :: forest
    type(MappedScalar3D) :: f
    integer :: bcids(1:6)
    integer :: i,j,k,iel,is
    real(prec) :: x,y,z,maxErr,baseVolume,volume
    integer,allocatable :: flag(:)

    r = 0
    bcids(1:6) = SELF_BC_PRESCRIBED

    call baseMesh%StructuredMesh(2,2,2,1,1,1,0.5_prec,0.5_prec,0.5_prec,bcids)

    ! Adapt: refine one corner element two levels, then 2:1-balance.
    call forest%Init(baseMesh)
    allocate(flag(1:forest%nLeaves))
    flag = OCTREE_KEEP
    flag(1) = OCTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    allocate(flag(1:forest%nLeaves))
    flag = OCTREE_KEEP
    flag(1:8) = OCTREE_REFINE
    call forest%AdaptFromFlags(flag)
    deallocate(flag)
    call forest%Balance2to1()

    call EmitMesh(forest,baseMesh,mesh)
    print*,"emitted nElem, nMortars :",mesh%nElem,mesh%nMortars
    if(mesh%nElem /= forest%nLeaves .or. mesh%nMortars <= 0) then
      print*,"FAIL: emitted mesh element/mortar count"
      r = 1
    endif

    ! Geometry on the base and emitted meshes.
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree, &
                     targetNodeType=UNIFORM)
    call baseGeom%Init(interp,baseMesh%nElem)
    call baseGeom%GenerateFromMesh(baseMesh)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    ! 1. Positive Jacobian and volume conservation.
    if(minval(geometry%J%interior) <= 0.0_prec) then
      print*,"FAIL: non-positive Jacobian on emitted mesh"
      r = 1
    endif
    baseVolume = integrateJ(baseGeom,interp)
    volume = integrateJ(geometry,interp)
    print*,"base volume, emitted volume :",baseVolume,volume
    if(abs(baseVolume-volume) > tolerance*abs(baseVolume)) then
      print*,"FAIL: emitted mesh does not preserve volume"
      r = 1
    endif

    ! 2. Mortar/conforming exchange of a globally linear field.
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
    call f%MortarExchange(mesh)
    call f%UpdateHost()

    ! On every interior face (conforming: sideInfo(3)>0, or mortar: sideInfo(1)>0) the
    ! external trace of a linear field must equal this element's own trace.
    maxErr = 0.0_prec
    do iel = 1,mesh%nElem
      do is = 1,6
        if(mesh%sideInfo(3,is,iel) > 0 .or. mesh%sideInfo(1,is,iel) > 0) then
          do j = 1,controlDegree+1
            do i = 1,controlDegree+1
              maxErr = max(maxErr,abs(f%extBoundary(i,j,is,iel,1)- &
                                      f%boundary(i,j,is,iel,1)))
            enddo
          enddo
        endif
      enddo
    enddo
    print*,"adaptive mortar/conforming exchange absmax error :",maxErr
    if(maxErr > tolerance) then
      print*,"FAIL: exchange did not reproduce the linear field on interior faces"
      r = 1
    endif

    ! 3. A uniformly refined forest is conforming, so EmitMesh must produce a mesh with
    ! no mortars (exercises the no-mortar emission path).
    block
      type(OctreeMesh3D) :: cForest
      type(Mesh3D) :: cMesh
      integer,allocatable :: cflag(:)
      call cForest%Init(baseMesh)
      allocate(cflag(1:cForest%nLeaves))
      cflag = OCTREE_REFINE
      call cForest%AdaptFromFlags(cflag)
      deallocate(cflag)
      call EmitMesh(cForest,baseMesh,cMesh)
      if(cMesh%nMortars /= 0) then
        print*,"FAIL: conforming (uniform) emit produced mortars :",cMesh%nMortars
        r = 1
      endif
      call cForest%Free()
      call cMesh%Free()
    endblock

    if(r == 0) print*,"ADAPTIVE MORTAR 3D CHECKS PASSED"

    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call baseGeom%Free()
    call interp%Free()
    call forest%Free()
    call mesh%Free()
    call baseMesh%Free()

  endfunction adaptive_mortar_3d

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
