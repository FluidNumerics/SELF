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

  exit_code = mesh2d_uniform_refine()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mesh2d_uniform_refine() result(r)
    !! Validates uniform h-refinement (AMR Stage 2) of a 2-D structured mesh:
    !!
    !!   1. Element count quadruples (each element -> 4 children).
    !!   2. The refined mesh represents the identical domain: the area integral (sum of J times
    !!      the quadrature weights over the geometry built from the refined mesh) equals the base
    !!      mesh's area to roundoff, and every Jacobian is strictly positive (no inverted
    !!      children).
    !!   3. The refined connectivity is self-consistent: for every interior side the neighbor's
    !!      neighbor points back with the matching side and flip.
    !!   4. The number of physical-boundary sides doubles (each boundary edge splits in two).
    !!
    !! Together these confirm the isoparametric child geometry and the inherited connectivity are
    !! correct without relying on element ordering.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MeshRefinement_2D

    implicit none

    integer,parameter :: controlDegree = 4
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-11_prec
#else
    real(prec),parameter :: tolerance = 1.0e-4_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: baseMesh
    type(Mesh2D),target :: fineMesh
    type(SEMQuad),target :: baseGeom
    type(SEMQuad),target :: fineGeom
    integer :: bcids(1:4)
    integer :: e,s,n,ns,fineBoundary,baseBoundary
    real(prec) :: baseArea,fineArea

    r = 0

    bcids(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    ! 3 x 2 element structured base mesh on a 1.5 x 1.0 domain.
    call baseMesh%StructuredMesh(3,2,1,1,0.5_prec,0.5_prec,bcids)

    call UniformRefineMesh(baseMesh,fineMesh)

    ! ---- 1. Element count quadruples ----
    if(fineMesh%nElem /= 4*baseMesh%nElem) then
      print*,"FAIL: element count",fineMesh%nElem,"expected",4*baseMesh%nElem
      r = 1
    endif

    ! ---- Build geometry on both meshes ----
    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree,targetNodeType=UNIFORM)
    call baseGeom%Init(interp,baseMesh%nElem)
    call baseGeom%GenerateFromMesh(baseMesh)
    call fineGeom%Init(interp,fineMesh%nElem)
    call fineGeom%GenerateFromMesh(fineMesh)

    ! ---- 2a. Area conservation ----
    baseArea = totalArea(baseGeom,interp)
    fineArea = totalArea(fineGeom,interp)
    print*,"base area, fine area :",baseArea,fineArea
    if(abs(baseArea-fineArea) > tolerance*abs(baseArea)) then
      print*,"FAIL: area not conserved under refinement"
      r = 1
    endif

    ! ---- 2b. All child Jacobians strictly positive ----
    if(minval(fineGeom%J%interior) <= 0.0_prec) then
      print*,"FAIL: non-positive Jacobian in refined mesh",minval(fineGeom%J%interior)
      r = 1
    endif

    ! ---- 3. Connectivity reciprocity ----
    do e = 1,fineMesh%nElem
      do s = 1,4
        n = fineMesh%sideInfo(3,s,e)
        if(n > 0) then
          ns = fineMesh%sideInfo(4,s,e)/10
          if(fineMesh%sideInfo(3,ns,n) /= e) then
            print*,"FAIL: neighbor not reciprocal at elem,side",e,s
            r = 1
          endif
          if(fineMesh%sideInfo(4,ns,n)/10 /= s) then
            print*,"FAIL: neighbor side not reciprocal at elem,side",e,s
            r = 1
          endif
          if(mod(fineMesh%sideInfo(4,ns,n),10) /= mod(fineMesh%sideInfo(4,s,e),10)) then
            print*,"FAIL: flip not reciprocal at elem,side",e,s
            r = 1
          endif
        endif
      enddo
    enddo

    ! ---- 4. Boundary side count doubles ----
    baseBoundary = 0
    do e = 1,baseMesh%nElem
      do s = 1,4
        if(baseMesh%sideInfo(3,s,e) == 0) baseBoundary = baseBoundary+1
      enddo
    enddo
    fineBoundary = 0
    do e = 1,fineMesh%nElem
      do s = 1,4
        if(fineMesh%sideInfo(3,s,e) == 0) fineBoundary = fineBoundary+1
      enddo
    enddo
    if(fineBoundary /= 2*baseBoundary) then
      print*,"FAIL: boundary side count",fineBoundary,"expected",2*baseBoundary
      r = 1
    endif

    if(r == 0) print*,"UNIFORM REFINE CHECKS PASSED"

    call baseGeom%Free()
    call fineGeom%Free()
    call interp%Free()
    call fineMesh%Free()
    call baseMesh%Free()

  endfunction mesh2d_uniform_refine

  function totalArea(geom,interp) result(area)
    !! Physical area of the domain = sum over elements of the reference-square integral of the
    !! Jacobian, evaluated with the interpolant's quadrature weights.
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

  endfunction totalArea

endprogram test
