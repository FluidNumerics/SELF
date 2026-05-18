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

  exit_code = points_dirac_delta_2d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function points_dirac_delta_2d() result(r)
    !! Verifies the 2D DiracDelta scatter:
    !!   (1) Conservation. For every located point p with elements(p) > 0,
    !!       <f_p, 1>_h = sum_{i,j} f_ij^{iEl(p),p} * w_i*w_j*J_{ij}^{iEl(p)}
    !!       must equal 1 (post-mass form encodes a unit-strength delta whose
    !!       discrete L2 inner product with the unit function recovers S=1).
    !!   (2) Other-element columns are zero. For iEl' /= elements(p), every
    !!       entry of scalar%interior(:,:,iEl',p) must be exactly zero.
    !!   (3) Node-coincident source. When a point lands exactly on the
    !!       interior LGL node (i*,j*) of element iEl*, the only non-zero
    !!       entry is scalar%interior(i*,j*,iEl*,p) = 1/(w_{i*}*w_{j*}*J_{i*,j*}).
    !!
    !! Uses the Block2D fixture (affine [0,1]^2 mesh) so J is constant inside
    !! every element and conservation is exact up to floating-point round-off.

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_Points

    implicit none

    integer,parameter :: controlDegree = 6
    integer,parameter :: targetDegree = 14
    integer,parameter :: nPoints = 5
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: conservationTol = 1.0e-11_prec
    real(prec),parameter :: nodeTol = 1.0e-11_prec
#else
    real(prec),parameter :: conservationTol = 1.0e-4_prec
    real(prec),parameter :: nodeTol = 1.0e-4_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    type(Points) :: pts
    character(LEN=255) :: WORKSPACE
    real(prec) :: xIn(nPoints,2)
    real(prec) :: integral,wi,wj,jLoc,maxOther,expected,nodeErr
    real(prec) :: f_atNode
    integer :: p,i,j,iEl,iElP,iStar,jStar,k

    r = 0

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS_LOBATTO, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5")

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    ! One variable per point.
    call f%Init(interp,nPoints,mesh%nelem)
    call f%AssociateGeometry(geometry)

    ! Points (1..nPoints-1): spread across the interior of [0,1]^2.
    ! Point nPoints: anchored exactly on an interior LGL node of element 1, so
    !                we can check the node-coincident expected value.
    do p = 1,nPoints-1
      xIn(p,1) = 0.07_prec+0.86_prec*real(modulo(7*p+1,nPoints),prec)/real(nPoints,prec)
      xIn(p,2) = 0.07_prec+0.86_prec*real(modulo(11*p+3,nPoints),prec)/real(nPoints,prec)
    enddo
    iStar = controlDegree/2+1  ! a strictly interior node
    jStar = controlDegree/2+1
    xIn(nPoints,1) = geometry%x%interior(iStar,jStar,1,1,1)
    xIn(nPoints,2) = geometry%x%interior(iStar,jStar,1,1,2)

    call pts%Init(nPoints,2)
    call pts%SetPoints(xIn)
    call pts%LocatePoints(geometry)

    call pts%DiracDelta(geometry,f)

    ! --- (1) + (2) Conservation and other-element zeros ----------------------
    do p = 1,nPoints
      iElP = pts%elements(p)
      if(iElP <= 0) then
        ! Point not located on this rank: column p must be all zero.
        if(maxval(abs(f%interior(:,:,:,p))) /= 0.0_prec) then
          print*,"unlocated point p=",p,"has nonzero column; max=", &
            maxval(abs(f%interior(:,:,:,p)))
          r = 1
          return
        endif
        cycle
      endif

      ! Other elements: must be exactly zero (DiracDelta zeroed the field
      ! before scattering into iElP only).
      maxOther = 0.0_prec
      do iEl = 1,mesh%nelem
        if(iEl == iElP) cycle
        maxOther = max(maxOther,maxval(abs(f%interior(:,:,iEl,p))))
      enddo
      if(maxOther /= 0.0_prec) then
        print*,"point p=",p,"leaked into non-owning element; max=",maxOther
        r = 1
        return
      endif

      ! Discrete inner product against 1: sum_{i,j} f * w_i*w_j*J_{ij}.
      integral = 0.0_prec
      do j = 1,interp%N+1
        wj = interp%qWeights(j)
        do i = 1,interp%N+1
          wi = interp%qWeights(i)
          jLoc = geometry%J%interior(i,j,iElP,1)
          integral = integral+f%interior(i,j,iElP,p)*wi*wj*jLoc
        enddo
      enddo
      if(abs(integral-1.0_prec) > conservationTol) then
        print*,"conservation failure p=",p,"  integral=",integral, &
          "  err=",abs(integral-1.0_prec),"  tol=",conservationTol
        r = 1
        return
      endif
    enddo

    ! --- (3) Node-coincident: only one nonzero, equal to 1/(w_i w_j J) ------
    p = nPoints
    iElP = pts%elements(p)
    if(iElP <= 0) then
      print*,"node-coincident point not located"
      r = 1
      return
    endif
    f_atNode = f%interior(iStar,jStar,iElP,p)
    expected = 1.0_prec/(interp%qWeights(iStar)*interp%qWeights(jStar)* &
                         geometry%J%interior(iStar,jStar,iElP,1))
    nodeErr = abs(f_atNode-expected)
    if(nodeErr > nodeTol) then
      print*,"node value mismatch: f=",f_atNode,"  expected=",expected, &
        "  err=",nodeErr,"  tol=",nodeTol
      r = 1
      return
    endif

    ! All other (i,j) entries in element iElP must be ~0.
    maxOther = 0.0_prec
    do j = 1,interp%N+1
      do i = 1,interp%N+1
        if(i == iStar .and. j == jStar) cycle
        maxOther = max(maxOther,abs(f%interior(i,j,iElP,p)))
      enddo
    enddo
    if(maxOther > nodeTol) then
      print*,"node-coincident off-node mass too large=",maxOther,"  tol=",nodeTol
      r = 1
      return
    endif

    print*,"DiracDelta 2D: conservation, isolation, and node-coincident checks passed"

    call pts%Free()
    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

    ! Silence the unused-variable warning for k under some compilers.
    k = 0

  endfunction points_dirac_delta_2d
endprogram test
