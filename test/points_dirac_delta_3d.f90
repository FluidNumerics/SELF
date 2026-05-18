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

  exit_code = points_dirac_delta_3d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function points_dirac_delta_3d() result(r)
    !! 3D analogue of points_dirac_delta_2d. Verifies:
    !!   (1) Conservation: sum_{i,j,k} f * w_i*w_j*w_k*J_{ijk} = 1 for every
    !!       located point's own column on its containing element.
    !!   (2) Other elements zero.
    !!   (3) Node-coincident source reproduces 1/(w_{i*}*w_{j*}*w_{k*}*J_node).

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D
    use SELF_Points

    implicit none

    integer,parameter :: controlDegree = 4
    integer,parameter :: targetDegree = 8
    integer,parameter :: nPoints = 4
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: conservationTol = 1.0e-11_prec
    real(prec),parameter :: nodeTol = 1.0e-11_prec
#else
    real(prec),parameter :: conservationTol = 1.0e-4_prec
    real(prec),parameter :: nodeTol = 1.0e-4_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: f
    type(Points) :: pts
    character(LEN=255) :: WORKSPACE
    real(prec) :: xIn(nPoints,3)
    real(prec) :: integral,wi,wj,wk,jLoc,maxOther,expected,nodeErr
    real(prec) :: f_atNode
    integer :: p,i,j,k,iEl,iElP,iStar,jStar,kStar

    r = 0

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS_LOBATTO, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D_mesh.h5")

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nPoints,mesh%nelem)
    call f%AssociateGeometry(geometry)

    ! Points 1..nPoints-1 spread across the interior; point nPoints anchored
    ! exactly on an interior LGL node of element 1.
    do p = 1,nPoints-1
      xIn(p,1) = 0.10_prec+0.80_prec*real(modulo(3*p+1,nPoints),prec)/real(nPoints,prec)
      xIn(p,2) = 0.10_prec+0.80_prec*real(modulo(5*p+2,nPoints),prec)/real(nPoints,prec)
      xIn(p,3) = 0.10_prec+0.80_prec*real(modulo(7*p+3,nPoints),prec)/real(nPoints,prec)
    enddo
    iStar = controlDegree/2+1
    jStar = controlDegree/2+1
    kStar = controlDegree/2+1
    xIn(nPoints,1) = geometry%x%interior(iStar,jStar,kStar,1,1,1)
    xIn(nPoints,2) = geometry%x%interior(iStar,jStar,kStar,1,1,2)
    xIn(nPoints,3) = geometry%x%interior(iStar,jStar,kStar,1,1,3)

    call pts%Init(nPoints,3)
    call pts%SetPoints(xIn)
    call pts%LocatePoints(geometry)

    call pts%DiracDelta(geometry,f)

    do p = 1,nPoints
      iElP = pts%elements(p)
      if(iElP <= 0) then
        if(maxval(abs(f%interior(:,:,:,:,p))) /= 0.0_prec) then
          print*,"unlocated point p=",p,"has nonzero column; max=", &
            maxval(abs(f%interior(:,:,:,:,p)))
          r = 1
          return
        endif
        cycle
      endif

      maxOther = 0.0_prec
      do iEl = 1,mesh%nelem
        if(iEl == iElP) cycle
        maxOther = max(maxOther,maxval(abs(f%interior(:,:,:,iEl,p))))
      enddo
      if(maxOther /= 0.0_prec) then
        print*,"point p=",p,"leaked into non-owning element; max=",maxOther
        r = 1
        return
      endif

      integral = 0.0_prec
      do k = 1,interp%N+1
        wk = interp%qWeights(k)
        do j = 1,interp%N+1
          wj = interp%qWeights(j)
          do i = 1,interp%N+1
            wi = interp%qWeights(i)
            jLoc = geometry%J%interior(i,j,k,iElP,1)
            integral = integral+f%interior(i,j,k,iElP,p)*wi*wj*wk*jLoc
          enddo
        enddo
      enddo
      if(abs(integral-1.0_prec) > conservationTol) then
        print*,"conservation failure p=",p,"  integral=",integral, &
          "  err=",abs(integral-1.0_prec),"  tol=",conservationTol
        r = 1
        return
      endif
    enddo

    p = nPoints
    iElP = pts%elements(p)
    if(iElP <= 0) then
      print*,"node-coincident point not located"
      r = 1
      return
    endif
    f_atNode = f%interior(iStar,jStar,kStar,iElP,p)
    expected = 1.0_prec/(interp%qWeights(iStar)*interp%qWeights(jStar)* &
                         interp%qWeights(kStar)* &
                         geometry%J%interior(iStar,jStar,kStar,iElP,1))
    nodeErr = abs(f_atNode-expected)
    if(nodeErr > nodeTol) then
      print*,"node value mismatch: f=",f_atNode,"  expected=",expected, &
        "  err=",nodeErr,"  tol=",nodeTol
      r = 1
      return
    endif

    maxOther = 0.0_prec
    do k = 1,interp%N+1
      do j = 1,interp%N+1
        do i = 1,interp%N+1
          if(i == iStar .and. j == jStar .and. k == kStar) cycle
          maxOther = max(maxOther,abs(f%interior(i,j,k,iElP,p)))
        enddo
      enddo
    enddo
    if(maxOther > nodeTol) then
      print*,"node-coincident off-node mass too large=",maxOther,"  tol=",nodeTol
      r = 1
      return
    endif

    print*,"DiracDelta 3D: conservation, isolation, and node-coincident checks passed"

    call pts%Free()
    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

  endfunction points_dirac_delta_3d
endprogram test
