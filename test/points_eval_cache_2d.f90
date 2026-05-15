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

  exit_code = points_eval_cache_2d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function points_eval_cache_2d() result(r)
    !! Verify that the cached and fallback (non-cached) code paths in
    !! EvaluateScalar produce equivalent results, and that both match the
    !! analytic reference.

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_Points

    implicit none

    integer,parameter :: controlDegree = 6
    integer,parameter :: targetDegree = 14
    integer,parameter :: nvar = 2
    integer,parameter :: nPoints = 12
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: agreementTol = 1.0e-12_prec
    real(prec),parameter :: analyticTol = 1.0e-4_prec
#else
    real(prec),parameter :: agreementTol = 1.0e-5_prec
    real(prec),parameter :: analyticTol = 1.0e-2_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    type(Points) :: pts
    character(LEN=255) :: WORKSPACE
    real(prec) :: xIn(nPoints,2)
    real(prec) :: vCached(nPoints,nvar),vFallback(nPoints,nvar)
    real(prec) :: fExact1,fExact2,errAgree,errAnalytic
    integer :: p

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5")

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)
    call f%SetEquation(1,'f = sin(3.14159265358979*x)*cos(3.14159265358979*y)')
    call f%SetEquation(2,'f = x*x + 2.0*y')
    call f%SetInteriorFromEquation(geometry,0.0_prec)

    ! Spread of interior points (avoid the 0.2-grid corner lines).
    do p = 1,nPoints
      xIn(p,1) = 0.07_prec+0.86_prec*real(modulo(7*p+1,nPoints),prec)/real(nPoints,prec)
      xIn(p,2) = 0.07_prec+0.86_prec*real(modulo(11*p+3,nPoints),prec)/real(nPoints,prec)
    enddo

    call pts%Init(nPoints,2)
    call pts%SetPoints(xIn)
    call pts%LocatePoints(geometry)

    ! Sanity: the cache must have been populated by LocatePoints.
    if(pts%nCached /= controlDegree) then
      print*,"cache not populated: nCached=",pts%nCached,"  expected=",controlDegree
      r = 1
      return
    endif
    if(.not. allocated(pts%lS_cache) .or. .not. allocated(pts%lT_cache)) then
      print*,"cache arrays not allocated after LocatePoints"
      r = 1
      return
    endif

    ! Path 1: cached.
    call pts%EvaluateScalar(f,vCached)

    ! Path 2: fallback. Invalidating nCached forces the on-the-fly path while
    ! leaving coordinates(:,:) and elements(:) intact.
    pts%nCached = 0
    call pts%EvaluateScalar(f,vFallback)

    ! Both paths must agree to within round-off (they evaluate the same math
    ! by two different orderings of operations).
    errAgree = 0.0_prec
    errAgree = max(errAgree,maxval(abs(vCached-vFallback)))
    print*,"max |cached - fallback| =",errAgree,"  tol=",agreementTol
    if(errAgree > agreementTol) then
      r = 1
      return
    endif

    ! Cross-check against analytic values for both variables.
    errAnalytic = 0.0_prec
    do p = 1,nPoints
      fExact1 = sin(3.14159265358979_prec*xIn(p,1))*cos(3.14159265358979_prec*xIn(p,2))
      fExact2 = xIn(p,1)*xIn(p,1)+2.0_prec*xIn(p,2)
      errAnalytic = max(errAnalytic, &
                        abs(vCached(p,1)-fExact1),abs(vFallback(p,1)-fExact1), &
                        abs(vCached(p,2)-fExact2),abs(vFallback(p,2)-fExact2))
    enddo
    print*,"max analytic error      =",errAnalytic,"  tol=",analyticTol
    if(errAnalytic > analyticTol) then
      r = 1
      return
    endif

    call pts%Free()
    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

    r = 0

  endfunction points_eval_cache_2d
endprogram test
