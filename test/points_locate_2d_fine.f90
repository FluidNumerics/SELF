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

! Regression test for https://github.com/FluidNumerics/SELF/issues/136
!
! LocatePoints / NewtonInverse_2D must resolve interior points on a refined mesh
! even in single precision. The Newton convergence test is on the reference-space
! step delta = J^{-1} * res, with J ~ h/2 for an element of physical size h, so the
! smallest achievable step is delta_floor ~ (2/h) * |x| * epsilon. This floor grows
! as h shrinks; with a fixed absolute tolerance of 1e-8 the iteration could never
! converge in real32 on fine meshes, and LocatePoints rejected interior points.
!
! Here we build a fine structured mesh (50 x 50 elements over [0,0.1]^2, so
! h = 2e-3) and require that every strictly-interior point is located and that the
! inverse map reconstructs the physical coordinate. In real32, delta_floor at
! |x| ~ 0.05 is ~6e-6 -- far above the old 1e-8 tolerance (which this test would
! have failed) but well below the precision-aware tolerance sqrt(epsilon) ~ 3.4e-4.
program test

  implicit none
  integer :: exit_code

  exit_code = points_locate_2d_fine()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function points_locate_2d_fine() result(r)

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_Points

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
    integer,parameter :: nPoints = 25
    ! 50 x 50 elements (two 25-element tiles per direction) of width 2e-3 over
    ! the domain [0, 0.1]^2, giving a small element size h = 2e-3.
    integer,parameter :: nElemPerTile = 25
    integer,parameter :: nTile = 2
    real(prec),parameter :: dxElem = 2.0e-3_prec
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: locateTol = 1.0e-6_prec
    real(prec),parameter :: evalTol = 1.0e-5_prec
#else
    real(prec),parameter :: locateTol = 1.0e-4_prec
    real(prec),parameter :: evalTol = 1.0e-3_prec
#endif
    real(prec),parameter :: lDomain = real(nElemPerTile*nTile,prec)*dxElem ! = 0.1
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    type(Points) :: pts
    real(prec) :: xIn(nPoints,2),fExact,xLoc,yLoc,err
    real(prec) :: fSampled(nPoints,nvar)
    real(prec) :: lS(0:controlDegree),lT(0:controlDegree)
    real(prec) :: xReconstructed,yReconstructed
    integer :: bcids(1:4)
    integer :: p,iEl,i,j,found

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Fine structured quad mesh: 50 x 50 elements over [0, 0.1]^2.
    bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                  SELF_BC_PRESCRIBED, & ! East
                  SELF_BC_PRESCRIBED, & ! North
                  SELF_BC_PRESCRIBED] ! West
    call mesh%StructuredMesh(nElemPerTile,nElemPerTile,nTile,nTile, &
                             dxElem,dxElem,bcids)

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)
    call f%SetEquation(1,'f = sin(3.14159265358979*x)*cos(3.14159265358979*y)')
    call f%SetInteriorFromEquation(geometry,0.0_prec)

    ! Pseudo-random sample points strictly inside the domain, kept away from the
    ! origin so |x| is order 0.05 and the single-precision convergence floor is
    ! firmly above the old 1e-8 tolerance.
    do p = 1,nPoints
      xIn(p,1) = lDomain*(0.1_prec+0.8_prec*real(modulo(7*p+1,nPoints),prec)/real(nPoints,prec))
      xIn(p,2) = lDomain*(0.1_prec+0.8_prec*real(modulo(11*p+3,nPoints),prec)/real(nPoints,prec))
    enddo

    call pts%Init(nPoints,2)
    call pts%SetPoints(xIn)
    call pts%LocatePoints(geometry)

    ! Every point must have been located (all are inside the mesh). Pre-fix, this
    ! count dropped well below nPoints in single precision.
    found = 0
    do p = 1,nPoints
      if(pts%elements(p) > 0) found = found+1
    enddo
    print*,"located ",found," of ",nPoints," points"
    if(found /= nPoints) then
      r = 1
      return
    endif

    ! Verify inverse map: reconstruct x from (iEl, s, t) and compare to xIn.
    err = 0.0_prec
    do p = 1,nPoints
      iEl = pts%elements(p)
      lS = interp%CalculateLagrangePolynomials(pts%coordinates(p,1))
      lT = interp%CalculateLagrangePolynomials(pts%coordinates(p,2))
      xReconstructed = 0.0_prec
      yReconstructed = 0.0_prec
      do j = 1,controlDegree+1
        do i = 1,controlDegree+1
          xReconstructed = xReconstructed+lS(i-1)*lT(j-1)*geometry%x%interior(i,j,iEl,1,1)
          yReconstructed = yReconstructed+lS(i-1)*lT(j-1)*geometry%x%interior(i,j,iEl,1,2)
        enddo
      enddo
      err = max(err,abs(xReconstructed-xIn(p,1)),abs(yReconstructed-xIn(p,2)))
    enddo
    print*,"max inverse-map residual: ",err
    if(err > locateTol) then
      r = 1
      return
    endif

    ! Verify field sampling against the analytic expression.
    call pts%EvaluateScalar(f,fSampled)
    err = 0.0_prec
    do p = 1,nPoints
      xLoc = xIn(p,1)
      yLoc = xIn(p,2)
      fExact = sin(3.14159265358979_prec*xLoc)*cos(3.14159265358979_prec*yLoc)
      err = max(err,abs(fSampled(p,1)-fExact))
    enddo
    print*,"max scalar-sample error: ",err,"  tol=",evalTol
    if(err > evalTol) then
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

  endfunction points_locate_2d_fine
endprogram test
