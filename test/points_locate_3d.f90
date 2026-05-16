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

  exit_code = points_locate_3d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function points_locate_3d() result(r)

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D
    use SELF_Points

    implicit none

    integer,parameter :: controlDegree = 5
    integer,parameter :: targetDegree = 12
    integer,parameter :: nvar = 1
    integer,parameter :: nPoints = 16
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: locateTol = 1.0e-6_prec
    real(prec),parameter :: evalTol = 1.0e-4_prec
#else
    real(prec),parameter :: locateTol = 1.0e-3_prec
    real(prec),parameter :: evalTol = 1.0e-2_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: f
    type(Points) :: pts
    character(LEN=255) :: WORKSPACE
    real(prec) :: xIn(nPoints,3),fExact,err
    real(prec) :: fSampled(nPoints,nvar)
    real(prec) :: lS(0:controlDegree),lT(0:controlDegree),lU(0:controlDegree)
    real(prec) :: xR,yR,zR
    integer :: p,iEl,i,j,k,found

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D_mesh.h5")

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)
    call f%SetEquation(1,'f = x + 2.0*y - 3.0*z')
    call f%SetInteriorFromEquation(geometry,0.0_prec)

    do p = 1,nPoints
      xIn(p,1) = 0.05_prec+0.9_prec*real(modulo(7*p+1,nPoints),prec)/real(nPoints,prec)
      xIn(p,2) = 0.05_prec+0.9_prec*real(modulo(11*p+3,nPoints),prec)/real(nPoints,prec)
      xIn(p,3) = 0.05_prec+0.9_prec*real(modulo(13*p+5,nPoints),prec)/real(nPoints,prec)
    enddo

    call pts%Init(nPoints,3)
    call pts%SetPoints(xIn)
    call pts%LocatePoints(geometry)

    found = 0
    do p = 1,nPoints
      if(pts%elements(p) > 0) found = found+1
    enddo
    print*,"located ",found," of ",nPoints," points"
    if(found /= nPoints) then
      r = 1
      return
    endif

    err = 0.0_prec
    do p = 1,nPoints
      iEl = pts%elements(p)
      lS = interp%CalculateLagrangePolynomials(pts%coordinates(p,1))
      lT = interp%CalculateLagrangePolynomials(pts%coordinates(p,2))
      lU = interp%CalculateLagrangePolynomials(pts%coordinates(p,3))
      xR = 0.0_prec
      yR = 0.0_prec
      zR = 0.0_prec
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            xR = xR+lS(i-1)*lT(j-1)*lU(k-1)*geometry%x%interior(i,j,k,iEl,1,1)
            yR = yR+lS(i-1)*lT(j-1)*lU(k-1)*geometry%x%interior(i,j,k,iEl,1,2)
            zR = zR+lS(i-1)*lT(j-1)*lU(k-1)*geometry%x%interior(i,j,k,iEl,1,3)
          enddo
        enddo
      enddo
      err = max(err,abs(xR-xIn(p,1)),abs(yR-xIn(p,2)),abs(zR-xIn(p,3)))
    enddo
    print*,"max inverse-map residual: ",err
    if(err > locateTol) then
      r = 1
      return
    endif

    call pts%EvaluateScalar(f,fSampled)
    err = 0.0_prec
    do p = 1,nPoints
      fExact = xIn(p,1)+2.0_prec*xIn(p,2)-3.0_prec*xIn(p,3)
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

  endfunction points_locate_3d
endprogram test
