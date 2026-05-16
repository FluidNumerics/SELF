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

  exit_code = points_locate_2d_outside()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function points_locate_2d_outside() result(r)

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_Points

    implicit none

    integer,parameter :: controlDegree = 5
    integer,parameter :: targetDegree = 12
    integer,parameter :: nvar = 1
    integer,parameter :: nInside = 4
    integer,parameter :: nOutside = 4
    integer,parameter :: nPoints = nInside+nOutside
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    type(Points) :: pts
    character(LEN=255) :: WORKSPACE
    real(prec) :: xIn(nPoints,2),fSampled(nPoints,nvar)
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
    call f%SetEquation(1,'f = 1.0')
    call f%SetInteriorFromEquation(geometry,0.0_prec)

    ! Four points strictly inside an element (Block2D is 5x5 elements over
    ! [0,1]^2, so we avoid the element-corner lines at x,y in {0.2,0.4,0.6,0.8}).
    xIn(1,1) = 0.13_prec; xIn(1,2) = 0.07_prec
    xIn(2,1) = 0.43_prec; xIn(2,2) = 0.57_prec
    xIn(3,1) = 0.83_prec; xIn(3,2) = 0.17_prec
    xIn(4,1) = 0.55_prec; xIn(4,2) = 0.85_prec
    xIn(5,1) = -2.0_prec; xIn(5,2) = 0.50_prec
    xIn(6,1) = 3.0_prec; xIn(6,2) = 0.50_prec
    xIn(7,1) = 0.50_prec; xIn(7,2) = -1.5_prec
    xIn(8,1) = 0.50_prec; xIn(8,2) = 2.7_prec

    call pts%Init(nPoints,2)
    call pts%SetPoints(xIn)
    call pts%LocatePoints(geometry)

    ! Check assignment.
    r = 0
    do p = 1,nInside
      if(pts%elements(p) <= 0) then
        print*,"inside point ",p," was not located (got elements=",pts%elements(p),")"
        r = 1
      endif
    enddo
    do p = nInside+1,nPoints
      if(pts%elements(p) /= 0) then
        print*,"outside point ",p," was incorrectly located (elements=",pts%elements(p),")"
        r = 1
      endif
    enddo

    ! Sampling: outside points should return 0; inside points should give f=1.
    call pts%EvaluateScalar(f,fSampled)
    do p = 1,nInside
      if(abs(fSampled(p,1)-1.0_prec) > 1.0e-4_prec) then
        print*,"inside point ",p," sampled value ",fSampled(p,1)," != 1"
        r = 1
      endif
    enddo
    do p = nInside+1,nPoints
      if(abs(fSampled(p,1)) > 0.0_prec) then
        print*,"outside point ",p," sampled value ",fSampled(p,1)," != 0"
        r = 1
      endif
    enddo

    call pts%Free()
    call f%DissociateGeometry()
    call f%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

  endfunction points_locate_2d_outside
endprogram test
