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

  exit_code = mappedscalarmortarexchange_2d_linear()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mappedscalarmortarexchange_2d_linear() result(r)
    !! Fills a MappedScalar2D with degree-N polynomial fields (distinct per variable)
    !! on the six-element DoubleMortarMesh and verifies that MortarExchange reproduces
    !! the exact trace in extBoundary on every side of both 2:1 mortar interfaces,
    !! including the interface whose second sub-edge has flip = 1. Restriction of the
    !! big side's trace to the sub-edges is exact for polynomials of the control
    !! degree, and the L2 projection of the (globally polynomial) small-side traces
    !! recovers the big side's trace, so extBoundary must match boundary to roundoff
    !! on all mortar sides.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 2
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    integer :: i,iel,is,ivar,m,k
    integer :: eB,sB,eS,sS,offset
    integer :: bcids(1:4)
    real(prec) :: maxErr

    ! The mesh is constructed before the interpolant : the mesh's domain
    ! decomposition assigns each MPI rank its GPU device, and the interpolant's
    ! device arrays must be allocated on that device.
    bcids(1:4) = [SELF_BC_PRESCRIBED, & ! South
                  SELF_BC_PRESCRIBED, & ! East
                  SELF_BC_PRESCRIBED, & ! North
                  SELF_BC_PRESCRIBED] ! West
    call mesh%DoubleMortarMesh(0.1_prec,bcids)

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    ! Distinct full-degree fields per variable guard against variable-index
    ! cross-contamination in the exchange; both lie in the degree-N tensor space,
    ! so the mortar operators are exact on their traces.
    call f%SetEquation(1,'f = (1.0 + 0.8*x + 0.5*y)^7')
    call f%SetEquation(2,'f = (1.0 - 0.6*x + 0.7*y)^7')
    call f%SetInteriorFromEquation(geometry,0.0_prec)

    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%MortarExchange(mesh)
    call f%UpdateHost()

    ! On every mortar side of a globally linear field, the external state computed
    ! through the mortar operators must match this element's own boundary trace.
    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)
    maxErr = 0.0_prec
    do m = 1,mesh%nMortars

      eB = mesh%mortarInfo(1,m)
      sB = mesh%mortarInfo(2,m)
      if(mesh%decomp%elemToRank(eB) == mesh%decomp%rankId) then
        iel = eB-offset
        do ivar = 1,nvar
          do i = 1,interp%N+1
            maxErr = max(maxErr,abs(f%extBoundary(i,sB,iel,ivar)- &
                                    f%boundary(i,sB,iel,ivar)))
          enddo
        enddo
      endif

      do k = 1,2
        eS = mesh%mortarInfo(2*k+1,m)
        sS = mesh%mortarInfo(2*k+2,m)/10
        if(mesh%decomp%elemToRank(eS) == mesh%decomp%rankId) then
          iel = eS-offset
          do ivar = 1,nvar
            do i = 1,interp%N+1
              maxErr = max(maxErr,abs(f%extBoundary(i,sS,iel,ivar)- &
                                      f%boundary(i,sS,iel,ivar)))
            enddo
          enddo
        endif
      enddo

    enddo

    print*,"rank ",mesh%decomp%rankId," mortar exchange absmax error :",maxErr
    if(maxErr <= tolerance) then
      r = 0
    else
      r = 1
    endif

    ! Also verify that conforming machinery left non-mortar interior sides intact
    ! (side 8 of the mesh joins the two small elements conformally); any NaN or
    ! unfilled extBoundary on those sides would show up in downstream tests.
    do iel = 1,mesh%nElem
      do is = 1,4
        if(mesh%sideInfo(3,is,iel) > 0) then
          do ivar = 1,nvar
            do i = 1,interp%N+1
              if(abs(f%extBoundary(i,is,iel,ivar)- &
                     f%boundary(i,is,iel,ivar)) > tolerance) then
                print*,"rank ",mesh%decomp%rankId, &
                  " conforming side mismatch at iel,is :",iel,is
                r = 1
              endif
            enddo
          enddo
        endif
      enddo
    enddo

    call f%DissociateGeometry()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()
    call f%Free()

  endfunction mappedscalarmortarexchange_2d_linear
endprogram test
