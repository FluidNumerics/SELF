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

  exit_code = mappedscalarmortarexchange_3d_linear()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mappedscalarmortarexchange_3d_linear() result(r)
    !! Fills a MappedScalar3D with degree-N polynomial fields (distinct per variable)
    !! on the five-element SimpleMortarMesh and verifies that MortarExchange
    !! reproduces the exact trace in extBoundary on every face of the 2:1 mortar
    !! interface. Restriction of the big face's trace to the sub-faces is exact for
    !! polynomials of the control degree, and the L2 projection of the (globally
    !! polynomial) small-face traces recovers the big face's trace, so extBoundary
    !! must match boundary to roundoff on all mortar faces.
    !!
    !! The check runs on three meshes: the axis-aligned mesh (all mortar flips 0) and
    !! two meshes with rotated small elements realizing flips 1,4,6,3 and 2,5,7,0, so
    !! every one of the eight SELF face flips is exercised through the mortar staging,
    !! restriction scatter, and projection paths. The realized flips recorded in
    !! mortarInfo are asserted against the requested ones, which pins the constructor's
    !! rotation tables to the flip conventions of SideExchange.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D

    implicit none

    integer :: flipsA(1:4),flipsB(1:4)
    integer :: keepAliveBcids(1:6)
    type(Mesh3D),target :: keepAlive

    ! SELF tears down MPI when the number of live meshes drops to zero (see
    ! test/domaindecomposition_two_meshes.f90). This test builds three meshes in
    ! sequence, so an extra mesh is held open across all of them.
    keepAliveBcids(1:6) = SELF_BC_PRESCRIBED
    call keepAlive%StructuredMesh(1,1,1,1,1,1,0.5_prec,0.5_prec,0.5_prec, &
                                  keepAliveBcids)

    r = 0

    ! Axis-aligned : all mortar flips are 0
    r = r+CheckMortarExchange()

    ! Rotated small elements : East-face donors realizing flips 1,4,6,3
    flipsA = [1,4,6,3]
    r = r+CheckMortarExchange(flipsA)

    ! Rotated small elements : West-face donors realizing flips 2,5,7 (and 0)
    flipsB = [2,5,7,0]
    r = r+CheckMortarExchange(flipsB)

    call keepAlive%Free()

  endfunction mappedscalarmortarexchange_3d_linear

  integer function CheckMortarExchange(flips) result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D

    implicit none
    integer,intent(in),optional :: flips(1:4)

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 2
    real(prec),parameter :: dx = 0.1_prec
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: f
    integer :: i,j,iel,is,ivar,m,q
    integer :: eB,sB,eS,sS,offset,flip
    integer :: bcids(1:6)
    real(prec) :: maxErr

    r = 0

    ! The mesh is constructed before the interpolant : the mesh's domain
    ! decomposition assigns each MPI rank its GPU device, and the interpolant's
    ! device arrays must be allocated on that device.
    bcids(1:6) = SELF_BC_PRESCRIBED
    if(present(flips)) then
      call mesh%SimpleMortarMesh(dx,bcids,flips=flips)
    else
      call mesh%SimpleMortarMesh(dx,bcids)
    endif

    ! The realized mortar flips must match the requested ones; this pins the
    ! constructor's rotation tables to the flip conventions used by SideExchange
    ! and MortarFaceMap.
    if(present(flips)) then
      do q = 1,4
        sS = mesh%mortarInfo(2*q+2,1)/10
        flip = mesh%mortarInfo(2*q+2,1)-10*sS
        if(flip /= flips(q)) then
          print*,"FAIL: realized mortar flip does not match request : quadrant, ", &
            "flip, requested :",q,flip,flips(q)
          r = 1
        endif
      enddo
    endif

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    ! Rotated elements must still carry a positive Jacobian
    if(minval(geometry%J%interior) <= 0.0_prec) then
      print*,"FAIL: non-positive Jacobian on the mortar mesh"
      r = 1
    endif

    call f%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    ! Distinct full-degree fields per variable guard against variable-index
    ! cross-contamination in the exchange; both lie in the degree-N tensor space,
    ! so the mortar operators are exact on their traces.
    call f%SetEquation(1,'f = (1.0 + 0.8*x + 0.5*y + 0.3*z)^7')
    call f%SetEquation(2,'f = (1.0 - 0.6*x + 0.7*y - 0.4*z)^7')
    call f%SetInteriorFromEquation(geometry,0.0_prec)

    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%MortarExchange(mesh)
    call f%UpdateHost()

    ! On every mortar face of a globally polynomial field, the external state
    ! computed through the mortar operators must match this element's own boundary
    ! trace.
    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)
    maxErr = 0.0_prec
    do m = 1,mesh%nMortars

      eB = mesh%mortarInfo(1,m)
      sB = mesh%mortarInfo(2,m)
      if(mesh%decomp%elemToRank(eB) == mesh%decomp%rankId) then
        iel = eB-offset
        do ivar = 1,nvar
          do j = 1,interp%N+1
            do i = 1,interp%N+1
              maxErr = max(maxErr,abs(f%extBoundary(i,j,sB,iel,ivar)- &
                                      f%boundary(i,j,sB,iel,ivar)))
            enddo
          enddo
        enddo
      endif

      do q = 1,4
        eS = mesh%mortarInfo(2*q+1,m)
        sS = mesh%mortarInfo(2*q+2,m)/10
        if(mesh%decomp%elemToRank(eS) == mesh%decomp%rankId) then
          iel = eS-offset
          do ivar = 1,nvar
            do j = 1,interp%N+1
              do i = 1,interp%N+1
                maxErr = max(maxErr,abs(f%extBoundary(i,j,sS,iel,ivar)- &
                                        f%boundary(i,j,sS,iel,ivar)))
              enddo
            enddo
          enddo
        endif
      enddo

    enddo

    print*,"rank ",mesh%decomp%rankId," mortar exchange absmax error :",maxErr
    if(maxErr > tolerance) then
      r = 1
    endif

    ! Also verify that conforming machinery filled non-mortar interior faces (the
    ! faces joining the small elements to each other, which carry nonzero flips on
    ! the rotated meshes); any NaN or unfilled extBoundary on those faces would show
    ! up in downstream tests.
    do iel = 1,mesh%nElem
      do is = 1,6
        if(mesh%sideInfo(3,is,iel) > 0) then
          do ivar = 1,nvar
            do j = 1,interp%N+1
              do i = 1,interp%N+1
                if(abs(f%extBoundary(i,j,is,iel,ivar)- &
                       f%boundary(i,j,is,iel,ivar)) > tolerance) then
                  print*,"rank ",mesh%decomp%rankId, &
                    " conforming face mismatch at iel,is :",iel,is
                  r = 1
                endif
              enddo
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

  endfunction CheckMortarExchange
endprogram test
