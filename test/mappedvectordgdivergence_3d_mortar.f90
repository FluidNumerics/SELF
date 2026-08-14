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

  exit_code = mappedvectordgdivergence_3d_mortar()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mappedvectordgdivergence_3d_mortar() result(r)
    !! Computes the weak-form (DG) divergence of f = (x,y,z) on the five-element
    !! SimpleMortarMesh. The surface-flux integrands on the mortar's small faces are
    !! built from face-averaged states (as in the DG models' BoundaryFlux), the big
    !! face's integrand is replaced by MortarFluxCollect, and the divergence must
    !! equal 3 to roundoff on all elements. The test also checks the discrete
    !! conservation identity on the mortar interface : the big face's surface
    !! integral must equal minus the sum of the small faces' surface integrals.
    !!
    !! The check runs on the axis-aligned mesh (all mortar flips 0) and on two meshes
    !! with rotated small elements realizing flips 1,4,6,3 and 2,5,7,0, so every SELF
    !! face flip is exercised through the flux staging and projection paths.
    use SELF_Constants
    use SELF_Mesh_3D

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

    r = r+CheckMortarDivergence()

    flipsA = [1,4,6,3]
    r = r+CheckMortarDivergence(flipsA)

    flipsB = [2,5,7,0]
    r = r+CheckMortarDivergence(flipsB)

    call keepAlive%Free()

  endfunction mappedvectordgdivergence_3d_mortar

  integer function CheckMortarDivergence(flips) result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D
    use SELF_MappedVector_3D
    use mpi

    implicit none
    integer,intent(in),optional :: flips(1:4)

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
    real(prec),parameter :: dx = 0.1_prec
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedVector3D) :: f
    type(MappedScalar3D) :: df
    integer :: i,j,is,iel,e2,ivar,m,q
    integer :: eB,sB,eS,sS,offset
    integer :: bcids(1:6)
    real(prec) :: nhat(1:3),nmag,fx,fy,fz
    real(prec) :: bigIntegral,smallIntegrals
    integer :: ierror

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

    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call df%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    do ivar = 1,nvar
      call f%SetEquation(1,ivar,'f = x') ! x-component
      call f%SetEquation(2,ivar,'f = y') ! y-component
      call f%SetEquation(3,ivar,'f = z') ! z-component
    enddo

    call f%SetInteriorFromEquation(geometry,0.0_prec)

    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%MortarExchange(mesh)
    call f%UpdateHost()

    ! Prolong the boundary attribute to the domain boundaries for external state
    do iel = 1,f%nElem
      do is = 1,6
        e2 = mesh%sideInfo(3,is,iel)
        if(e2 == 0 .and. mesh%sideInfo(1,is,iel) == 0) then
          do j = 1,f%interp%N+1
            do i = 1,f%interp%N+1
              f%extBoundary(i,j,is,iel,1:nvar,1:3) = f%boundary(i,j,is,iel,1:nvar,1:3)
            enddo
          enddo
        endif
      enddo
    enddo

    ! Surface-flux integrand from face-averaged states, as in BoundaryFlux
    do ivar = 1,nvar
      do iEl = 1,f%nElem
        do is = 1,6
          do j = 1,f%interp%N+1
            do i = 1,f%interp%N+1
              nhat(1:3) = geometry%nHat%boundary(i,j,is,iEl,1,1:3)
              nmag = geometry%nScale%boundary(i,j,is,iEl,1)
              fx = 0.5_prec*(f%boundary(i,j,is,iEl,ivar,1)+ &
                             f%extboundary(i,j,is,iEl,ivar,1))
              fy = 0.5_prec*(f%boundary(i,j,is,iEl,ivar,2)+ &
                             f%extboundary(i,j,is,iEl,ivar,2))
              fz = 0.5_prec*(f%boundary(i,j,is,iEl,ivar,3)+ &
                             f%extboundary(i,j,is,iEl,ivar,3))
              f%boundaryNormal(i,j,is,iEl,ivar) = (fx*nhat(1)+fy*nhat(2)+ &
                                                   fz*nhat(3))*nmag
            enddo
          enddo
        enddo
      enddo
    enddo

    call f%UpdateDevice()
    call f%MortarFluxCollect(mesh)
    call f%UpdateHost()

    ! Discrete conservation across the mortar interface : the big face's surface
    ! integral balances the small faces' surface integrals (opposite normals)
    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)
    do m = 1,mesh%nMortars
      bigIntegral = 0.0_prec
      smallIntegrals = 0.0_prec
      eB = mesh%mortarInfo(1,m)
      sB = mesh%mortarInfo(2,m)
      if(mesh%decomp%elemToRank(eB) == mesh%decomp%rankId) then
        do j = 1,interp%N+1
          do i = 1,interp%N+1
            bigIntegral = bigIntegral+interp%qWeights(i)*interp%qWeights(j)* &
                          f%boundaryNormal(i,j,sB,eB-offset,1)
          enddo
        enddo
      endif
      do q = 1,4
        eS = mesh%mortarInfo(2*q+1,m)
        sS = mesh%mortarInfo(2*q+2,m)/10
        if(mesh%decomp%elemToRank(eS) == mesh%decomp%rankId) then
          do j = 1,interp%N+1
            do i = 1,interp%N+1
              smallIntegrals = smallIntegrals+interp%qWeights(i)*interp%qWeights(j)* &
                               f%boundaryNormal(i,j,sS,eS-offset,1)
            enddo
          enddo
        endif
      enddo
      if(mesh%decomp%mpiEnabled) then
        call mpi_allreduce(MPI_IN_PLACE,bigIntegral,1,mesh%decomp%mpiPrec, &
                           MPI_SUM,mesh%decomp%mpiComm,ierror)
        call mpi_allreduce(MPI_IN_PLACE,smallIntegrals,1,mesh%decomp%mpiPrec, &
                           MPI_SUM,mesh%decomp%mpiComm,ierror)
      endif
      print*,"rank ",mesh%decomp%rankId," mortar conservation defect :", &
        abs(bigIntegral+smallIntegrals)
      if(abs(bigIntegral+smallIntegrals) > tolerance) then
        print*,"mortar interface is not conservative"
        r = 1
      endif
    enddo

#ifdef ENABLE_GPU
    call f%MappedDGDivergence(df%interior_gpu)
#else
    call f%MappedDGDivergence(df%interior)
#endif
    call df%UpdateHost()

    ! Calculate diff from exact
    df%interior = abs(df%interior-3.0_prec)

    print*,"rank ",mesh%decomp%rankId," absmax error :",maxval(df%interior)
    if(maxval(df%interior) > tolerance) then
      r = 1
    endif

    call f%DissociateGeometry()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()
    call f%Free()
    call df%Free()

  endfunction CheckMortarDivergence
endprogram test
