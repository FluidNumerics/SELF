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

  exit_code = mappedvectordgdivergence_2d_mortar_mpi()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mappedvectordgdivergence_2d_mortar_mpi() result(r)
    !! Computes the weak-form (DG) divergence of f = (x,y) on the six-element
    !! DoubleMortarMesh (two mortar interfaces, one with a flip = 1 sub-edge). The surface-flux integrands on the mortar's small sides are
    !! built from side-averaged states (as in the DG models' BoundaryFlux), the big
    !! side's integrand is replaced by MortarFluxCollect, and the divergence must
    !! equal 2 to roundoff on all elements. The test also checks the discrete
    !! conservation identity on each mortar interface : the big side's surface
    !! integral must equal minus the sum of the small sides' surface integrals.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_MappedVector_2D
    use mpi

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedVector2D) :: f
    type(MappedScalar2D) :: df
    integer :: i,j,iel,e2,ivar,m,k
    integer :: eB,sB,eS,sS,offset
    integer :: bcids(1:4)
    real(prec) :: nhat(1:2),nmag,fx,fy
    real(prec) :: bigIntegral,smallIntegrals
    integer :: ierror

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
    call df%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    do ivar = 1,nvar
      call f%SetEquation(1,ivar,'f = x') ! x-component
      call f%SetEquation(2,ivar,'f = y') ! y-component
    enddo

    call f%SetInteriorFromEquation(geometry,0.0_prec)

    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%MortarExchange(mesh)
    call f%UpdateHost()

    ! Prolong the boundary attribute to the domain boundaries for external state
    do iel = 1,f%nElem
      do j = 1,4
        e2 = mesh%sideInfo(3,j,iel)
        if(e2 == 0 .and. mesh%sideInfo(1,j,iel) == 0) then
          do i = 1,f%interp%N+1
            f%extBoundary(i,j,iel,1:nvar,1:2) = f%boundary(i,j,iel,1:nvar,1:2)
          enddo
        endif
      enddo
    enddo

    ! Surface-flux integrand from side-averaged states, as in BoundaryFlux
    do ivar = 1,nvar
      do iEl = 1,f%nElem
        do j = 1,4
          do i = 1,f%interp%N+1
            nhat(1:2) = geometry%nHat%boundary(i,j,iEl,1,1:2)
            nmag = geometry%nScale%boundary(i,j,iEl,1)
            fx = 0.5_prec*(f%boundary(i,j,iEl,ivar,1)+f%extboundary(i,j,iEl,ivar,1))
            fy = 0.5_prec*(f%boundary(i,j,iEl,ivar,2)+f%extboundary(i,j,iEl,ivar,2))
            f%boundaryNormal(i,j,iEl,ivar) = (fx*nhat(1)+fy*nhat(2))*nmag
          enddo
        enddo
      enddo
    enddo

    call f%UpdateDevice()
    call f%MortarFluxCollect(mesh)
    call f%UpdateHost()

    ! Discrete conservation across the mortar interface : the big side's surface
    ! integral balances the small sides' surface integrals (opposite normals)
    offset = mesh%decomp%offsetElem(mesh%decomp%rankId+1)
    do m = 1,mesh%nMortars
      bigIntegral = 0.0_prec
      smallIntegrals = 0.0_prec
      eB = mesh%mortarInfo(1,m)
      sB = mesh%mortarInfo(2,m)
      if(mesh%decomp%elemToRank(eB) == mesh%decomp%rankId) then
        do i = 1,interp%N+1
          bigIntegral = bigIntegral+interp%qWeights(i)* &
                        f%boundaryNormal(i,sB,eB-offset,1)
        enddo
      endif
      do k = 1,2
        eS = mesh%mortarInfo(2*k+1,m)
        sS = mesh%mortarInfo(2*k+2,m)/10
        if(mesh%decomp%elemToRank(eS) == mesh%decomp%rankId) then
          do i = 1,interp%N+1
            smallIntegrals = smallIntegrals+interp%qWeights(i)* &
                             f%boundaryNormal(i,sS,eS-offset,1)
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
        return
      endif
    enddo

#ifdef ENABLE_GPU
    call f%MappedDGDivergence(df%interior_gpu)
#else
    call f%MappedDGDivergence(df%interior)
#endif
    call df%UpdateHost()

    ! Calculate diff from exact
    df%interior = abs(df%interior-2.0_prec)

    print*,"rank ",mesh%decomp%rankId," absmax error :",maxval(df%interior)
    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    call f%DissociateGeometry()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()
    call f%Free()
    call df%Free()

  endfunction mappedvectordgdivergence_2d_mortar_mpi
endprogram test
