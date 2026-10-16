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
!! Spectral accuracy of the DG (Bassi-Rebay) gradient on per-axis periodic structured
!! meshes built with Mesh2D%PeriodicStructuredMesh.
!!
!! The scalar is f = sin(2 pi x/Lx + 0.3) cos(2 pi y/Ly), which is periodic in both
!! directions, on [0,Lx] x [0,Ly] = [0,1] x [0,0.8]. Two meshes are tested: doubly
!! periodic, which has no domain boundary at all, and a channel that is periodic in x
!! with walls in y, where the exterior state is set to the interior trace.
!!
!! Across a periodic wrap the exterior state comes only from SideExchange. The test
!! asserts, for every interior side including the wrapped ones, that the neighbor
!! trace (extBoundary) agrees with the local trace to the interpolation error of f,
!! which a side wired to the wrong element (or left unset) cannot. It then asserts
!! that the gradient converges spectrally: the maximum nodal error at N = 10 is below
!! `tolerance` and at least 100 times smaller than at N = 6.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Geometry_2D
  use SELF_MappedScalar_2D
  use SELF_MappedVector_2D
  use mpi

  implicit none
  integer :: exit_code

  exit_code = mappedscalarbrgradient_2d_periodic()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  subroutine GradientError(mesh,N,traceError,gradError)
    !! Builds the geometry on `mesh` with a degree N interpolant, takes the DG gradient
    !! of f, and returns the largest neighbor-trace mismatch on interior sides and the
    !! largest nodal error of the gradient, both reduced over all ranks.
    implicit none
    type(Mesh2D),intent(inout) :: mesh
    integer,intent(in) :: N
    real(prec),intent(out) :: traceError
    real(prec),intent(out) :: gradError
    ! Local
    real(prec),parameter :: twopi = 2.0_prec*pi
    real(prec),parameter :: Lx = 1.0_prec,Ly = 0.8_prec
    type(Lagrange),target :: interp
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    type(MappedVector2D) :: df
    integer :: iel,iside,i,j,e2,ierror
    real(prec) :: x,y,dfdx,dfdy,localErr(1:2),globalErr(1:2)

    call interp%Init(N=N, &
                     controlNodeType=GAUSS, &
                     M=N+2, &
                     targetNodeType=UNIFORM)

    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,1,mesh%nElem)
    call df%Init(interp,1,mesh%nElem)
    call f%AssociateGeometry(geometry)

    call f%SetEquation(1,'f = sin(6.283185307179586*x + 0.3)*cos(7.853981633974483*y)')
    call f%SetInteriorFromEquation(geometry,0.0_prec)
    call f%BoundaryInterp()
    call f%SideExchange(mesh)
    call f%UpdateHost()

    localErr = 0.0_prec
    do iel = 1,mesh%nElem
      do iside = 1,4
        e2 = mesh%sideInfo(3,iside,iel)
        if(e2 == 0) then
          ! Wall of the channel: prolong the interior trace
          do i = 1,N+1
            f%extBoundary(i,iside,iel,1) = f%boundary(i,iside,iel,1)
          enddo
        else
          do i = 1,N+1
            localErr(1) = max(localErr(1),abs(f%extBoundary(i,iside,iel,1)-f%boundary(i,iside,iel,1)))
          enddo
        endif
      enddo
    enddo

    call f%UpdateDevice()
    call f%AverageSides()
#ifdef ENABLE_GPU
    call f%MappedDGGradient(df%interior_gpu)
#else
    call f%MappedDGGradient(df%interior)
#endif
    call df%UpdateHost()

    do iel = 1,mesh%nElem
      do j = 1,N+1
        do i = 1,N+1
          x = geometry%x%interior(i,j,iel,1,1)
          y = geometry%x%interior(i,j,iel,1,2)
          dfdx = (twopi/Lx)*cos(twopi*x/Lx+0.3_prec)*cos(twopi*y/Ly)
          dfdy = -(twopi/Ly)*sin(twopi*x/Lx+0.3_prec)*sin(twopi*y/Ly)
          localErr(2) = max(localErr(2),abs(df%interior(i,j,iel,1,1)-dfdx))
          localErr(2) = max(localErr(2),abs(df%interior(i,j,iel,1,2)-dfdy))
        enddo
      enddo
    enddo

    call MPI_Allreduce(localErr,globalErr,2,mesh%decomp%mpiPrec,MPI_MAX,mesh%decomp%mpiComm,ierror)
    traceError = globalErr(1)
    gradError = globalErr(2)

    call f%DissociateGeometry()
    call f%Free()
    call df%Free()
    call geometry%Free()
    call interp%Free()

  endsubroutine GradientError

  integer function CheckMesh(periodicY) result(r)
    implicit none
    logical,intent(in) :: periodicY
    ! Local
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
    real(prec),parameter :: traceTolerance = 10.0_prec**(-6)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-2)
    real(prec),parameter :: traceTolerance = 10.0_prec**(-2)
#endif
    type(Mesh2D),target :: mesh
    real(prec) :: traceErr6,gradErr6,traceErr10,gradErr10
    integer :: bcids(1:4)

    bcids = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    ! 8 x 8 elements of 0.125 x 0.1 on [0,1] x [0,0.8], four tiles
    call mesh%PeriodicStructuredMesh(4,4,2,2,0.125_prec,0.1_prec,.true.,periodicY,bcids)

    call GradientError(mesh,6,traceErr6,gradErr6)
    call GradientError(mesh,10,traceErr10,gradErr10)

    if(mesh%decomp%rankId == 0) then
      print*,"periodicY = ",periodicY
      print*,"  N = 6 : max trace mismatch ",traceErr6,", max gradient error ",gradErr6
      print*,"  N = 10: max trace mismatch ",traceErr10,", max gradient error ",gradErr10
    endif

    r = 0
    if(traceErr6 > traceTolerance .or. traceErr10 > traceTolerance) then
      if(mesh%decomp%rankId == 0) print*,"FAIL: neighbor trace mismatch exceeds ",traceTolerance
      r = 1
    endif
    if(gradErr10 > tolerance) then
      if(mesh%decomp%rankId == 0) print*,"FAIL: gradient error at N = 10 exceeds ",tolerance
      r = 1
    endif
#ifdef DOUBLE_PRECISION
    if(gradErr10 > 0.01_prec*gradErr6) then
      if(mesh%decomp%rankId == 0) print*,"FAIL: gradient error did not drop by 100 from N = 6 to N = 10"
      r = 1
    endif
#endif

    call mesh%Free()

  endfunction CheckMesh

  integer function mappedscalarbrgradient_2d_periodic() result(r)
    implicit none
    ! Local
    type(Mesh2D) :: keepalive

    ! Keeps one decomposition alive so that freeing a mesh under test does not
    ! finalize MPI before the next one is built.
    call keepalive%StructuredMesh(2,1,1,1,1.0_prec,1.0_prec,[1,1,1,1])

    r = 0
    r = max(r,CheckMesh(.true.)) ! doubly periodic
    r = max(r,CheckMesh(.false.)) ! channel: periodic in x, walls in y

    call keepalive%Free()

  endfunction mappedscalarbrgradient_2d_periodic

endprogram test
