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

  exit_code = mappedscalarbrgradient_3d_linear()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mappedscalarbrgradient_3d_linear() result(r)

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MappedScalar_3D
    use SELF_MappedVector_3D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-2)
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: f
    type(MappedVector3D) :: df
    integer :: iel
    integer :: iside
    integer :: i
    integer :: j
    integer :: k
    integer :: e2,s2,bcid
    integer :: bcids(1:6)

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Create a uniform block mesh
    bcids(1:6) = [SELF_BC_PRESCRIBED, & ! Bottom
                  SELF_BC_PRESCRIBED, & ! South
                  SELF_BC_PRESCRIBED, & ! East
                  SELF_BC_PRESCRIBED, & ! North
                  SELF_BC_PRESCRIBED, & ! West
                  SELF_BC_PRESCRIBED] ! Top

    call mesh%StructuredMesh(2,2,2, &
                             2,2,2, &
                             0.1_prec,0.1_prec,0.1_prec, &
                             bcids)

    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call df%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    call f%SetName(1,"f")
    call df%SetName(1,"df")

    call f%SetEquation(1,'f = x*y*z')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    print*,"min, max (interior)",minval(f%interior),maxval(f%interior)

    call f%BoundaryInterp()
    call f%UpdateHost()
    print*,"min, max (boundary)",minval(f%boundary),maxval(f%boundary)

    call f%SideExchange(mesh)
    call f%UpdateHost()

    ! Set boundary conditions by prolonging the "boundary" attribute to the domain boundaries
    do iel = 1,f%nElem
      do iside = 1,6
        e2 = mesh%sideInfo(3,iside,iel) ! Neighboring Element ID
        if(e2 == 0) then
          do j = 1,f%interp%N+1
            do i = 1,f%interp%N+1
              f%extBoundary(i,j,iside,iel,1) = f%boundary(i,j,iside,iel,1)
            enddo
          enddo
        endif
      enddo
    enddo

    print*,"min, max (extboundary)",minval(f%extBoundary),maxval(f%extBoundary)

    call f%UpdateDevice()
    call f%AverageSides()
    call f%UpdateHost()
    print*,"min, max (avgboundary)",minval(f%avgBoundary),maxval(f%avgBoundary)

#ifdef ENABLE_GPU
    call f%MappedDGGradient(df%interior_gpu)
#else
    call f%MappedDGGradient(df%interior)
#endif
    call df%UpdateHost()

    ! Calculate diff from exact
    do iel = 1,mesh%nelem
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            df%interior(i,j,k,iel,1,1) = abs(df%interior(i,j,k,iel,1,1)- &
                                             geometry%x%interior(i,j,k,iel,1,2)* &
                                             geometry%x%interior(i,j,k,iel,1,3)) ! df/dx = y*z
            df%interior(i,j,k,iel,1,2) = abs(df%interior(i,j,k,iel,1,2)- &
                                             geometry%x%interior(i,j,k,iel,1,1)* &
                                             geometry%x%interior(i,j,k,iel,1,3)) ! df/dy = x*z
            df%interior(i,j,k,iel,1,3) = abs(df%interior(i,j,k,iel,1,3)- &
                                             geometry%x%interior(i,j,k,iel,1,1)* &
                                             geometry%x%interior(i,j,k,iel,1,2)) ! df/dy = x*y
          enddo
        enddo
      enddo
    enddo
    print*,"maxval(df_error)",maxval(df%interior(:,:,:,:,:,1)),tolerance
    print*,"maxval(df_error)",maxval(df%interior(:,:,:,:,:,2)),tolerance
    print*,"maxval(df_error)",maxval(df%interior(:,:,:,:,:,3)),tolerance

    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    ! Clean up
    call f%DissociateGeometry()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()
    call f%free()
    call df%free()

  endfunction mappedscalarbrgradient_3d_linear
endprogram test
