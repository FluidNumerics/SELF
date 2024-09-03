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

  exit_code = mappedvectordgdivergence_2d_cpu_linear()
  stop exit_code

contains
  integer function mappedvectordgdivergence_2d_cpu_linear() result(r)

    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_MappedVector_2D

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
    type(MPILayer),target :: decomp
    character(LEN=255) :: WORKSPACE
    integer :: i,j,iel
    real(prec) :: nhat(1:2),nmag

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    call decomp%Init(enableMPI=.false.)

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Create a uniform block mesh
    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5",decomp)

    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call df%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    call f%SetEquation(1,1,'f = x') ! x-component
    call f%SetEquation(2,1,'f = y') ! y-component

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    print*,"min, max (interior)",minval(f%interior),maxval(f%interior)

    call f%boundaryInterp()
    call f%UpdateHost()

    do iEl = 1,f%nElem
      do j = 1,4
        do i = 1,f%interp%N+1

          ! Get the boundary normals on cell edges from the mesh geometry
          nhat(1:2) = geometry%nHat%boundary(i,j,iEl,1,1:2)
          nmag = geometry%nScale%boundary(i,j,iEl,1)

          f%boundaryNormal(i,j,iEl,1) = (f%boundary(i,j,iEl,1,1)*nhat(1)+ &
                                         f%boundary(i,j,iEl,1,2)*nhat(2))*nmag

        enddo
      enddo
    enddo

    call f%UpdateDevice()

#ifdef ENABLE_GPU
    call f%MappedDGDivergence(df%interior_gpu)
#else
    call f%MappedDGDivergence(df%interior)
#endif
    call df%UpdateHost()

    ! Calculate diff from exact
    df%interior = abs(df%interior-2.0_prec)

    print*,"absmax error :",maxval(df%interior)
    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    ! Clean up
    call f%DissociateGeometry()
    call decomp%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()
    call f%free()
    call df%free()

  endfunction mappedvectordgdivergence_2d_cpu_linear
endprogram test