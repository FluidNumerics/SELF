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

  exit_code = mappedscalarderivative_1d_cpu_constant()
  stop exit_code

contains
  integer function mappedscalarderivative_1d_cpu_constant() result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_MappedScalar_1D
    use SELF_Mesh_1D
    use SELF_Geometry_1D
    use SELF_MPI

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
    integer,parameter :: nelem = 100
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(MappedScalar1D) :: f
    type(MappedScalar1D) :: df
    type(Lagrange),target :: interp
    type(Mesh1D),target :: mesh
    type(Geometry1D),target :: geometry
    type(MPILayer),target :: decomp

    call mesh%UniformBlockMesh(nGeo=1, &
                               nElem=nelem, &
                               x=(/0.0_prec,10.0_prec/))

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    ! Initialize scalars
    call f%Init(interp,nvar,nelem)
    call df%Init(interp,nvar,nelem)
    call f%AssociateGeometry(geometry)
    call df%AssociateGeometry(geometry)

    call f%SetEquation(1,'f = 1.0')

    call f%SetInteriorFromEquation(0.0_prec)

#ifdef ENABLE_GPU
    call f%MappedDerivative(df%interior_gpu)
#else
    call f%MappedDerivative(df%interior)
#endif

    call df%UpdateHost()

    ! Calculate diff from exact
    df%interior = abs(df%interior-0.0_prec)

    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    ! Clean up
    call f%DissociateGeometry()
    call df%DissociateGeometry()

    call mesh%Free()
    call geometry%Free()
    call interp%free()
    call f%free()
    call df%free()

  endfunction mappedscalarderivative_1d_cpu_constant
endprogram test
