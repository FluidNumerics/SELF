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

program self_external_consumer
  !! Builds against an installed SELF through find_package(SELF) and checks that
  !! the package hands this program what it needs.
  !!
  !! A 4x4x4 structured box of side 1 m is built on one rank and the quadrature
  !! volume
  !!
  !!   V = sum_e sum_ijk w_i w_j w_k J_ijk^e
  !!
  !! is accumulated. It must equal 1 m^3 to within 100 machine epsilon, and the
  !! mesh must hold 64 elements. The program also checks that the precision
  !! selected by the DOUBLE_PRECISION definition it was compiled with, which
  !! comes from SELF::self, is the precision of SELF's real kind prec.

  use iso_c_binding,only:c_double,c_float
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_3D
  use SELF_Geometry_3D

  implicit none

  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  real(prec),parameter :: dx = 0.25_prec ! element width [m]
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  integer :: i,j,k,iel
  real(prec) :: volume

  ! The definition must agree with the library, or prec seen here and the
  ! kind SELF was compiled with would differ.
#ifdef DOUBLE_PRECISION
  if(prec /= c_double) then
    print*,"Error: DOUBLE_PRECISION is defined but SELF was built in single precision"
    stop 1
  endif
#else
  if(prec /= c_float) then
    print*,"Error: DOUBLE_PRECISION is not defined but SELF was built in double precision"
    stop 1
  endif
#endif

  bcids(1:6) = 0
  call mesh%StructuredMesh(2,2,2,2,2,2,dx,dx,dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  volume = 0.0_prec
  do iel = 1,mesh%nElem
    do k = 1,controlDegree+1
      do j = 1,controlDegree+1
        do i = 1,controlDegree+1
          volume = volume+interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)* &
                   geometry%J%interior(i,j,k,iel,1)
        enddo
      enddo
    enddo
  enddo

  if(mesh%nElem /= 64) then
    print*,"Error: expected 64 elements, saw ",mesh%nElem
    stop 1
  endif
  if(abs(volume-1.0_prec) > 100.0_prec*epsilon(1.0_prec)) then
    print*,"Error: expected unit volume, saw ",volume
    stop 1
  endif

  print*,"self_external_consumer : volume = ",volume

  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram self_external_consumer
