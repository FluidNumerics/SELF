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

   exit_code = mappedvectordgdivergence_3d_linear()
   if (exit_code /= 0) then
      stop exit_code
   end if

contains
   integer function mappedvectordgdivergence_3d_linear() result(r)

      use SELF_Constants
      use SELF_Lagrange
      use SELF_Mesh_3D
      use SELF_Geometry_3D
      use SELF_MappedScalar_3D
      use SELF_MappedVector_3D

      implicit none

      integer, parameter :: controlDegree = 7
      integer, parameter :: targetDegree = 16
      integer, parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
      real(prec), parameter :: tolerance = 10.0_prec**(-7)
#else
      real(prec), parameter :: tolerance = 10.0_prec**(-3)
#endif
      type(Lagrange), target :: interp
      type(Mesh3D), target :: mesh
      type(SEMHex), target :: geometry
      type(MappedVector3D) :: f
      type(MappedScalar3D) :: df
      integer :: i, j, k, iel, e2
      real(prec) :: nhat(1:3), nmag, fx, fy, fz
      integer :: bcids(1:6)

      ! Create a uniform block mesh
      bcids(1:6) = [SELF_BC_PRESCRIBED, & ! Bottom
                    SELF_BC_PRESCRIBED, & ! South
                    SELF_BC_PRESCRIBED, & ! East
                    SELF_BC_PRESCRIBED, & ! North
                    SELF_BC_PRESCRIBED, & ! West
                    SELF_BC_PRESCRIBED] ! Top

      call mesh%StructuredMesh(5, 5, 5, &
                               2, 2, 2, &
                               0.1_prec, 0.1_prec, 0.1_prec, &
                               bcids)

      ! Create an interpolant
      call interp%Init(N=controlDegree, &
                       controlNodeType=GAUSS, &
                       M=targetDegree, &
                       targetNodeType=UNIFORM)

      ! Generate geometry (metric terms) from the mesh elements
      call geometry%Init(interp, mesh%nElem)
      call geometry%GenerateFromMesh(mesh)

      call f%Init(interp, nvar, mesh%nelem)
      call df%Init(interp, nvar, mesh%nelem)
      call f%AssociateGeometry(geometry)

      call f%SetEquation(1, 1, 'f = x') ! x-component
      call f%SetEquation(2, 1, 'f = y') ! y-component
      call f%SetEquation(3, 1, 'f = 0') ! z-component

      call f%SetInteriorFromEquation(geometry, 0.0_prec)
      print *, "min, max (interior)", minval(f%interior), maxval(f%interior)
      call f%boundaryInterp()

      print *, "Exchanging data on element faces"

      call f%SideExchange(mesh)
      call f%UpdateHost()

      print *, "Setting boundary conditions"
      ! Set boundary conditions
      do iEl = 1, f%nElem
         do k = 1, 6
            e2 = mesh%sideInfo(3, k, iel) ! Neighbor Element (global id)

            if (e2 == 0) then ! Exterior edge
               do j = 1, f%interp%N + 1
                  do i = 1, f%interp%N + 1
                     f%extboundary(i, j, k, iEl, 1, 1) = f%boundary(i, j, k, iEl, 1, 1)
                     f%extboundary(i, j, k, iEl, 1, 2) = f%boundary(i, j, k, iEl, 1, 2)
                     f%extboundary(i, j, k, iEl, 1, 3) = f%boundary(i, j, k, iEl, 1, 3)
                  end do
               end do
            end if
         end do
      end do

      print *, "Calculating boundary normal flux"
      do iEl = 1, f%nElem
         do k = 1, 6
            do j = 1, f%interp%N + 1
               do i = 1, f%interp%N + 1

                  ! Get the boundary normals on cell edges from the mesh geometry
                  nhat(1:3) = geometry%nHat%boundary(i, j, k, iEl, 1, 1:3)
                  nmag = geometry%nScale%boundary(i, j, k, iEl, 1)
                  fx = 0.5*(f%boundary(i, j, k, iEl, 1, 1) + f%extboundary(i, j, k, iEl, 1, 1))
                  fy = 0.5*(f%boundary(i, j, k, iEl, 1, 2) + f%extboundary(i, j, k, iEl, 1, 2))
                  fz = 0.5*(f%boundary(i, j, k, iEl, 1, 3) + f%extboundary(i, j, k, iEl, 1, 3))

                  f%boundaryNormal(i, j, k, iEl, 1) = (fx*nhat(1) + fy*nhat(2) + fz*nhat(3))*nmag
               end do
            end do
         end do
      end do
      call f%UpdateDevice()

#ifdef ENABLE_GPU
      call f%MappedDGDivergence(df%interior_gpu)
#else
      call f%MappedDGDivergence(df%interior)
#endif
      call df%UpdateHost()

      ! Calculate diff from exact
      df%interior = abs(df%interior - 2.0_prec)

      if (maxval(df%interior) <= tolerance) then
         r = 0
      else
         print *, "max error (tolerance)", maxval(df%interior), tolerance
         r = 1
      end if

      ! Clean up
      call f%DissociateGeometry()
      call geometry%Free()
      call mesh%Free()
      call interp%Free()
      call f%free()
      call df%free()

   end function mappedvectordgdivergence_3d_linear
end program test
