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
!! Cross-checks Mesh3D%ExtrudedMesh against HOHQMesh's own extrusion over a bottom
!! topography, using the fixtures in share/mesh/Seamount (generated with HOHQMesh
!! 1.5.6 from the .control files stored beside them).
!!
!! Seamount2D.mesh is the quad mesh of a circle of radius 2 (21 elements, order 4).
!! Seamount3D.mesh sweeps the same quad mesh upward from the bottom z = f(x,y),
!! f = 0.5 exp(-((x-0.5)^2 + (y-0.3)^2)), to the flat top z = H = 3 in nz = 3
!! subdivisions. HOHQMesh (Source/3DSource/SimpleSweep.f90) spreads the topography
!! over the whole column: level j sits at z_j = f (1 - j/nz) + H j/nz, and elements are
!! written layer by layer, element (k-1)*nElem2D + e2 being layer k above quad e2. In
!! SELF's convention, z' = z - H and h = H - f, those are the sigma levels
!! z'_j = sigma_j h with sigma_j = j/nz - 1.
!!
!! The test extrudes Seamount2D.mesh with those levels and compares, element by
!! element, with Seamount3D.mesh read by Read_HOHQMesh and shifted down by H:
!!
!!   1. the corner nodes of every element agree to round-off (asserted). This checks
!!      the level formula, the column ordering against HOHQMesh's layer ordering, and
!!      that both meshes orient the hexahedra alike.
!!   2. the full sets of geometry nodes are compared and the differences reported per
!!      layer. They are not asserted: HOHQMesh's ISM file carries curved face data only
!!      for faces on the topography or on a curved boundary, so the faces of the upper
!!      layers that lie inside the domain are rebuilt from their corners, whereas the
!!      extruded mesh samples the terrain-following levels at every node.
!!   3. on the bottom faces of the bottom layer, which HOHQMesh does write as curved
!!      topography data, all nodes agree to `bottomTolerance` (asserted).
!!   4. the horizontal coordinates of all nodes agree to round-off (asserted), so the
!!      differences reported in 2 are purely vertical.

  use SELF_Constants
  use SELF_Mesh_2D
  use SELF_Mesh_3D

  implicit none
  integer :: exit_code

  exit_code = mesh3d_extruded_hohqmesh()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function mesh3d_extruded_hohqmesh() result(r)
    implicit none
    ! Local
    integer,parameter :: nz = 3
    real(prec),parameter :: topHeight = 3.0_prec ! H, the height of the flat top
    real(prec),parameter :: cornerTolerance = 1.0e-12_prec
    real(prec),parameter :: bottomTolerance = 1.0e-12_prec
    type(Mesh2D) :: surf
    type(Mesh3D) :: mesh,hohq
    character(LEN=255) :: WORKSPACE
    real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
    real(prec) :: sigma(0:nz),x,y,d(1:3)
    real(prec) :: errCorner,errBottom,errLayer(1:nz),errHorizontal
    integer :: nGeo,e2,k,i,j,l,e3,eH,n

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call surf%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Seamount/Seamount2D.mesh")
    call hohq%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Seamount/Seamount3D.mesh")

    r = 1
    nGeo = surf%nGeo
    n = nGeo+1
    if(hohq%nElem /= nz*surf%nElem .or. hohq%nGeo /= nGeo) then
      print*,"FAIL: the fixtures do not match: ",hohq%nElem,surf%nElem,hohq%nGeo,nGeo
      return
    endif

    allocate(h(1:n,1:n,1:surf%nElem),zLevels(1:n,1:n,0:nz,1:surf%nElem))
    do e2 = 1,surf%nElem
      do j = 1,n
        do i = 1,n
          x = surf%nodeCoords(1,i,j,e2)
          y = surf%nodeCoords(2,i,j,e2)
          h(i,j,e2) = topHeight-0.5_prec*exp(-((x-0.5_prec)**2+(y-0.3_prec)**2))
        enddo
      enddo
    enddo
    do l = 0,nz
      sigma(l) = real(l,prec)/real(nz,prec)-1.0_prec
    enddo
    call ExtrudedSigmaLevels(h,sigma,zLevels)
    call mesh%ExtrudedMesh(surf,zLevels,1,2)

    errCorner = 0.0_prec
    errBottom = 0.0_prec
    errLayer = 0.0_prec
    errHorizontal = 0.0_prec
    do e2 = 1,surf%nElem
      do k = 1,nz
        e3 = (e2-1)*nz+k
        eH = (k-1)*surf%nElem+e2
        do l = 1,n
          do j = 1,n
            do i = 1,n
              d = hohq%nodeCoords(1:3,i,j,l,eH)-mesh%nodeCoords(1:3,i,j,l,e3)
              d(3) = d(3)-topHeight
              errLayer(k) = max(errLayer(k),maxval(abs(d)))
              errHorizontal = max(errHorizontal,maxval(abs(d(1:2))))
              if((i == 1 .or. i == n) .and. (j == 1 .or. j == n) .and. (l == 1 .or. l == n)) then
                errCorner = max(errCorner,maxval(abs(d)))
              endif
              if(k == 1 .and. l == 1) errBottom = max(errBottom,maxval(abs(d)))
            enddo
          enddo
        enddo
      enddo
    enddo

    print*,"max corner-node difference         ",errCorner," (tolerance ",cornerTolerance,")"
    print*,"max bottom-face difference, layer 1 ",errBottom," (tolerance ",bottomTolerance,")"
    print*,"max horizontal node difference     ",errHorizontal," (tolerance ",cornerTolerance,")"
    do k = 1,nz
      print*,"max node difference, layer ",k," : ",errLayer(k)," (reported)"
    enddo

    r = 0
    if(errCorner > cornerTolerance) then
      print*,"FAIL: corner nodes differ from HOHQMesh's extrusion"
      r = 1
    endif
    if(errHorizontal > cornerTolerance) then
      print*,"FAIL: horizontal node coordinates differ from HOHQMesh's"
      r = 1
    endif
    if(errBottom > bottomTolerance) then
      print*,"FAIL: the bottom topography differs from HOHQMesh's"
      r = 1
    endif

    call mesh%Free()
    call hohq%Free()
    call surf%Free()
    deallocate(h,zLevels)

  endfunction mesh3d_extruded_hohqmesh

endprogram test
