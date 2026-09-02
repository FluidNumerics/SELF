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

program DGModel2D_Unmapped_Boundary_Mortar
!! Guard for the unmapped-boundary scan added for issue #180 on a nonconforming mesh.
!!
!! Mortar sides carry sideInfo(3) = 0 exactly like a physical boundary, and sideInfo(5) = 0,
!! so a scan that keyed only on sideInfo(3) would report every mortar side as an unmapped
!! boundary. The scan excludes them via sideInfo(1), which holds the mortar index.
!!
!! ResetBoundaryConditionType carries the same exclusion, for the same reason, so this test
!! also round-trips through it: a reset that tagged the mortar sides would put interior faces
!! into a boundary condition's element/side list, and that condition would then overwrite the
!! exterior state the mortar exchange had just written.
!!
!! Assertions:
!!  (a) the mesh really does carry mortars, so the guard is exercised;
!!  (b) after a ResetBoundaryConditionType, no registered boundary condition claims a
!!      mortar side;
!!  (c) with every domain edge tagged to a registered bcid, nothing is reported unmapped.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d
  use SELF_BoundaryConditions

  implicit none

  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  real(prec),parameter :: dx = 0.1_prec

  type(LinearEuler2D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  ! Radiation is registered by AdditionalInit_LinearEuler2D_t
  bcids(1:4) = [SELF_BC_RADIATION, & ! south
                SELF_BC_RADIATION, & ! east
                SELF_BC_RADIATION, & ! north
                SELF_BC_RADIATION] ! west

  ! Three-element mesh with a single 2:1 mortar interface at x = 2*dx
  call mesh%SimpleMortarMesh(dx,bcids)

  ! (a) the guard is only meaningful on a mesh that has mortar sides
  if(mesh%nMortars <= 0) then
    print*,"Error: SimpleMortarMesh produced no mortars; the guard is not exercised."
    stop 1
  endif

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)

  ! Re-tag through ResetBoundaryConditionType and re-map. This is the model-level
  ! consequence of the mortar guard in that routine: if it wrote a bcid onto the mortar
  ! sides, MapBoundaryConditions would hand them to the radiation boundary condition, whose
  ! element/side list would then include interior faces.
  ! MapBoundaryConditions re-maps on the host only, so on a GPU build the device BC
  ! arrays are left stale. Nothing steps after this point.
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)
  call modelobj%MapBoundaryConditions()
  call CheckMortarSidesUnclaimed(mesh,modelobj)

  ! (b) mortar sides must not be mistaken for unmapped physical boundaries
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: ",modelobj%nUnmappedBoundaries, &
      " boundary edges reported unmapped on a fully tagged mortar mesh; bcid ", &
      modelobj%unmappedBoundaryID
    stop 1
  endif

  print*,"dgmodel2d_unmapped_boundary_mortar : all assertions passed"

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

contains

  subroutine CheckMortarSidesUnclaimed(mesh,modelobj)
    !! Assert that no registered boundary condition claims a mortar side.
    implicit none
    type(Mesh2D),intent(in) :: mesh
    type(LinearEuler2D),intent(in) :: modelobj
    ! Local
    type(BoundaryCondition),pointer :: bc
    integer :: m,k,n,e,s

    bc => modelobj%hyperbolicBCs%head
    do while(associated(bc))
      do n = 1,bc%nBoundaries
        do m = 1,mesh%nMortars
          ! big side, then the two small sides; element ids are global, and this test is
          ! serial, so they are already local
          if(bc%elements(n) == mesh%mortarInfo(1,m) .and. &
             bc%sides(n) == mesh%mortarInfo(2,m)) then
            print*,"Error: boundary condition ",trim(bc%bcname)," claims mortar ",m, &
              " big side (element ",bc%elements(n),", side ",bc%sides(n),")"
            stop 1
          endif
          do k = 1,2
            e = mesh%mortarInfo(2*k+1,m)
            s = mesh%mortarInfo(2*k+2,m)/10
            if(bc%elements(n) == e .and. bc%sides(n) == s) then
              print*,"Error: boundary condition ",trim(bc%bcname)," claims mortar ",m, &
                " small side ",k," (element ",e,", side ",s,")"
              stop 1
            endif
          enddo
        enddo
      enddo
      bc => bc%next
    enddo

  endsubroutine CheckMortarSidesUnclaimed

endprogram DGModel2D_Unmapped_Boundary_Mortar
