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

module dgmodel3d_unmapped_boundary_mortar_model
!! A LinearEuler3D extension registering bcid 0, for the same reason as the 2-D twin: zero is
!! a legal boundary condition id under the documented contract, and it is also what a mortar
!! face carries in sideInfo(5). A forward mapping pass keyed on sideInfo(5) alone would hand
!! every mortar face to this condition, and SetBoundaryCondition would then overwrite what the
!! mortar exchange had just written. MapBoundaryConditions is a separate routine per dimension,
!! so the 2-D test cannot cover the 3-D six-face loops.

  use self_lineareuler3d
  use SELF_BoundaryConditions

  implicit none

  integer,parameter :: MORTARTEST3D_BC_ZERO = 0

  type,extends(LinearEuler3D) :: lineareuler3d_zerobc

  contains

    procedure :: AdditionalInit => AdditionalInit_lineareuler3d_zerobc

  endtype lineareuler3d_zerobc

contains

  subroutine AdditionalInit_lineareuler3d_zerobc(this)
    implicit none
    class(lineareuler3d_zerobc),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    ! Keep the parent's registrations, through the parent COMPONENT: the CPU LinearEuler3D
    ! inherits AdditionalInit_LinearEuler3D_t while the GPU one overrides it with device
    ! kernels, and naming either routine directly would be wrong on the other backend.
    call this%LinearEuler3D%AdditionalInit()

    bcfunc => hbc3d_Zero_lineareuler3d_zerobc
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      MORTARTEST3D_BC_ZERO,"zero",bcfunc)

  endsubroutine AdditionalInit_lineareuler3d_zerobc

  subroutine hbc3d_Zero_lineareuler3d_zerobc(bc,mymodel)
    !! Never expected to run on this mesh: nothing is tagged 0, and mortar faces must be
    !! excluded from the mapping. The assertion is on bc%nBoundaries, not on any state this
    !! would write.
    implicit none
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    if(.false.) print*,bc%nBoundaries,mymodel%nvar ! suppress unused-dummy-argument warnings

  endsubroutine hbc3d_Zero_lineareuler3d_zerobc

endmodule dgmodel3d_unmapped_boundary_mortar_model

program DGModel3D_Unmapped_Boundary_Mortar
!! Model-level cover for the 3-D mortar exclusion in MapBoundaryConditions.
!!
!! Assertions:
!!  (a) the mesh really does carry mortars, so the guard is exercised;
!!  (b) the bcid-0 condition owns no faces, so the forward passes are not matching mortar
!!      faces on their sideInfo(5) = 0;
!!  (c) no registered condition claims a mortar face after a ResetBoundaryConditionType;
!!  (d) the reverse scan reports no false positives on a fully tagged mortar mesh.

  use self_data
  use self_lineareuler3d
  use self_mesh_3d
  use SELF_BoundaryConditions
  use dgmodel3d_unmapped_boundary_mortar_model

  implicit none

  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  real(prec),parameter :: dx = 0.1_prec

  type(lineareuler3d_zerobc) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)

  ! Radiation is registered by AdditionalInit_LinearEuler3D_t
  bcids(1:6) = SELF_BC_RADIATION

  call mesh%SimpleMortarMesh(dx,bcids)

  ! (a) the guard is only meaningful on a mesh that has mortar faces
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

  call CheckZeroBCUnclaimed(modelobj)
  call CheckMortarFacesUnclaimed(mesh,modelobj)

  ! (d) mortar faces must not be mistaken for unmapped physical boundaries
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: ",modelobj%nUnmappedBoundaries, &
      " boundary faces reported unmapped on a fully tagged mortar mesh; bcid ", &
      modelobj%unmappedBoundaryID
    stop 1
  endif

  ! Re-tag and re-map: ResetBoundaryConditionType carries the same exclusion, and the mapping
  ! must still leave the mortar faces alone. MapBoundaryConditions re-maps on the host only,
  ! so on a GPU build the device BC arrays are left stale; nothing steps after this point.
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION)
  call modelobj%MapBoundaryConditions()
  call CheckZeroBCUnclaimed(modelobj)
  call CheckMortarFacesUnclaimed(mesh,modelobj)

  print*,"dgmodel3d_unmapped_boundary_mortar : all assertions passed"

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

contains

  subroutine CheckZeroBCUnclaimed(modelobj)
    !! (b) the bcid-0 condition owns no faces on a mesh where nothing is tagged 0.
    implicit none
    type(lineareuler3d_zerobc),intent(in) :: modelobj
    ! Local
    type(BoundaryCondition),pointer :: bc

    bc => modelobj%hyperbolicBCs%GetBCForID(MORTARTEST3D_BC_ZERO)
    if(.not. associated(bc)) then
      print*,"Error: the bcid-0 boundary condition was not registered; (b) is vacuous."
      stop 1
    endif
    if(bc%nBoundaries /= 0) then
      print*,"Error: the bcid-0 boundary condition claimed ",bc%nBoundaries," faces on a mesh"// &
        " where nothing is tagged 0 - the forward mapping passes are matching mortar faces."
      stop 1
    endif

  endsubroutine CheckZeroBCUnclaimed

  subroutine CheckMortarFacesUnclaimed(mesh,modelobj)
    !! (c) no registered condition claims a mortar face. The big face is mortarInfo(1:2,m) and
    !! the four small faces are (2*q+1, 2*q+2/10) for q = 1..4; element ids there are global,
    !! and this test is serial, so they are already local.
    implicit none
    type(Mesh3D),intent(in) :: mesh
    type(lineareuler3d_zerobc),intent(in) :: modelobj
    ! Local
    type(BoundaryCondition),pointer :: bc
    integer :: m,q,n,e,s

    bc => modelobj%hyperbolicBCs%head
    do while(associated(bc))
      do n = 1,bc%nBoundaries
        do m = 1,mesh%nMortars
          if(bc%elements(n) == mesh%mortarInfo(1,m) .and. &
             bc%sides(n) == mesh%mortarInfo(2,m)) then
            print*,"Error: boundary condition ",trim(bc%bcname)," claims mortar ",m, &
              " big face (element ",bc%elements(n),", face ",bc%sides(n),")"
            stop 1
          endif
          do q = 1,4
            e = mesh%mortarInfo(2*q+1,m)
            s = mesh%mortarInfo(2*q+2,m)/10
            if(bc%elements(n) == e .and. bc%sides(n) == s) then
              print*,"Error: boundary condition ",trim(bc%bcname)," claims mortar ",m, &
                " small face ",q," (element ",e,", face ",s,")"
              stop 1
            endif
          enddo
        enddo
      enddo
      bc => bc%next
    enddo

  endsubroutine CheckMortarFacesUnclaimed

endprogram DGModel3D_Unmapped_Boundary_Mortar
