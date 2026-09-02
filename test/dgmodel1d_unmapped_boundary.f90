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

module dgmodel1d_unmapped_boundary_model
!! A minimal burgers1D extension that registers one boundary condition against a bcid of its
!! own choosing. No built-in 1-D model registers any boundary condition, so the test needs a
!! model that does in order to tell "registered" from "unregistered".
!!
!! The bcid is a plain module parameter, not one of the SELF_BC_* values: a bcid is any
!! integer the model author picks, and the SELF_BC_* constants are only the values the
!! built-in models happen to use.

  use self_Burgers1D
  use SELF_BoundaryConditions

  implicit none

  integer,parameter :: MYMODEL_BC_WALL = 4242
  !! Registered on the PARABOLIC list only. SetBoundaryCondition dispatches the
  !! hyperbolic list alone, so this id does not write solution%extBoundary and an
  !! endpoint carrying it is still unhandled.
  integer,parameter :: MYMODEL_BC_STRESSONLY = 4243

  type,extends(burgers1D) :: burgers1d_walled

  contains

    procedure :: AdditionalInit => AdditionalInit_burgers1d_walled

  endtype burgers1d_walled

contains

  subroutine AdditionalInit_burgers1d_walled(this)
    implicit none
    class(burgers1d_walled),intent(inout) :: this
    ! Local
    procedure(SELF_bcMethod),pointer :: bcfunc

    bcfunc => hbc1d_Wall_burgers1d_walled
    call this%hyperbolicBCs%RegisterBoundaryCondition( &
      MYMODEL_BC_WALL,"wall",bcfunc)

    ! Deliberately parabolic-only, to prove the scan does not accept it as handling the
    ! solution trace.
    bcfunc => pbc1d_StressOnly_burgers1d_walled
    call this%parabolicBCs%RegisterBoundaryCondition( &
      MYMODEL_BC_STRESSONLY,"stress_only",bcfunc)

  endsubroutine AdditionalInit_burgers1d_walled

  subroutine hbc1d_Wall_burgers1d_walled(bc,mymodel)
    !! Copies the interior trace to the exterior. This exists so that the bcid it is
    !! registered against is a REGISTERED one - the test asserts on the mapping counters and
    !! never steps the model. It is not a wall: for Burgers, u_ext = u_int is extrapolation
    !! and leaves the physical flux u^2/2 at the endpoint rather than imposing zero flux.
    implicit none
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel
    ! Local
    integer :: n,iEl,s

    select type(m => mymodel)
    class is(burgers1d_walled)
      do n = 1,bc%nBoundaries
        iEl = bc%elements(n)
        s = bc%sides(n)
        m%solution%extBoundary(s,iEl,1) = m%solution%boundary(s,iEl,1)
      enddo
    endselect

  endsubroutine hbc1d_Wall_burgers1d_walled

  subroutine pbc1d_StressOnly_burgers1d_walled(bc,mymodel)
    !! Parabolic-list counterpart, registered so that MYMODEL_BC_STRESSONLY is a real
    !! registration on that list and nothing else. Never invoked here: the test does not step.
    implicit none
    class(BoundaryCondition),intent(in) :: bc
    class(Model),intent(inout) :: mymodel

    if(.false.) print*,bc%nBoundaries,mymodel%nvar ! suppress unused-dummy-argument warnings

  endsubroutine pbc1d_StressOnly_burgers1d_walled

endmodule dgmodel1d_unmapped_boundary_model

program DGModel1D_Unmapped_Boundary
!! Regression test for issue #180 on the 1-D path.
!!
!! In 1-D an endpoint whose bcid matches no registration keeps the periodic default that
!! SetBoundaryCondition seeds, rather than a zero exterior state, but it is still not the
!! condition the bcid was meant to select. A bcid of 0 IS the periodic default and must not
!! be reported.
!!
!! Assertions:
!!  (a) an unregistered bcid on both endpoints counts both, and names the offending id;
!!  (b) the model's own registered bcid counts nothing;
!!  (c) the default bcid of 0 (periodic) counts nothing.

  use self_data
  use SELF_BoundaryConditions
  use dgmodel1d_unmapped_boundary_model

  implicit none

  integer,parameter :: nelem = 10
  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  integer,parameter :: bogusBCID = 777

  type(burgers1d_walled) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh1D),target :: mesh
  type(Geometry1D),target :: geometry
  type(BoundaryCondition),pointer :: bcnode

  call mesh%StructuredMesh(nElem=nelem,x=(/0.0_prec,1.0_prec/))
  call mesh%ResetBoundaryConditionType(bogusBCID,bogusBCID)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%tecplot_enabled = .false.

  ! (a) both endpoints are unmapped
  if(modelobj%nUnmappedBoundaries /= 2) then
    print*,"Error: expected 2 unmapped endpoints, counted ",modelobj%nUnmappedBoundaries
    stop 1
  endif
  if(modelobj%unmappedBoundaryID /= bogusBCID) then
    print*,"Error: expected the unregistered bcid ",bogusBCID," to be reported, got ", &
      modelobj%unmappedBoundaryID
    stop 1
  endif

  ! (b) the model's own bcid is registered, so nothing is unmapped
  call mesh%ResetBoundaryConditionType(MYMODEL_BC_WALL,MYMODEL_BC_WALL)
  call modelobj%MapBoundaryConditions()
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: the registered bcid ",MYMODEL_BC_WALL," was counted as unmapped ", &
      modelobj%nUnmappedBoundaries," times."
    stop 1
  endif

  ! (c) bcid 0 is the deliberate periodic default and is never a fault
  call mesh%ResetBoundaryConditionType(0,0)
  call modelobj%MapBoundaryConditions()
  if(modelobj%nUnmappedBoundaries /= 0) then
    print*,"Error: the periodic default (bcid 0) was counted as unmapped ", &
      modelobj%nUnmappedBoundaries," times."
    stop 1
  endif

  ! (d) a parabolic-only registration does NOT handle the solution trace, so both endpoints
  ! are still unmapped. SetBoundaryCondition dispatches the hyperbolic list alone.
  call mesh%ResetBoundaryConditionType(MYMODEL_BC_STRESSONLY,MYMODEL_BC_STRESSONLY)
  call modelobj%MapBoundaryConditions()
  if(modelobj%nUnmappedBoundaries /= 2) then
    print*,"Error: a parabolic-only bcid was accepted as handling the solution trace; ", &
      "expected 2 unmapped endpoints, counted ",modelobj%nUnmappedBoundaries
    stop 1
  endif

  ! (e) re-tagging away from a mapped bcid must clear that condition's mapping. Otherwise
  ! SetBoundaryCondition keeps invoking it on endpoints it no longer owns, and the counters
  ! above would describe something other than what runs. MYMODEL_BC_WALL owned both endpoints
  ! at (b); after (d) it must own none.
  bcnode => modelobj%hyperbolicBCs%GetBCForID(MYMODEL_BC_WALL)
  if(.not. associated(bcnode)) then
    print*,"Error: MYMODEL_BC_WALL is not registered; (e) is vacuous."
    stop 1
  endif
  if(bcnode%nBoundaries /= 0) then
    print*,"Error: MYMODEL_BC_WALL still claims ",bcnode%nBoundaries," endpoints after the "// &
      "mesh was re-tagged away from it - a stale mapping survived the re-map."
    stop 1
  endif

  print*,"dgmodel1d_unmapped_boundary : all assertions passed"

  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram DGModel1D_Unmapped_Boundary
