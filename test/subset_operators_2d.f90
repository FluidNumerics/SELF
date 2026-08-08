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

program Subset_Operators_2D
!! Bitwise regression test for the element/mortar subsetting added for local time stepping
!! (AMR Stage 7).
!!
!! Every tensor-product operator in the 2-D pipeline now traverses an index list instead of
!! 1:nElem, defaulting to the identity list when no subset is supplied. That default must be
!! indistinguishable from the original code, and a genuine subset must produce, on the
!! elements it covers, exactly what the full evaluation produced there. Both are checked here
!! to the LAST BIT: the whole point of routing the loops through an index list is that the
!! arithmetic and its ordering are untouched, so anything short of bitwise equality means the
!! restriction changed the discretization.
!!
!! Two meshes are exercised:
!!
!!   (a) a conforming structured mesh - covers BoundaryInterp, SideExchange, BoundaryFlux,
!!       FluxMethod, SourceMethod, MappedDGDivergence and the dSdt loop;
!!   (b) the three-element SimpleMortarMesh - additionally covers MortarExchange and
!!       MortarFluxCollect, whose subset is a list of mortars rather than elements.
!!
!! The subsets used are the odd- and even-numbered elements. Their union is the whole mesh,
!! and neither is closed under "is a neighbour of", so the test also pins down the claim that
!! restricting BoundaryInterp to the active elements is safe: the neighbour traces a subset
!! evaluation reads were written by an earlier call and must still be the right ones.

  use self_data
  use self_lineareuler2d
  use self_mesh_2d

  implicit none

  integer,parameter :: controlDegree = 4
  integer,parameter :: targetDegree = 7
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: c0 = 1.0_prec
  real(prec),parameter :: rho0 = 1.0_prec
  real(prec),parameter :: amp = 1.0e-4_prec
  real(prec),parameter :: Lr = 0.2_prec

  ! SELF owns the MPI lifecycle and releases it when the last live mesh is freed, so a
  ! keepalive mesh spans both cases; without it, freeing the first case's mesh would
  ! finalize MPI and the second case could not start (see domaindecomposition_two_meshes).
  type(Mesh2D),target :: keepalive
  integer :: keepaliveBCs(1:4)

  keepaliveBCs(1:4) = [SELF_BC_RADIATION,SELF_BC_RADIATION, &
                       SELF_BC_RADIATION,SELF_BC_RADIATION]
  call keepalive%StructuredMesh(2,2,1,1,dx,dx,keepaliveBCs)

  if(.not. RunCase("structured")) stop 1
  if(.not. RunCase("mortar")) stop 1

  call keepalive%Free()

  print*,"subset operators: full and subset tendencies agree bitwise on both meshes"

contains

  logical function RunCase(meshKind) result(ok)
    !! Build a model on the named mesh, evaluate the tendency over the whole mesh, then
    !! re-evaluate it over the odd and then the even elements, and require the restricted
    !! results to reproduce the full one bitwise on the elements each covers.
    implicit none
    character(*),intent(in) :: meshKind
    ! Local
    type(LinearEuler2D) :: modelobj
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    integer :: bcids(1:4)
    integer :: nEl,nVar,Np,i,nOdd,nEven
    ! Subset lists are pointers because the subset-aware operators take pointer dummies;
    ! see ResolveIndexList in SELF_Data for why that is the safe association.
    integer,pointer,contiguous :: oddElem(:),evenElem(:),allMortar(:)
    real(prec),allocatable :: dSdtFull(:,:,:,:)

    ok = .false.

    bcids(1:4) = [SELF_BC_RADIATION,SELF_BC_RADIATION, &
                  SELF_BC_RADIATION,SELF_BC_RADIATION]

    if(meshKind == "mortar") then
      call mesh%SimpleMortarMesh(dx,bcids)
    else
      call mesh%StructuredMesh(5,4,2,2,dx,dx,bcids)
    endif

    call interp%Init(N=controlDegree,controlNodeType=GAUSS, &
                     M=targetDegree,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call modelobj%Init(mesh,geometry)
    modelobj%prescribed_bcs_enabled = .false.
    modelobj%tecplot_enabled = .false.
    modelobj%rho0 = rho0

    ! An off-centre pulse, so the tendency is non-trivial on essentially every element.
    call modelobj%SphericalSoundWave(amp,Lr,2.0_prec*dx,1.5_prec*dx,c0)

    nEl = mesh%nElem
    nVar = modelobj%solution%nVar
    Np = controlDegree+1

    ! ---- Reference : the full-mesh tendency, with no subset set ----
    call modelobj%CalculateTendency()
    allocate(dSdtFull(1:Np,1:Np,1:nEl,1:nVar))
    dSdtFull = modelobj%dSdt%interior(1:Np,1:Np,1:nEl,1:nVar)

    if(all(dSdtFull == 0.0_prec)) then
      print*,"Error: reference tendency is identically zero on the ",meshKind," mesh; ", &
        "the test would pass vacuously."
      return
    endif

    ! ---- Build the odd/even element subsets ----
    nOdd = (nEl+1)/2
    nEven = nEl-nOdd
    allocate(oddElem(1:nOdd),evenElem(1:max(nEven,1)))
    nOdd = 0
    nEven = 0
    do i = 1,nEl
      if(mod(i,2) == 1) then
        nOdd = nOdd+1
        oddElem(nOdd) = i
      else
        nEven = nEven+1
        evenElem(nEven) = i
      endif
    enddo

    ! Every mortar is active in both passes: a mortar couples two elements that the odd/even
    ! split may place on opposite sides, so restricting the mortar list here would drop a
    ! coupling that the reference evaluation included.
    allocate(allMortar(1:max(mesh%nMortars,1)))
    do i = 1,mesh%nMortars
      allMortar(i) = i
    enddo

    ! ---- Restricted evaluations ----
    ! Poison dSdt first, so that an operator which silently skips an element it should have
    ! written is caught rather than reading a stale correct value.
    modelobj%dSdt%interior = -huge(1.0_prec)

    if(mesh%nMortars > 0) modelobj%activeMortar => allMortar

    modelobj%activeElem => oddElem
    call modelobj%CalculateTendency()
    if(.not. Matches(modelobj,dSdtFull,oddElem(1:nOdd),Np,nVar,meshKind,"odd")) return

    modelobj%activeElem => evenElem(1:nEven)
    call modelobj%CalculateTendency()
    if(.not. Matches(modelobj,dSdtFull,evenElem(1:nEven),Np,nVar,meshKind,"even")) return

    ! ---- Back to the unrestricted path : it must still reproduce the reference exactly ----
    modelobj%activeElem => null()
    modelobj%activeMortar => null()
    modelobj%dSdt%interior = -huge(1.0_prec)
    call modelobj%CalculateTendency()
    if(any(modelobj%dSdt%interior(1:Np,1:Np,1:nEl,1:nVar) /= dSdtFull)) then
      print*,"Error: repeating the unrestricted tendency on the ",meshKind, &
        " mesh did not reproduce the reference bitwise."
      return
    endif

    deallocate(dSdtFull,oddElem,evenElem,allMortar)
    call modelobj%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()

    ok = .true.

  endfunction RunCase

  logical function Matches(modelobj,dSdtFull,elems,Np,nVar,meshKind,label) result(ok)
    !! Require the restricted tendency to equal the reference bitwise on every element of
    !! elems. Elements outside elems are deliberately not checked: the restricted evaluation
    !! is expected to leave them alone (here, still poisoned).
    implicit none
    type(LinearEuler2D),intent(in) :: modelobj
    real(prec),intent(in) :: dSdtFull(:,:,:,:)
    integer,intent(in) :: elems(:)
    integer,intent(in) :: Np,nVar
    character(*),intent(in) :: meshKind,label
    ! Local
    integer :: k,iel

    ok = .true.
    do k = 1,size(elems)
      iel = elems(k)
      if(any(modelobj%dSdt%interior(1:Np,1:Np,iel,1:nVar) /= &
             dSdtFull(1:Np,1:Np,iel,1:nVar))) then
        print*,"Error: ",label," subset tendency differs from the full tendency on the ", &
          meshKind," mesh, element ",iel
        print*,"  max |difference| : ", &
          maxval(abs(modelobj%dSdt%interior(1:Np,1:Np,iel,1:nVar)- &
                     dSdtFull(1:Np,1:Np,iel,1:nVar)))
        ok = .false.
        return
      endif
    enddo

  endfunction Matches

endprogram Subset_Operators_2D
