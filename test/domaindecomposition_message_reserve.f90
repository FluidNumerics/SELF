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
!! Checks the three routines an exchange uses to size its request scratch:
!! CountRemoteSides, CountRemoteMortarSides and ReserveMessages.
!!
!! Every caller of these lives behind `if(mesh%decomp%mpiEnabled)`, so on one
!! rank nothing reaches them and a wrong count is invisible -- and an
!! undercount is exactly what overruns the scratch, silently, on the runs that
!! do reach them. This test calls them directly instead.
!!
!! The partition is simulated by writing decomp%elemToRank, which is what the
!! counters and the exchanges both read to decide whether a neighbor is
!! rank-local. Marking one interior element as owned by another rank gives a
!! count that can be written down by hand: every side of that element's
!! neighbors that faces it, and nothing else.

  use SELF_Constants
  use SELF_Mesh_2D
  use SELF_Mesh_3D

  implicit none
  integer :: exit_code

  exit_code = domaindecomposition_message_reserve()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function CheckCount(got,expected,label) result(r)
    !! Reports one count against the value worked out by hand.
    implicit none
    integer,intent(in) :: got
    integer,intent(in) :: expected
    character(*),intent(in) :: label

    if(got /= expected) then
      print*,"FAIL ("//label//"): expected ",expected," remote sides, got ",got
      r = 1
    else
      print*,"PASS ("//label//"): ",got," remote sides"
      r = 0
    endif

  endfunction CheckCount

  integer function domaindecomposition_message_reserve() result(r)
    !! The meshes are held simultaneously and freed at the end: SELF owns the MPI
    !! lifecycle here, and freeing the last live mesh finalizes MPI, which would
    !! leave any mesh created afterwards without one.
    implicit none
    ! The centre element of a 3x3 mesh (i + 3*(j-1), i = j = 2) and of a 3x3x3
    ! mesh (i + 3*(j-1) + 9*(k-1), i = j = k = 2). Both are interior, so every
    ! side has a neighbor.
    integer,parameter :: centre2D = 5
    integer,parameter :: centre3D = 14
    type(Mesh2D),target :: mesh2d,mortar2d
    type(Mesh3D),target :: mesh3d,mortar3d
    integer :: bcids2d(1:4),bcids3d(1:6)
    integer :: ri,maxMsg0

    bcids2d(1:4) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]
    bcids3d(1:6) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED, &
                    SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    r = 0

    call mesh2d%StructuredMesh(3,3,1,1,0.1_prec,0.1_prec,bcids2d)
    call mesh3d%StructuredMesh(3,3,3,1,1,1,0.1_prec,0.1_prec,0.1_prec,bcids3d)
    call mortar2d%SimpleMortarMesh(0.1_prec,bcids2d)
    call mortar3d%SimpleMortarMesh(0.1_prec,bcids3d)

    ! Nothing is remote on one rank, whatever the mesh.
    ri = CheckCount(mesh2d%decomp%CountRemoteSides(mesh2d%sideInfo,mesh2d%nElem,4), &
                    0,"2-D structured, single rank")
    r = max(r,ri)
    ri = CheckCount(mesh3d%decomp%CountRemoteSides(mesh3d%sideInfo,mesh3d%nElem,6), &
                    0,"3-D structured, single rank")
    r = max(r,ri)

    ! Hand the centre element to another rank. Its four (2-D) or six (3-D)
    ! neighbors each then have exactly one side facing a remote element, and the
    ! centre element's own sides face rank-local neighbors, so they do not count.
    mesh2d%decomp%elemToRank(centre2D) = mesh2d%decomp%rankId+1
    ri = CheckCount(mesh2d%decomp%CountRemoteSides(mesh2d%sideInfo,mesh2d%nElem,4), &
                    4,"2-D structured, centre element remote")
    r = max(r,ri)

    mesh3d%decomp%elemToRank(centre3D) = mesh3d%decomp%rankId+1
    ri = CheckCount(mesh3d%decomp%CountRemoteSides(mesh3d%sideInfo,mesh3d%nElem,6), &
                    6,"3-D structured, centre element remote")
    r = max(r,ri)

    ! A mesh with no mortars never allocates the mortar table, so the count has
    ! to read the null pointer as "nothing to exchange" rather than dereference it.
    ri = CheckCount(mesh2d%decomp%CountRemoteMortarSides(mesh2d%mortarInfo,mesh2d%nMortars,2), &
                    0,"2-D structured, no mortar table")
    r = max(r,ri)

    ! A mortar entirely on one rank is served from memory.
    ri = CheckCount(mortar2d%decomp%CountRemoteMortarSides(mortar2d%mortarInfo, &
                                                           mortar2d%nMortars,2), &
                    0,"2-D mortar, single rank")
    r = max(r,ri)
    ri = CheckCount(mortar3d%decomp%CountRemoteMortarSides(mortar3d%mortarInfo, &
                                                           mortar3d%nMortars,4), &
                    0,"3-D mortar, single rank")
    r = max(r,ri)

    ! Hand the big element to another rank: every sub-side of the interface then
    ! crosses the rank boundary. SimpleMortarMesh builds one mortar, with two
    ! sub-edges in 2-D and four sub-faces in 3-D.
    mortar2d%decomp%elemToRank(mortar2d%mortarInfo(1,1)) = mortar2d%decomp%rankId+1
    ri = CheckCount(mortar2d%decomp%CountRemoteMortarSides(mortar2d%mortarInfo, &
                                                           mortar2d%nMortars,2), &
                    2,"2-D mortar, big element remote")
    r = max(r,ri)

    mortar3d%decomp%elemToRank(mortar3d%mortarInfo(1,1)) = mortar3d%decomp%rankId+1
    ri = CheckCount(mortar3d%decomp%CountRemoteMortarSides(mortar3d%mortarInfo, &
                                                           mortar3d%nMortars,4), &
                    4,"3-D mortar, big element remote")
    r = max(r,ri)

    ! ReserveMessages only ever grows: an exchange that needs less than the last
    ! one must not hand back memory another exchange is still sized for.
    maxMsg0 = mesh2d%decomp%maxMsg
    call mesh2d%decomp%ReserveMessages(maxMsg0)
    if(mesh2d%decomp%maxMsg /= maxMsg0 .or. size(mesh2d%decomp%requests) /= maxMsg0) then
      print*,"FAIL: reserving what is already there resized the scratch to ", &
        mesh2d%decomp%maxMsg,size(mesh2d%decomp%requests)
      r = 1
    endif

    call mesh2d%decomp%ReserveMessages(maxMsg0+7)
    if(mesh2d%decomp%maxMsg /= maxMsg0+7) then
      print*,"FAIL: maxMsg is ",mesh2d%decomp%maxMsg," after reserving ",maxMsg0+7
      r = 1
    endif
    if(size(mesh2d%decomp%requests) /= maxMsg0+7) then
      print*,"FAIL: the request scratch holds ",size(mesh2d%decomp%requests), &
        " after reserving ",maxMsg0+7
      r = 1
    endif
    if(size(mesh2d%decomp%stats,2) /= maxMsg0+7) then
      print*,"FAIL: the status scratch holds ",size(mesh2d%decomp%stats,2), &
        " after reserving ",maxMsg0+7
      r = 1
    endif

    call mesh2d%decomp%ReserveMessages(1)
    if(mesh2d%decomp%maxMsg /= maxMsg0+7) then
      print*,"FAIL: a smaller reservation shrank the scratch to ",mesh2d%decomp%maxMsg
      r = 1
    endif

    if(r == 0) print*,"PASS: ReserveMessages grows to the requested size and never shrinks"

    call mesh2d%Free()
    call mesh3d%Free()
    call mortar2d%Free()
    call mortar3d%Free()

  endfunction domaindecomposition_message_reserve

endprogram test
