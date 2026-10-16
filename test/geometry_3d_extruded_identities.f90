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
!! Checks the metric identities of extruded elements (Mesh3D%ExtrudedMesh) in the
!! discrete metric terms SEMHex computes: the covariant basis by differentiating the
!! interpolated coordinates, the Jacobian J as its determinant, and the contravariant
!! basis J a^i in Kopriva's curl-invariant form.
!!
!! For the map x = (X(xi,eta), Y(xi,eta), Z(xi,eta,zeta)) the exact metric terms are
!!
!!   J a^xi   = Z_zeta (Y_eta, -X_eta, 0)
!!   J a^eta  = Z_zeta (-Y_xi, X_xi, 0)
!!   J a^zeta = (Y_xi Z_eta - Z_xi Y_eta, Z_xi X_eta - X_xi Z_eta, J_2)
!!
!! with J_2 = X_xi Y_eta - X_eta Y_xi, and the test measures, at every node,
!!
!!   (G1) the vertical components of J a^xi and J a^eta, and the difference of their
!!        horizontal components from Z_zeta times the 2-D contravariant vectors that
!!        SEMQuad computes on the surface mesh, relative to max |J a^xi|, |J a^eta|;
!!   (G2) J a^zeta_z minus the 2-D Jacobian J_2 from SEMQuad, relative to max |J_2|;
!!   (I0) sum_i J a^i D_i Z - J zhat, relative to max |J|, which is the statement that
!!        the discrete gradient of the coordinate z is the unit vector zhat.
!!
!! Each term the curl form differentiates is a product of two interpolated factors.
!! The products are degree 2 nGeo in the horizontal reference coordinates, so the
!! discrete product rule, and with it all three identities, is exact once
!! N >= 2 nGeo. On an element whose levels are flat (Z depends on zeta only), G1 and
!! the horizontal components of I0 hold for every N, because there the factor Z_zeta
!! is constant over the element. G2 and the vertical component of I0 are a 2-D
!! property: J a^zeta_z is the curl form of J_2, which equals the determinant J_2 only
!! when the 2-D products are resolved (N >= 2 nGeo) or the surface element is affine.
!!
!! The column sweeps rely on three weaker properties, which hold to round-off for every
!! N because each one only involves products that are linear in zeta or constant along
!! a vertical grid line:
!!
!!   (a) the vertical components of J a^xi and J a^eta vanish;
!!   (b) J a^zeta_z is constant along every vertical grid line, inside an element and
!!       from layer to layer of a column (it depends on X and Y only);
!!   (c) the slope components J a^zeta_x and J a^zeta_y are continuous across the faces
!!       between layers (they are linear in zeta inside a layer, and the two layers
!!       evaluate them on the same level surface).
!!
!! The test also reports (d), J - J a^zeta_z Z_zeta, with J the determinant SEMHex
!! stores; it is the vertical part of I0.
!!
!! The test asserts every identity to round-off, 100 N**4 eps relative, where that
!! analysis says it is exact, and prints the measured truncation error where it is not. Cases: a structured
!! (affine, nGeo = 1) surface and the curved HOHQMesh Circle surface (nGeo = 7), each
!! with flat levels, terrain-following sigma levels over a sloping bathymetry, and
!! hybrid levels (one flat layer over terrain-following ones), at N = 2 nGeo and at an
!! N below it.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Mesh_3D
  use SELF_Geometry_2D
  use SELF_Geometry_3D

  implicit none
  integer :: exit_code

  exit_code = geometry_3d_extruded_identities()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  subroutine BuildLevels(surf,levelType,zLevels)
    !! Level heights for one of the three level types: 1 = flat (constant depth 1,
    !! hybrid with one flat interface at depth 0.3), 2 = sigma levels over the sloping
    !! depth h = 1 + 0.3 x + 0.2 y^2, 3 = hybrid levels over the same sloping depth.
    implicit none
    type(Mesh2D),intent(in) :: surf
    integer,intent(in) :: levelType
    real(prec),allocatable,intent(out) :: zLevels(:,:,:,:)
    ! Local
    real(prec),allocatable :: h(:,:,:)
    integer :: nGeo,e,i,j
    real(prec) :: x,y

    nGeo = surf%nGeo
    allocate(h(1:nGeo+1,1:nGeo+1,1:surf%nElem))
    do e = 1,surf%nElem
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          x = surf%nodeCoords(1,i,j,e)
          y = surf%nodeCoords(2,i,j,e)
          if(levelType == 1) then
            h(i,j,e) = 1.0_prec
          else
            h(i,j,e) = 1.0_prec+0.3_prec*x+0.2_prec*y*y
          endif
        enddo
      enddo
    enddo

    if(levelType == 2) then
      allocate(zLevels(1:nGeo+1,1:nGeo+1,0:2,1:surf%nElem))
      call ExtrudedSigmaLevels(h,[-1.0_prec,-0.45_prec,0.0_prec],zLevels)
    else
      allocate(zLevels(1:nGeo+1,1:nGeo+1,0:3,1:surf%nElem))
      call ExtrudedHybridLevels(h,[0.3_prec],[-1.0_prec,-0.4_prec,0.0_prec],zLevels)
    endif
    deallocate(h)

  endsubroutine BuildLevels

  integer function CheckCase(surf,label,levelType,N) result(r)
    implicit none
    type(Mesh2D),intent(in) :: surf
    character(*),intent(in) :: label
    integer,intent(in) :: levelType
    integer,intent(in) :: N
    ! Local
    real(prec) :: tolerance
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geom3
    type(SEMQuad),target :: geom2
    real(prec),allocatable :: zLevels(:,:,:,:)
    integer :: nz,e3,e2,k,i,j,l,kk,idir,nGeo
    integer :: nElems(1:2) ! number of flat and of sloping elements
    logical :: flat,resolved
    real(prec) :: scaleJa,scaleJ2,scaleJ,zzeta,res
    ! Errors, (1) on flat elements and (2) on sloping elements
    real(prec) :: errG1(1:2),errG2(1:2),errI0h(1:2),errI0z(1:2)
    ! Column-sweep properties (a), (b), (c) and the Jacobian consistency (d), absolute
    real(prec) :: errA,errB,errC,errD,scaleS,top,bot
    integer :: ii,c

    r = 0
    nGeo = surf%nGeo
    resolved = N >= 2*nGeo
    ! Round-off threshold. The identities differentiate products of differentiated
    ! coordinates, so round-off grows like N**4 eps (the derivative matrix entries are
    ! O(N**2)); measured, the exact cases sit 3 to 30 times below 100 N**4 eps.
    tolerance = 100.0_prec*real(N,prec)**4*epsilon(1.0_prec)

    call BuildLevels(surf,levelType,zLevels)
    nz = size(zLevels,3)-1
    call mesh%ExtrudedMesh(surf,zLevels,SELF_BC_NONORMALFLOW,SELF_BC_PRESCRIBED)

    call interp%Init(N=N,controlNodeType=GAUSS,M=N,targetNodeType=UNIFORM)
    call geom3%Init(interp,mesh%nElem)
    call geom3%GenerateFromMesh(mesh)
    call geom2%Init(interp,surf%nElem)
    call geom2%GenerateFromMesh(surf)

    scaleJa = maxval(abs(geom3%dsdx%interior(:,:,:,:,1,:,1:2)))
    scaleJ2 = maxval(abs(geom2%J%interior))
    scaleJ = maxval(abs(geom3%J%interior))

    errG1 = 0.0_prec
    errG2 = 0.0_prec
    errI0h = 0.0_prec
    errI0z = 0.0_prec
    nElems = 0

    do e2 = 1,surf%nElem
      do k = 1,nz
        e3 = (e2-1)*nz+k
        ! An element is flat when both of its bounding levels are constant over it
        flat = maxval(zLevels(:,:,k-1,e2)) == minval(zLevels(:,:,k-1,e2)) .and. &
               maxval(zLevels(:,:,k,e2)) == minval(zLevels(:,:,k,e2))
        l = 2
        if(flat) l = 1
        nElems(l) = nElems(l)+1
        do kk = 1,N+1
          do j = 1,N+1
            do i = 1,N+1
              ! dxds(...,1,c,d) = d x_c / d xi_d ; dsdx(...,1,c,d) = J a^d, component c
              zzeta = geom3%dxds%interior(i,j,kk,e3,1,3,3)

              ! (G1)
              do idir = 1,2
                errG1(l) = max(errG1(l),abs(geom3%dsdx%interior(i,j,kk,e3,1,3,idir)))
                errG1(l) = max(errG1(l),abs(geom3%dsdx%interior(i,j,kk,e3,1,1,idir)- &
                                            zzeta*geom2%dsdx%interior(i,j,e2,1,1,idir)))
                errG1(l) = max(errG1(l),abs(geom3%dsdx%interior(i,j,kk,e3,1,2,idir)- &
                                            zzeta*geom2%dsdx%interior(i,j,e2,1,2,idir)))
              enddo

              ! (G2)
              errG2(l) = max(errG2(l),abs(geom3%dsdx%interior(i,j,kk,e3,1,3,3)-geom2%J%interior(i,j,e2,1)))

              ! (I0): sum_i J a^i_c D_i Z - J delta_{c3}
              do idir = 1,3
                res = geom3%dsdx%interior(i,j,kk,e3,1,idir,1)*geom3%dxds%interior(i,j,kk,e3,1,3,1)+ &
                      geom3%dsdx%interior(i,j,kk,e3,1,idir,2)*geom3%dxds%interior(i,j,kk,e3,1,3,2)+ &
                      geom3%dsdx%interior(i,j,kk,e3,1,idir,3)*geom3%dxds%interior(i,j,kk,e3,1,3,3)
                if(idir == 3) then
                  errI0z(l) = max(errI0z(l),abs(res-geom3%J%interior(i,j,kk,e3,1)))
                else
                  errI0h(l) = max(errI0h(l),abs(res))
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    ! Properties the column sweeps rely on, measured on every element for every N:
    !  (a) the vertical components of J a^xi and J a^eta vanish;
    !  (b) J a^zeta_z, as stored, does not vary along a vertical grid line, neither
    !      with the zeta node inside an element nor from layer to layer of a column;
    !  (c) the slope components J a^zeta_x, J a^zeta_y are continuous across the
    !      faces between layers. Both traces are extrapolated here from the interior
    !      with the interpolant's boundary matrix, before SEMHex flips the sign of the
    !      stored bottom-face trace;
    !  (d) J, which SEMHex stores as the determinant of the covariant basis, against
    !      J a^zeta_z Z_zeta (reported only; it is the vertical part of I0).
    errA = 0.0_prec
    errB = 0.0_prec
    errC = 0.0_prec
    errD = 0.0_prec
    scaleS = maxval(abs(geom3%dsdx%interior(:,:,:,:,1,1:2,3)))
    do e2 = 1,surf%nElem
      do k = 1,nz
        e3 = (e2-1)*nz+k
        do kk = 1,N+1
          do j = 1,N+1
            do i = 1,N+1
              errA = max(errA,abs(geom3%dsdx%interior(i,j,kk,e3,1,3,1)), &
                         abs(geom3%dsdx%interior(i,j,kk,e3,1,3,2)))
              errB = max(errB,abs(geom3%dsdx%interior(i,j,kk,e3,1,3,3)- &
                                  geom3%dsdx%interior(i,j,1,(e2-1)*nz+1,1,3,3)))
              errD = max(errD,abs(geom3%J%interior(i,j,kk,e3,1)- &
                                  geom3%dsdx%interior(i,j,kk,e3,1,3,3)*geom3%dxds%interior(i,j,kk,e3,1,3,3)))
            enddo
          enddo
        enddo
        if(k < nz) then
          do c = 1,2
            do j = 1,N+1
              do i = 1,N+1
                top = 0.0_prec
                bot = 0.0_prec
                do ii = 1,N+1
                  top = top+interp%bMatrix(ii,2)*geom3%dsdx%interior(i,j,ii,e3,1,c,3)
                  bot = bot+interp%bMatrix(ii,1)*geom3%dsdx%interior(i,j,ii,e3+1,1,c,3)
                enddo
                errC = max(errC,abs(top-bot))
              enddo
            enddo
          enddo
        endif
      enddo
    enddo

    errG1 = errG1/scaleJa
    errG2 = errG2/scaleJ2
    errI0h = errI0h/scaleJ
    errI0z = errI0z/scaleJ

    print'(A,A,A,I0,A,I0,A,L1,A,I0,A,I0,A,ES10.3)',"Case ",label,": N = ",N,", nGeo = ",nGeo,", N >= 2 nGeo : ",resolved, &
      ", flat elements ",nElems(1),", sloping elements ",nElems(2),", round-off threshold ",tolerance
    print'(A,ES12.4,A,ES12.4)',"    (a) |Ja^xi_z|,|Ja^eta_z|      max abs ",errA,", relative ",errA/scaleJa
    print'(A,ES12.4,A,ES12.4)',"    (b) Ja^zeta_z along a column   max abs ",errB,", relative ",errB/scaleJ2
    ! (c) is relative to the slope magnitude, or to J_2 on flat levels where there is none
    scaleS = max(scaleS,scaleJ2)
    print'(A,ES12.4,A,ES12.4)',"    (c) Ja^zeta_x,y jump at levels max abs ",errC,", relative ",errC/scaleS
    print'(A,ES12.4,A,ES12.4)',"    (d) J - Ja^zeta_z Z_zeta       max abs ",errD,", relative ",errD/scaleJ
    call Report(r,tolerance,"(a)             ",errA/scaleJa,.true.)
    call Report(r,tolerance,"(b)             ",errB/scaleJ2,.true.)
    call Report(r,tolerance,"(c)             ",errC/scaleS,.true.)
    call Report(r,tolerance,"(d)             ",errD/scaleJ,resolved)
    if(nElems(1) > 0) then
      call Report(r,tolerance,"flat    G1      ",errG1(1),.true.)
      call Report(r,tolerance,"flat    I0 x,y  ",errI0h(1),.true.)
      call Report(r,tolerance,"flat    G2      ",errG2(1),resolved)
      call Report(r,tolerance,"flat    I0 z    ",errI0z(1),resolved)
    endif
    if(nElems(2) > 0) then
      call Report(r,tolerance,"sloping G1      ",errG1(2),resolved)
      call Report(r,tolerance,"sloping I0 x,y  ",errI0h(2),resolved)
      call Report(r,tolerance,"sloping G2      ",errG2(2),resolved)
      call Report(r,tolerance,"sloping I0 z    ",errI0z(2),resolved)
    endif

    call geom3%Free()
    call geom2%Free()
    call interp%Free()
    call mesh%Free()
    deallocate(zLevels)

  endfunction CheckCase

  subroutine Report(r,tolerance,name,err,asserted)
    !! Prints one measured identity error; when the identity is expected to hold to
    !! round-off, sets r = 1 if it does not.
    implicit none
    integer,intent(inout) :: r
    real(prec),intent(in) :: tolerance
    character(*),intent(in) :: name
    real(prec),intent(in) :: err
    logical,intent(in) :: asserted
    ! Local
    character(len=8) :: verdict

    if(asserted) then
      if(err <= tolerance) then
        verdict = "exact"
      else
        verdict = "FAIL"
        r = 1
      endif
    else
      verdict = "reported"
    endif
    print'(A,A,ES12.4,A,A)',"    ",name,err,"  ",trim(verdict)

  endsubroutine Report

  integer function geometry_3d_extruded_identities() result(r)
    implicit none
    ! Local
    type(Mesh2D) :: structured,circle
    character(LEN=255) :: WORKSPACE
    integer :: levelType
    character(len=16),parameter :: levelName(1:3) = ["flat levels     ", &
                                                     "sigma levels    ", &
                                                     "hybrid levels   "]

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call structured%StructuredMesh(3,2,2,2,0.3_prec,0.25_prec,[1,1,1,1])
    call circle%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Circle/Circle.mesh")

    r = 0
    do levelType = 1,3
      r = max(r,CheckCase(structured,"structured, "//trim(levelName(levelType)),levelType,2))
      r = max(r,CheckCase(structured,"structured, "//trim(levelName(levelType)),levelType,5))
      r = max(r,CheckCase(circle,"Circle, "//trim(levelName(levelType)),levelType,7))
      r = max(r,CheckCase(circle,"Circle, "//trim(levelName(levelType)),levelType,14))
    enddo

    call structured%Free()
    call circle%Free()

  endfunction geometry_3d_extruded_identities

endprogram test
