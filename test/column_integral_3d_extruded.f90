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
!! Checks ColumnIntegral (SELF_ColumnOperators_3D) on extruded meshes,
!!
!!   (C q)_ij = sum_{k=1..nz} sum_{l=0..N} w_l (Z_zeta)^(k)_ijl q^(k)_ijl ,
!!
!! against exact column integrals.
!!
!!   1. Flat levels over a structured surface, total depth H = 2.5 in three uneven
!!      layers. Z_zeta is constant in every layer, so the quadrature integrates z^p
!!      exactly up to p = 2N-1 on Gauss-Lobatto nodes and p = 2N+1 on Gauss nodes. The
!!      test integrates z^p for p = 0, 1, N and that maximum degree and compares with
!!      -(-H)^(p+1)/(p+1), relative error asserted to round-off.
!!   2. Terrain-following levels (sigma levels, and hybrid levels with one flat layer
!!      on top) over a sloping bathymetry. Inside a layer z is linear in zeta at every
!!      solution node and Z_zeta is constant in zeta, so on every mesh and for every N
!!
!!        C 1 = z_nz - z_0 = h_d ,   C z = (z_nz^2 - z_0^2)/2 = -h_d^2/2 ,
!!
!!      where h_d = -z_0 is the discrete depth at the solution node: the bathymetry the
!!      mesh carries, i.e. its degree nGeo interpolant evaluated at the solution node.
!!      Both are asserted to round-off. Against the analytic depth h(x,y) at the node
!!      they hold only where the interpolant reproduces h: the structured surface
!!      (nGeo = 1) with the linear bathymetry h = 1 + 0.3 x + 0.2 y, which is asserted.
!!      On the curved surfaces the bathymetry h = 1 + 0.3 x + 0.2 y^2 composed with the
!!      curved map is not a polynomial of degree nGeo, and the difference is reported.
!!
!! Surfaces: StructuredMesh and the HOPr Circle (both decomposed across ranks), and,
!! on one rank only, the HOHQMesh Circle and Seamount surfaces (the HOHQMesh reader
!! does not decompose the mesh). The test runs on one and on two ranks; the column
!! integral is rank local, and a checksum of the results keyed by global node is
!! printed so that the two runs can be compared bitwise.
!!
!! Tolerance: 100 N^2 eps relative. Z_zeta is a derivative of the interpolated
!! coordinate, which carries O(N^2 eps) round-off; the sum adds nz (N+1) terms.

  use iso_fortran_env,only:int64
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Mesh_3D
  use SELF_Geometry_3D
  use SELF_Scalar_2D
  use SELF_MappedScalar_3D
  use SELF_ColumnOperators_3D
  use mpi

  implicit none
  integer :: exit_code

  exit_code = column_integral_3d_extruded()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  real(prec) function Bathymetry(x,y,linear) result(h)
    !! Water depth in metres: 1 + 0.3 x + 0.2 y (linear) or 1 + 0.3 x + 0.2 y^2
    implicit none
    real(prec),intent(in) :: x,y
    logical,intent(in) :: linear

    if(linear) then
      h = 1.0_prec+0.3_prec*x+0.2_prec*y
    else
      h = 1.0_prec+0.3_prec*x+0.2_prec*y*y
    endif

  endfunction Bathymetry

  integer(int64) function NodeBits(v) result(b)
    !! The bit pattern of a real, for the checksum
    implicit none
    real(prec),intent(in) :: v

#ifdef DOUBLE_PRECISION
    b = transfer(v,0_int64)
#else
    b = int(transfer(v,0),int64)
#endif

  endfunction NodeBits

  subroutine Checksum(surf,qbar,label)
    !! Prints an order independent checksum (bitwise XOR across ranks) of qbar, each
    !! value rotated by an amount that depends on its global node and variable index
    implicit none
    type(Mesh2D),intent(in) :: surf
    type(Scalar2D),intent(in) :: qbar
    character(*),intent(in) :: label
    ! Local
    integer(int64) :: h,hGlobal,key
    integer :: i,j,e2,ivar,np,ierror

    np = qbar%N+1
    h = 0_int64
    do ivar = 1,qbar%nVar
      do e2 = 1,qbar%nElem
        do j = 1,np
          do i = 1,np
            key = ((int(surf%decomp%offsetElem(surf%decomp%rankId+1)+e2-1,int64)*np+(j-1))*np+(i-1))*qbar%nVar+ivar
            h = ieor(h,ishftc(NodeBits(qbar%interior(i,j,e2,ivar)),int(mod(key,61_int64))+1))
          enddo
        enddo
      enddo
    enddo
    call MPI_Allreduce(h,hGlobal,1,MPI_INTEGER8,MPI_BXOR,surf%decomp%mpiComm,ierror)
    if(surf%decomp%rankId == 0) print'(A,A,A,Z16.16)',"    checksum (",label,") ",hGlobal

  endsubroutine Checksum

  subroutine Report(r,tolerance,name,err,asserted,surf)
    !! Reduces err over the ranks, prints it, and sets r = 1 if it is asserted and
    !! exceeds the tolerance
    implicit none
    integer,intent(inout) :: r
    real(prec),intent(in) :: tolerance
    character(*),intent(in) :: name
    real(prec),intent(in) :: err
    logical,intent(in) :: asserted
    type(Mesh2D),intent(in) :: surf
    ! Local
    real(prec) :: errGlobal
    integer :: ierror,rankId
    character(len=8) :: verdict

    call MPI_Allreduce(err,errGlobal,1,surf%decomp%mpiPrec,MPI_MAX,surf%decomp%mpiComm,ierror)
    rankId = surf%decomp%rankId
    if(asserted) then
      if(errGlobal <= tolerance) then
        verdict = "exact"
      else
        verdict = "FAIL"
        r = 1
      endif
    else
      verdict = "reported"
    endif
    if(rankId == 0) print'(A,A,ES12.4,A,ES10.3,A,A)',"    ",name,errGlobal,"  (tolerance ",tolerance,")  ",trim(verdict)

  endsubroutine Report

  integer function CheckFlat(surf,label,N,nodeType) result(r)
    !! Case 1: flat levels, integrals of z^p
    implicit none
    type(Mesh2D),intent(in) :: surf
    character(*),intent(in) :: label
    integer,intent(in) :: N
    integer,intent(in) :: nodeType
    ! Local
    integer,parameter :: nvar = 4
    real(prec),parameter :: depth = 2.5_prec
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: q
    type(Scalar2D) :: qbar
    type(ColumnOperators3D) :: cols
    real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
    integer :: p(1:nvar),nz,i,j,l,iel,e2,ivar
    real(prec) :: z,exact,err,tolerance

    r = 0
    tolerance = 100.0_prec*real(N,prec)**2*epsilon(1.0_prec)
    p(1:3) = [0,1,N]
    if(nodeType == GAUSS_LOBATTO) then
      p(4) = 2*N-1
    else
      p(4) = 2*N+1
    endif

    allocate(h(1:surf%nGeo+1,1:surf%nGeo+1,1:surf%nElem))
    allocate(zLevels(1:surf%nGeo+1,1:surf%nGeo+1,0:3,1:surf%nElem))
    h = depth
    call ExtrudedSigmaLevels(h,[-1.0_prec,-0.6_prec,-0.25_prec,0.0_prec],zLevels)
    nz = 3
    call mesh%ExtrudedMesh(surf,zLevels,SELF_BC_NONORMALFLOW,SELF_BC_PRESCRIBED)

    call interp%Init(N=N,controlNodeType=nodeType,M=N,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)
    call cols%Init(interp,mesh,geometry)

    call q%Init(interp,nvar,mesh%nElem)
    call qbar%Init(interp,nvar,surf%nElem)
    do ivar = 1,nvar
      do iel = 1,mesh%nElem
        do l = 1,N+1
          do j = 1,N+1
            do i = 1,N+1
              z = geometry%x%interior(i,j,l,iel,1,3)
              q%interior(i,j,l,iel,ivar) = z**p(ivar)
            enddo
          enddo
        enddo
      enddo
    enddo
    call q%UpdateDevice()

    call cols%ColumnIntegral(q,qbar)
    call qbar%UpdateHost()

    if(surf%decomp%rankId == 0) print'(A,A,A,I0,A,I0,A,L1)',"Case ",label,", flat levels: N = ",N, &
      ", nz = ",nz,", Gauss-Lobatto ",nodeType == GAUSS_LOBATTO
    do ivar = 1,nvar
      exact = -(-depth)**(p(ivar)+1)/real(p(ivar)+1,prec)
      err = 0.0_prec
      do e2 = 1,surf%nElem
        do j = 1,N+1
          do i = 1,N+1
            err = max(err,abs(qbar%interior(i,j,e2,ivar)-exact)/abs(exact))
          enddo
        enddo
      enddo
      call Report(r,tolerance,"C z^"//trim(Str(p(ivar)))//" vs analytic       ",err,.true.,surf)
    enddo
    call Checksum(surf,qbar,label//", flat")

    call cols%Free()
    call qbar%Free()
    call q%Free()
    call geometry%Free()
    call interp%Free()
    call mesh%Free()
    deallocate(h,zLevels)

  endfunction CheckFlat

  integer function CheckSloping(surf,label,N,hybrid,linear) result(r)
    !! Case 2: terrain-following levels, C 1 = h and C z = -h^2/2
    implicit none
    type(Mesh2D),intent(in) :: surf
    character(*),intent(in) :: label
    integer,intent(in) :: N
    logical,intent(in) :: hybrid
    !! hybrid levels (one flat layer of depth 0.3 over two sigma layers), else 3 sigma layers
    logical,intent(in) :: linear
    !! linear bathymetry, which a structured (nGeo = 1) surface represents exactly
    ! Local
    integer,parameter :: nvar = 2
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: q
    type(Scalar2D) :: qbar
    type(ColumnOperators3D) :: cols
    real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
    integer :: nGeo,nz,i,j,l,iel,e2,ierror
    real(prec) :: hd,ha,x,y,tolerance,scale
    real(prec) :: err1d,err2d,err1a,err2a

    r = 0
    tolerance = 100.0_prec*real(N,prec)**2*epsilon(1.0_prec)
    nGeo = surf%nGeo
    allocate(h(1:nGeo+1,1:nGeo+1,1:surf%nElem))
    do e2 = 1,surf%nElem
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          h(i,j,e2) = Bathymetry(surf%nodeCoords(1,i,j,e2),surf%nodeCoords(2,i,j,e2),linear)
        enddo
      enddo
    enddo
    if(hybrid) then
      allocate(zLevels(1:nGeo+1,1:nGeo+1,0:3,1:surf%nElem))
      call ExtrudedHybridLevels(h,[0.3_prec],[-1.0_prec,-0.4_prec,0.0_prec],zLevels)
    else
      allocate(zLevels(1:nGeo+1,1:nGeo+1,0:3,1:surf%nElem))
      call ExtrudedSigmaLevels(h,[-1.0_prec,-0.55_prec,-0.2_prec,0.0_prec],zLevels)
    endif
    nz = 3
    call mesh%ExtrudedMesh(surf,zLevels,SELF_BC_NONORMALFLOW,SELF_BC_PRESCRIBED)

    call interp%Init(N=N,controlNodeType=GAUSS_LOBATTO,M=N,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)
    call cols%Init(interp,mesh,geometry)

    call q%Init(interp,nvar,mesh%nElem)
    call qbar%Init(interp,nvar,surf%nElem)
    do iel = 1,mesh%nElem
      do l = 1,N+1
        do j = 1,N+1
          do i = 1,N+1
            q%interior(i,j,l,iel,1) = 1.0_prec
            q%interior(i,j,l,iel,2) = geometry%x%interior(i,j,l,iel,1,3)
          enddo
        enddo
      enddo
    enddo
    call q%UpdateDevice()

    call cols%ColumnIntegral(q,qbar)
    call qbar%UpdateHost()

    ! The discrete depth is minus the height of the bottom node of the bottom layer
    err1d = 0.0_prec
    err2d = 0.0_prec
    err1a = 0.0_prec
    err2a = 0.0_prec
    scale = 0.0_prec
    do e2 = 1,surf%nElem
      iel = (e2-1)*nz+1
      do j = 1,N+1
        do i = 1,N+1
          hd = -geometry%x%interior(i,j,1,iel,1,3)
          x = geometry%x%interior(i,j,1,iel,1,1)
          y = geometry%x%interior(i,j,1,iel,1,2)
          ha = Bathymetry(x,y,linear)
          scale = max(scale,hd)
          err1d = max(err1d,abs(qbar%interior(i,j,e2,1)-hd))
          err2d = max(err2d,abs(qbar%interior(i,j,e2,2)+0.5_prec*hd*hd))
          err1a = max(err1a,abs(qbar%interior(i,j,e2,1)-ha))
          err2a = max(err2a,abs(qbar%interior(i,j,e2,2)+0.5_prec*ha*ha))
        enddo
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,scale,1,surf%decomp%mpiPrec,MPI_MAX,surf%decomp%mpiComm,ierror)

    if(surf%decomp%rankId == 0) then
      if(hybrid) then
        print'(A,A,A,I0,A,I0,A,I0)',"Case ",label,", hybrid levels: N = ",N,", nGeo = ",nGeo,", nz = ",nz
      else
        print'(A,A,A,I0,A,I0,A,I0)',"Case ",label,", sigma levels: N = ",N,", nGeo = ",nGeo,", nz = ",nz
      endif
    endif
    call Report(r,tolerance,"C 1 - h_d,      relative        ",err1d/scale,.true.,surf)
    call Report(r,tolerance,"C z + h_d^2/2,  relative        ",err2d/(0.5_prec*scale*scale),.true.,surf)
    call Report(r,tolerance,"C 1 - h(x,y),   relative        ",err1a/scale,linear,surf)
    call Report(r,tolerance,"C z + h(x,y)^2/2, relative      ",err2a/(0.5_prec*scale*scale),linear,surf)
    if(hybrid) then
      call Checksum(surf,qbar,label//", hybrid")
    else
      call Checksum(surf,qbar,label//", sigma")
    endif

    call cols%Free()
    call qbar%Free()
    call q%Free()
    call geometry%Free()
    call interp%Free()
    call mesh%Free()
    deallocate(h,zLevels)

  endfunction CheckSloping

  function Str(i) result(s)
    implicit none
    integer,intent(in) :: i
    character(len=8) :: s

    write(s,'(I0)') i

  endfunction Str

  integer function column_integral_3d_extruded() result(r)
    implicit none
    ! Local
    type(Mesh2D) :: structured,circleHOPr,circle,seamount
    character(LEN=255) :: WORKSPACE
    logical :: serial

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call structured%StructuredMesh(3,2,2,2,0.3_prec,0.25_prec, &
                                   [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED])
    call circleHOPr%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Circle/Circle_mesh.h5")
    serial = structured%decomp%nRanks == 1

    r = 0
    r = max(r,CheckFlat(structured,"structured",4,GAUSS_LOBATTO))
    r = max(r,CheckFlat(structured,"structured",7,GAUSS_LOBATTO))
    r = max(r,CheckFlat(structured,"structured",7,GAUSS))
    r = max(r,CheckFlat(circleHOPr,"HOPr Circle",7,GAUSS_LOBATTO))
    r = max(r,CheckSloping(structured,"structured",4,.false.,.true.))
    r = max(r,CheckSloping(structured,"structured",7,.true.,.true.))
    r = max(r,CheckSloping(structured,"structured",7,.false.,.false.))
    r = max(r,CheckSloping(circleHOPr,"HOPr Circle",4,.false.,.false.))
    r = max(r,CheckSloping(circleHOPr,"HOPr Circle",7,.true.,.false.))

    if(serial) then
      ! The HOHQMesh reader holds every element on every rank, so these run on one rank only
      call circle%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Circle/Circle.mesh")
      call seamount%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Seamount/Seamount2D.mesh")
      r = max(r,CheckSloping(circle,"HOHQMesh Circle",7,.false.,.false.))
      r = max(r,CheckSloping(circle,"HOHQMesh Circle",14,.true.,.false.))
      r = max(r,CheckSloping(seamount,"HOHQMesh Seamount",5,.false.,.false.))
      r = max(r,CheckSloping(seamount,"HOHQMesh Seamount",8,.true.,.false.))
      call circle%Free()
      call seamount%Free()
    endif

    call circleHOPr%Free()
    call structured%Free()

  endfunction column_integral_3d_extruded

endprogram test
