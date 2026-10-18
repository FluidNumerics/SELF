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
!! Checks the identity (I1) of SELF-HPE, that Broadcast (B) is the adjoint of
!! ColumnIntegral (C) on extruded meshes,
!!
!!   <B A, q>_3 = <A, C q>_2 ,
!!
!!   <a,c>_3 = sum_e3 sum_ijl J3 w_i w_j w_l a c ,   <A,C>_2 = sum_e2 sum_ij J2 w_i w_j A C ,
!!
!! for random surface fields A and volume fields q. Expanding both sides, the identity
!! holds term by term when J3 = J2 Z_zeta at every node. It is measured with two pairs
!! of Jacobians:
!!
!!   curl form    J2 = J2c = (J a^zeta)_z read from the top node of the top layer of the
!!                column, and J3 = Jc = (J a^zeta)_z Z_zeta at each volume node. These are
!!                the Jacobians SELF-HPE uses. (J a^zeta)_z is constant along a vertical
!!                grid line to round-off (geometry_3d_extruded_identities, property b), so
!!                the identity holds to round-off on every extruded mesh and for every N.
!!   determinant  J2 = the determinant Jacobian SEMQuad stores for the surface and J3 = the
!!                determinant Jacobian SEMHex stores for the volume. Since the covariant
!!                basis of an extruded element has dX/dzeta = dY/dzeta = 0 (to round-off),
!!                det J3 = Z_zeta det J2 exactly, so this pair is expected to be exact too.
!!                The test measures J3 - Z_zeta J2 node by node and asserts both.
!!
!! The error is |<B A, q>_3 - <A, C q>_2| relative to sum |J3 w_i w_j w_l A q|, reduced over
!! the ranks. Cases: structured and curved (HOPr Circle, and on one rank the HOHQMesh
!! Circle) surfaces, sigma and hybrid (sloping) levels, N = 4 and 7, Gauss-Lobatto nodes,
!! and Gauss nodes at N = 7. Tolerance 100 N^4 eps, as in geometry_3d_extruded_identities.

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh_2D
  use SELF_Mesh_3D
  use SELF_Geometry_2D
  use SELF_Geometry_3D
  use SELF_MappedScalar_2D
  use SELF_MappedScalar_3D
  use SELF_ColumnOperators_3D
  use mpi

  implicit none
  integer :: exit_code

  exit_code = column_broadcast_adjoint_3d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  subroutine Report(r,tolerance,name,err,surf)
    !! Prints an error already reduced over the ranks; sets r = 1 if it exceeds the tolerance
    implicit none
    integer,intent(inout) :: r
    real(prec),intent(in) :: tolerance
    character(*),intent(in) :: name
    real(prec),intent(in) :: err
    type(Mesh2D),intent(in) :: surf
    ! Local
    character(len=8) :: verdict

    if(err <= tolerance) then
      verdict = "exact"
    else
      verdict = "FAIL"
      r = 1
    endif
    if(surf%decomp%rankId == 0) print'(A,A,ES12.4,A,ES10.3,A,A)',"    ",name,err,"  (tolerance ",tolerance,")  ",trim(verdict)

  endsubroutine Report

  integer function CheckCase(surf,label,N,nodeType,hybrid) result(r)
    implicit none
    type(Mesh2D),intent(in) :: surf
    character(*),intent(in) :: label
    integer,intent(in) :: N
    integer,intent(in) :: nodeType
    logical,intent(in) :: hybrid
    ! Local
    integer,parameter :: nvar = 3
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geom3
    type(SEMQuad),target :: geom2
    type(MappedScalar3D) :: q,bA
    type(MappedScalar2D) :: A,cq
    type(ColumnOperators3D) :: cols
    real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
    integer :: nGeo,nz,i,j,l,iel,e2,ivar,ierror
    real(prec) :: tolerance,w3,w2,j3c,j2c,zz
    ! sums(1:3) curl form <BA,q>, <A,Cq>, scale; sums(4:6) the same with the determinant
    real(prec) :: sums(1:6),errJ,scaleJ

    r = 0
    tolerance = 100.0_prec*real(N,prec)**4*epsilon(1.0_prec)
    nGeo = surf%nGeo
    allocate(h(1:nGeo+1,1:nGeo+1,1:surf%nElem))
    do e2 = 1,surf%nElem
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          h(i,j,e2) = 1.0_prec+0.3_prec*surf%nodeCoords(1,i,j,e2)+0.2_prec*surf%nodeCoords(2,i,j,e2)**2
        enddo
      enddo
    enddo
    allocate(zLevels(1:nGeo+1,1:nGeo+1,0:3,1:surf%nElem))
    if(hybrid) then
      call ExtrudedHybridLevels(h,[0.3_prec],[-1.0_prec,-0.4_prec,0.0_prec],zLevels)
    else
      call ExtrudedSigmaLevels(h,[-1.0_prec,-0.55_prec,-0.2_prec,0.0_prec],zLevels)
    endif
    nz = 3
    call mesh%ExtrudedMesh(surf,zLevels,SELF_BC_NONORMALFLOW,SELF_BC_PRESCRIBED)

    call interp%Init(N=N,controlNodeType=nodeType,M=N,targetNodeType=UNIFORM)
    call geom3%Init(interp,mesh%nElem)
    call geom3%GenerateFromMesh(mesh)
    call geom2%Init(interp,surf%nElem)
    call geom2%GenerateFromMesh(surf)
    call cols%Init(interp,mesh,geom3)

    call q%Init(interp,nvar,mesh%nElem)
    call bA%Init(interp,nvar,mesh%nElem)
    call A%Init(interp,nvar,surf%nElem)
    call cq%Init(interp,nvar,surf%nElem)
    call q%AssociateGeometry(geom3)
    call A%AssociateGeometry(geom2)
    call random_number(q%interior)
    call random_number(A%interior)
    q%interior = q%interior-0.5_prec
    A%interior = A%interior-0.5_prec
    call q%UpdateDevice()
    call A%UpdateDevice()

    call cols%Broadcast(A,bA)
    call cols%ColumnIntegral(q,cq)
    call bA%UpdateHost()
    call cq%UpdateHost()

    sums = 0.0_prec
    errJ = 0.0_prec
    scaleJ = 0.0_prec
    do ivar = 1,nvar
      do iel = 1,mesh%nElem
        e2 = (iel-1)/nz+1
        do l = 1,N+1
          do j = 1,N+1
            do i = 1,N+1
              w3 = interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(l)
              zz = geom3%dxds%interior(i,j,l,iel,1,3,3)
              j3c = geom3%dsdx%interior(i,j,l,iel,1,3,3)*zz
              sums(1) = sums(1)+j3c*w3*bA%interior(i,j,l,iel,ivar)*q%interior(i,j,l,iel,ivar)
              sums(3) = sums(3)+abs(j3c*w3*bA%interior(i,j,l,iel,ivar)*q%interior(i,j,l,iel,ivar))
              sums(4) = sums(4)+geom3%J%interior(i,j,l,iel,1)*w3*bA%interior(i,j,l,iel,ivar)*q%interior(i,j,l,iel,ivar)
              sums(6) = sums(6)+abs(geom3%J%interior(i,j,l,iel,1)*w3*bA%interior(i,j,l,iel,ivar)*q%interior(i,j,l,iel,ivar))
              if(ivar == 1) then
                errJ = max(errJ,abs(geom3%J%interior(i,j,l,iel,1)-zz*geom2%J%interior(i,j,e2,1)))
                scaleJ = max(scaleJ,abs(geom3%J%interior(i,j,l,iel,1)))
              endif
            enddo
          enddo
        enddo
      enddo
      do e2 = 1,surf%nElem
        do j = 1,N+1
          do i = 1,N+1
            w2 = interp%qWeights(i)*interp%qWeights(j)
            j2c = geom3%dsdx%interior(i,j,N+1,e2*nz,1,3,3)
            sums(2) = sums(2)+j2c*w2*A%interior(i,j,e2,ivar)*cq%interior(i,j,e2,ivar)
            sums(5) = sums(5)+geom2%J%interior(i,j,e2,1)*w2*A%interior(i,j,e2,ivar)*cq%interior(i,j,e2,ivar)
          enddo
        enddo
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,sums,6,surf%decomp%mpiPrec,MPI_SUM,surf%decomp%mpiComm,ierror)
    call MPI_Allreduce(MPI_IN_PLACE,errJ,1,surf%decomp%mpiPrec,MPI_MAX,surf%decomp%mpiComm,ierror)
    call MPI_Allreduce(MPI_IN_PLACE,scaleJ,1,surf%decomp%mpiPrec,MPI_MAX,surf%decomp%mpiComm,ierror)

    if(surf%decomp%rankId == 0) then
      print'(A,A,A,I0,A,I0,A,L1,A,L1)',"Case ",label,": N = ",N,", nGeo = ",nGeo, &
        ", Gauss-Lobatto ",nodeType == GAUSS_LOBATTO,", hybrid levels ",hybrid
      print'(A,ES22.14,A,ES22.14)',"    curl form   <BA,q>_3 = ",sums(1),"  <A,Cq>_2 = ",sums(2)
      print'(A,ES22.14,A,ES22.14)',"    determinant <BA,q>_3 = ",sums(4),"  <A,Cq>_2 = ",sums(5)
    endif
    call Report(r,tolerance,"(I1) curl-form Jacobian,   relative ",abs(sums(1)-sums(2))/sums(3),surf)
    call Report(r,tolerance,"(I1) determinant Jacobian, relative ",abs(sums(4)-sums(5))/sums(6),surf)
    call Report(r,tolerance,"det J3 - Z_zeta det J2,    relative ",errJ/scaleJ,surf)

    call q%DissociateGeometry()
    call A%DissociateGeometry()
    call cols%Free()
    call cq%Free()
    call A%Free()
    call bA%Free()
    call q%Free()
    call geom2%Free()
    call geom3%Free()
    call interp%Free()
    call mesh%Free()
    deallocate(h,zLevels)

  endfunction CheckCase

  integer function column_broadcast_adjoint_3d() result(r)
    implicit none
    ! Local
    type(Mesh2D) :: structured,circleHOPr,circle
    character(LEN=255) :: WORKSPACE
    integer :: n,seedSize
    integer,allocatable :: seed(:)

    call get_environment_variable("WORKSPACE",WORKSPACE)
    call structured%StructuredMesh(3,2,2,2,0.3_prec,0.25_prec, &
                                   [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED])
    call circleHOPr%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Circle/Circle_mesh.h5")

    ! A fixed seed per rank, so that a failure can be reproduced
    call random_seed(size=seedSize)
    allocate(seed(1:seedSize))
    do n = 1,seedSize
      seed(n) = 7919*n+104729*structured%decomp%rankId
    enddo
    call random_seed(put=seed)
    deallocate(seed)

    r = 0
    r = max(r,CheckCase(structured,"structured",4,GAUSS_LOBATTO,.false.))
    r = max(r,CheckCase(structured,"structured",7,GAUSS_LOBATTO,.true.))
    r = max(r,CheckCase(structured,"structured",7,GAUSS,.false.))
    r = max(r,CheckCase(circleHOPr,"HOPr Circle",4,GAUSS_LOBATTO,.false.))
    r = max(r,CheckCase(circleHOPr,"HOPr Circle",7,GAUSS_LOBATTO,.true.))
    r = max(r,CheckCase(circleHOPr,"HOPr Circle",7,GAUSS,.false.))
    if(structured%decomp%nRanks == 1) then
      ! The HOHQMesh reader holds every element on every rank, so this runs on one rank only
      call circle%Read_HOHQMesh(trim(WORKSPACE)//"/share/mesh/Circle/Circle.mesh")
      r = max(r,CheckCase(circle,"HOHQMesh Circle",4,GAUSS_LOBATTO,.false.))
      r = max(r,CheckCase(circle,"HOHQMesh Circle",7,GAUSS_LOBATTO,.true.))
      call circle%Free()
    endif

    call circleHOPr%Free()
    call structured%Free()

  endfunction column_broadcast_adjoint_3d

endprogram test
