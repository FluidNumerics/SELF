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
!! Checks VerticalSweepUp and VerticalSweepDown (SELF_ColumnOperators_3D), the upwind
!! DG solves along the vertical lines of an extruded mesh,
!!
!!   up:    (D + M^{-1} e_0 e_0^T) q^(k) = r^(k) + M^{-1} e_0 q_in^(k) ,  q_in^(1) = qBottom ,  q_in^(k) = q^(k-1)_N
!!   down:  (D - M^{-1} e_N e_N^T) q^(k) = r^(k) - M^{-1} e_N q_in^(k) ,  q_in^(nz) = qTop ,   q_in^(k) = q^(k+1)_0
!!
!! on Gauss-Lobatto nodes. Three properties are asserted.
!!
!!   1. Inverse accuracy. For N = 1..12 the precomputed inverses satisfy
!!      max |A A^{-1} - I| and max |A^{-1} A - I| <= 10 cond_inf(A) eps, with A built
!!      from interp%dMatrix in prec and the products formed in double precision.
!!   2. Exactness. On every vertical line, q^(k)(zeta) = f(k-1 + (1+zeta)/2) for a
!!      polynomial f of degree N of the layer coordinate s in [0,nz], with coefficients
!!      that vary from line to line and variable to variable. Each q^(k) is a
!!      polynomial of degree N in zeta and q is continuous across layers. With
!!      r^(k) = D q^(k), the upward sweep from f(0) and the downward sweep from f(nz)
!!      both return q, to 100 N^2 nz eps relative to max |q|.
!!   3. Adjointness (I2). With w from the upward sweep and p from the downward sweep, for
!!      random r_w, r_p and random inflows, on every vertical line,
!!
!!        sum_k [ p^(k)T M r_w^(k) + w^(k)T M r_p^(k) ] = p_in^(nz) w^(nz)_N - p^(1)_0 w_in^(1) ,
!!
!!      to 100 N^2 eps relative to the sum of the absolute values of all the terms.
!!      It follows from M D + D^T M = e_N e_N^T - e_0 e_0^T, exact on Gauss-Lobatto nodes:
!!      per element p^T M A w = -(At p)^T M w, and the inflow terms telescope.
!!
!! Cases 2 and 3 run on structured and HOPr Circle surfaces with sloping levels at
!! N = 4 and 7 (the sweeps do not use the geometry, but the mesh fixes the layout), on
!! one and on two ranks. A checksum of the exactness results keyed by global node is
!! printed so that the two runs can be compared bitwise.

  use iso_fortran_env,only:int64,real64
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

  exit_code = vertical_sweep_3d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

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

  subroutine Report(r,tolerance,name,err,surf)
    !! Reduces err over the ranks, prints it, and sets r = 1 if it exceeds the tolerance
    implicit none
    integer,intent(inout) :: r
    real(prec),intent(in) :: tolerance
    character(*),intent(in) :: name
    real(prec),intent(in) :: err
    type(Mesh2D),intent(in) :: surf
    ! Local
    real(prec) :: errGlobal
    integer :: ierror
    character(len=8) :: verdict

    call MPI_Allreduce(err,errGlobal,1,surf%decomp%mpiPrec,MPI_MAX,surf%decomp%mpiComm,ierror)
    if(errGlobal <= tolerance) then
      verdict = "exact"
    else
      verdict = "FAIL"
      r = 1
    endif
    if(surf%decomp%rankId == 0) print'(A,A,ES12.4,A,ES10.3,A,A)',"    ",name,errGlobal,"  (tolerance ",tolerance,")  ", &
      trim(verdict)

  endsubroutine Report

  subroutine BuildMesh(surf,nz,mesh)
    !! Extrudes surf with nz sigma levels over the sloping depth 1 + 0.3 x + 0.2 y^2
    implicit none
    type(Mesh2D),intent(in) :: surf
    integer,intent(in) :: nz
    type(Mesh3D),intent(out) :: mesh
    ! Local
    real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
    real(prec) :: sigma(0:nz)
    integer :: nGeo,i,j,e2,k

    nGeo = surf%nGeo
    allocate(h(1:nGeo+1,1:nGeo+1,1:surf%nElem))
    allocate(zLevels(1:nGeo+1,1:nGeo+1,0:nz,1:surf%nElem))
    do e2 = 1,surf%nElem
      do j = 1,nGeo+1
        do i = 1,nGeo+1
          h(i,j,e2) = 1.0_prec+0.3_prec*surf%nodeCoords(1,i,j,e2)+0.2_prec*surf%nodeCoords(2,i,j,e2)**2
        enddo
      enddo
    enddo
    ! Uneven stretching, sigma_k = -(1 - k/nz)^1.5
    do k = 0,nz
      sigma(k) = -(1.0_prec-real(k,prec)/real(nz,prec))**1.5_prec
    enddo
    sigma(0) = -1.0_prec
    sigma(nz) = 0.0_prec
    call ExtrudedSigmaLevels(h,sigma,zLevels)
    call mesh%ExtrudedMesh(surf,zLevels,SELF_BC_NONORMALFLOW,SELF_BC_PRESCRIBED)
    deallocate(h,zLevels)

  endsubroutine BuildMesh

  integer function CheckInverses(surf) result(r)
    !! Case 1, for N = 1..12
    implicit none
    type(Mesh2D),intent(in) :: surf
    ! Local
    integer,parameter :: nMax = 12
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(ColumnOperators3D) :: cols
    integer :: N,np,i,m,c,dir
    real(real64) :: a(1:nMax+1,1:nMax+1),aInv(1:nMax+1,1:nMax+1)
    real(real64) :: s1,s2,err,normA,normAinv
    real(prec) :: tolerance

    r = 0
    call BuildMesh(surf,2,mesh)
    if(surf%decomp%rankId == 0) print'(A)',"Case inverse accuracy, Gauss-Lobatto nodes"
    do N = 1,nMax
      np = N+1
      call interp%Init(N=N,controlNodeType=GAUSS_LOBATTO,M=N,targetNodeType=UNIFORM)
      call geometry%Init(interp,mesh%nElem)
      call geometry%GenerateFromMesh(mesh)
      call cols%Init(interp,mesh,geometry)
      do dir = 1,2
        ! a(i,m) = A(i,m) in row, column order; aInv(i,m) likewise
        do m = 1,np
          do i = 1,np
            a(i,m) = real(interp%dMatrix(m,i),real64)
            if(dir == 1) then
              aInv(i,m) = real(cols%upInverse(m,i),real64)
            else
              aInv(i,m) = real(cols%downInverse(m,i),real64)
            endif
          enddo
        enddo
        if(dir == 1) then
          a(1,1) = a(1,1)+1.0_real64/real(interp%qWeights(1),real64)
        else
          a(np,np) = a(np,np)-1.0_real64/real(interp%qWeights(np),real64)
        endif
        err = 0.0_real64
        do m = 1,np
          do i = 1,np
            s1 = 0.0_real64
            s2 = 0.0_real64
            do c = 1,np
              s1 = s1+a(i,c)*aInv(c,m)
              s2 = s2+aInv(i,c)*a(c,m)
            enddo
            if(i == m) then
              s1 = s1-1.0_real64
              s2 = s2-1.0_real64
            endif
            err = max(err,abs(s1),abs(s2))
          enddo
        enddo
        normA = 0.0_real64
        normAinv = 0.0_real64
        do i = 1,np
          normA = max(normA,sum(abs(a(i,1:np))))
          normAinv = max(normAinv,sum(abs(aInv(i,1:np))))
        enddo
        tolerance = real(10.0_real64*normA*normAinv*epsilon(1.0_prec),prec)
        if(dir == 1) then
          call Report(r,tolerance,"N = "//trim(Str(N))//", |A Ainv - I|   cond "//trim(StrE(normA*normAinv)),real(err,prec),surf)
        else
          call Report(r,tolerance,"N = "//trim(Str(N))//", |At Atinv - I| cond "//trim(StrE(normA*normAinv)),real(err,prec),surf)
        endif
      enddo
      call cols%Free()
      call geometry%Free()
      call interp%Free()
    enddo
    call mesh%Free()

  endfunction CheckInverses

  integer function CheckSweeps(surf,label,N,nz) result(r)
    !! Cases 2 and 3
    implicit none
    type(Mesh2D),intent(in) :: surf
    character(*),intent(in) :: label
    integer,intent(in) :: N
    integer,intent(in) :: nz
    ! Local
    integer,parameter :: nvar = 3
    type(Lagrange),target :: interp
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: qExact,r1,r2,w,p
    type(Scalar2D) :: inBottom,inTop
    type(ColumnOperators3D) :: cols
    integer :: i,j,l,m,k,e2,iel,ivar,ip,np,ierror
    real(prec) :: coef(0:N),s,f,x,y,dq,lhs,rhs,scale
    real(prec) :: errUp,errDown,errI2,qMax,tolExact,tolI2
    integer(int64) :: h,hGlobal,key

    r = 0
    np = N+1
    tolExact = 100.0_prec*real(N,prec)**2*real(nz,prec)*epsilon(1.0_prec)
    tolI2 = 100.0_prec*real(N,prec)**2*epsilon(1.0_prec)
    call BuildMesh(surf,nz,mesh)
    call interp%Init(N=N,controlNodeType=GAUSS_LOBATTO,M=N,targetNodeType=UNIFORM)
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)
    call cols%Init(interp,mesh,geometry)

    call qExact%Init(interp,nvar,mesh%nElem)
    call r1%Init(interp,nvar,mesh%nElem)
    call r2%Init(interp,nvar,mesh%nElem)
    call w%Init(interp,nvar,mesh%nElem)
    call p%Init(interp,nvar,mesh%nElem)
    call inBottom%Init(interp,nvar,surf%nElem)
    call inTop%Init(interp,nvar,surf%nElem)

    ! ---- 2. Exactness ----
    ! f(s) = sum_p coef_p (s/nz)^p, with coefficients that depend on the position of the
    ! vertical line and on the variable, so that the data are the same on any rank count
    do ivar = 1,nvar
      do e2 = 1,surf%nElem
        do j = 1,np
          do i = 1,np
            x = geometry%x%interior(i,j,1,(e2-1)*nz+1,1,1)
            y = geometry%x%interior(i,j,1,(e2-1)*nz+1,1,2)
            do ip = 0,N
              coef(ip) = cos(real(ip+1,prec)*(x+0.7_prec*y)+real(ivar,prec))
            enddo
            do k = 1,nz
              iel = (e2-1)*nz+k
              do l = 1,np
                s = (real(k-1,prec)+0.5_prec*(1.0_prec+interp%controlPoints(l)))/real(nz,prec)
                f = 0.0_prec
                do ip = N,0,-1
                  f = f*s+coef(ip)
                enddo
                qExact%interior(i,j,l,iel,ivar) = f
              enddo
            enddo
            inBottom%interior(i,j,e2,ivar) = coef(0)
            inTop%interior(i,j,e2,ivar) = sum(coef(0:N))
          enddo
        enddo
      enddo
    enddo
    ! r = D q in each element, D(l,m) = dMatrix(m,l)
    do ivar = 1,nvar
      do iel = 1,mesh%nElem
        do l = 1,np
          do j = 1,np
            do i = 1,np
              dq = 0.0_prec
              do m = 1,np
                dq = dq+interp%dMatrix(m,l)*qExact%interior(i,j,m,iel,ivar)
              enddo
              r1%interior(i,j,l,iel,ivar) = dq
            enddo
          enddo
        enddo
      enddo
    enddo
    call r1%UpdateDevice()
    call inBottom%UpdateDevice()
    call inTop%UpdateDevice()
    call cols%VerticalSweepUp(r1,inBottom,w)
    call cols%VerticalSweepDown(r1,inTop,p)
    call w%UpdateHost()
    call p%UpdateHost()

    qMax = maxval(abs(qExact%interior))
    call MPI_Allreduce(MPI_IN_PLACE,qMax,1,surf%decomp%mpiPrec,MPI_MAX,surf%decomp%mpiComm,ierror)
    errUp = maxval(abs(w%interior-qExact%interior))/qMax
    errDown = maxval(abs(p%interior-qExact%interior))/qMax

    ! Checksum of both results, keyed by global node
    h = 0_int64
    do ivar = 1,nvar
      do iel = 1,mesh%nElem
        do l = 1,np
          do j = 1,np
            do i = 1,np
              key = (((int(mesh%decomp%offsetElem(mesh%decomp%rankId+1)+iel-1,int64)*np+(l-1))*np+(j-1))*np+(i-1))*nvar+ivar
              h = ieor(h,ishftc(NodeBits(w%interior(i,j,l,iel,ivar)),int(mod(key,61_int64))+1))
              h = ieor(h,ishftc(NodeBits(p%interior(i,j,l,iel,ivar)),int(mod(key+29_int64,61_int64))+1))
            enddo
          enddo
        enddo
      enddo
    enddo
    call MPI_Allreduce(h,hGlobal,1,MPI_INTEGER8,MPI_BXOR,surf%decomp%mpiComm,ierror)

    if(surf%decomp%rankId == 0) then
      print'(A,A,A,I0,A,I0,A,I0)',"Case ",label,": N = ",N,", nz = ",nz,", columns on rank 0 ",surf%nElem
    endif
    call Report(r,tolExact,"exactness, upward sweep,   relative ",errUp,surf)
    call Report(r,tolExact,"exactness, downward sweep, relative ",errDown,surf)
    if(surf%decomp%rankId == 0) print'(A,Z16.16)',"    checksum (exactness results) ",hGlobal

    ! ---- 3. Adjointness (I2), random data ----
    call random_number(r1%interior)
    call random_number(r2%interior)
    call random_number(inBottom%interior)
    call random_number(inTop%interior)
    r1%interior = r1%interior-0.5_prec
    r2%interior = r2%interior-0.5_prec
    inBottom%interior = inBottom%interior-0.5_prec
    inTop%interior = inTop%interior-0.5_prec
    call r1%UpdateDevice()
    call r2%UpdateDevice()
    call inBottom%UpdateDevice()
    call inTop%UpdateDevice()
    call cols%VerticalSweepUp(r1,inBottom,w)
    call cols%VerticalSweepDown(r2,inTop,p)
    call w%UpdateHost()
    call p%UpdateHost()

    errI2 = 0.0_prec
    do ivar = 1,nvar
      do e2 = 1,surf%nElem
        do j = 1,np
          do i = 1,np
            lhs = 0.0_prec
            scale = 0.0_prec
            do k = 1,nz
              iel = (e2-1)*nz+k
              do l = 1,np
                lhs = lhs+interp%qWeights(l)*(p%interior(i,j,l,iel,ivar)*r1%interior(i,j,l,iel,ivar)+ &
                                              w%interior(i,j,l,iel,ivar)*r2%interior(i,j,l,iel,ivar))
                scale = scale+interp%qWeights(l)*(abs(p%interior(i,j,l,iel,ivar)*r1%interior(i,j,l,iel,ivar))+ &
                                                  abs(w%interior(i,j,l,iel,ivar)*r2%interior(i,j,l,iel,ivar)))
              enddo
            enddo
            rhs = inTop%interior(i,j,e2,ivar)*w%interior(i,j,np,e2*nz,ivar)- &
                  p%interior(i,j,1,(e2-1)*nz+1,ivar)*inBottom%interior(i,j,e2,ivar)
            scale = scale+abs(inTop%interior(i,j,e2,ivar)*w%interior(i,j,np,e2*nz,ivar))+ &
                    abs(p%interior(i,j,1,(e2-1)*nz+1,ivar)*inBottom%interior(i,j,e2,ivar))
            errI2 = max(errI2,abs(lhs-rhs)/scale)
          enddo
        enddo
      enddo
    enddo
    call Report(r,tolI2,"(I2) per vertical line,   relative ",errI2,surf)

    call cols%Free()
    call inTop%Free()
    call inBottom%Free()
    call p%Free()
    call w%Free()
    call r2%Free()
    call r1%Free()
    call qExact%Free()
    call geometry%Free()
    call interp%Free()
    call mesh%Free()

  endfunction CheckSweeps

  function Str(i) result(s)
    implicit none
    integer,intent(in) :: i
    character(len=8) :: s

    write(s,'(I0)') i

  endfunction Str

  function StrE(v) result(s)
    implicit none
    real(real64),intent(in) :: v
    character(len=12) :: s

    write(s,'(ES10.3)') v

  endfunction StrE

  integer function vertical_sweep_3d() result(r)
    implicit none
    ! Local
    type(Mesh2D) :: structured,circleHOPr
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
      seed(n) = 6151*n+15485863*structured%decomp%rankId
    enddo
    call random_seed(put=seed)
    deallocate(seed)

    r = 0
    r = max(r,CheckInverses(structured))
    r = max(r,CheckSweeps(structured,"structured",4,5))
    r = max(r,CheckSweeps(structured,"structured",7,4))
    r = max(r,CheckSweeps(circleHOPr,"HOPr Circle",4,3))
    r = max(r,CheckSweeps(circleHOPr,"HOPr Circle",7,4))

    call circleHOPr%Free()
    call structured%Free()

  endfunction vertical_sweep_3d

endprogram test
