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

module SELF_ColumnOperators_3D_t
!! Column operators on extruded meshes (Mesh3D%ExtrudedMesh).
!!
!! An extruded mesh stacks nz = mesh%nLayers hexahedra above every surface element
!! e2 of a Mesh2D. Local element e3 = (e2-1)*nz + k is layer k of that column, k = 1 at
!! the bottom, and a column is never split across MPI ranks. Every element maps
!!
!!   x(xi,eta,zeta) = ( X(xi,eta), Y(xi,eta), Z(xi,eta,zeta) ),
!!
!! so the lines of constant (xi,eta), the third index of SELF's interior arrays, are
!! vertical. The operators below act along those lines. They are local to a column,
!! so none of them communicates.
!!
!! Notation. N is the polynomial degree, w_l (l = 0..N) the quadrature weights,
!! M = diag(w), D the nodal derivative matrix, D(row=i,col=m) = interp%dMatrix(m,i),
!! e_0 and e_N the unit vectors of the bottom (zeta = -1) and the top (zeta = 1) node,
!! and Z_zeta = dz/dzeta, the (z,zeta) entry dxds%interior(...,1,3,3) of the covariant
!! basis SEMHex stores at the solution nodes. A superscript (k) is the layer.
!!
!! ColumnIntegral, volume to surface, in the unit of q times metres:
!!
!!   (C q)_ij = sum_{k=1..nz} sum_{l=0..N} w_l (Z_zeta)^(k)_ijl q^(k)_ijl .
!!
!! Broadcast, surface to volume, unit unchanged:
!!
!!   (B Q)^(k)_ijl = Q_ij .
!!
!! Because the Jacobian of an extruded element factors as J = J_2 Z_zeta at every node,
!! with J_2 the surface Jacobian, B is the adjoint of C,
!!
!!   <B A, q>_3 = <A, C q>_2 ,  <a,c>_3 = sum J w_i w_j w_l a c ,  <A,C>_2 = sum J_2 w_i w_j A C .
!!
!! VerticalSweepUp and VerticalSweepDown solve, on every vertical line of every element,
!! the one dimensional upwind DG problems of SELF-HPE's continuity equation and hydrostatic
!! relation. Upward, from a given bottom trace,
!!
!!   A q^(k) = r^(k) + M^{-1} e_0 q_in^(k) ,  A = D + M^{-1} e_0 e_0^T ,
!!   q_in^(1) = qBottom ,  q_in^(k) = q^(k-1)_N ;
!!
!! downward, from a given top trace,
!!
!!   At q^(k) = r^(k) - M^{-1} e_N q_in^(k) ,  At = D - M^{-1} e_N e_N^T ,
!!   q_in^(nz) = qTop ,  q_in^(k) = q^(k+1)_0 .
!!
!! The right hand side r is in reference scaling: the sweeps solve dq/dzeta = r in the
!! reference coordinate of each element, and they apply no metric factor. A caller that
!! solves a physical equation such as J2c dw/dzeta = -(horizontal terms) folds the metric
!! factors into r before the call. Since A 1 = M^{-1} e_0 and At 1 = -M^{-1} e_N, the
!! solution is the inflow value plus an element-local particular solution,
!!
!!   q^(k) = q_in^(k) 1 + A^{-1} r^(k)    (upward),
!!   q^(k) = q_in^(k) 1 + At^{-1} r^(k)   (downward),
!!
!! with A^{-1} and At^{-1} precomputed once in Init. The layers of a column are visited
!! in order, because each takes its inflow from the one before; columns, the nodes (i,j)
!! of a column and the variables are independent. By the summation-by-parts property
!! M D + D^T M = e_N e_N^T - e_0 e_0^T of Gauss-Lobatto nodes the two sweeps are adjoint:
!! with w from the upward sweep and p from the downward sweep, on every vertical line,
!!
!!   sum_k [ p^(k)T M r_w^(k) + w^(k)T M r_p^(k) ] = p_in^(nz) w^(nz)_N - p^(1)_0 w_in^(1) .
!!
!! The sweeps need nodes at zeta = -1 and zeta = 1, so they require Gauss-Lobatto control
!! nodes. The column integral and the broadcast work with any quadrature.
!!
!! Fields. Volume fields are Scalar3D, or any extension such as MappedScalar3D, on the
!! extruded mesh with nElem = nz*nColumns. Surface fields are Scalar2D, or any extension
!! such as MappedScalar2D, with nElem = nColumns, the same interpolant degree and the same
!! number of variables; surface element e2 is the column above the caller's surface
!! element e2. Only nodal values are touched, so neither field needs a geometry. On GPU
!! builds every operator works on the device copies (interior_gpu) and leaves the host
!! copies untouched.

  use iso_fortran_env,only:real64
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Scalar_2D
  use SELF_Scalar_3D
  use SELF_Mesh_3D
  use SELF_Geometry_3D

  implicit none

  type,public :: ColumnOperators3D_t
    integer :: N = 0
    !! Polynomial degree of the solution
    integer :: nLayers = 0
    !! Layers per column, nz
    integer :: nColumns = 0
    !! Rank-local columns, which is the rank-local element count of the surface mesh
    integer :: nElem = 0
    !! Rank-local 3-D elements, nLayers*nColumns
    logical :: lobatto = .false.
    !! True when the control nodes are Gauss-Lobatto, which the sweeps require
    type(Lagrange),pointer :: interp => null()
    type(SEMHex),pointer :: geometry => null()
    real(prec),pointer,contiguous,dimension(:,:) :: upInverse => null()
    !! A^{-1} in SELF's matrix convention: (A^{-1} r)_l = sum_m upInverse(m,l) r_m
    real(prec),pointer,contiguous,dimension(:,:) :: downInverse => null()
    !! At^{-1} in the same convention

  contains

    procedure,public :: Init => Init_ColumnOperators3D_t
    procedure,public :: Free => Free_ColumnOperators3D_t
    procedure,public :: ColumnIntegral => ColumnIntegral_ColumnOperators3D_t
    procedure,public :: Broadcast => Broadcast_ColumnOperators3D_t
    procedure,public :: VerticalSweepUp => VerticalSweepUp_ColumnOperators3D_t
    procedure,public :: VerticalSweepDown => VerticalSweepDown_ColumnOperators3D_t

  endtype ColumnOperators3D_t

  public :: Init_ColumnOperators3D_t
  public :: Free_ColumnOperators3D_t
  public :: CheckFields_ColumnOperators3D_t
  public :: CheckSweep_ColumnOperators3D_t

contains

  subroutine Init_ColumnOperators3D_t(this,interp,mesh,geometry)
    !! Binds the column operators to an extruded mesh and its geometry and, on
    !! Gauss-Lobatto nodes, precomputes the inverses of the two sweep matrices
    !!
    !!   A = D + M^{-1} e_0 e_0^T ,  At = D - M^{-1} e_N e_N^T .
    !!
    !! Both are (N+1) x (N+1). They are built from interp%dMatrix and inverted by
    !! Gauss-Jordan elimination with partial pivoting in double precision, then rounded
    !! to prec. The geometry must stay alive while the operators are used.
    !!
    !!  Input
    !!    - interp : the interpolant of the fields and of the geometry, any N >= 1
    !!    - mesh : a mesh built by Mesh3D%ExtrudedMesh (nLayers >= 1)
    !!    - geometry : SEMHex generated from mesh with interp
    implicit none
    class(ColumnOperators3D_t),intent(out) :: this
    type(Lagrange),target,intent(in) :: interp
    type(Mesh3D),intent(in) :: mesh
    type(SEMHex),target,intent(in) :: geometry
    ! Local
    integer :: i,m,np
    real(real64),allocatable :: a(:,:),aInv(:,:)

    if(mesh%nLayers <= 0) then
      error stop "ColumnOperators3D : the mesh is not extruded (nLayers = 0); build it with Mesh3D%ExtrudedMesh"
    endif
    if(mod(mesh%nElem,mesh%nLayers) /= 0) then
      error stop "ColumnOperators3D : the element count of the mesh is not a multiple of its layer count"
    endif
    if(geometry%nElem /= mesh%nElem) then
      error stop "ColumnOperators3D : the geometry does not have the element count of the mesh"
    endif
    if(geometry%dxds%N /= interp%N) then
      error stop "ColumnOperators3D : the geometry was not generated with this interpolant degree"
    endif

    this%interp => interp
    this%geometry => geometry
    this%N = interp%N
    this%nLayers = mesh%nLayers
    this%nElem = mesh%nElem
    this%nColumns = mesh%nElem/mesh%nLayers
    this%lobatto = interp%controlNodeType == GAUSS_LOBATTO

    if(this%lobatto) then
      np = interp%N+1
      allocate(this%upInverse(1:np,1:np),this%downInverse(1:np,1:np))
      allocate(a(1:np,1:np),aInv(1:np,1:np))

      ! A(i,m) = D(i,m) + delta_{i,0} delta_{m,0}/w_0, in the usual row, column order
      do m = 1,np
        do i = 1,np
          a(i,m) = real(interp%dMatrix(m,i),real64)
        enddo
      enddo
      a(1,1) = a(1,1)+1.0_real64/real(interp%qWeights(1),real64)
      call InvertMatrix(a,aInv,np)
      do m = 1,np
        do i = 1,np
          this%upInverse(m,i) = real(aInv(i,m),prec)
        enddo
      enddo

      ! At(i,m) = D(i,m) - delta_{i,N} delta_{m,N}/w_N
      do m = 1,np
        do i = 1,np
          a(i,m) = real(interp%dMatrix(m,i),real64)
        enddo
      enddo
      a(np,np) = a(np,np)-1.0_real64/real(interp%qWeights(np),real64)
      call InvertMatrix(a,aInv,np)
      do m = 1,np
        do i = 1,np
          this%downInverse(m,i) = real(aInv(i,m),prec)
        enddo
      enddo

      deallocate(a,aInv)
    endif

  endsubroutine Init_ColumnOperators3D_t

  subroutine Free_ColumnOperators3D_t(this)
    implicit none
    class(ColumnOperators3D_t),intent(inout) :: this

    if(associated(this%upInverse)) deallocate(this%upInverse)
    if(associated(this%downInverse)) deallocate(this%downInverse)
    this%upInverse => null()
    this%downInverse => null()
    this%interp => null()
    this%geometry => null()
    this%N = 0
    this%nLayers = 0
    this%nColumns = 0
    this%nElem = 0
    this%lobatto = .false.

  endsubroutine Free_ColumnOperators3D_t

  subroutine CheckFields_ColumnOperators3D_t(this,volume,surface,caller)
    !! Stops when a volume and a surface field do not match the extruded mesh the
    !! operators were initialized for, or each other.
    implicit none
    class(ColumnOperators3D_t),intent(in) :: this
    class(Scalar3D),intent(in) :: volume
    class(Scalar2D),intent(in) :: surface
    character(*),intent(in) :: caller

    if(this%nLayers <= 0) then
      print*,"ColumnOperators3D%",caller," : Error : the column operators are not initialized"
      error stop 1
    endif
    if(volume%N /= this%N .or. surface%N /= this%N) then
      print*,"ColumnOperators3D%",caller," : Error : field degree differs from the operators' degree ",this%N
      error stop 1
    endif
    if(volume%nElem /= this%nElem) then
      print*,"ColumnOperators3D%",caller," : Error : volume field has ",volume%nElem," elements, the mesh ",this%nElem
      error stop 1
    endif
    if(surface%nElem /= this%nColumns) then
      print*,"ColumnOperators3D%",caller," : Error : surface field has ",surface%nElem," elements, the mesh ", &
        this%nColumns," columns"
      error stop 1
    endif
    if(volume%nVar /= surface%nVar) then
      print*,"ColumnOperators3D%",caller," : Error : volume field has ",volume%nVar," variables, the surface field ", &
        surface%nVar
      error stop 1
    endif

  endsubroutine CheckFields_ColumnOperators3D_t

  subroutine CheckSweep_ColumnOperators3D_t(this,r,inflow,q,caller)
    !! The checks of CheckFields for a vertical sweep, plus the requirement of
    !! Gauss-Lobatto nodes and of a right hand side and solution of the same shape.
    implicit none
    class(ColumnOperators3D_t),intent(in) :: this
    class(Scalar3D),intent(in) :: r
    class(Scalar2D),intent(in) :: inflow
    class(Scalar3D),intent(in) :: q
    character(*),intent(in) :: caller

    call CheckFields_ColumnOperators3D_t(this,r,inflow,caller)
    if(.not. this%lobatto) then
      print*,"ColumnOperators3D%",caller," : Error : the vertical sweeps need Gauss-Lobatto control nodes,", &
        " since the inflow is imposed at the end nodes zeta = -1 and zeta = 1"
      error stop 1
    endif
    if(q%N /= r%N .or. q%nElem /= r%nElem .or. q%nVar /= r%nVar) then
      print*,"ColumnOperators3D%",caller," : Error : the solution and the right hand side differ in shape"
      error stop 1
    endif

  endsubroutine CheckSweep_ColumnOperators3D_t

  subroutine ColumnIntegral_ColumnOperators3D_t(this,q,qbar)
    !! Column integral of a volume field,
    !!
    !!   qbar_ij = sum_{k=1..nz} sum_{l=0..N} w_l (Z_zeta)^(k)_ijl q^(k)_ijl ,
    !!
    !! for every column and variable. With Z_zeta in metres, qbar has the unit of q
    !! times metres; the column integral of 1 is the local depth. On Gauss-Lobatto nodes
    !! the quadrature is exact when Z_zeta q is a polynomial of degree 2N-1 in zeta in
    !! every layer, on Gauss nodes when it is of degree 2N+1.
    !!
    !!  Input
    !!    - q : volume field on the extruded mesh
    !!  Output
    !!    - qbar : surface field, overwritten
    implicit none
    class(ColumnOperators3D_t),intent(in) :: this
    class(Scalar3D),intent(in) :: q
    class(Scalar2D),intent(inout) :: qbar
    ! Local
    integer :: i,j,l,k,e2,iel,ivar,nz

    call CheckFields_ColumnOperators3D_t(this,q,qbar,"ColumnIntegral")
    nz = this%nLayers

    ! The sum runs over the layers from the bottom up and, inside a layer, over the
    ! nodes from the bottom up, the same order as the device kernel. The innermost loop
    ! is the contiguous index i.
    do ivar = 1,q%nVar
      do e2 = 1,this%nColumns
        do j = 1,this%N+1
          do i = 1,this%N+1
            qbar%interior(i,j,e2,ivar) = 0.0_prec
          enddo
        enddo
        do k = 1,nz
          iel = (e2-1)*nz+k
          do l = 1,this%N+1
            do j = 1,this%N+1
              do i = 1,this%N+1
                qbar%interior(i,j,e2,ivar) = qbar%interior(i,j,e2,ivar)+ &
                                             this%interp%qWeights(l)* &
                                             this%geometry%dxds%interior(i,j,l,iel,1,3,3)* &
                                             q%interior(i,j,l,iel,ivar)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

  endsubroutine ColumnIntegral_ColumnOperators3D_t

  subroutine Broadcast_ColumnOperators3D_t(this,qbar,q)
    !! Broadcast of a surface field to every node of its column,
    !!
    !!   q^(k)_ijl = qbar_ij ,  k = 1..nz, l = 0..N ,
    !!
    !! the adjoint of ColumnIntegral in the inner products described in the module
    !! header. The unit is unchanged.
    !!
    !!  Input
    !!    - qbar : surface field
    !!  Output
    !!    - q : volume field on the extruded mesh, overwritten
    implicit none
    class(ColumnOperators3D_t),intent(in) :: this
    class(Scalar2D),intent(in) :: qbar
    class(Scalar3D),intent(inout) :: q
    ! Local
    integer :: i,j,l,iel,ivar

    call CheckFields_ColumnOperators3D_t(this,q,qbar,"Broadcast")

    do concurrent(i=1:this%N+1,j=1:this%N+1,l=1:this%N+1, &
                  iel=1:this%nElem,ivar=1:q%nVar)
      q%interior(i,j,l,iel,ivar) = qbar%interior(i,j,(iel-1)/this%nLayers+1,ivar)
    enddo

  endsubroutine Broadcast_ColumnOperators3D_t

  subroutine VerticalSweepUp_ColumnOperators3D_t(this,r,qBottom,q)
    !! Upward sweep: solves, layer by layer from the bottom of every column,
    !!
    !!   (D + M^{-1} e_0 e_0^T) q^(k) = r^(k) + M^{-1} e_0 q_in^(k) ,
    !!   q_in^(1) = qBottom ,  q_in^(k) = q^(k-1)_N ,
    !!
    !! on every vertical line (i,j) and for every variable, as q^(k) = q_in^(k) 1 + A^{-1} r^(k).
    !! This is the upwind DG discretization of dq/dzeta = r with the inflow at the bottom
    !! of each element. r is in reference scaling: its unit is that of q, per unit zeta,
    !! and no metric factor is applied (see the module header). Requires Gauss-Lobatto
    !! control nodes.
    !!
    !!  Input
    !!    - r : right hand side, volume field
    !!    - qBottom : inflow at the bottom of every column, surface field
    !!  Output
    !!    - q : solution, volume field, overwritten. It must be a different object from r.
    implicit none
    class(ColumnOperators3D_t),intent(in) :: this
    class(Scalar3D),intent(in) :: r
    class(Scalar2D),intent(in) :: qBottom
    class(Scalar3D),intent(inout) :: q
    ! Local
    integer :: i,j,l,m,k,e2,iel,ivar,nz,np

    call CheckSweep_ColumnOperators3D_t(this,r,qBottom,q,"VerticalSweepUp")
    nz = this%nLayers
    np = this%N+1

    ! For each node l the particular solution sum_m upInverse(m,l) r_m is accumulated in
    ! q in the order m = 1..N+1, and the inflow is added last, the same operations as
    ! the device kernel. The innermost loop is the contiguous index i.
    do ivar = 1,r%nVar
      do e2 = 1,this%nColumns
        do k = 1,nz
          iel = (e2-1)*nz+k
          do l = 1,np
            do j = 1,np
              do i = 1,np
                q%interior(i,j,l,iel,ivar) = 0.0_prec
              enddo
            enddo
            do m = 1,np
              do j = 1,np
                do i = 1,np
                  q%interior(i,j,l,iel,ivar) = q%interior(i,j,l,iel,ivar)+ &
                                               this%upInverse(m,l)*r%interior(i,j,m,iel,ivar)
                enddo
              enddo
            enddo
          enddo
          ! Inflow: the bottom trace for the bottom layer, else the top node of the layer below
          if(k == 1) then
            do l = 1,np
              do j = 1,np
                do i = 1,np
                  q%interior(i,j,l,iel,ivar) = qBottom%interior(i,j,e2,ivar)+q%interior(i,j,l,iel,ivar)
                enddo
              enddo
            enddo
          else
            do l = 1,np
              do j = 1,np
                do i = 1,np
                  q%interior(i,j,l,iel,ivar) = q%interior(i,j,np,iel-1,ivar)+q%interior(i,j,l,iel,ivar)
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
    enddo

  endsubroutine VerticalSweepUp_ColumnOperators3D_t

  subroutine VerticalSweepDown_ColumnOperators3D_t(this,r,qTop,q)
    !! Downward sweep: solves, layer by layer from the top of every column,
    !!
    !!   (D - M^{-1} e_N e_N^T) q^(k) = r^(k) - M^{-1} e_N q_in^(k) ,
    !!   q_in^(nz) = qTop ,  q_in^(k) = q^(k+1)_0 ,
    !!
    !! on every vertical line (i,j) and for every variable, as q^(k) = q_in^(k) 1 + At^{-1} r^(k).
    !! This is the upwind DG discretization of dq/dzeta = r with the inflow at the top of
    !! each element, the adjoint of VerticalSweepUp. r is in reference scaling, as for
    !! VerticalSweepUp. Requires Gauss-Lobatto control nodes.
    !!
    !!  Input
    !!    - r : right hand side, volume field
    !!    - qTop : inflow at the top of every column, surface field
    !!  Output
    !!    - q : solution, volume field, overwritten. It must be a different object from r.
    implicit none
    class(ColumnOperators3D_t),intent(in) :: this
    class(Scalar3D),intent(in) :: r
    class(Scalar2D),intent(in) :: qTop
    class(Scalar3D),intent(inout) :: q
    ! Local
    integer :: i,j,l,m,k,e2,iel,ivar,nz,np

    call CheckSweep_ColumnOperators3D_t(this,r,qTop,q,"VerticalSweepDown")
    nz = this%nLayers
    np = this%N+1

    do ivar = 1,r%nVar
      do e2 = 1,this%nColumns
        do k = nz,1,-1
          iel = (e2-1)*nz+k
          do l = 1,np
            do j = 1,np
              do i = 1,np
                q%interior(i,j,l,iel,ivar) = 0.0_prec
              enddo
            enddo
            do m = 1,np
              do j = 1,np
                do i = 1,np
                  q%interior(i,j,l,iel,ivar) = q%interior(i,j,l,iel,ivar)+ &
                                               this%downInverse(m,l)*r%interior(i,j,m,iel,ivar)
                enddo
              enddo
            enddo
          enddo
          ! Inflow: the top trace for the top layer, else the bottom node of the layer above
          if(k == nz) then
            do l = 1,np
              do j = 1,np
                do i = 1,np
                  q%interior(i,j,l,iel,ivar) = qTop%interior(i,j,e2,ivar)+q%interior(i,j,l,iel,ivar)
                enddo
              enddo
            enddo
          else
            do l = 1,np
              do j = 1,np
                do i = 1,np
                  q%interior(i,j,l,iel,ivar) = q%interior(i,j,1,iel+1,ivar)+q%interior(i,j,l,iel,ivar)
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
    enddo

  endsubroutine VerticalSweepDown_ColumnOperators3D_t

  subroutine InvertMatrix(a,aInv,n)
    !! Inverts the n x n matrix a by Gauss-Jordan elimination with partial pivoting, in
    !! double precision. Used once in Init for the (N+1) x (N+1) sweep matrices, which
    !! are small and well enough conditioned for a direct solve; this follows the
    !! routine of the same name in SELF_RefinementIndicator_3D_t.
    implicit none
    integer,intent(in) :: n
    real(real64),intent(in) :: a(1:n,1:n)
    real(real64),intent(out) :: aInv(1:n,1:n)
    ! Local
    real(real64) :: b(1:n,1:n)
    integer :: i,j,k,piv
    real(real64) :: pmax,factor,tmp

    b = a
    aInv = 0.0_real64
    do i = 1,n
      aInv(i,i) = 1.0_real64
    enddo

    do k = 1,n
      ! Partial pivot: the largest entry of column k at or below the diagonal
      piv = k
      pmax = abs(b(k,k))
      do i = k+1,n
        if(abs(b(i,k)) > pmax) then
          pmax = abs(b(i,k))
          piv = i
        endif
      enddo
      if(pmax <= 0.0_real64) then
        error stop "ColumnOperators3D : singular sweep matrix"
      endif
      if(piv /= k) then
        do j = 1,n
          tmp = b(k,j); b(k,j) = b(piv,j); b(piv,j) = tmp
          tmp = aInv(k,j); aInv(k,j) = aInv(piv,j); aInv(piv,j) = tmp
        enddo
      endif

      ! Normalize the pivot row
      factor = b(k,k)
      do j = 1,n
        b(k,j) = b(k,j)/factor
        aInv(k,j) = aInv(k,j)/factor
      enddo

      ! Eliminate column k from every other row
      do i = 1,n
        if(i /= k) then
          factor = b(i,k)
          do j = 1,n
            b(i,j) = b(i,j)-factor*b(k,j)
            aInv(i,j) = aInv(i,j)-factor*aInv(k,j)
          enddo
        endif
      enddo
    enddo

  endsubroutine InvertMatrix

endmodule SELF_ColumnOperators_3D_t
