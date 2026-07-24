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

module SELF_Lagrange_t

  use iso_fortran_env
  use iso_c_binding

  use SELF_Constants
  use SELF_SupportRoutines
  use SELF_Quadrature
  use SELF_HDF5
  use HDF5

  use iso_c_binding

  implicit none

  type,public :: Lagrange_t
    !! A data structure for working with Lagrange Interpolating Polynomials in one, two, and three dimensions.
    !! The Lagrange data-structure stores the information necessary to interpolate between two
    !! sets of grid-points and to estimate the derivative of data at native grid points. Routines for
    !! multidimensional interpolation are based on the tensor product of 1-D interpolants. It is
    !! assumed that the polynomial degree (and the interpolation nodes) are the same in each direction.
    !! This assumption permits the storage of only one array of interpolation nodes and barycentric
    !! weights and is what allows this data structure to be flexible.

    integer :: N
      !! The number of control points.

    integer :: controlNodeType

    integer :: M
      !! The number of target points.

    integer :: targetNodeType

    real(prec),pointer,contiguous,dimension(:) :: controlPoints
      !! The set of nodes in one dimension where data is known.
      !! To create higher dimension interpolation and differentiation operators, structured grids in two and three
      !! dimensions are created by tensor products of the controlPoints. This design decision implies that all
      !! spectral element methods supported by the Lagrange class have the same polynomial degree in each
      !! computational/spatial dimension. In practice, the controlPoints are the Legendre-Gauss, Legendre-Gauss-Lobatto,
      !! Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, or Chebyshev-Gauss-Radau quadrature points over
      !! the domain [-1,1] (computational space). The Init routine for this class restricts controlPoints to one of
      !! these quadrature types or uniform points on [-1,1].

    real(prec),pointer,contiguous,dimension(:) :: targetPoints
      !! The set of nodes in one dimension where data is to be interpolated to. To create higher dimension interpolation
      !! and differentiation operators, structured grids in two and three dimensions are created by tensor products of
      !! the targetPoints. In practice, the targetPoints are set to a uniformly distributed set of points between [-1,1]
      !! (computational space) to allow for interpolation from unevenly spaced quadrature points to a plotting grid.

    real(prec),pointer,contiguous,dimension(:) :: bWeights
      !! The barycentric weights that are calculated from the controlPoints and used for interpolation.

    real(prec),pointer,contiguous,dimension(:) :: qWeights
      !! The quadrature weights for discrete integration. The quadradture weights depend on the type of controlPoints
      !! provided; one of Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss,
      !! Chebyshev-Gauss-Lobatto, Chebyshev-Gauss Radau, or Uniform. If Uniform, the quadrature weights are constant
      !! $$dx = \frac{2.0}{N+1}$$.

    real(prec),pointer,contiguous,dimension(:,:) :: iMatrix
      !! The interpolation matrix (transpose) for mapping data from the control grid to the target grid.

    real(prec),pointer,contiguous,dimension(:,:) :: dMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The
      !! dMatrix is based on a strong form of the derivative.

    real(prec),pointer,contiguous,dimension(:,:) :: dgMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The dgMatrix is based
      !! on a weak form of the derivative. It must be used with bMatrix to account for boundary contributions in the weak form.

    real(prec),pointer,contiguous,dimension(:,:) :: dSplitMatrix
      !! The split-form derivative matrix D_split = D - 0.5*M^{-1}*B, where B is the SBP boundary operator.
      !!
      !! D_split is skew-symmetric under the M inner product: M*D_split + D_split^T*M = 0. Unlike D, it is
      !! NOT an SBP operator (the weighted symmetric part is zero, not the boundary term B). This property
      !! makes D_split ideal for the EC-DGSEM split-form volume integral:
      !!
      !!   du/dt = -(2/J) sum_n D_split[n,i] * F_EC(u_i, u_n)  +  (1/J) M^{-1} B^T f_Riemann
      !!
      !! The skew-symmetry of the volume term guarantees it contributes zero to the entropy rate, so all
      !! entropy change passes through the surface term. An entropy-dissipative Riemann solver then makes
      !! the full scheme entropy-stable.
      !!
      !! Using D_split in the volume combined with plain f_Riemann on the surface is algebraically identical
      !! to using D in the volume with the penalty (f_Riemann - f_local) on the surface (Trixi.jl convention).
      !!
      !! In SELF index convention (dSplitMatrix(ii,i) = D_split[i-1, ii-1]):
      !!   dSplitMatrix(ii,i) = dMatrix(ii,i)
      !!                       - 0.5*(bMatrix(i,2)*bMatrix(ii,2)
      !!                              - bMatrix(i,1)*bMatrix(ii,1)) / qWeights(i)

    real(prec),pointer,contiguous,dimension(:,:) :: bMatrix
      !! The boundary interpolation matrix that is used to map a grid of nodal values at the control points to the element boundaries.

    real(prec),pointer,contiguous,dimension(:,:,:) :: mortarR
      !! Mortar restriction matrices for 2:1 nonconforming (mortar) interfaces.
      !!
      !! An element edge trace u(xi) = sum_i u_i l_i(xi) on [-1,1] ("big" side) is restricted to
      !! the two half-edges ("small" sides) of its 2:1 neighbors. Sub-edge k=1 occupies
      !! xi in [-1,0] and k=2 occupies xi in [0,1], with the sub-edge coordinate s in [-1,1]
      !! related to the big-edge coordinate through xi_1(s) = (s-1)/2 and xi_2(s) = (s+1)/2.
      !!
      !! mortarR(ii,i,k) = l_ii( xi_k(s_i) ), where s_i are the control points; following the
      !! SELF matrix convention, the restricted trace is u_k(i) = sum_ii mortarR(ii,i,k)*u(ii).
      !! Restriction of a degree N polynomial is exact.

    real(prec),pointer,contiguous,dimension(:,:,:) :: mortarP
      !! Mortar L2-projection matrices, the adjoints of mortarR under the exact L2 inner
      !! product on [-1,1]: given traces g_k(s) on the two sub-edges, the projected big-side
      !! trace is u(i) = sum_k sum_ii mortarP(ii,i,k)*g_k(ii).
      !!
      !! The matrices are built as P_k = M^{-1} B_k with the exact (dense) 1-D mass matrix
      !! M(m,i) = int l_m l_i dxi and B_k(m,j) = (1/2) int l_m(xi_k(s)) l_j(s) ds, evaluated
      !! with an internal Gauss rule that is exact for the degree 2N integrands. Consequently
      !!   sum_k matmul(P_k,R_k) = I   (a polynomial split across both sub-edges is recovered)
      !!   sum_m w_m (P_k g)_m = (1/2) sum_j w_j g_k(j)   (discrete conservation)
      !! for any control node type. mortarP carries the 1/2 sub-edge Jacobian appropriate for
      !! projecting *solution* traces; when projecting surface-flux integrands (which carry the
      !! small side's nScale = nScale_big/2), scale by 2 (see MortarFluxCollect in the 2-D
      !! mapped data classes).

  contains

    procedure,public :: Init => Init_Lagrange_t
    procedure,public :: Free => Free_Lagrange_t

    procedure,public :: WriteHDF5 => WriteHDF5_Lagrange_t

    procedure,public :: CalculateBarycentricWeights
    procedure,public :: CalculateInterpolationMatrix
    procedure,public :: CalculateDerivativeMatrix
    procedure,public :: CalculateLagrangePolynomials
    procedure,public :: CalculateMortarMatrices

  endtype Lagrange_t

contains

  subroutine Init_Lagrange_t(this,N,controlNodeType,M,targetNodeType)
    !! Initialize an instance of the Lagrange_t class
    !! On output, all of the attributes for the Lagrange_t class are allocated and values are initialized according to the number of
    !! control points, number of target points, and the types for the control and target nodes.
    !! If a GPU is available, device pointers for the Lagrange_t attributes are allocated and initialized.
    implicit none
    class(Lagrange_t),intent(out) :: this
    !! Lagrange_t class instance
    integer,intent(in)          :: N
    !! The number of control points for interpolant
    integer,intent(in)          :: M
    !! The number of target points for the interpolant
    integer,intent(in)          :: controlNodeType
    !! The integer code specifying the type of control points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    integer,intent(in)          :: targetNodeType
    !! The integer code specifying the type of target points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    ! -------!
    ! Local
    real(prec) :: q(0:M)

    this%N = N
    this%M = M
    this%controlNodeType = controlNodeType
    this%targetNodeType = targetNodeType
    allocate(this%controlPoints(1:N+1), &
             this%targetPoints(1:M+1), &
             this%bWeights(1:N+1), &
             this%qWeights(1:N+1), &
             this%iMatrix(1:N+1,1:M+1), &
             this%dMatrix(1:N+1,1:N+1), &
             this%dgMatrix(1:N+1,1:N+1), &
             this%dSplitMatrix(1:N+1,1:N+1), &
             this%bMatrix(1:N+1,1:2), &
             this%mortarR(1:N+1,1:N+1,1:2), &
             this%mortarP(1:N+1,1:N+1,1:2))

    if(controlNodeType == GAUSS .or. controlNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(N, &
                              this%controlPoints, &
                              this%qWeights, &
                              controlNodeType)

    elseif(controlNodeType == CHEBYSHEV_GAUSS .or. controlNodeType == CHEBYSHEV_GAUSS_LOBATTO) then

      call ChebyshevQuadrature(N, &
                               this%controlPoints, &
                               this%qWeights, &
                               controlNodeType)

    elseif(controlNodeType == UNIFORM) then

      this%controlPoints = UniformPoints(-1.0_prec,1.0_prec,0,N)
      this%qWeights = 2.0_prec/real(N,prec)

    endif

    ! Target Points
    if(targetNodeType == GAUSS .or. targetNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(M, &
                              this%targetPoints, &
                              q, &
                              targetNodeType)

    elseif(targetNodeType == UNIFORM) then

      this%targetPoints = UniformPoints(-1.0_prec,1.0_prec,0,M)

    endif

    call this%CalculateBarycentricWeights()
    call this%CalculateInterpolationMatrix()
    call this%CalculateDerivativeMatrix()
    this%bMatrix(1:N+1,1) = this%CalculateLagrangePolynomials(-1.0_prec)
    this%bMatrix(1:N+1,2) = this%CalculateLagrangePolynomials(1.0_prec)

    ! D_split = D - 0.5*M^{-1}*B  (skew-symmetric under the M inner product)
    ! In SELF index convention: dSplitMatrix(ii,i) = dMatrix(ii,i)
    !   - 0.5*(bMatrix(i,2)*bMatrix(ii,2) - bMatrix(i,1)*bMatrix(ii,1)) / qWeights(i)
    block
      integer :: i,ii
      do i = 1,N+1
        do ii = 1,N+1
          this%dSplitMatrix(ii,i) = this%dMatrix(ii,i)- &
                                    0.5_prec*(this%bMatrix(i,2)*this%bMatrix(ii,2)- &
                                              this%bMatrix(i,1)*this%bMatrix(ii,1))/ &
                                    this%qWeights(i)
        enddo
      enddo
    endblock

    call this%CalculateMortarMatrices()

  endsubroutine Init_Lagrange_t

  subroutine Free_Lagrange_t(this)
    !! Frees all memory (host and device) associated with an instance of the Lagrange_t class
    implicit none
    class(Lagrange_t),intent(inout) :: this
    !! Lagrange_t class instance

    deallocate(this%controlPoints)
    deallocate(this%targetPoints)
    deallocate(this%bWeights)
    deallocate(this%qWeights)
    deallocate(this%iMatrix)
    deallocate(this%dMatrix)
    deallocate(this%dgMatrix)
    deallocate(this%dSplitMatrix)
    deallocate(this%bMatrix)
    deallocate(this%mortarR)
    deallocate(this%mortarP)

  endsubroutine Free_Lagrange_t

! ================================================================================================ !
!
! CalculateBarycentricWeights (PRIVATE)
!
!   A PRIVATE routine that calculates and stores the barycentric weights for the Lagrange_t
!   data-structure.
!
!   This routine is from Alg. 30 on pg. 75 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateBarycentricWeights(this)
    implicit none
    class(Lagrange_t),intent(inout) :: this
    ! Local
    integer :: i,j
    real(real64) :: bWeights(0:this%N)
    real(real64) :: controlPoints(0:this%N)

    do i = 0,this%N
      bWeights(i) = 1.0_real64
      controlPoints(i) = real(this%controlPoints(i+1),real64)
    enddo

    ! Computes the product w_k = w_k*(s_k - s_j), k /= j
    do j = 1,this%N
      do i = 0,j-1

        bWeights(i) = bWeights(i)*(controlPoints(i)-controlPoints(j))
        bWeights(j) = bWeights(j)*(controlPoints(j)-controlPoints(i))

      enddo
    enddo

    do j = 0,this%N
      bWeights(j) = 1.0_prec/bWeights(j)
      this%bWeights(j+1) = real(bWeights(j),prec)
    enddo

  endsubroutine CalculateBarycentricWeights

! ================================================================================================ !
!
! CalculateInterpolationMatrix (PRIVATE)
!
!   A PRIVATE routine that fills in the interpolation matrix for the Lagrange_t data structure.
!
!   This function is from Alg. 32 on pg. 76 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateInterpolationMatrix(this)
    implicit none
    class(Lagrange_t),intent(inout) :: this
    ! Local
    integer    :: row,col
    logical    :: rowHasMatch
    real(real64) :: temp1,temp2
    real(real64) :: iMatrix(0:this%M,0:this%N)
    real(real64) :: bWeights(0:this%N)
    real(real64) :: controlPoints(0:this%N)
    real(real64) :: targetPoints(0:this%M)

    do col = 0,this%N
      controlPoints(col) = real(this%controlPoints(col+1),real64)
      bWeights(col) = real(this%bWeights(col+1),real64)
    enddo
    do row = 0,this%M
      targetPoints(row) = real(this%targetPoints(row+1),real64)
    enddo

    do row = 0,this%M

      rowHasMatch = .false.

      do col = 0,this%N

        iMatrix(row,col) = 0.0_real64

        if(AlmostEqual(targetPoints(row),controlPoints(col))) then
          rowHasMatch = .true.
          iMatrix(row,col) = 1.0_real64
        endif

      enddo

      if(.not.(rowHasMatch)) then

        temp1 = 0.0_real64

        do col = 0,this%N
          temp2 = bWeights(col)/ &
                  (targetPoints(row)- &
                   controlPoints(col))
          iMatrix(row,col) = temp2
          temp1 = temp1+temp2
        enddo

        do col = 0,this%N
          iMatrix(row,col) = iMatrix(row,col)/temp1
        enddo

      endif

    enddo

    do row = 0,this%M
      do col = 0,this%N
        this%iMatrix(col+1,row+1) = real(iMatrix(row,col),prec)
      enddo
    enddo

  endsubroutine CalculateInterpolationMatrix

! ================================================================================================ !
!
! CalculateDerivativeMatrix (PRIVATE)
!
!   Calculates and stores the derivative matrix and its transpose.
!   Generates a matrix that can be used to approximate derivatives at the interpolation nodes.
!
!   This function is from Alg. 37 on pg. 82 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateDerivativeMatrix(this)
    implicit none
    class(Lagrange_t),intent(inout) :: this
    ! Local
    integer      :: row,col
    real(real64) :: dmat(0:this%N,0:this%N)
    real(real64) :: dgmat(0:this%N,0:this%N)
    real(real64) :: bWeights(0:this%N)
    real(real64) :: qWeights(0:this%N)
    real(real64) :: controlPoints(0:this%N)

    do row = 0,this%N
      bWeights(row) = real(this%bWeights(row+1),real64)
      qWeights(row) = real(this%qWeights(row+1),real64)
      controlPoints(row) = real(this%controlPoints(row+1),real64)
    enddo

    do row = 0,this%N

      dmat(row,row) = 0.0_prec

      do col = 0,this%N

        if(.not.(col == row)) then

          dmat(row,col) = bWeights(col)/ &
                          (bWeights(row)* &
                           (controlPoints(row)- &
                            controlPoints(col)))

          dmat(row,row) = dmat(row,row)-dmat(row,col)

        endif

      enddo

    enddo

    do row = 0,this%N
      do col = 0,this%N
        dgmat(row,col) = -dmat(col,row)* &
                         qWeights(col)/ &
                         qWeights(row)
      enddo
    enddo

    do row = 0,this%N
      do col = 0,this%N
        this%dMatrix(row+1,col+1) = real(dmat(col,row),prec)
        this%dgMatrix(row+1,col+1) = real(dgmat(col,row),prec)
      enddo
    enddo

  endsubroutine CalculateDerivativeMatrix

! ================================================================================================ !
!
! CalculateLagrangePolynomials
!
!   Evaluates each of the 1-D Lagrange interpolating polynomials at a specified point.
!
!   This function is from Alg. 34 on pg. 77 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  function CalculateLagrangePolynomials(this,sE) result(lAtS)
    implicit none
    class(Lagrange_t) :: this
    real(prec)      :: sE
    real(prec)      :: lAtS(0:this%N)
    ! Local
    integer    :: j
    logical    :: xMatchesNode
    real(real64) :: temp1,temp2
    real(real64) :: sELocal
    real(real64) :: controlPoints(0:this%N)
    real(real64) :: bWeights(0:this%N)
    real(real64) :: lS(0:this%N)

    sELocal = real(sE,real64)
    do j = 0,this%N
      controlPoints(j) = real(this%controlPoints(j+1),real64)
      bWeights(j) = real(this%bWeights(j+1),real64)
    enddo

    xMatchesNode = .false.

    do j = 0,this%N

      lS(j) = 0.0_real64
      if(AlmostEqual(sELocal,controlPoints(j))) then
        lS(j) = 1.0_real64
        xMatchesNode = .true.
      endif

    enddo

    if(xMatchesNode) then
      do j = 0,this%N
        lAtS(j) = real(lS(j),prec)
      enddo
      return
    endif

    temp1 = 0.0_real64

    do j = 0,this%N
      temp2 = bWeights(j)/(sE-controlPoints(j))
      lS(j) = temp2
      temp1 = temp1+temp2
    enddo

    lS = lS/temp1

    do j = 0,this%N
      lAtS(j) = real(lS(j),prec)
    enddo

  endfunction CalculateLagrangePolynomials

! ================================================================================================ !
!
! CalculateMortarMatrices
!
!   Builds the 2:1 mortar restriction (mortarR) and L2 projection (mortarP) matrices used at
!   nonconforming element interfaces (see the type-bound attribute documentation above for the
!   mathematical definitions and discrete identities).
!
!   The projection uses the exact 1-D mass matrix M(m,i) = int l_m l_i dxi and the mixed mass
!   matrices B_k(m,j) = (1/2) int l_m(xi_k(s)) l_j(s) ds. All integrals are computed with an
!   internal (N+2)-point Legendre-Gauss rule, which is exact for the degree 2N integrands, so
!   the identities sum_k P_k R_k = I and the discrete conservation property hold to roundoff
!   for any control node type. The linear systems M P_k = B_k are solved in float64 with
!   partially-pivoted Gaussian elimination (the systems are (N+1)x(N+1)).
!
! ================================================================================================ !

  subroutine CalculateMortarMatrices(this)
    implicit none
    class(Lagrange_t),intent(inout) :: this
    ! Local
    integer :: i,ii,k,q,nq
    real(prec) :: gpts(0:this%N+1)
    real(prec) :: gwts(0:this%N+1)
    real(prec) :: lbig(0:this%N)
    real(prec) :: lsub(0:this%N)
    real(prec) :: xik
    real(real64) :: mmat(1:this%N+1,1:this%N+1)
    real(real64) :: mwork(1:this%N+1,1:this%N+1)
    real(real64) :: bmat(1:this%N+1,1:this%N+1,1:2)
    real(real64) :: pmat(1:this%N+1,1:this%N+1)

    nq = this%N+1 ! Internal quadrature rule has nq+1 = N+2 points; exact for degree 2N+3

    call LegendreQuadrature(nq,gpts,gwts,GAUSS)

    mmat = 0.0_real64
    bmat = 0.0_real64
    do q = 0,nq

      ! Values of the control-point Lagrange basis at the internal quadrature point; used
      ! both for the sub-edge (small side) basis in B_k and for the big-side basis in M.
      lsub = this%CalculateLagrangePolynomials(gpts(q))

      do i = 1,this%N+1
        do ii = 1,this%N+1
          mmat(ii,i) = mmat(ii,i)+ &
                       real(gwts(q),real64)*real(lsub(ii-1),real64)*real(lsub(i-1),real64)
        enddo
      enddo

      do k = 1,2
        ! Map the sub-edge coordinate s = gpts(q) to the big-edge coordinate
        if(k == 1) then
          xik = 0.5_prec*(gpts(q)-1.0_prec) ! xi in [-1,0]
        else
          xik = 0.5_prec*(gpts(q)+1.0_prec) ! xi in [0,1]
        endif
        lbig = this%CalculateLagrangePolynomials(xik)

        ! B_k(m,j) = (1/2) int l_m(xi_k(s)) l_j(s) ds ; the 1/2 is the sub-edge Jacobian
        do i = 1,this%N+1 ! j index (sub-edge basis)
          do ii = 1,this%N+1 ! m index (big-edge basis)
            bmat(ii,i,k) = bmat(ii,i,k)+ &
                           0.5_real64*real(gwts(q),real64)* &
                           real(lbig(ii-1),real64)*real(lsub(i-1),real64)
          enddo
        enddo
      enddo

    enddo

    do k = 1,2

      ! Solve M P_k = B_k for the projection matrix; the solver overwrites its
      ! matrix argument with the factorization, so factor a fresh copy of M
      mwork = mmat
      pmat = bmat(1:this%N+1,1:this%N+1,k)
      call SolveLinearSystems(this%N+1,mwork,pmat)

      ! Store the transpose, following the SELF operator convention
      ! result(i) = sum_ii matrix(ii,i)*input(ii)
      do i = 1,this%N+1
        do ii = 1,this%N+1
          this%mortarP(ii,i,k) = real(pmat(i,ii),prec)
        enddo
      enddo

      ! Restriction: evaluate the big-side basis at the sub-edge images of the control points
      do i = 1,this%N+1
        if(k == 1) then
          xik = 0.5_prec*(this%controlPoints(i)-1.0_prec)
        else
          xik = 0.5_prec*(this%controlPoints(i)+1.0_prec)
        endif
        lbig = this%CalculateLagrangePolynomials(xik)
        do ii = 1,this%N+1
          this%mortarR(ii,i,k) = lbig(ii-1)
        enddo
      enddo

    enddo

  endsubroutine CalculateMortarMatrices

! ================================================================================================ !
!
! SolveLinearSystems (PRIVATE)
!
!   Solves A X = B for X, where A is n x n and B holds n right-hand sides, using Gaussian
!   elimination with partial pivoting in float64. On exit, B is overwritten with X and A is
!   overwritten with its factorization. Intended for the small (N+1)x(N+1) operator systems
!   assembled at initialization; not for use in hot loops.
!
! ================================================================================================ !

  subroutine SolveLinearSystems(n,a,b)
    implicit none
    integer,intent(in) :: n
    real(real64),intent(inout) :: a(1:n,1:n)
    real(real64),intent(inout) :: b(1:n,1:n)
    ! Local
    integer :: i,j,p,ipiv
    real(real64) :: pivmax,factor,swap(1:n)

    do j = 1,n-1

      ! Partial pivoting
      ipiv = j
      pivmax = abs(a(j,j))
      do p = j+1,n
        if(abs(a(p,j)) > pivmax) then
          ipiv = p
          pivmax = abs(a(p,j))
        endif
      enddo
      if(ipiv /= j) then
        swap = a(j,1:n)
        a(j,1:n) = a(ipiv,1:n)
        a(ipiv,1:n) = swap
        swap = b(j,1:n)
        b(j,1:n) = b(ipiv,1:n)
        b(ipiv,1:n) = swap
      endif

      ! Elimination
      do i = j+1,n
        factor = a(i,j)/a(j,j)
        a(i,j:n) = a(i,j:n)-factor*a(j,j:n)
        b(i,1:n) = b(i,1:n)-factor*b(j,1:n)
      enddo

    enddo

    ! Back substitution
    do j = n,1,-1
      do i = j+1,n
        b(j,1:n) = b(j,1:n)-a(j,i)*b(i,1:n)
      enddo
      b(j,1:n) = b(j,1:n)/a(j,j)
    enddo

  endsubroutine SolveLinearSystems

  subroutine WriteHDF5_Lagrange_t(this,fileId)
    implicit none
    class(Lagrange_t),intent(in) :: this
    integer(HID_T),intent(in) :: fileId

    call CreateGroup_HDF5(fileId,'/interp')

    call WriteArray_HDF5(fileId,'/interp/controlpoints', &
                         this%controlPoints)

    call WriteArray_HDF5(fileId,'/interp/qweights', &
                         this%qWeights)

    call WriteArray_HDF5(fileId,'/interp/dgmatrix', &
                         this%dgMatrix)

    call WriteArray_HDF5(fileId,'/interp/dmatrix', &
                         this%dMatrix)

    call WriteArray_HDF5(fileId,'/interp/bmatrix', &
                         this%bMatrix)

    call WriteArray_HDF5(fileId,'/interp/imatrix', &
                         this%iMatrix)

  endsubroutine WriteHDF5_Lagrange_t

endmodule SELF_Lagrange_t
