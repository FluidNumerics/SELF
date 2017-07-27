! Lagrange_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
! Copyright 2016 Joseph Schoonover <jschoonover@lanl.gov>, Los Alamos National Laboratory
!
! The SELF and accompanying documentation were produced in part under the 
! support of Florida State University and the National Science Foundation 
! through Grant OCE-1049131 during 2015 and in part  the support of the 
! Center for Nonlinear Studies and the Department of Energy through the 
! LANL/LDRD program in 2016.
!
! Lagrange_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
! 
! Licensed under the Apache License, Version 2.0 (the "License"); 
! You may obtain a copy of the License at 
!
! http://www.apache.org/licenses/LICENSE-2.0 
!
! Unless required by applicable law or agreed to in writing, software 
! distributed under the License is distributed on an "AS IS" BASIS, 
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and  
! limitations under the License.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file Lagrange_Class.f90
!! Contains the \ref Lagrange_Class module, and <BR>
!! defines the \ref Lagrange data-structure.


!> \defgroup Lagrange_Class Lagrange_Class 
!! This module defines the Lagrange data-structure and its associated routines.

MODULE Lagrange_Class

!src/common
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
! src/interp
USE InterpolationSupportRoutines

IMPLICIT NONE

!> \addtogroup Lagrange_Class 
!! @{

!> \struct Lagrange
!! A data-structure for handling Lagrange interpolation in one, two, or three dimensions
!!
!! The Lagrange data-structure stores the information necessary to interpolate between two
!! sets of grid-points and to estimate the derivative of data at native grid points. Routines for
!! multidimensional interpolation are based on the tensor product on two 1-D interpolants. It is 
!! assumed that the polynomial degree (and the interpolation nodes) are the same in each direction.
!! This assumption permits the storage of only one array of interpolation nodes and barycentric 
!! weights and is what allows this data structure to be flexible.
!!
!! <H2> Lagrange </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Number of native grid points
!!       <tr> <th> M <td> INTEGER <td> Number of target grid points
!!       <tr> <th> s(0:N) <td> REAL(prec) <td> Locations where we have observations (native grid points)
!!       <tr> <th> bWs(0:N) <td> REAL(prec) <td> Barycentric interpolation weights 
!!       <tr> <th> so(0:M) <td> REAL(prec) <td> Locations where we want observations (target grid points)
!!       <tr> <th> Ts(0:M,0:N) <td> REAL(prec) <td> Interpolation matrix to help map an 
!!                                    array of data given at the native nodes to the target nodes.
!!       <tr> <th> Tp(0:M,0:N) <td> REAL(prec) <td> Transpose of the interpolation matrix to help
!!                                     map an array of data given at the native nodes to the 
!!                                     target nodes.
!!       <tr> <th> D(0:N,0:N) <td> REAL(prec) <td> Derivative matrix to estimate the 
!!                                    derivative of an interpolant at the native grid points in the 
!!                                    first computational direction.
!!       <tr> <th> DTr(0:N,0:N) <td> REAL(prec) <td> Derivative matrix to estimate the 
!!                                    derivative of an interpolant at the native grid points in the 
!!                                    second computational direction.
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Lagrange_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_Lagrange
!!       <tr> <th> Trash <td> Trash_Lagrange
!!       <tr> <th> CalculateLagrangePolynomials <td> CalculateLagrangePolynomials_Lagrange
!!       <tr> <th> Interpolate_1D <td> Interpolate_1D_Lagrange
!!       <tr> <th> Interpolate_2D <td> Interpolate_2D_Lagrange
!!       <tr> <th> Interpolate_3D <td> Interpolate_3D_Lagrange
!!       <tr> <th> Differentiate_1D <td> Differentiate_1D_Lagrange
!!       <tr> <th> Differentiate_2D <td> Differentiate_2D_Lagrange
!!       <tr> <th> Differentiate_3D <td> Differentiate_3D_Lagrange
!!       <tr> <th> ApplyInterpolationMatrix_1D <td> ApplyInterpolationMatrix_1D_Lagrange
!!       <tr> <th> ApplyInterpolationMatrix_2D <td> ApplyInterpolationMatrix_2D_Lagrange
!!       <tr> <th> ApplyInterpolationMatrix_3D <td> ApplyInterpolationMatrix_3D_Lagrange
!!       <tr> <th> ApplyDerivativeMatrix_1D <td> ApplyDerivativeMatrix_1D_Lagrange
!!       <tr> <th> ApplyDerivativeMatrix_2D <td> ApplyDerivativeMatrix_2D_Lagrange
!!       <tr> <th> ApplyDerivativeMatrix_3D <td> ApplyDerivativeMatrix_3D_Lagrange
!!       <tr> <th> WriteTecplot_1D <td> WriteTecplot_1D_Lagrange
!!       <tr> <th> WriteTecplot_2D <td> WriteTecplot_2D_Lagrange
!!       <tr> <th> WriteTecplot_3D <td> WriteTecplot_3D_Lagrange
!!    </table>
!!

!>@}

   TYPE, PUBLIC :: Lagrange
      INTEGER                 :: N     ! number of nodal points where we have observations
      INTEGER                 :: Nc    ! Number of columns needed to collapse a 3-D array to a 2-D array (for differentiation)
      INTEGER                 :: M     ! number of nodal points where we want observations
      INTEGER                 :: Mc    ! Number of columns needed to collapse a 3-D array to a 2-D array (for interpolation)
      REAL(prec), ALLOCATABLE :: s(:)    ! locations where we have obervations
      REAL(prec), ALLOCATABLE :: bWs(:)  ! barycentric weights
      REAL(prec), ALLOCATABLE :: so(:)   ! Locations where we want observations
      REAL(prec), ALLOCATABLE :: Ts(:,:) ! Interpolation matrix to get us from what we have to what
                                         ! we want, in the "s" direction.
      REAL(prec), ALLOCATABLE :: Tp(:,:) ! Transpose of "Ts". Interpolation from one set of points
                                         ! to another may occur frequently in application. To avoid
                                         ! computing the transpose repeatedly, we store it here 
      REAL(prec), ALLOCATABLE :: D(:,:) ! Derivative matrix to calculate the derivative of an 
                                         ! interpolant at the interpolation nodes ("s"). 
                                         ! This derivative matrix is used to calculate the 
                                         ! derivative in the "s" computational direction.
      REAL(prec), ALLOCATABLE :: DTr(:,:) ! Derivative matrix to calculate the derivative of an 
                                         ! interpolant at the interpolation nodes ("s").
                                         ! This derivative matrix is used to calculate the 
                                         ! derivative in the "p" computational direction. We store
                                         ! a separate matrix for the "p" computational direction,
                                         ! b/c it is the transpose of "D" and in practice,
                                         ! the derivative operation will be called frequently. Due
                                         ! to the expected frequency of calls to compute the 
                                         ! derivative, storing "DTr" avoids having to compute the
                                         ! transpose repeatedly.
      CONTAINS
      
      !-------------!
      ! Constructors/Destructors
      PROCEDURE :: Build => Build_Lagrange
      PROCEDURE :: Trash => Trash_Lagrange

      ! Type-Specific
      PROCEDURE, PRIVATE :: CalculateBarycentricWeights  => CalculateBarycentricWeights_Lagrange
      PROCEDURE, PRIVATE :: CalculateInterpolationMatrix => CalculateInterpolationMatrix_Lagrange
      PROCEDURE, PRIVATE :: CalculateDerivativeMatrix    => CalculateDerivativeMatrix_Lagrange
      
      PROCEDURE :: CalculateLagrangePolynomials => CalculateLagrangePolynomials_Lagrange
      PROCEDURE :: Interpolate_1D               => Interpolate_1D_Lagrange
      PROCEDURE :: Interpolate_2D               => Interpolate_2D_Lagrange
      PROCEDURE :: Interpolate_3D               => Interpolate_3D_Lagrange
      PROCEDURE :: Differentiate_1D             => Differentiate_1D_Lagrange
      PROCEDURE :: Differentiate_2D             => Differentiate_2D_Lagrange
      PROCEDURE :: Differentiate_3D             => Differentiate_3D_Lagrange
      PROCEDURE :: ApplyInterpolationMatrix_1D  => ApplyInterpolationMatrix_1D_Lagrange
      PROCEDURE :: ApplyInterpolationMatrix_2D  => ApplyInterpolationMatrix_2D_Lagrange
      PROCEDURE :: ApplyInterpolationMatrix_3D  => ApplyInterpolationMatrix_3D_Lagrange
      PROCEDURE :: ApplyInterpolationMatrices_3D  => ApplyInterpolationMatrices_3D_Lagrange
      PROCEDURE :: ApplyDerivativeMatrix_1D     => ApplyDerivativeMatrix_1D_Lagrange
      PROCEDURE :: ApplyDerivativeMatrix_2D     => ApplyDerivativeMatrix_2D_Lagrange
      PROCEDURE :: ApplyDerivativeMatrix_3D     => ApplyDerivativeMatrix_3D_Lagrange
      PROCEDURE :: ApplyDerivativeMatrices_3D   => ApplyDerivativeMatrices_3D_Lagrange
      
      ! File I/O
      PROCEDURE :: WriteTecplot_1D => WriteTecplot_1D_Lagrange
      PROCEDURE :: WriteTecplot_2D => WriteTecplot_2D_Lagrange
      PROCEDURE :: WriteTecplot_3D => WriteTecplot_3D_Lagrange
      
    END TYPE Lagrange

 INTEGER, PRIVATE :: Ncim = 2
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R Build
!
!>  \fn Build_Lagrange
!!  A manual constructor for the Lagrange class that allocates memory and fills in data 
!!  for the attributes of the Lagrange class.
!! 
!!  The Build subroutine allocates memory for the native and non-native grid points, barycentric
!!  weights interpolation matrix, and derivative matrix. The native and non-native grid points are
!!  filled in using the REAL(prec) input arrays "s" and "so". The barycentric weights are then 
!!  calculated and stored. Once the barycentric weights are calculated, the interpolation and
!!  derivative matrices are calculated and stored.
!!
!!  <H2> Usage : </H2>
!!     <B>TYPE</B>(Lagrange) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( N, M, s, so ) <BR>
!!
!!  <table>
!!       <tr> <td> in/out <th> myPoly <td> TYPE(Lagrange) <td> The Lagrange data structure to 
!!                                                                be constructed
!!       <tr> <td> in <th> N <td> INTEGER <td> The number of native grid points
!!       <tr> <td> in <th> M <td> INTEGER <td> The number of target grid points
!!       <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> The native grid points
!!       <tr> <td> in <th> so(0:N) <td> REAL(prec) <td> The target grid points
!!  </table>
!!
! =============================================================================================== !
!>@}
 SUBROUTINE Build_Lagrange( myPoly, N, M, s, so )

   IMPLICIT NONE
   CLASS(Lagrange), INTENT(out)   :: myPoly
   INTEGER, INTENT(in)            :: N, M
   REAL(prec), INTENT(in)         :: s(0:N), so(0:M)
   
      ! Set the number of observations (those we have and those we want)
      myPoly % N  = N
      myPoly % Nc = N*(N+2)
      myPoly % M  = M
      myPoly % Mc = M*(M+2)
      
      ! Allocate storage
      ALLOCATE( myPoly % s(0:N), myPoly % bWs(0:N) )
      ALLOCATE( myPoly % so(0:M), myPoly % Ts(0:M,0:N), myPoly % Tp(0:N,0:M) )
      ALLOCATE( myPoly % D(0:N,0:N), myPoly % DTr(0:N,0:N) )
      myPoly % s   = 0.0_prec
      myPoly % bWs = 0.0_prec
      myPoly % so  = 0.0_prec
      myPoly % Ts  = 0.0_prec
      myPoly % Tp  = 0.0_prec
      myPoly % D   = 0.0_prec
      myPoly % DTr = 0.0_prec
      
      ! Fill in the nodal locations
      myPoly % s(0:N)  = s(0:N)
      myPoly % so(0:M) = so(0:M)

      ! and calculate the barycentric weights for quick interpolation.
      CALL myPoly % CalculateBarycentricWeights( )

      ! Using the two nodal locations, we can construct the interpolation matrix. The interpolation
      ! matrix enables quick interpolation.
      CALL myPoly % CalculateInterpolationMatrix( )

      CALL myPoly % CalculateDerivativeMatrix( )
 
 END SUBROUTINE Build_Lagrange
!
!
!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_Lagrange
!! A manual destructor for the Lagrange class that deallocates the memory held by its 
!!  attributes. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!!  
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myPoly <td> TYPE(Lagrange) <td> 
!!                       On <B>input</B>, the Lagrange data structure with attributes filled in. <BR>
!!                       On <B>output</B>,the memory associated with this data-structure is freed.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
SUBROUTINE Trash_Lagrange(myPoly)

   IMPLICIT NONE
   CLASS(Lagrange), INTENT(inout) :: myPoly

      DEALLOCATE( myPoly % s, myPoly % bWs )
      DEALLOCATE( myPoly % so, myPoly % Ts, myPoly % Tp )
      DEALLOCATE( myPoly % D, myPoly % DTr )

 END SUBROUTINE Trash_Lagrange
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateBarycentricWeights
! 
!> \fn CalculateBarycentricWeights_Lagrange
!!  A PRIVATE routine that calculates and stores the barycentric weights for the Lagrange 
!!  data-structure.
!! 
!!  Calculates the barycentric weights from the interpolation nodes and stores them in the "bWs"
!!  attribute. This routine should be called after the native interpolation nodes have been 
!!  assigned.
!!
!!  This subroutine depends on
!!       Module \ref InterpolatioNupportRoutiness.f90 : Function \ref barycentricweights
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateBarycentricWeights( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myPoly <td> TYPE(Lagrange) <td>
!!           On <B>input</B>, myPoly is the Lagrange data structure is sent in with the <BR>
!!           native interpolation nodes already filled in. <BR>
!!           On <B>output</B>, myPoly has the barycentric weights filled in.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE CalculateBarycentricWeights_Lagrange( myPoly )

   IMPLICIT NONE
   CLASS(Lagrange), INTENT(inout) :: myPoly
   ! Local
   REAL(prec) :: s(0:myPoly % N)
   INTEGER    :: N

      N = myPoly % N
      s = myPoly % s

      myPoly % bWs = BarycentricWeights( s, N )
 
 END SUBROUTINE CalculateBarycentricWeights_Lagrange
!
!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateLagrangePolynomials 
! 
!> \fn CalculateLagrangePolynomials_Lagrange  
!! Evaluates each of the (1-D, yes 1-D!) Lagrange interpolating polynomials at a specified point. 
!! 
!! The Lagrange interpolating polynomials are given by 
!!   \f[ l_j(\xi) = \prod_{i=0,\neq j}^N \left( \frac{\xi-\xi_i}{\xi_j-\xi_i} \right)  
!!   \f] 
!! where the \f$ \xi_i \f$ are the native nodes (the "s" attribute). Given an input value (\f$\xi\f$),
!! this function returns each of the Lagrange interpolating polynomials \f$ l_j(\xi) \f$ evaluated
!! at this point. This is useful if you have multiple arrays of data that are given at the same
!! native nodes and require interpolation onto a single point.
!!
!!  This function depends on
!!       Module \ref InterpolatioNupportRoutiness.f90 : Function \ref lagrangepolynomials
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)        :: lAtS(0:this % N,0:this % N) <BR>
!! <B>REAL</B>(prec)        :: sE, pE <BR>
!!         .... <BR>
!!     lAtS = this % CalculateLagrangePolynomials( sE ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange structure. <BR>
!!                     The interpolation nodes and barycentric weights are required to produce 
!!                     sensible output.
!!   <tr> <td> in <th> sE <td> REAL(prec) <td> 
!!                     Location to evaluate the Lagrange interpolating polyomials
!!   <tr> <td> out <th> lAtS(0:myPoly % N) <td> REAL(prec) <td>
!!                      Array containing the value of each Lagrange interpolating polynomial 
!!                      at sE. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION CalculateLagrangePolynomials_Lagrange( myPoly, sE ) RESULT( lAtS )  

   IMPLICIT NONE
   CLASS(Lagrange) :: myPoly
   REAL(prec)      :: sE
   REAL(prec)      :: lAtS(0:myPoly % N)
   ! LOCAL
   REAL(prec) :: s(0:myPoly % N), w(0:myPoly % N)
   INTEGER    :: N

      N = myPoly % N
      s = myPoly % s
      w = myPoly % bWs

      lAtS = LagrangePolynomials( s, w, sE, N )


 END FUNCTION CalculateLagrangePolynomials_Lagrange
!
!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateInterpolationMatrix (PRIVATE) 
! 
!> \fn CalculateInterpolationMatrix_Lagrange 
!! A PRIVATE routine that fills in the interpolation matrix for the Lagrange data structure.
!! 
!! We can write the operations of interpolating data from one set of points to another (in this 
!! case from "s" to "so") as
!!            \f[ f_{m,n} = \sum_{i,j=0}^N f_{i,j} l_i(\xi^1_m)l_j(\xi^2_n) , \hspace{2mm} m,n=0,1,2,...,M
!!            \f]
!! where \f$ l_i(\xi^1) \f$ and \f$ l_j(\xi^2) \f$ are the Lagrange interpolating polynomials at the 
!! \f$ \lbrace ( \xi^1_i, \xi^2_j ) \rbrace_{i,j=0}^N \f$ nodes. Evaluation of the Lagrange 
!! interpolating polynomials at each of the new points 
!! \f$ \lbrace (\xi^1_m, \xi^2_n )\rbrace_{m,n=0}^M \f$ in each computational directions can be
!! written as two matrices, where 
!!       \f[ T^{(1)}_{m,i} = l_i(\xi^1_m).
!!        \f]
!! and 
!!       \f[ T^{(2)}_{n,j} = l_j(\xi^2_n).
!!        \f]
!! For simplicity, the native interpolation nodes are assumed identical in each computational
!! direction. Similarly, the target nodes are assumed identical in each direction. This assumption
!! implies that \f$ T^{(1)} \f$ and \f$ T^{(2)} \f$ are identical.
!!
!! In the SELF, interpolation onto the target grid (in 2-D) is executed via two matrix-matrix 
!! multiplications. The native data, \f$ f_{i,j} \f$ can be viewed as an \f$ N+1 \times N+1 \f$ 
!! matrix. Matrix multiplication (on the left) by \f$ T \f$ maps the native data to the target nodes
!! in the \f$ \xi^1 \f$ direction,
!!
!! \f[
!!      \tilde{f} = T f 
!! \f]
!! where \f$ \tilde{f} \f$ is now viewed as an \f$ M+1 \times N+1 \f$ matrix. Multiplication on the
!! right by the transpose of the interpolation matrix completes the 2-D interpolation onto the 
!! target nodes
!! \f[
!!      f_{target} = \tilde{f} T^{T} 
!! \f]
!!
!! Because of this implementation, this routine fills in the "Ts" attribute with the interpolation
!! matrix, and the "Tp" attribute with the transpose of the interpolation matrix.
!!
!! This subroutine depends on <BR>
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref interpolationmatrix
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateInterpolationMatrix(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myPoly <td> TYPE <td> 
!!             A previously constructed Lagrange data-structure. <BR>
!!             On <B>input</B>, the native interpolation nodes, target interpolation nodes,
!!             and the barycentric weights must be defined. <BR>
!!             On <B>output</B>, the interpolation matrix and its transpose are filled in.
!!  </table>
!!
! ================================================================================================ ! 
!>@}
 SUBROUTINE CalculateInterpolationMatrix_Lagrange( myPoly )  

   IMPLICIT NONE
   CLASS(Lagrange), INTENT(inout) :: myPoly
   ! LOCAL
   REAL(prec) :: s(0:myPoly % N), w(0:myPoly % N)
   REAL(prec) :: so(0:myPoly % M), T(0:myPoly % M, 0:myPoly % N)
   INTEGER    :: N, Nnew

      N    = myPoly % N
      nNew = myPoly % M
      s    = myPoly % s
      w    = myPoly % bWs
      so   = myPoly % so

      T = InterpolationMatrix( s, w, so, N, nNew )

      myPoly % Ts = T
      myPoly % Tp = TRANSPOSE( T )

 END SUBROUTINE CalculateInterpolationMatrix_Lagrange
!
!> \addtogroup Lagrange_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateDerivativeMatrix 
! 
!> \fn CalculateDerivativeMatrix_Lagrange  
!! Calculates and stores the derivative matrix and its transpose to estimate the derivative of a 
!! function ( in both computational directions ) at the native interpolation nodes.
!! 
!! Given nodal values of an interpolant, the derivative can be estimated at the interpolation 
!! nodes using the summations
!!      \f[ \lbrace \frac{\partial f}{\partial \xi^1} \rbrace_{m,n=0}^N = \sum_{i=0}^N( f_{i,n} l'_i(\xi^1_m) )
!!      \f]
!! and
!!      \f[ \lbrace \frac{\partial f}{\partial \xi^2} \rbrace_{m,n=0}^N = \sum_{j=0}^N( f_{m,j} l'_j(\xi^2_n) )
!!      \f]
!!
!! The native interpolation nodes are assumed identical in each direction so that the derivative
!! matrices in each direction are identical. We can write the derivatives as matrix-matrix products
!!      \f[ \lbrace \frac{\partial f}{\partial \xi^1} \rbrace_{m,n=0}^N = D f
!!      \f]
!! and
!!      \f[ \lbrace \frac{\partial f}{\partial \xi^2} \rbrace_{m,n=0}^N = f D^T
!!      \f] 
!! where 
!!      \f[ D_{j,i} = l'_i(\xi_j) 
!!      \f]
!! and \f$ f \f$ is the 2-D array of the native data that is viewed as a matrix.
!!
!! This subroutine calculates the derivative matrix and its transpose and stores them in the
!!  data-structure attributes "D" and "DTr" respectively.
!! 
!! This subroutine depends on <BR>
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref derivativematrix
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CaclulateDerivativeMatrix(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myPoly <td> TYPE <td> 
!!             A previously constructed Lagrange data-structure. <BR>
!!             On <B>input</B>, the native interpolation nodes
!!             and the barycentric weights must be defined. <BR>
!!             On <B>output</B>, the derivative matrix is filled in.
!!  </table>
! ================================================================================================ ! 
!>@}
 SUBROUTINE CalculateDerivativeMatrix_Lagrange( myPoly )  

   IMPLICIT NONE
   CLASS(Lagrange), INTENT(inout) :: myPoly
   ! LOCAL
   REAL(prec) :: s(0:myPoly % N), w(0:myPoly % N)
   REAL(prec) :: dMat(0:myPoly % N, 0:myPoly % N)
   INTEGER    :: N

      N = myPoly % N
      s = myPoly % s
      w = myPoly % bWs
      
      dMat = DerivativeMatrix( w, s, N )
  
      myPoly % D = dMat
      myPoly % DTr = TRANSPOSE( dMat )

 END SUBROUTINE CalculateDerivativeMatrix_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! Function Interpolate 
! 
!> \fn Interpolate_1D_Lagrange
!!  Interpolates an array of data onto a desired location using Lagrange interpolation.
!!
!!  This function depends on <BR>
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref interpolation
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: f(0:this % N) <BR>
!! <B>REAL</B>(prec)     :: sE <BR>
!! <B>REAL</B>(prec)     :: fAtSE <BR>
!!         .... <BR>
!!     fAtSE = this % Interpolate_1D( f, sE ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N) <td> REAL(prec) <td> 
!!                     An array of function nodal values located at the native interpolation nodes.
!!   <tr> <td> in <th> sE <td> REAL(prec) <td> 
!!                     The location where you want to interpolate to.
!!   <tr> <td> out <th> fAtSE <td> REAL(prec) <td> The interpolant evaluated at sE.
!!  </table>  
!!
! ================================================================================================ ! 
!>@}
 FUNCTION Interpolate_1D_Lagrange( myPoly, f, sE ) RESULT( interpF )  

   IMPLICIT NONE
   CLASS(Lagrange) :: myPoly
   REAL(prec)      :: sE
   REAL(prec)      :: f(0:myPoly % N)
   REAL(prec)      :: interpF
   ! LOCAL
   REAL(prec) :: w(0:myPoly % N), s(0:myPoly % N)
   INTEGER    :: N 
   
     N = myPoly % N
     w = myPoly % bWs
     s = myPoly % s
     
     interpF = Interpolation( s, w, f, N, sE )

 END FUNCTION Interpolate_1D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! Function Interpolate_2D 
! 
!> \fn Interpolate_2D_Lagrange 
!!  Interpolates an array of data onto a desired location using Lagrange interpolation.
!!
!!  This function depends on <BR>
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref interpolation
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: f(0:this % N, 0:this % N) <BR>
!! <B>REAL</B>(prec)     :: sE(1:2) <BR>
!! <B>REAL</B>(prec)     :: fAtSE <BR>
!!         .... <BR>
!!     fAtSE = this % Interpolate( f, sE ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N) <td> REAL(prec) <td> 
!!                     An array of function nodal values located at the native interpolation nodes.
!!   <tr> <td> in <th> sE(1:2) <td> REAL(prec) <td> 
!!                     The (2-D) location where you want to interpolate to.
!!   <tr> <td> out <th> fAtSE <td> REAL(prec) <td> The interpolant evaluated at sE.
!!  </table>  
!!
! ================================================================================================ ! 
!>@}
  FUNCTION Interpolate_2D_Lagrange( myPoly, f, sE ) RESULT( interpF )  
  
   IMPLICIT NONE
   CLASS(Lagrange) :: myPoly
   REAL(prec)      :: sE(1:2)
   REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N)
   REAL(prec)      :: interpF
   ! LOCAL
   REAL(prec) :: w(0:myPoly % N), s(0:myPoly % N), ls(0:myPoly % N), lp(0:myPoly % N)
   REAL(prec) :: f1d(0:myPoly % N),fT(0:myPoly % N, 0:myPoly % N)
   INTEGER    :: N

      N = myPoly % N
      s = myPoly % s
      w = myPoly % bWs
     
      ls = LagrangePolynomials( s, w, sE(1), N ) 
      lp = LagrangePolynomials( s, w, sE(2), N ) 
      fT = TRANSPOSE( f )
      f1d = MATMUL( fT, ls )
      interpF = DOT_PRODUCT( f1d, lp )

 END FUNCTION Interpolate_2D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! Function Interpolate_3D 
! 
!> \fn Interpolate_3D_Lagrange 
!!  Interpolates an array of data onto a desired location using Lagrange interpolation.
!!
!!  This function depends on <BR>
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref interpolation
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: f(0:this % N, 0:this % N, 0:this % N) <BR>
!! <B>REAL</B>(prec)     :: sE(1:3) <BR>
!! <B>REAL</B>(prec)     :: fAtSE <BR>
!!         .... <BR>
!!     fAtSE = this % Interpolate( f, sE ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange_1D) <td> 
!!                     A previously constructed Lagrange_1D data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N,0:myPoly % N) <td> REAL(prec) <td> 
!!                     An array of function nodal values located at the native interpolation nodes.
!!   <tr> <td> in <th> sE(1:3) <td> REAL(prec) <td> 
!!                     The (3-D) location where you want to interpolate to.
!!   <tr> <td> out <th> interpF <td> REAL(prec) <td> The interpolant evaluated at sE.
!!  </table>  
!!
! ================================================================================================ ! 
!>@}
  FUNCTION Interpolate_3D_Lagrange( myPoly, f, sE ) RESULT( interpF )  
  
   IMPLICIT NONE
   CLASS(Lagrange) :: myPoly
   REAL(prec)      :: sE(1:3)
   REAL(prec)      :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N)
   REAL(prec)      :: interpF
   ! LOCAL
   REAL(prec) :: fjk, fk, locF
   REAL(prec) :: ls(0:myPoly % N)
   REAL(prec) :: lp(0:myPoly % N)
   REAL(prec) :: lq(0:myPoly % N)
   INTEGER    ::  i, j, k, e

      ls = LagrangePolynomials( myPoly % s, myPoly % bWs, sE(1), myPoly % N ) 
      lp = LagrangePolynomials( myPoly % s, myPoly % bWs, sE(2), myPoly % N )
      lq = LagrangePolynomials( myPoly % s, myPoly % bWs, sE(3), myPoly % N ) 
      
      locF = 0.0_prec
      DO k = 0, myPoly % N
      
         fk = 0.0_prec
         DO j = 0, myPoly % N
         
            fjk = 0.0_prec
            DO i = 0, myPoly % N
               fjk = fjk + f(i,j,k)*ls(i)
            ENDDO
            
            fk = fk + fjk*lp(j)
         ENDDO
         
         locF = locF + fk*lq(k)
      ENDDO
      
      interpF = locF
      
 END FUNCTION Interpolate_3D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! Function Differentiate 
! 
!> \fn Differentiate_1D_Lagrange 
!!  Differentiates the Lagrange interpolant at a specified point.
!!
!!  This function depends on <BR>
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref differentiation
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: f(0:this % N) <BR>
!! <B>REAL</B>(prec)     :: sE <BR>
!! <B>REAL</B>(prec)     :: dfds <BR>
!!         .... <BR>
!!     dfds = this % Differentiate_1D( f, sE ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N) <td> REAL(prec) <td> 
!!                     An array of function nodal values located at the native interpolation nodes.
!!   <tr> <td> in <th> sE <td> REAL(prec) <td> 
!!                     The location where you want to interpolate to.
!!   <tr> <td> out <th> dInFdx <td> REAL(prec) <td> The derivative of the interpolant at sE.
!!  </table>  
!!
! ================================================================================================ ! 
!>@}
 FUNCTION Differentiate_1D_Lagrange( myPoly, f, sE ) RESULT( dInFdx )  

   IMPLICIT NONE
   CLASS(Lagrange), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)      :: sE
   REAL(prec), INTENT(in)      :: f(0:myPoly % N)
   REAL(prec)                  :: dInFdx
   ! LOCAL
   REAL(prec) :: s(0:myPoly % N), w(0:myPoly % N)
   INTEGER    :: N

      N = myPoly % N
      s = myPoly % s
      w = myPoly % bWs
     
      dInFdx = Differentiation( s, w, f, N, sE )

 END FUNCTION Differentiate_1D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! Function Differentiate_2D 
! 
!> \fn Differentiate_2D_Lagrange 
!!  Differentiates the 2-D Lagrange interpolant at a specified point.
!!
!!  This function depends on <BR>
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref differentiation
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: f(0:this % N, 0:this % N) <BR>
!! <B>REAL</B>(prec)     :: sE(1:2) <BR>
!! <B>REAL</B>(prec)     :: dfds <BR>
!!         .... <BR>
!!     dfds = this % Differentiate_2D( f, sE ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N,0:myPoly % N) <td> REAL(prec) <td> 
!!                     An array of function nodal values located at the native interpolation nodes.
!!   <tr> <td> in <th> sE(1:2) <td> REAL(prec) <td> 
!!                     The (2-D) location where you want to interpolate to.
!!   <tr> <td> out <th> dInFdx(1:2) <td> REAL(prec) <td> The derivative of the interpolant at sE.
!!  </table>  
!!
! ================================================================================================ ! 
!>@}
 FUNCTION Differentiate_2D_Lagrange( myPoly, f, sE ) RESULT( dInFdx )  

   IMPLICIT NONE
   CLASS(Lagrange), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)      :: sE(1:2)
   REAL(prec), INTENT(in)      :: f(0:myPoly % N, 0:myPoly % N)
   REAL(prec)                  :: dInFdx(1:2)
   ! LOCAL
   REAL(prec) :: s(0:myPoly % N), w(0:myPoly % N), ls(0:myPoly % N), lp(0:myPoly % N)
   REAL(prec) :: f1d(0:myPoly % N)
   REAL(prec) :: dfds, dfdp
   INTEGER    :: i, N

      N = myPoly % N
      s = myPoly % s
      w = myPoly % bWs
     
      ls = LagrangePolynomials( s, w, sE(1), N ) 
      lp = LagrangePolynomials( s, w, sE(2), N ) 

      dfds = ZERO
      dfdp = ZERO

      DO i = 0, N ! Loop over p-points
         f1d = f(0:N,i)
         ! Evaluate s-derivative
         dfds = dfds + ( Differentiation( s, w, f1d, N, sE(1) ) )*lp(i)
      ENDDO
      DO i = 0, N ! Loop over p-points
         f1d = f(i,0:N)
         ! Evaluate s-derivative
         dfdp = dfdp + ( Differentiation( s, w, f1d, N, sE(2) ) )*ls(i)
      ENDDO
     
      dInFdx(1) = dfds
      dInFdx(2) = dfdp

 END FUNCTION Differentiate_2D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! Function Differentiate_3D 
! 
!> \fn Differentiate_3D_Lagrange 
!!  Differentiates the 3-D Lagrange interpolant at a specified point.
!!
!!  This function depends on <BR>
!!   Module \ref InterpolatioNupportRoutines.f90 : \ref differentiation
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: f(0:this % N, 0:this % N, 0:this % N) <BR>
!! <B>REAL</B>(prec)     :: sE(1:3) <BR>
!! <B>REAL</B>(prec)     :: dfds <BR>
!!         .... <BR>
!!     dfds = this % Differentiate_3D( f, sE ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:this % N, 0:this % N, 0:this % N) <td> REAL(prec) <td> 
!!                     An array of function nodal values located at the native interpolation nodes.
!!   <tr> <td> in <th> sE(1:3) <td> REAL(prec) <td> 
!!                     The (2-D) location where you want to interpolate to.
!!   <tr> <td> out <th> dInFdx(1:3) <td> REAL(prec) <td> The derivative of the interpolant at sE.
!!  </table>  
!!
! ================================================================================================ ! 
!>@}
 FUNCTION Differentiate_3D_Lagrange( myPoly, f, sE ) RESULT( dInFdx )  

   IMPLICIT NONE
   CLASS(Lagrange), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)      :: sE(1:3)
   REAL(prec), INTENT(in)      :: f(0:myPoly % N, 0:myPoly % N, 0:myPoly % N)
   REAL(prec)                  :: dInFdx(1:3)
   ! LOCAL
   REAL(prec) :: s(0:myPoly % N), w(0:myPoly % N)
   REAL(prec) :: ls(0:myPoly % N), lp(0:myPoly % N), lq(0:myPoly % N)
   REAL(prec) :: f1d(0:myPoly % N)
   REAL(prec) :: dfds, dfdp, dfdq
   INTEGER    :: i, j, N

      N = myPoly % N
      s = myPoly % s
      w = myPoly % bWs
     
      ls = LagrangePolynomials( s, w, sE(1), N ) 
      lp = LagrangePolynomials( s, w, sE(2), N ) 
      lq = LagrangePolynomials( s, w, sE(3), N ) 

      dfds = ZERO
      dfdp = ZERO
      dfdq = ZERO
      
      DO j = 0, N
         DO i = 0, N 
            f1d = f(0:N,i,j)
            ! Evaluate s-derivative
            dfds = dfds + ( Differentiation( s, w, f1d, N, sE(1) ) )*lp(i)*lq(j)
         ENDDO
      ENDDO
      DO j = 0, N
         DO i = 0, N 
            f1d = f(i,0:N,j)
            ! Evaluate p-derivative
            dfdp = dfdp + ( Differentiation( s, w, f1d, N, sE(2) ) )*ls(i)*lq(j)
         ENDDO
      ENDDO
      DO j = 0, N
         DO i = 0, N 
            f1d = f(i,j,0:N)
            ! Evaluate q-derivative
            dfdq = dfdq + ( Differentiation( s, w, f1d, N, sE(2) ) )*ls(i)*lp(j)
         ENDDO
      ENDDO
     
      dInFdx(1) = dfds
      dInFdx(2) = dfdp
      dInFdx(3) = dfdq

 END FUNCTION Differentiate_3D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyInterpolationMatrix_1D 
! 
!> \fn ApplyInterpolationMatrix_1D_Lagrange 
!! Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!! 
!! As described in calculateinterpolationmatrix_lagrange_1d, 
!! we can write the operations of interpolating data from one set of points to another (in this 
!! case from "s" to "so") as
!!            \f[ f_j = \sum_{i=0}^N f_i l_i(\xi_j), \hspace{2mm} j=0,1,2,...,M
!!            \f]
!! where \f$ l_i(\xi) \f$ are the Lagrange interpolating polynomials at the 
!! \f$ \lbrace \xi_i \rbrace_{i=0}^N \f$ nodes. This routine performs the matrix-multiply that
!! maps an array of nodal values from the native interpolation nodes to the target nodes 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N) <BR>
!! <B>REAL</B>(prec)     :: ftarget(0:this % M) <BR>
!!         .... <BR>
!!     ftarget = this % ApplyInterpolationMatrix_1D( fnative ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N) <td> REAL(prec) <td>
!!                     Array of function nodal values at the native interpolation nodes.
!!   <tr> <td> out <th> fNew(0:myPoly % M) <td> REAL(prec) <td> 
!!                     Array of function nodal values at the target interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ApplyInterpolationMatrix_1D_Lagrange( myPoly, f ) RESULT( fNew )  

  IMPLICIT NONE
  CLASS(Lagrange) :: myPoly
  REAL(prec)      :: f(0:myPoly % N)
  REAL(prec)      :: fNew(0:myPoly % M)

     fNew = MATMUL( myPoly % Ts, f )

 END FUNCTION ApplyInterpolationMatrix_1D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyInterpolationMatrix_2D 
! 
!> \fn ApplyInterpolationMatrix_2D_Lagrange  
!! Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!! 
!! As described in \ref calculateinterpolationmatrix_Lagrange, 
!! we can write the operations of interpolating data from one set of points to another (in this 
!! case from "s" to "so") as a sequence of matrix-matrix multiplications
!! \f[ 
!!       \tilde{f} = Tf 
!! \f]
!! \f[ 
!!       f_{target} = \tilde{f} T^T 
!! \f]
!!
!! This routine performs the matrix-multiplications that map an array of nodal values from the 
!! native interpolation nodes to the target nodes. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N,0:this % N) <BR>
!! <B>REAL</B>(prec)     :: ftarget(0:this % M,0:this % M) <BR>
!!         .... <BR>
!!     ftarget = this % ApplyInterpolationMatrix_2D( fnative ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N,0:this % N) <td> REAL(prec) <td>
!!                     2-D Array of function nodal values at the native interpolation nodes.
!!   <tr> <td> out <th> fNew(0:myPoly % M,0:this % M) <td> REAL(prec) <td> 
!!                     2-D Array of function nodal values at the target interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ApplyInterpolationMatrix_2D_Lagrange( myPoly, f ) RESULT( fNew )  

  IMPLICIT NONE
  CLASS(Lagrange) :: myPoly
  REAL(prec)         :: f(0:myPoly % N,0:myPoly % N)
  REAL(prec)         :: fNew(0:myPoly % M,0:myPoly % M)
  ! Local
  REAL(prec) :: fInt(0:myPoly % M, 0:myPoly % N)

     fInt = MATMUL( myPoly % Ts, f )
     fNew = MATMUL( fInt, myPoly % Tp )

 END FUNCTION ApplyInterpolationMatrix_2D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyInterpolationMatrix_3D 
! 
!> \fn ApplyInterpolationMatrix_3D_Lagrange  
!! Interpolates an array of nodal values at the native nodes to nodal values at the target nodes.
!! 
!! As described in \ref calculateinterpolationmatrix_lagrange_3D, 
!! we can write the operations of interpolating data from one set of points to another (in this 
!! case from "s" to "so") as a sequence of matrix-matrix multiplications. First a sequence of 2-D
!! interpolations are performed,
!! \f[ 
!!       \tilde{f}_k = Tf_kT ^T
!! \f]
!! Then, by collapsing the first two dimensions of the data, a single matrix-matrix multiplication
!! followed by an unpacking of the first two dimensions results in the data interpolated onto the 
!! target nodes
!! \f[ 
!!       f_{target} = UNPACK( PACK( \tilde{f} T^T ) ) 
!! \f]
!!
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N,0:this % N,0:this % N) <BR>
!! <B>REAL</B>(prec)     :: ftarget(0:this % M,0:this % M,0:this % M) <BR>
!!         .... <BR>
!!     ftarget = this % ApplyInterpolationMatrix_3D( fnative ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N,0:this % N,0:this % N) <td> REAL(prec) <td>
!!                     3-D Array of function nodal values at the native interpolation nodes.
!!   <tr> <td> out <th> fNew(0:myPoly % M,0:this % M,0:this % M) <td> REAL(prec) <td> 
!!                     3-D Array of function nodal values at the target interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ApplyInterpolationMatrices_3D_Lagrange( myPoly, f, nElems ) RESULT( fNew )  

   IMPLICIT NONE
   CLASS(Lagrange) :: myPoly
   INTEGER         :: nElems
   REAL(prec)      :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nElems)
   REAL(prec)      :: fNew(0:myPoly % M,0:myPoly % M,0:myPoly % M,1:nElems)
   ! Local
   INTEGER :: i, j, k, m, n, p, iEl
   REAL(prec) :: floc(0:myPoly % N,0:myPoly % N,0:myPoly % N)
   REAL(prec) :: fm, fmn, fmnp
   
      DO iEl = 1, nElems
         
         floc = f(:,:,:,iEl)
         
         DO p = 0, myPoly % M
            DO n = 0, myPoly % M
               DO m = 0, myPoly % M
               
                  fmnp = 0.0_prec
                  DO k = 0, myPoly % N
                  
                     fmn = 0.0_prec
                     DO j = 0, myPoly % N
                     
                        fm = 0.0_prec
                        DO i = 0, myPoly % N
                           fm = fm + floc(i,j,k)*myPoly % Tp(i,m)
                        ENDDO
                        
                        fmn = fmn + fm*myPoly % Tp(j,n)
                     ENDDO
                     
                     fmnp = fmnp + fmn*myPoly % Tp(k,p)
                  ENDDO
                  
                  fnew(m,n,p,iEl) = fmnp
               ENDDO
            ENDDO
         ENDDO
         
      ENDDO

 END FUNCTION ApplyInterpolationMatrices_3D_Lagrange
!
 FUNCTION ApplyInterpolationMatrix_3D_Lagrange( myPoly, f ) RESULT( fNew )  

   IMPLICIT NONE
   CLASS(Lagrange) :: myPoly
   REAL(prec)      :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N)
   REAL(prec)      :: fNew(0:myPoly % M,0:myPoly % M,0:myPoly % M)
   ! Local
   INTEGER :: i, j, k, m, n, p
   REAL(prec) :: fm, fmn, fmnp
   
   !   DO iEl = 1, nElems
         DO p = 0, myPoly % M
            DO n = 0, myPoly % M
               DO m = 0, myPoly % M
               
                  fmnp = 0.0_prec
                  DO k = 0, myPoly % N
                  
                     fmn = 0.0_prec
                     DO j = 0, myPoly % N
                     
                        fm = 0.0_prec
                        DO i = 0, myPoly % N
                           fm = fm + f(i,j,k)*myPoly % Tp(i,m)
                        ENDDO
                        
                        fmn = fmn + fm*myPoly % Tp(j,n)
                     ENDDO
                     
                     fmnp = fmnp + fmn*myPoly % Tp(k,p)
                  ENDDO
                  
                  fnew(m,n,p) = fmnp
               ENDDO
            ENDDO
         ENDDO
         
     ! ENDDO
!
 END FUNCTION ApplyInterpolationMatrix_3D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyDerivativeMatrix_1D 
! 
!> \fn ApplyDerivativeMatrix_1D_Lagrange 
!! Calculates the derivative of the Lagrange interpolant given a set of nodal function values at
!! the native interpolation nodes
!! 
!! As described in calculatederivativematrix_lagrange_1d, 
!! given nodal values of an interpolant, the derivative can be estimated at the interpolation 
!! nodes using the summation
!!      \f[ \lbrace f' \rbrace_{j=0}^N = \sum_{i}^N( f_i l'_i(\xi_j) )
!!      \f]
!! Where \f$ l'_i(\xi_j) \f$ is the derivative of the \f$ i^{th} \f$ Lagrange interpolating 
!! polynomial at the \f$ j^{th} \f$ native interpolation node. This routine performs the 
!! matrix-multiply that results in an estimate of the derivative of a function whose nodal values
!! are the  \f$ f_i \f$ at the native interpolation nodes.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N) <BR>
!! <B>REAL</B>(prec)     :: dfds(0:this % N) <BR>
!!         .... <BR>
!!     dfds = this % ApplyDerivativeMatrix_1D( fnative ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N) <td> REAL(prec) <td>
!!                     Array of function nodal values at the native interpolation nodes.
!!   <tr> <td> out <th> derF(0:myPoly % M) <td> REAL(prec) <td> 
!!                     Array of estimated derivative values at the native interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ApplyDerivativeMatrix_1D_Lagrange( myPoly, f ) RESULT( derF )  

  IMPLICIT NONE
  CLASS(Lagrange) :: myPoly
  REAL(prec)      :: f(0:myPoly % N)
  REAL(prec)      :: derF(0:myPoly % N)

     derF = MATMUL( myPoly % D, f )

 END FUNCTION ApplyDerivativeMatrix_1D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyDerivativeMatrix_2D 
! 
!> \fn ApplyDerivativeMatrix_2D_Lagrange  
!! Calculates the derivative of the Lagrange interpolant, in each computational direction, given a 
!! set of nodal function values at the native interpolation nodes.
!!
!! As described in \ref calculatederivativematrix_Lagrange, 
!! we can write the derivative calculations as a set of two matrix-matrix products
!! \f[ 
!!       \frac{\partial f}{\partial s} = D f  
!! \f]
!! \f[ 
!!       \frac{\partial f}{\partial p} = f D^T 
!! \f]
!! This routine performs the matrix-multiplications that result in the derivative of the interpolant
!! at the native interpolation nodes. This serves as an estimate of the derivative of the underlying
!! function.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N,0:this % N) <BR>
!! <B>REAL</B>(prec)     :: derF(0:this % N,0:this % N,1:2) <BR>
!!         .... <BR>
!!     derF = this % ApplyDerivativeMatrix_2D( fnative ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N,0:this % N) <td> REAL(prec) <td>
!!                     2-D Array of function nodal values at the native Derivative nodes.
!!   <tr> <td> out <th> derF(0:myPoly % N,0:this % N,1:2) <td> REAL(prec) <td> 
!!                     3-D Array containing the derivative (in each direction) of the interpolant at
!!                     the native interpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ApplyDerivativeMatrix_2D_Lagrange( myPoly, f ) RESULT( derF )  

  IMPLICIT NONE
  CLASS(Lagrange) :: myPoly
  REAL(prec)      :: f(0:myPoly % N,0:myPoly % N)
  REAL(prec)      :: derF(0:myPoly % N,0:myPoly % N, 1:2)

     derF(:,:,1) = MATMUL( myPoly % D, f )
     derF(:,:,2) = MATMUL( f, myPoly % DTr )

 END FUNCTION ApplyDerivativeMatrix_2D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R ApplyDerivativeMatrix_3D 
! 
!> \fn ApplyDerivativeMatrix_3D_Lagrange  
!! Calculates the derivative of the Lagrange interpolant, in each computational direction, given a 
!! set of nodal function values at the native interpolation nodes.
!!
!! As described in \ref calculatederivativematrix_lagrange_3D, to compute the derivative in a single
!! computational direction, the other array dimensions can be collapsed enabling a matrix-matrix 
!! product. After computing the product, the result can be unpacked into a 3-D array.
!! See the <a href="http://www.numeriphy.com/SELF/doc/SELF-TechnicalDocumentation.pdf"> SELF Technical Documentation </a> for more details.
!!
!! This routine performs the matrix-multiplications that result in the derivative of the interpolant
!! at the native interpolation nodes. This serves as an estimate of the derivative of the underlying
!! function.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N,0:this % N,0:this % N) <BR>
!! <B>REAL</B>(prec)     :: derF(0:this % N,0:this % N,0:this % N,1:3) <BR>
!!         .... <BR>
!!     derF = this % ApplyDerivativeMatrix_3D( fnative ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> f(0:myPoly % N,0:this % N,0:this % N) <td> REAL(prec) <td>
!!                     3-D Array of function nodal values at the native Derivative nodes.
!!   <tr> <td> out <th> derF(0:myPoly % N,0:this % N,1:2) <td> REAL(prec) <td> 
!!                     4-D Array containing the derivative of the interpolant (in each direction) at
!!                     the nativeinterpolation nodes.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ApplyDerivativeMatrices_3D_Lagrange( myPoly, f, nElems ) RESULT( derF )  

   IMPLICIT NONE
   CLASS(Lagrange) :: myPoly
   INTEGER         :: nElems
   REAL(prec)      :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nElems)
   REAL(prec)      :: derF(1:3,0:myPoly % N,0:myPoly % N,0:myPoly % N,1:nElems)
  ! Local
  INTEGER :: i, j, k, iEl, ii
  
      ! Typically, N will vary between 4 and 10
      ! and nElems can range up to 150,000 or more
      
      DO iEl = 1, nElems  ! << Largest Loop O( 100,000 )
         DO k = 0, myPoly % N ! << O( 10 )
            DO j = 0, myPoly % N ! << O( 10 )
               DO i = 0, myPoly % N ! << O( 10 )
               
                  derf(1:3,i,j,k,iEl) = 0.0_prec
                  DO ii = 0, myPoly % N ! << O( 10 ), Reduction loop ( vector-vector operation 
                                        ! for a single i,j,k )
                     derf(1,i,j,k,iEl) = derf(1,i,j,k,iEl) + myPoly % DTr(ii,i)*f(ii,j,k,iEl)
                     derf(2,i,j,k,iEl) = derf(2,i,j,k,iEl) + myPoly % DTr(ii,j)*f(i,ii,k,iEl)
                     derf(3,i,j,k,iEl) = derf(3,i,j,k,iEl) + myPoly % DTr(ii,k)*f(i,j,ii,iEl)
                  ENDDO
               ENDDO
               
            ENDDO
         ENDDO
      ENDDO
      
 END FUNCTION ApplyDerivativeMatrices_3D_Lagrange
!
 FUNCTION ApplyDerivativeMatrix_3D_Lagrange( myPoly, f ) RESULT( derF )  

   IMPLICIT NONE
   CLASS(Lagrange) :: myPoly
   REAL(prec)      :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N)
   REAL(prec)      :: derF(0:myPoly % N,0:myPoly % N,0:myPoly % N,1:3)
   ! Local
   INTEGER :: i, j, k,  ii
  
         DO k = 0, myPoly % N ! << O( 10 )
            DO j = 0, myPoly % N ! << O( 10 )
               DO i = 0, myPoly % N ! << O( 10 )
               
                  derf(i,j,k,1:3) = 0.0_prec
                  DO ii = 0, myPoly % N ! << O( 10 ), Reduction loop ( vector-vector operation 
                                        ! for a single i,j,k )
                     derf(i,j,k,1) = derf(i,j,k,1) + myPoly % DTr(ii,i)*f(ii,j,k)
                     derf(i,j,k,2) = derf(i,j,k,2) + myPoly % DTr(ii,j)*f(i,ii,k)
                     derf(i,j,k,3) = derf(i,j,k,3) + myPoly % DTr(ii,k)*f(i,j,ii)
                  ENDDO
               ENDDO
               
            ENDDO
         ENDDO
      
 END FUNCTION ApplyDerivativeMatrix_3D_Lagrange
!
!
!==================================================================================================!
!--------------------------------- File I/O Routines ----------------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R WriteTecplot_1D 
! 
!> \fn WriteTecplot_1D_Lagrange 
!! Writes an ASCII tecplot file for 1-D data given a set of function nodal values at the native
!! interpolation nodes.
!! 
!! Passing in filename ="example" results in a file called "example.curve" that conforms to the 
!! 1-D tecplot file format. This file can be viewed in any data-visualization software that can 
!! read tecplot files.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N) <BR>
!! <B>CHARACTER</B>(len) :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % WriteTecplot_1D( fnative, filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange_1D data-structure.
!!   <tr> <td> in <th> fnative(0:myPoly % N) <td> REAL(prec) <td>
!!                     An array of function nodal values at the native interpolation nodes
!!   <tr> <td> in <th> filename <td> CHARACTER(*) <td>
!!                     Name of the file where the native interpolation nodes and function <BR>
!!                     nodal values will be written.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteTecplot_1D_Lagrange( myPoly, f, filename )

   IMPLICIT NONE
   CLASS( Lagrange ), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)        :: f(0:myPoly % N)
   CHARACTER(*), INTENT(in)      :: filename
   ! Local
   INTEGER :: fUnit, iS

      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = TRIM(filename)//'.curve', &
            FORM   = 'FORMATTED', &
            STATUS = 'REPLACE' )

      WRITE( fUnit, * )'#f'
      DO iS = 0, myPoly % N
         WRITE( fUnit, * ) myPoly % s(iS), f(iS)
      ENDDO

      CLOSE( fUnit )

 END SUBROUTINE WriteTecplot_1D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R WriteTecplot 
! 
!> \fn WriteTecplot_2D_Lagrange  
!! Writes an ASCII tecplot file for 2-D data given a set of function nodal values at the native
!! interpolation nodes.
!! 
!! Passing in filename ="example" results in a file called "example.tec" that conforms to the 
!! 2-D tecplot FEM file format. This file can be viewed in any data-visualization software that can 
!! read tecplot files.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N,0:this % N) <BR>
!! <B>CHARACTER</B>(len) :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % WriteTecplot_2D( fnative, filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange data-structure.
!!   <tr> <td> in <th> fnative(0:myPoly % N, 0:myPoly % N)<td> REAL(prec) <td>
!!                     An array of function nodal values at the native interpolation nodes
!!   <tr> <td> in <th> filename <td> CHARACTER(*) <td>
!!                     Name of the file where the native interpolation nodes and function <BR>
!!                     nodal values will be written.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteTecplot_2D_Lagrange( myPoly, f, filename )

   IMPLICIT NONE
   CLASS( Lagrange ), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)        :: f(0:myPoly % N,0:myPoly % N)
   CHARACTER(*), INTENT(in)      :: filename
   ! Local
   INTEGER    :: fUnit, i, j, N
   REAL(prec) :: s(0:myPoly % N)
   
      N = myPoly % N
      s = myPoly % N

      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = TRIM(filename)//'.tec', &
            FORM   = 'FORMATTED', &
            STATUS = 'REPLACE' )

      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "f" '
      WRITE(fUnit,*) 'ZONE T="el00", I=',N+1,', J=', N+1,',F=POINT'
      
      DO j = 0, N
         DO i = 0, N
            WRITE( fUnit, * ) s(i), s(j), f(i,j)
         ENDDO
      ENDDO

      CLOSE( fUnit )

 END SUBROUTINE WriteTecplot_2D_Lagrange
!
!> \addtogroup Lagrange_Class
!! @{ 
! ================================================================================================ !
! S/R WriteTecplot 
! 
!> \fn WriteTecplot_Lagrange  
!! Writes an ASCII tecplot file for 3-D data given a set of function nodal values at the native
!! interpolation nodes.
!! 
!! Passing in filename ="example" results in a file called "example.tec" that conforms to the 
!! 3-D tecplot FEM file format. This file can be viewed in any data-visualization software that can 
!! read tecplot files.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Lagrange) :: this <BR>
!! <B>REAL</B>(prec)     :: fnative(0:this % N,0:this % N,0:this % N) <BR>
!! <B>CHARACTER</B>(len) :: filename <BR>
!!         .... <BR>
!!     <B>CALL</B> this % WriteTecplot_3D( fnative, filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myPoly <td> TYPE(Lagrange) <td> 
!!                     A previously constructed Lagrange_1D data-structure.
!!   <tr> <td> in <th> fnative(0:myPoly % N, 0:myPoly % N,0:this % N)<td> REAL(prec) <td>
!!                     An array of function nodal values at the native interpolation nodes
!!   <tr> <td> in <th> filename <td> CHARACTER(*) <td>
!!                     Name of the file where the native interpolation nodes and function <BR>
!!                     nodal values will be written.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteTecplot_3D_Lagrange( myPoly, f, filename )

   IMPLICIT NONE
   CLASS( Lagrange ), INTENT(in) :: myPoly
   REAL(prec), INTENT(in)        :: f(0:myPoly % N,0:myPoly % N,0:myPoly % N)
   CHARACTER(*), INTENT(in)      :: filename
   ! Local
   INTEGER    :: fUnit, i, j, k, N
   REAL(prec) :: s(0:myPoly % N)
   
      N = myPoly % N
      s = myPoly % s

      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = TRIM(filename)//'.tec', &
            FORM   = 'FORMATTED', &
            STATUS = 'REPLACE' )

      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "f" '
      WRITE(fUnit,*) 'ZONE T="el00", I=',N+1,', J=', N+1,', K=', N+1,',F=POINT'
      
      DO k = 0, N
         DO j = 0, N
            DO i = 0, N
               WRITE( fUnit, * ) myPoly % s(i), myPoly % s(j), myPoly % s(k), f(i,j,k)
            ENDDO
         ENDDO
      ENDDO

      CLOSE( fUnit )

 END SUBROUTINE WriteTecplot_3D_Lagrange
!
!
!
END MODULE Lagrange_Class
