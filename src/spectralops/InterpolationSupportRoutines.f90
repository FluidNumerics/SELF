! InterpolationSupportRoutines.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file InterpolationSupportRoutines.f90
!! Contains the \ref InterpolationSupportRoutines module.

!> \defgroup InterpolationSupportRoutines InterpolationSupportRoutines
!! Contains routines from D.A. Kopriva, 2009, "Implementing Spectral Methods for Partial 
!! Differential Equations: Algorithms for Scientists and Engineers", Springer for assisting in the 
!! implementation of Lagrange interpolation, the basis of polynomial Spectral Element Methods.
 
MODULE InterpolationSupportRoutines

! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines

 IMPLICIT NONE

 CONTAINS

!> \addtogroup InterpolationSupportRoutines 
!! @{ 
! ================================================================================================ !
! Function BarycentricWeights
! 
!> \fn BarycentricWeights  
!! Given an array of interpolation nodes, this function computes and returns the barycentric 
!! interpolation weights associate with Lagrange Interpolation. 
!! 
!! From a set of interpolation nodes \f$ \lbrace \xi_j \rbrace_{j=0}^N \f$, the Lagrange 
!! interpolating polynomials are given as
!!     \f[
!!           l_j = \prod_{i=0,\neq j}^N \frac{\xi-\xi_i}{\xi_j-\xi_i}
!!     \f]
!! For efficient interpolation with favorable round-off error, the "barycentric weights" are usually 
!! stored when performing Lagrange interpolation. The barycentric weights are
!!     \f[
!!           w_j = \prod_{i=0,\neq j}^N \frac{1}{\xi_j-\xi_i}
!!     \f]
!! 
!!   This function is from Alg. 30 on pg. 75 of D.A. Kopriva, 2009.
!!
!! <H2> Usage : </H2> 
!! <B>INTEGER</B> :: N <BR>
!! <B>REAL</B>(prec) :: s(0:N), w(0:N) <BR>
!!         .... <BR>
!!     w = BarycentricWeights( s, N ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> Array of interpolation nodes.
!!   <tr> <td> in <th> N <td> INTEGER <td> Degree of Lagrange interpolating polynomials that pass 
!!                                         through the interpolation nodes.
!!   <tr> <td> out <th> w(0:N) <td> REAL(prec) <td> Array of barycentric weights.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION BarycentricWeights( s, N ) RESULT( w )

   IMPLICIT NONE
   INTEGER    :: N
   REAL(prec) :: s(0:N)
   REAL(prec) :: w(0:N)
   ! Local (to this function)
   INTEGER :: i, j
   REAL(prec) :: wlocal(0:N)
   
       wlocal(0:N) = ONE ! initialize the weights to 1

      ! Computes the product w_k = w_k*(s_k - s_j), k /= j
      DO j = 1,N ! loop over the interpolation nodes
         DO i = 0, j-1 ! loop to perform multiplication for weights

            wlocal(i) = wlocal(i)*( s(i) - s(j) )
            wlocal(j) = wlocal(j)*( s(j) - s(i) )

         ENDDO 
      ENDDO 
 
     ! Invert
      DO j = 0, N
        wlocal(j) = ONE/wlocal(j)
      ENDDO 

      w = wlocal

 END FUNCTION BarycentricWeights
!
!> \addtogroup InterpolationSupportRoutines 
!! @{ 
! ================================================================================================ !
! Function LagrangePolynomials
! 
!> \fn LagrangePolynomials  
!! This function returns the value of each Lagrange interpolating polynomial associated with a 
!! set of interpolation nodes and barycentric weights at a specified point.
!! 
!! From a set of interpolation nodes \f$ \lbrace \xi_j \rbrace_{j=0}^N \f$, the Lagrange 
!! interpolating polynomials are given as
!!     \f[
!!           l_j = \prod_{i=0,\neq j}^N \frac{\xi-\xi_i}{\xi_j-\xi_i}
!!     \f] 
!! This function evaluates the Lagrange interpolating polynomials using the "Barycentric Formulation"
!! (Eq. 3.36 of Kopriva (2009), pg. 74 (with \f$ f_j = \delta_{i,j} \f$))
!!     \f[
!!            l_j = \frac{\frac{w_j}{\xi-\xi_i}}{\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f] 
!!
!!   This function is from Alg. 34 on pg. 77 of D.A. Kopriva, 2009.
!!
!! <H2> Usage : </H2> 
!! <B>INTEGER</B>    :: N <BR>
!! <B>REAL</B>(prec) :: s(0:N), w(0:N), f(0:N), lAtS(0:N), sE <BR>
!!         .... <BR>
!!     lAtS = LagrangePolynomials( s, w, f, N, sE )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> Array of interpolation nodes
!!   <tr> <td> in <th> w(0:N) <td> REAL(prec) <td> Array of barycentric weights
!!   <tr> <td> in <th> sE <td> REAL(prec) <td> Location where we want to evaluate the Lagrange 
!!                                             interpolating polynomials
!!   <tr> <td> in <th> N <td> INTEGER <td> Degree of the Lagrange interpolating polynomials.
!!   <tr> <td> out <th> lAtS(0:N) <td> REAL(prec) <td> Array of the Lagrange interpolating 
!!                                                     polynomials evaluated at sE.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
FUNCTION LagrangePolynomials( s, w, sE, N ) RESULT( lAtS )  

   IMPLICIT NONE
   INTEGER    :: N
   REAL(prec) :: s(0:N)
   REAL(prec) :: w(0:N)
   REAL(prec) :: sE
   REAL(prec) :: lAtS(0:N)
   ! LOCAL
   REAL(prec) :: temp1, temp2
   INTEGER    :: j
   LOGICAL    :: xMatchesNode

      xMatchesNode = .FALSE.

      DO j = 0, N
        
         lAtS(j) = ZERO

         IF( AlmostEqual(sE, s(j)) ) THEN
            lAtS(j) = ONE
            xMatchesNode = .TRUE.
         ENDIF 

      ENDDO

      IF( xMatchesNode )THEN 
         RETURN
      ENDIF

      temp1 = ZERO
     
      DO j = 0, N 
         temp2 = w(j)/(sE - s(j))
         lAtS(j) = temp2
         temp1 = temp1 + temp2
      ENDDO 
  
      lAtS = lAtS/temp1 
     

 END FUNCTION LagrangePolynomials
!
!> \addtogroup InterpolationSupportRoutines 
!! @{ 
! ================================================================================================ !
! Function Interpolation
! 
!> \fn Interpolation 
!! Interpolates a discrete array of nodal function values to a specified point using Lagrange interpolation 
!! 
!! From a set of interpolation nodes \f$ \lbrace \xi_j \rbrace_{j=0}^N \f$,  and nodal function
!! values \f$ \lbrace f_j \rbrace_{j=0}^N \f$, Lagrange interpolation in the "barycentric 
!! formulation" (Eq. 3.36 of Kopriva (2009), pg. 74  )
!!     \f[
!!            I_N(f) = \frac{\sum_{i=0}^N f_i\frac{w_i}{\xi-\xi_i}}{\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f] 
!!
!!   This function is from Alg. 31 on pg. 75 of D.A. Kopriva, 2009.
!!
!! <H2> Usage : </H2> 
!! <B>INTEGER</B>    :: N <BR>
!! <B>REAL</B>(prec) :: s(0:N), w(0:N), f(0:N), interpF, sE <BR>
!!         .... <BR>
!!     interpF = Interpolation( s, w, f, N, sE )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> Array of interpolation nodes
!!   <tr> <td> in <th> w(0:N) <td> REAL(prec) <td> Array of barycentric weights
!!   <tr> <td> in <th> sE <td> REAL(prec) <td> Location where we want to evaluate the Lagrange 
!!                                             interpolating polynomials
!!   <tr> <td> in <th> N <td> INTEGER <td> Degree of the Lagrange interpolating polynomials.
!!   <tr> <td> out <th> interpF <td> REAL(prec) <td> Value of the Lagrange interpolant at sE.
!!  </table>   
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION Interpolation( s, w, f, N, sE ) RESULT( interpF )  

   IMPLICIT NONE
   INTEGER    :: N
   REAL(prec) :: s(0:N)
   REAL(prec) :: w(0:N)
   REAL(prec) :: f(0:N)
   REAL(prec) :: sE
   REAL(prec) :: interpF
   ! LOCAL
   REAL(prec) :: num, den, t
   INTEGER    :: j
 
      num = ZERO
      den = ZERO

      DO j = 0, N

         IF( AlmostEqual(sE, s(j)) ) THEN 
           interpF = f(j)
           RETURN
         ELSE 
           
            t = w(j)/(sE - s(j))
            num = num + t*f(j)
            den = den + t

         ENDIF
        
      ENDDO

      interpF = num/den

 END FUNCTION Interpolation
!
!> \addtogroup InterpolationSupportRoutines 
!! @{ 
! ================================================================================================ !
! Function Differentiation 
! 
!> \fn Differentiation  
!!  Evaluates the derivative of the Lagrange interpolating polynomial associated with a set of 
!!  interpolation nodes, barycentric weights, and function values at a specified point.
!! 
!! Differentiation of the Lagrange interpolating polynomial
!!     \f[
!!            I_N(f) = \frac{\sum_{i=0}^N f_i\frac{w_i}{\xi-\xi_i}}{\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f]
!! can be rearranged to give (Eq. 3.45 of Kopriva (2009), pg. 80  ) 
!!     \f[
!!            I'_N(f) = \frac{\sum_{i=0}^N f_i\frac{w_i}{\xi-\xi_i}\frac{I_N(f)-f_i}{\xi-\xi_i}}
!!                           {\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f]
!!
!!   This function is from Alg. 36 on pg. 80 of D.A. Kopriva, 2009.
!!
!! <H2> Usage : </H2> 
!! <B>INTEGER</B>    :: N <BR>
!! <B>REAL</B>(prec) :: s(0:N), w(0:N), f(0:N), dInFdx, sE <BR>
!!         .... <BR>
!!     dInFdx = Differentiation( s, w, f, N, sE )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> Array of interpolation nodes
!!   <tr> <td> in <th> w(0:N) <td> REAL(prec) <td> Array of barycentric weights
!!   <tr> <td> in <th> sE <td> REAL(prec) <td> Location where we want to evaluate the Lagrange 
!!                                             interpolating polynomials
!!   <tr> <td> in <th> N <td> INTEGER <td> Degree of the Lagrange interpolating polynomials.
!!   <tr> <td> out <th> dInFdx <td> REAL(prec) <td> Derivative of the Lagrange interpolant at sE.
!!  </table>   
!!   
! ================================================================================================ ! 
!>@}
FUNCTION Differentiation( s, w, f, N, sE ) RESULT( dInFdx )  

   IMPLICIT NONE
   INTEGER    :: N
   REAL(prec) :: s(0:N)
   REAL(prec) :: w(0:N)
   REAL(prec) :: f(0:N)
   REAL(prec) :: sE
   REAL(prec) :: dInFdx
   ! LOCAL
   REAL(prec) :: num, den, t, p
   INTEGER    :: j, i
   LOGICAL    :: atNode 

      num = ZERO
      atNode = .FALSE.

      DO j = 0, N  

         IF( AlmostEqual(sE, s(j)) ) THEN 
            atNode = .TRUE.
            p = f(j)           
            den = -w(j)
            i = j
         ENDIF  
        
      ENDDO 

      IF( atNode ) THEN 

         DO j = 0, N 
            IF( .NOT.(j == i) ) THEN 
               num = num + w(j)*(p - f(j))/(sE - s(j))
            ENDIF
         ENDDO 

      ELSE !

         den = ZERO
         p = Interpolation( s, w, f, N, sE ) 
        
         DO j = 0, N ! loop over the interpolation s
            t = w(j)/(sE - s(j))
            num = num + t*(p - f(j))/(sE - s(j))
            den = den + t
         ENDDO ! j, loop over the interpolation s

      ENDIF ! conditional, IF we're on an interpolating node

      dInFdx = num/den

 END FUNCTION Differentiation
!
!> \addtogroup InterpolationSupportRoutines 
!! @{ 
! ================================================================================================ !
! Function InterpolationMatrix 
! 
!> \fn InterpolationMatrix  
!! Generates a matrix that maps function nodal values from one set of interpolation nodes to another.
!!
!! Given the interpolation formula
!!  \f[
!!       I_N(f) = \sum_{i=0}^N f_i l_i(\xi)
!!  \f] 
!! we can map \f$ \lbrace f_i \rbrace_{i=0}^N\f$ to \f$ \lbrace \tilde{f}_j \rbrace_{j=0}^M\f$ by
!! computing
!!  \f[
!!       \tilde{f}_j = \sum_{i=0}^N f_i l_i(\xi_j)
!!  \f]
!! Row j, column i of the "interpolation matrix" is 
!!  \f[
!!     T_{j,i} = l_i(\xi_j)
!!  \f]
!!
!!   This function is from Alg. 32 on pg. 76 of D.A. Kopriva, 2009.
!!
!! <H2> Usage : </H2> 
!! <B>INTEGER</B>    :: N, nNew <BR>
!! <B>REAL</B>(prec) :: s(0:N), so(0:nNew), w(0:N), T(0:nNew,0:N) <BR>
!!         .... <BR>
!!     T = InterpolationMatrix( s, w, so, N, nNew)
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> Array of native interpolation nodes
!!   <tr> <td> in <th> w(0:N) <td> REAL(prec) <td> Array of barycentric weights
!!   <tr> <td> in <th> so(0:nNew) <td> REAL(prec) <td> Target interpolation nodes
!!   <tr> <td> in <th> N <td> INTEGER <td> Degree of the Lagrange interpolating polynomials 
!!                                         associated with the native interpolation nodes
!!   <tr> <td> in <th> nNew <td> INTEGER <td> Degree of the Lagrange interpolating polynomials 
!!                                           associated with the target interpolation nodes
!!   <tr> <td> out <th> T(0:nNew,0:N) <td> REAL(prec) <td> Interpolation matrix
!!  </table> 
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION InterpolationMatrix( s, w, so, N, nNew ) RESULT( T )

   IMPLICIT NONE
   INTEGER    :: N, nNew
   REAL(prec) :: s(0:N)
   REAL(prec) :: w(0:N)
   REAL(prec) :: so(0:nNew)
   REAL(prec) :: T(0:nNew,0:N)
   ! LOCAL
   REAL(prec) :: temp1, temp2
   INTEGER    :: row, col
   LOGICAL    :: rowHasMatch 


      DO row = 0, nNew ! loop over the new interpolation nodes ("so")

         rowHasMatch = .FALSE.
       
         DO col = 0, N ! loop over the old interpolation nodes ("s")

            T(row,col) = ZERO
           
            IF( AlmostEqual( so(row), s(col) ) )THEN
               rowHasMatch = .TRUE.
               T(row,col) = ONE
            ENDIF

         ENDDO 

         IF( .NOT.(rowHasMatch) )THEN 

            temp1 = ZERO

            DO col = 0, N ! loop over the old interpolation nodes ("s")         
               temp2 = w(col)/( so(row) - s(col) )
               T(row,col) = temp2
               temp1 = temp1 + temp2
            ENDDO 

            DO col = 0, N 
               T(row,col) = T(row,col)/temp1
            ENDDO

         ENDIF 

      ENDDO


 END FUNCTION InterpolationMatrix
!
!> \addtogroup InterpolationSupportRoutines 
!! @{ 
! ================================================================================================ !
! Function DerivativeMatrix 
! 
!> \fn DerivativeMatrix  
!! Generates a matrix that can be used to approximate derivatives at the interpolation nodes.
!!
!! Differentiation of the Lagrange interpolating polynomial
!!     \f[
!!            I_N(f) = \frac{\sum_{i=0}^N f_i\frac{w_i}{\xi-\xi_i}}{\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f]
!! can be rearranged to give (Eq. 3.46 of Kopriva (2009), pg. 80  ) 
!!     \f[
!!            I'_N(f)|_{\xi_j} = \frac{\sum_{i=0}^N f_i\frac{w_i}{\xi-\xi_i}\frac{f_j-f_i}{\xi-\xi_i}}
!!                           {\sum_{i=0}^N \frac{w_i}{\xi-\xi_i}}
!!     \f]
!! when evaluated at each interpolation node
!!
!!   This function is from Alg. 37 on pg. 82 of D.A. Kopriva, 2009.
!!
!! <H2> Usage : </H2> 
!! <B>INTEGER</B>    :: N  <BR>
!! <B>REAL</B>(prec) :: s(0:N), w(0:N), dMat(0:nNew,0:N) <BR>
!!         .... <BR>
!!     dMat = DerivativeMatrix( s, w, N )
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> Array of interpolation nodes
!!   <tr> <td> in <th> w(0:N) <td> REAL(prec) <td> Array of barycentric weights
!!   <tr> <td> in <th> N <td> INTEGER <td> Degree of the Lagrange interpolating polynomials 
!!                                         associated with the  interpolation nodes
!!   <tr> <td> out <th> dMat(0:N,0:N) <td> REAL(prec) <td> Interpolation matrix
!!  </table> 
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION DerivativeMatrix( w, s, N ) RESULT( dMat )  

  IMPLICIT NONE
  INTEGER    :: N
  REAL(prec) :: w(0:N)
  REAL(prec) :: s(0:N)
  REAL(prec) :: dMat(0:N,0:N)
  ! LOCAL
  INTEGER    :: row, col

      DO row = 0, N ! loop over the interpolation s
         dMat(row,row) = ZERO
         DO col = 0, N ! loop over the interpolation s (again)
           
            IF( .NOT. (col == row) )THEN
               dMat(row,col) = w(col)/( w(row)*( s(row) - s(col) ) )
               dMat(row,row) = dMat(row,row) - dMat(row,col)
            ENDIF
        
         ENDDO 
      ENDDO 


 END FUNCTION DerivativeMatrix

END MODULE InterpolationSupportRoutines
