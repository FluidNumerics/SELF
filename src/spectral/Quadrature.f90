! Quadrature.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file Quadrature.f90
!! Contains the \ref Quadrature module.


!> \defgroup Quadrature Quadrature 
!! Contains routines from D.A. Kopriva, 2009, "Implementing Spectral Methods for Partial 
!! Differential Equations: Algorithms for Scientists and Engineers", Springer.
!!
!! Routines are defined for computing Legendre and Chebyshev Gauss and Gauss-Lobatto
!! quadrature nodes and weights.
!!
!! Only PUBLIC routines are documented here. However, references are given below for the PRIVATE 
!! routines.
!!  <table> 
!!   <tr> <th> ChebyshevGauss <td> Alg. 26 on pg. 67 of D.A. Kopriva, 2009. 
!!   <tr> <th> ChebyshevGaussLobatto <td> Alg. 27 on pg. 68 of D.A. Kopriva, 2009. 
!!   <tr> <th> LegendreGauss <td> Alg. 23 on pg. 64 of D.A. Kopriva, 2009. 
!!   <tr> <th> LegendreGaussLobatto <td> Alg. 25 on pg. 66 of D.A. Kopriva, 2009. 
!!   <tr> <th> LegendreQandL <td> Alg. 24 on pg. 65 of D.A. Kopriva, 2009.
!!   <tr> <th> LegendrePolynomial <td> Alg. 22 on pg. 63 of D.A. Kopriva, 2009. 
!!  </table>

 
MODULE Quadrature

! src/common/
USE ModelPrecision
USE ConstantsDictionary


IMPLICIT NONE

 PUBLIC  :: ChebyshevQuadrature, LegendreQuadrature
 PRIVATE :: ChebyshevGauss, ChebyshevGaussLobatto, LegendreGauss, LegendreGaussLobatto, &
            LegendreGaussRadau, LegendreQandL!, LegendrePolynomial

 CONTAINS
!
! ================================================================================================ !
! -------------------------------------- Chebyshev ----------------------------------------------- !
! ================================================================================================ !
!
!> \addtogroup Quadrature 
!! @{ 
! ================================================================================================ !
! S/R ChebyshevQuadrature 
! 
!> \fn ChebyshevQuadrature  
!! Returns the specified Chebyshev quadrature nodes and integration weights. 
!! 
!! Given a polynomial degree, and quadrature type (Gauss or Gauss Lobatto), this subroutine manages
!! the calls to underlying private routines to generate the desired Chebyshev quadrature. 
!! 
!! <H2> Usage : </H2> 
!! <B>INTEGER</B> :: N, quadType <BR>
!! <B>REAL</B>(prec) :: nodes(0:N), weights(0:N) <BR>
!!         .... <BR>
!!     <B>CALL</B> ChebyshevQuadrature( N, quadType, nodes, weights ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> N <td> INTEGER <td> Degree of the quadrature
!!   <tr> <td> in <th> quadType <td> INTEGER <td> Flag specifying the quadrature type. Can be set
!!                                                to GAUSS or GAUSS_LOBATTO (See \ref ModelFlags.f90 )
!!   <tr> <td> out <th> nodes(0:N) <td> REAL(prec) <td> Array of quadrature nodes
!!   <tr> <td> out <th> weights(0:N) <td> REAL(prec) <td> Array of quadrature weights
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ChebyshevQuadrature( N, quadType, nodes, weights )

   IMPLICIT NONE
   INTEGER, INTENT(in)     :: N
   REAL(prec), INTENT(out) :: nodes(0:N)
   REAL(prec), INTENT(out) :: weights(0:N)
   INTEGER, INTENT(in)     :: QuadType
   
   
      ! Compute the quadrature nodes and weights
      IF( QuadType  == GAUSS_LOBATTO )then ! Gauss Lobatto quadrature
         CALL ChebyshevGaussLobatto( N, nodes, weights )
      ELSEIF( QuadType == GAUSS )then  ! Gauss Quadrature
         CALL ChebyshevGauss( N, nodes, weights )
      ELSE
         PRINT*,'Module Chebyshev.f90 : S/R GenerateChebyshevQuadrature: Invalid form. Stopping'
         STOP
      ENDIF

 END SUBROUTINE ChebyshevQuadrature
!
!
!
 SUBROUTINE ChebyshevGauss(N, nodes, weights)  
 ! S/R ChebyshevGauss
 !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 67
 !   Algorithm 26
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER       :: N
   REAL(prec)    :: nodes(0:N)
   REAL(prec)    :: weights(0:N)
   ! LOCAL
   INTEGER    :: j ! Loop counter 
   REAL(prec) :: nReal, jReal, den


      nReal = REAL(N, prec)
      den = nReal + ONE

      DO j = 0, N
         jReal = REAL(j, prec)
         weights(j) = pi/den
         nodes(j) = -cos( 0.5_prec*(2.0_prec*jReal + ONE)*weights(j) )
      ENDDO

 END SUBROUTINE ChebyshevGauss
!
!
!
 SUBROUTINE ChebyshevGaussLobatto(N, nodes, weights)  
 ! S/R ChebyshevGaussLobatto
 !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 68
 !   Algorithm 27
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER       :: N
   REAL(prec)    :: nodes(0:N)
   REAL(prec)    :: weights(0:N)
   ! LOCAL
   INTEGER    :: j ! Loop counter 
   REAL(prec) :: nReal, jReal

      nReal = REAL(N, prec)

      DO j = 0, N
         jReal = REAL(j, prec)
         weights(j) = pi/nReal
         nodes(j) = -cos( jReal*weights(j) )
      ENDDO

      weights(0)  = weights(0)*HALF
      weights(N) = weights(N)*HALF 

 END SUBROUTINE ChebyshevGaussLobatto
!
! ================================================================================================ !
! ---------------------------------------- Legendre ---------------------------------------------- !
! ================================================================================================ !
!
!> \addtogroup Quadrature 
!! @{ 
! ================================================================================================ !
! S/R LegendreQuadrature 
! 
!> \fn LegendreQuadrature  
!! Returns the specified Legendre quadrature nodes and integration weights. 
!! 
!! Given a polynomial degree, and quadrature type (Gauss or Gauss Lobatto), this subroutine manages
!! the calls to underlying private routines to generate the desired Legendre quadrature. 
!! 
!! <H2> Usage : </H2> 
!! <B>INTEGER</B> :: N, quadType <BR>
!! <B>REAL</B>(prec) :: nodes(0:N), weights(0:N) <BR>
!!         .... <BR>
!!     <B>CALL</B> LegendreQuadrature( N, quadType, nodes, weights ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> N <td> INTEGER <td> Degree of the quadrature
!!   <tr> <td> in <th> quadType <td> INTEGER <td> Flag specifying the quadrature type. Can be set
!!                                                to GAUSS or GAUSS_LOBATTO (See \ref ModelFlags.f90 )
!!   <tr> <td> out <th> nodes(0:N) <td> REAL(prec) <td> Array of quadrature nodes
!!   <tr> <td> out <th> weights(0:N) <td> REAL(prec) <td> Array of quadrature weights
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE LegendreQuadrature( N, nodes, weights, QuadType )

   IMPLICIT NONE
   INTEGER, INTENT(in)     :: N
   REAL(prec), INTENT(out) :: nodes(0:N)
   REAL(prec), INTENT(out) :: weights(0:N)
   INTEGER, INTENT(in)     :: QuadType
   
   
      ! Compute the quadrature nodes and weights
      IF( QuadType  == GAUSS_LOBATTO )then ! Gauss Lobatto quadrature
         CALL LegendreGaussLobatto( N, nodes, weights )
      ELSEIF( QuadType == GAUSS )then  ! Gauss Quadrature
         CALL LegendreGauss( N, nodes, weights )
      ELSE
         PRINT*,'Module Legendre.f90 : S/R GenerateLegendreQuadrature: Invalid form. Stopping'
         STOP
      ENDIF

 END SUBROUTINE LegendreQuadrature
!
!
!
 SUBROUTINE LegendreGauss(N, nodes, weights)  
 ! S/R LegendreGauss
 !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 64
 !   Algorithm 23
 ! =============================================================================================== !
  IMPLICIT NONE
  INTEGER :: N
  REAL(prec)    :: nodes(0:N)
  REAL(prec)    :: weights(0:N)
  ! LOCAL
  REAL(prec)    :: lN1, dlN1  ! Legendre polynomial and derivative
  REAL(prec)    :: delta
  INTEGER :: j, kIt ! Loop counter 
 
      IF( N == 0 ) then

         nodes(0) = ZERO
         weights(0) = TWO

      ELSEIF( N == 1 ) then
     
         nodes(0) = -sqrt(ONE/3.0_prec)
         weights(0) = ONE
         nodes(1) = -nodes(0)
         weights(1) = weights(0)

      ELSE ! use Newton's method
     
         DO j = 0, ( (N+1)/2 ) ! Loop over the roots

            nodes(j) = -cos( (TWO*REAL(j,prec) + ONE)*pi/(TWO*REAL(N,prec) + ONE) )

            DO kIt = 1, newtonMax ! Loop over the Newton's iterations

               CALL LegendrePolynomial(N+1, nodes(j), lN1, dlN1)
               delta = -lN1/dlN1
               nodes(j) = nodes(j) + delta
               IF( abs(delta) <= TOL*nodes(j) ) EXIT

            ENDDO ! kIt, loop over the Newton's iterations

            CALL LegendrePolynomial(N+1, nodes(j), lN1, dlN1)
            weights(j) = TWO/( (ONE - nodes(j)*nodes(j))*dlN1*dlN1 )
            weights(N - j) = weights(j) ! uses symmetry to assign weights
            nodes(N - j) = -nodes(j)

         ENDDO ! j, loop over all of the roots
 
      ENDIF ! conditional on whether to use newton's method

      IF( mod(REAL(N,prec),TWO) == ZERO)then ! odd number of roots - get the weight for xRoot=0.0
         
         CALL LegendrePolynomial(N+1, ZERO, lN1, dlN1)
         nodes(N/2) = ZERO
         weights(N/2) = 2.0/(dlN1*dlN1)

      ENDIF ! IF we are looking for an odd number of roots.

 END SUBROUTINE LegendreGauss
!
!
!
 SUBROUTINE LegendreGaussLobatto(N, nodes, weights)  
 ! S/R LegendreGaussLobatto
 !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 66
 !   Algorithm 25
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER    :: N
   REAL(prec) :: nodes(0:N)
   REAL(prec) :: weights(0:N)
   ! LOCAL
   REAL(prec) :: delta, q, qprime, lN
   INTEGER    :: j, kIt  
 

      IF( N == 1 ) then

         nodes(0) = -ONE
         weights(0) = ONE
         nodes(1) = ONE
         weights(1) = ONE

      ELSE ! use Newton's method

         nodes(0) = -ONE
         weights(0) = TWO/(REAL(N,prec)*(REAL(N,prec) + ONE) )

         nodes(N) = ONE
         weights(N) = weights(0)

         DO j = 1, ( (N+1)/2 -1 ) ! Loop over the roots

            nodes(j) = -cos( (REAL(j,prec) + 0.25_prec)*pi/REAL(N,prec) - &
                               3.0_prec/(8.0_prec*REAL(N,prec)*pi*(REAL(j,prec) + 0.25_prec) ) )

            DO kIt = 1, newtonMax ! Loop over the Newton's iterations

              CALL LegendreQandL(N, nodes(j), q, qprime, lN)
              delta = -q/qprime
              nodes(j) = nodes(j) + delta
              IF( abs(delta) <= TOL*nodes(j) ) EXIT

            ENDDO ! kIt, loop over the Newton's iterations
            
            CALL LegendreQandL(N, nodes(j), q, qprime, lN)
            weights(j) = TWO/( REAL(N,prec)*(REAL(N,prec) + ONE)*lN*lN )
            weights(N - j) = weights(j) ! uses symmetry to assign weights
            nodes(N - j) = -nodes(j)

         ENDDO ! j, loop over all of the roots

      ENDIF ! conditional on whether to use newton's method

     
      IF( mod(REAL(N,prec),TWO) == ZERO)then ! odd number of roots - get the weight for xRoot=0.0
         
         CALL LegendreQandL(N, ZERO, q, qprime, lN)
         nodes(N/2) = ZERO
         weights(N/2) = TWO/( REAL(N,prec)*(REAL(N,prec) + ONE)*lN*lN )

      ENDIF ! IF we are looking for an odd number of roots.

 END SUBROUTINE LegendreGaussLobatto
!
!
!
 SUBROUTINE LegendreGaussRadau(N, nodes, weights)  
 ! S/R LegendreGaussRadau
 !  
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER    :: N
   REAL(prec) :: nodes(0:N)
   REAL(prec) :: weights(0:N)
   ! LOCAL
   REAL(prec) :: delta, dLdx, lN, q, qprime
   INTEGER    :: j, kIt ! Loop counter 


      nodes(0) = ONE
      weights(0) = TWO/(REAL(N,prec)*(REAL(N,prec)) )

      DO j = 1, N ! Loop over the roots

        nodes(j) =  -cos( HALF*pi*(TWO*REAL(j+1,prec)-ONE)/(TWO*REAL(N,prec)-ONE) )
        
         DO kIt = 1, newtonMax ! Loop over the Newton's iterations

            CALL LegendrePolynomial(N+1, nodes(j), lN, dLdx)
            q = lN
            qprime = dLdx*(ONE + nodes(j)) - lN
            CALL LegendrePolynomial(N, nodes(j), lN, dLdx)
            q = (q + lN)
            qprime = (qprime + (ONE + nodes(j))*dLdx - lN)/( ONE + nodes(j) )
            delta = -q/qprime 
            nodes(j) = nodes(j) + delta
 
            IF( abs(delta) <= TOL*nodes(j) ) EXIT

         ENDDO ! kIt, loop over the Newton's iterations

         CALL LegendrePolynomial(N, nodes(j), lN, dLdx)
         weights(j) = ONE/( (ONE-nodes(j))*dLdx*dLdx )
         nodes(j) = -nodes(j)

      ENDDO ! j, loop over all of the roots

 END SUBROUTINE LegendreGaussRadau
!
!
!
 SUBROUTINE LegendrePolynomial(N, x, lAtX, dLdxAtX)
 ! S/R LegendrePolynomial
 !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 63
 !   Algorithm 22
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER       :: N
   REAL(prec)    :: x
   REAL(prec)    :: lAtX, dLdxAtX
   ! LOCAL
   REAL(prec) :: lNm1, lNm2, dlNm1, dlNm2
   INTEGER :: i

      IF( N == 0 )then
 
         lAtX = ONE    ! Legendre Polynomial
         dLdxAtX = ZERO ! Derivative
      
      ELSEIF( N == 1)then

         lAtX = x       ! Legendre Polynomial
         dLdxAtX = ONE  ! Derivative

      ELSE  ! Then we turn to the recursive relation for higher order Legendre Polynomials
  
         lnM2 = ONE     ! Initializing the recursive loop
         lnM1 = x
         dlnM2 = ZERO
         dlnM1 = ONE

         DO i = 2,N ! Recursive relation for the legendre polynomials
        
            lAtX = ((TWO*REAL(i,prec) - ONE)*x*lnM1 -&
                    (REAL(i,prec) - 1.0)*lnM2)/(REAL(i,prec))

            dldxAtX = dlnM2 + (TWO*REAL(i,prec)-ONE)*lnM1
            lnM2 = lnM1
            lnM1 = lAtX
            dlnM2 = dlnM1
            dlnM1 = dldxAtX
 
         ENDDO ! i, loop over the legendre polynomial degrees.

      ENDIF

 END SUBROUTINE LegendrePolynomial
!
!
!
 SUBROUTINE LegendreQandL(N, x, q, qprime, lN)
 ! S/R LegendreQandL
 !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 65
 !   Algorithm 24
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER       :: N
   REAL(prec)    :: x
   REAL(prec)    :: lN, q, qprime
   ! LOCAL
   REAL(prec) :: lNm1, lNm2, dlNm1, dlNm2, dlN, lN1, dlN1
   INTEGER :: i

      lNm2 = ONE
      lNm1 = x
      dlNm2 = ZERO
      dlNm1 = ONE

      DO i = 2, N
         lN = (TWO*i - ONE)/(REAL(i,prec))*x*lNm1 - (REAL(i,prec) - ONE)/(REAL(i,prec))*lNm2
         dlN = dlNm2 + (TWO*REAL(i,prec) - ONE)*lNm1
         lNm2 = lNm1
         lNm1 = lN
         dlNm2 = dlNm1
         dlNm1 = dlN
      ENDDO

      i = N + 1
      lN1 = (TWO*i - ONE)/(REAL(i,prec))*x*lN - (REAL(i,prec) - ONE)/(REAL(i,prec))*lNm2
      dlN1 = dlNm2 + (TWO*REAL(i,prec) - ONE)*lNm1
      q = lN1 - lNm2
      qprime = dlN1 - dlNm2

 END SUBROUTINE LegendreQandL

END MODULE Quadrature
