! Quadrature.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Contains routines from D.A. Kopriva, 2009, "Implementing Spectral Methods for Partial 
! Differential Equations: Algorithms for Scientists and Engineers", Springer.
!
! Routines are defined for computing Legendre and Chebyshev Gauss and Gauss-Lobatto
! quadrature nodes and weights.
 
MODULE Quadrature

USE ModelPrecision
USE ConstantsDictionary


IMPLICIT NONE

  PUBLIC  :: ChebyshevQuadrature, LegendreQuadrature
  PRIVATE :: ChebyshevGauss, ChebyshevGaussLobatto, &
             LegendreGauss, LegendreGaussLobatto, &
             LegendreQandL

  CONTAINS

! =============================================================================================== !
! LegendreQuadrature  
!   Returns the specified Legendre quadrature nodes and integration weights. 
! 
!   Given a polynomial degree, and quadrature type (Gauss or Gauss Lobatto), this subroutine manages
!   the calls to underlying private routines to generate the desired Legendre quadrature. 
! 
!   Usage :
!
!     INTEGER    :: N, quadType 
!     REAL(prec) :: nodes(0:N), weights(0:N) 
!
!       CALL LegendreQuadrature( N, quadType, nodes, weights ) 
! 
!   Parameters :
!
!     N (in)
!       Degree of the quadrature
!
!     quadType (in) 
!       Flag specifying the quadrature type. Can be set to GAUSS or GAUSS_LOBATTO
!
!     nodes(0:N) (out)
!       Array of quadrature nodes
!
!     weights(0:N) (out)
!       Array of quadrature weights
!   
! =============================================================================================== !

  SUBROUTINE LegendreQuadrature( N, nodes, weights, QuadType )
    IMPLICIT NONE
    INTEGER, INTENT(in)     :: N
    REAL(prec), INTENT(out) :: nodes(0:N)
    REAL(prec), INTENT(out) :: weights(0:N)
    INTEGER, INTENT(in)     :: QuadType
   
      IF( QuadType  == GAUSS_LOBATTO )THEN

        CALL LegendreGaussLobatto( N, nodes, weights )

      ELSEIF( QuadType == GAUSS )THEN

        CALL LegendreGauss( N, nodes, weights )

      ENDIF

  END SUBROUTINE LegendreQuadrature

! =============================================================================================== !
! ChebyshevQuadrature  
!
!   Returns the specified Chebyshev quadrature nodes and integration weights. 
! 
!   Given a polynomial degree, and quadrature type (Gauss or Gauss Lobatto), this subroutine manages
!   the calls to underlying private routines to generate the desired Chebyshev quadrature. 
! 
!   Usage :
!
!     INTEGER    :: N, quadType 
!     REAL(prec) :: nodes(0:N), weights(0:N) 
!
!       CALL ChebyshevQuadrature( N, quadType, nodes, weights ) 
! 
!   Input/Output : 
!
!     N (in) 
!       Degree of the quadrature
!
!     quadType (in)
!       Flag specifying the quadrature type. Can be set to GAUSS or GAUSS_LOBATTO
!
!     nodes(0:N) (out)
!       Array of quadrature nodes
!
!     weights(0:N) (out)
!       Array of quadrature weights
!   
! ================================================================================================ ! 

  SUBROUTINE ChebyshevQuadrature( N, quadType, nodes, weights )
    IMPLICIT NONE
    INTEGER, INTENT(in)     :: N
    REAL(prec), INTENT(out) :: nodes(0:N)
    REAL(prec), INTENT(out) :: weights(0:N)
    INTEGER, INTENT(in)     :: QuadType
   
   
      IF( QuadType  == GAUSS_LOBATTO )then

        CALL ChebyshevGaussLobatto( N, nodes, weights )

      ELSEIF( QuadType == GAUSS )then 

        CALL ChebyshevGauss( N, nodes, weights )

      ENDIF

  END SUBROUTINE ChebyshevQuadrature

! =============================================================================================== !
! S/R ChebyshevGauss
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 67
!   Algorithm 26
! =============================================================================================== !

  SUBROUTINE ChebyshevGauss(N, nodes, weights)  
    IMPLICIT NONE
    INTEGER       :: N
    REAL(prec)    :: nodes(0:N)
    REAL(prec)    :: weights(0:N)
    ! Local
    INTEGER    :: j

      DO j = 0, N

        weights(j) = pi/( REAL(N, prec) + 1.0_prec )
        nodes(j) = -cos( 0.5_prec*(2.0_prec*REAL(j, prec) + 1.0_prec)*weights(j) )

      ENDDO

  END SUBROUTINE ChebyshevGauss

! =============================================================================================== !
! S/R ChebyshevGaussLobatto
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 68
!   Algorithm 27
! =============================================================================================== !

  SUBROUTINE ChebyshevGaussLobatto(N, nodes, weights)  
    IMPLICIT NONE
    INTEGER       :: N
    REAL(prec)    :: nodes(0:N)
    REAL(prec)    :: weights(0:N)
    ! LOCAL
    INTEGER    :: j 

      DO j = 0, N

        weights(j) = pi/REAL(N, prec)
        nodes(j) = -cos( REAL(j, prec)*weights(j) )

      ENDDO

      weights(0)  = weights(0)*0.5_prec
      weights(N) = weights(N)*0.5_prec 

  END SUBROUTINE ChebyshevGaussLobatto


! =============================================================================================== !
! S/R LegendreGauss
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 64
!   Algorithm 23
! =============================================================================================== !

  SUBROUTINE LegendreGauss( N, nodes, weights )  
    IMPLICIT NONE
    INTEGER    :: N
    REAL(prec) :: nodes(0:N)
    REAL(prec) :: weights(0:N)
    ! Local
    REAL(prec) :: lN1, dlN1
    REAL(prec) :: delta
    INTEGER    :: j, kIt
 
      IF( N == 0 ) then

        nodes(0) = 0.0_prec
        weights(0) = 2.0_prec

      ELSEIF( N == 1 ) then
     
        nodes(0) = -SQRT(1.0_prec/3.0_prec)
        weights(0) = 1.0_prec
        nodes(1) = -nodes(0)
        weights(1) = weights(0)

      ELSE
     
        DO j = 0, ( (N+1)/2 )

          nodes(j) = -cos( (2.0_prec*REAL(j,prec) + 1.0_prec)*pi/(2.0_prec*REAL(N,prec) + 1.0_prec) )

          DO kIt = 1, newtonMax

            CALL LegendrePolynomial(N+1, nodes(j), lN1, dlN1)
            delta = -lN1/dlN1
            nodes(j) = nodes(j) + delta
            IF( abs(delta) <= TOL*nodes(j) ) EXIT

          ENDDO

          CALL LegendrePolynomial(N+1, nodes(j), lN1, dlN1)
          weights(j) = 2.0_prec/( (1.0_prec - nodes(j)*nodes(j))*dlN1*dlN1 )
          weights(N - j) = weights(j)
          nodes(N - j) = -nodes(j)

        ENDDO
 
      ENDIF

      IF( MOD(REAL(N,prec),2.0_prec) == 0.0_prec)then 
         
        CALL LegendrePolynomial(N+1, 0.0_prec, lN1, dlN1)
        nodes(N/2) = 0.0_prec
        weights(N/2) = 2.0/(dlN1*dlN1)

      ENDIF

  END SUBROUTINE LegendreGauss

  ! =============================================================================================== !
  ! S/R LegendreGaussLobatto
  !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 66
  !   Algorithm 25
  ! =============================================================================================== !

  SUBROUTINE LegendreGaussLobatto(N, nodes, weights)  
    IMPLICIT NONE
    INTEGER    :: N
    REAL(prec) :: nodes(0:N)
    REAL(prec) :: weights(0:N)
    ! Local
    REAL(prec) :: delta, q, qprime, lN
    INTEGER    :: j, kIt  
 
      IF( N == 1 ) then

        nodes(0) = -1.0_prec
        weights(0) = 1.0_prec
        nodes(1) = 1.0_prec
        weights(1) = 1.0_prec

      ELSE

        nodes(0) = -1.0_prec
        weights(0) = 2.0_prec/(REAL(N,prec)*(REAL(N,prec) + 1.0_prec) )
        nodes(N) = 1.0_prec
        weights(N) = weights(0)

        DO j = 1, ( (N+1)/2 -1 )

          nodes(j) = -COS( (REAL(j,prec) + 0.25_prec)*pi/REAL(N,prec) - &
                             3.0_prec/(8.0_prec*REAL(N,prec)*pi*(REAL(j,prec) + 0.25_prec) ) )

          DO kIt = 1, newtonMax 

            CALL LegendreQandL(N, nodes(j), q, qprime, lN)

            delta = -q/qprime
            nodes(j) = nodes(j) + delta
            IF( ABS(delta) <= TOL*nodes(j) ) EXIT

          ENDDO
           
          CALL LegendreQandL(N, nodes(j), q, qprime, lN)

          weights(j) = 2.0_prec/( REAL(N,prec)*(REAL(N,prec) + 1.0_prec)*lN*lN )
          weights(N - j) = weights(j)
          nodes(N - j) = -nodes(j)

        ENDDO

      ENDIF

     
      IF( MOD(REAL(N,prec),2.0_prec) == 0.0_prec )THEN
         
        CALL LegendreQandL(N, 0.0_prec, q, qprime, lN)

        nodes(N/2) = 0.0_prec
        weights(N/2) = 2.0_prec/( REAL(N,prec)*(REAL(N,prec) + 1.0_prec)*lN*lN )

      ENDIF

  END SUBROUTINE LegendreGaussLobatto

! =============================================================================================== !
! S/R LegendrePolynomial
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 63
!   Algorithm 22
! =============================================================================================== !

  SUBROUTINE LegendrePolynomial(N, x, lAtX, dLdxAtX)
    IMPLICIT NONE
    INTEGER       :: N
    REAL(prec)    :: x
    REAL(prec)    :: lAtX, dLdxAtX
    ! Local
    REAL(prec) :: lNm1, lNm2, dlNm1, dlNm2
    INTEGER :: i

      IF( N == 0 )then
 
        lAtX = 1.0_prec    
        dLdxAtX = 0.0_prec 
      
      ELSEIF( N == 1)then

        lAtX = x            
        dLdxAtX = 1.0_prec  

      ELSE 
  
        lnM2 = 1.0_prec 
        lnM1 = x
        dlnM2 = 0.0_prec
        dlnM1 = 1.0_prec

        DO i = 2,N
        
          lAtX = ((2.0_prec*REAL(i,prec) - 1.0_prec)*x*lnM1 -&
                  (REAL(i,prec) - 1.0_prec)*lnM2)/(REAL(i,prec))

          dldxAtX = dlnM2 + (2.0_prec*REAL(i,prec)-1.0_prec)*lnM1
          lnM2 = lnM1
          lnM1 = lAtX
          dlnM2 = dlnM1
          dlnM1 = dldxAtX
 
        ENDDO

      ENDIF

  END SUBROUTINE LegendrePolynomial

! =============================================================================================== !
! S/R LegendreQandL
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 65
!   Algorithm 24
! =============================================================================================== !

  SUBROUTINE LegendreQandL(N, x, q, qprime, lN)
    IMPLICIT NONE
    INTEGER    :: N
    REAL(prec) :: x
    REAL(prec) :: lN, q, qprime
    ! Local
    REAL(prec) :: lNm1, lNm2, dlNm1, dlNm2, dlN, lN1, dlN1
    INTEGER    :: i

      lNm2 = 1.0_prec
      lNm1 = x
      dlNm2 = 0.0_prec
      dlNm1 = 1.0_prec

      DO i = 2, N

        lN = (2.0_prec*i - 1.0_prec)/(REAL(i,prec))*x*lNm1 - (REAL(i,prec) - 1.0_prec)/(REAL(i,prec))*lNm2
        dlN = dlNm2 + (2.0_prec*REAL(i,prec) - 1.0_prec)*lNm1
        lNm2 = lNm1
        lNm1 = lN
        dlNm2 = dlNm1
        dlNm1 = dlN

      ENDDO

      i = N + 1
      lN1 = (2.0_prec*i - 1.0_prec)/(REAL(i,prec))*x*lN - (REAL(i,prec) - 1.0_prec)/(REAL(i,prec))*lNm2
      dlN1 = dlNm2 + (2.0_prec*REAL(i,prec) - 1.0_prec)*lNm1
      q = lN1 - lNm2
      qprime = dlN1 - dlNm2

  END SUBROUTINE LegendreQandL

END MODULE Quadrature
