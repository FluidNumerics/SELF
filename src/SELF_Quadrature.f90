! SELF_Quadrature.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Contains routines from D.A. Kopriva, 2009, "Implementing Spectral Methods for Partial
! Differential Equations: Algorithms for Scientists and Engineers", Springer.
!
! Routines are defined for computing Legendre and Chebyshev Gauss and Gauss-Lobatto
! quadrature nodes and weights.

MODULE SELF_Quadrature

  USE ISO_FORTRAN_ENV
  USE SELF_Constants

  IMPLICIT NONE

  PUBLIC  :: ChebyshevQuadrature,LegendreQuadrature
  PRIVATE :: ChebyshevGauss,ChebyshevGaussLobatto, &
             LegendreGauss,LegendreGaussLobatto, &
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
!     REAL(real64) :: nodes(0:N), weights(0:N)
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

  SUBROUTINE LegendreQuadrature(N,nodes,weights,QuadType)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: N
    REAL(prec),INTENT(out) :: nodes(0:N)
    REAL(prec),INTENT(out) :: weights(0:N)
    INTEGER,INTENT(in)     :: QuadType
    ! Local
    REAL(real64) :: nodesLocal(0:N)
    REAL(real64) :: weightsLocal(0:N)
    INTEGER :: i


    IF (QuadType == GAUSS_LOBATTO) THEN

      CALL LegendreGaussLobatto(N,nodesLocal,weightsLocal)

    ELSEIF (QuadType == GAUSS) THEN

      CALL LegendreGauss(N,nodesLocal,weightsLocal)

    END IF

    DO i = 0, N
      nodes(i) = REAL(nodesLocal(i),prec)
      weights(i) = REAL(weightsLocal(i),prec)
    ENDDO

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
!     REAL(real64) :: nodes(0:N), weights(0:N)
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

  SUBROUTINE ChebyshevQuadrature(N,quadType,nodes,weights)
    IMPLICIT NONE
    INTEGER,INTENT(in)     :: N
    REAL(real64),INTENT(out) :: nodes(0:N)
    REAL(real64),INTENT(out) :: weights(0:N)
    INTEGER,INTENT(in)     :: QuadType
    ! Local
    REAL(real64) :: nodesLocal(0:N)
    REAL(real64) :: weightsLocal(0:N)
    INTEGER :: i

    IF (QuadType == GAUSS_LOBATTO) then

      CALL ChebyshevGaussLobatto(N,nodesLocal,weightsLocal)

    ELSEIF (QuadType == GAUSS) then

      CALL ChebyshevGauss(N,nodesLocal,weightsLocal)

    END IF

    DO i = 0, N
      nodes(i) = REAL(nodesLocal(i),prec)
      weights(i) = REAL(weightsLocal(i),prec)
    ENDDO

  END SUBROUTINE ChebyshevQuadrature

! =============================================================================================== !
! S/R ChebyshevGauss
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 67
!   Algorithm 26
! =============================================================================================== !

  SUBROUTINE ChebyshevGauss(N,nodes,weights)
    IMPLICIT NONE
    INTEGER       :: N
    REAL(real64)    :: nodes(0:N)
    REAL(real64)    :: weights(0:N)
    ! Local
    INTEGER    :: j

    DO j = 0,N

      weights(j) = pi/(REAL(N,real64) + 1.0_real64)
      nodes(j) = -cos(0.5_real64*(2.0_real64*REAL(j,real64) + 1.0_real64)*weights(j))

    END DO

  END SUBROUTINE ChebyshevGauss

! =============================================================================================== !
! S/R ChebyshevGaussLobatto
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 68
!   Algorithm 27
! =============================================================================================== !

  SUBROUTINE ChebyshevGaussLobatto(N,nodes,weights)
    IMPLICIT NONE
    INTEGER       :: N
    REAL(real64)    :: nodes(0:N)
    REAL(real64)    :: weights(0:N)
    ! LOCAL
    INTEGER    :: j

    DO j = 0,N

      weights(j) = pi/REAL(N,real64)
      nodes(j) = -cos(REAL(j,real64)*weights(j))

    END DO

    weights(0) = weights(0)*0.5_real64
    weights(N) = weights(N)*0.5_real64

  END SUBROUTINE ChebyshevGaussLobatto

! =============================================================================================== !
! S/R LegendreGauss
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 64
!   Algorithm 23
! =============================================================================================== !

  SUBROUTINE LegendreGauss(N,nodes,weights)
    IMPLICIT NONE
    INTEGER    :: N
    REAL(real64) :: nodes(0:N)
    REAL(real64) :: weights(0:N)
    ! Local
    REAL(real64) :: nodes_local(0:N)
    REAL(real64) :: weights_local(0:N)
    REAL(real64) :: lN1,dlN1
    REAL(real64) :: delta
    INTEGER  :: j,kIt

    IF (N == 0) then

      nodes_local(0) = 0.0_real64
      weights_local(0) = 2.0_real64

    ELSEIF (N == 1) then

      nodes_local(0) = -SQRT(1.0_real64/3.0_real64)
      weights_local(0) = 1.0_real64
      nodes_local(1) = -nodes(0)
      weights_local(1) = weights(0)

    ELSE

      DO j = 0, ((N + 1)/2)

        nodes_local(j) = -cos((2.0_real64*REAL(j,real64) + 1.0_real64)*pi/(2.0_real64*REAL(N,real64) + 1.0_real64))

        DO kIt = 1,newtonMax

          CALL LegendrePolynomial(N + 1,nodes_local(j),lN1,dlN1)
          delta = -lN1/dlN1
          nodes_local(j) = nodes_local(j) + delta
          IF (abs(delta) <= TOL*nodes_local(j)) EXIT

        END DO

        CALL LegendrePolynomial(N + 1,nodes_local(j),lN1,dlN1)
        weights_local(j) = 2.0_real64/((1.0_real64 - nodes_local(j)*nodes_local(j))*dlN1*dlN1)
        weights_local(N - j) = weights_local(j)
        nodes_local(N - j) = -nodes_local(j)

      END DO

    END IF

    IF (MOD(REAL(N,real64),2.0_real64) == 0.0_real64) then

      CALL LegendrePolynomial(N + 1,0.0_real64,lN1,dlN1)
      nodes_local(N/2) = 0.0_real64
      weights_local(N/2) = 2.0/(dlN1*dlN1)

    END IF

    DO j = 0,N
      nodes(j) = REAL(nodes_local(j),real64)
      weights(j) = REAL(weights_local(j),real64)
    END DO

  END SUBROUTINE LegendreGauss

  ! =============================================================================================== !
  ! S/R LegendreGaussLobatto
  !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 66
  !   Algorithm 25
  ! =============================================================================================== !

  SUBROUTINE LegendreGaussLobatto(N,nodes,weights)
    IMPLICIT NONE
    INTEGER    :: N
    REAL(real64) :: nodes(0:N)
    REAL(real64) :: weights(0:N)
    ! Local
    REAL(real64) :: nodes_local(0:N)
    REAL(real64) :: weights_local(0:N)
    REAL(real64) :: delta,q,qprime,lN
    INTEGER  :: j,kIt

    IF (N == 1) then

      nodes_local(0) = -1.0_real64
      weights_local(0) = 1.0_real64
      nodes_local(1) = 1.0_real64
      weights_local(1) = 1.0_real64

    ELSE

      nodes_local(0) = -1.0_real64
      weights_local(0) = 2.0_real64/(REAL(N,real64)*(REAL(N,real64) + 1.0_real64))
      nodes_local(N) = 1.0_real64
      weights_local(N) = weights_local(0)

      DO j = 1, ((N + 1)/2 - 1)

        nodes_local(j) = -COS((REAL(j,real64) + 0.25_real64)*pi/REAL(N,real64) - &
                              3.0_real64/(8.0_real64*REAL(N,real64)*pi*(REAL(j,real64) + 0.25_real64)))

        DO kIt = 1,newtonMax

          CALL LegendreQandL(N,nodes_local(j),q,qprime,lN)

          delta = -q/qprime
          nodes_local(j) = nodes_local(j) + delta
          IF (ABS(delta) <= TOL*nodes_local(j)) EXIT

        END DO

        CALL LegendreQandL(N,nodes_local(j),q,qprime,lN)

        weights_local(j) = 2.0_real64/(REAL(N,real64)*(REAL(N,real64) + 1.0_real64)*lN*lN)
        weights_local(N - j) = weights_local(j)
        nodes_local(N - j) = -nodes_local(j)

      END DO

    END IF

    IF (MOD(REAL(N,real64),2.0_real64) == 0.0_real64) THEN

      CALL LegendreQandL(N,0.0_real64,q,qprime,lN)

      nodes_local(N/2) = 0.0_real64
      weights_local(N/2) = 2.0_real64/(REAL(N,real64)*(REAL(N,real64) + 1.0_real64)*lN*lN)

    END IF

    DO j = 0,N
      nodes(j) = REAL(nodes_local(j),real64)
      weights(j) = REAL(weights_local(j),real64)
    END DO

  END SUBROUTINE LegendreGaussLobatto

! =============================================================================================== !
! S/R LegendrePolynomial
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 63
!   Algorithm 22
! =============================================================================================== !

  SUBROUTINE LegendrePolynomial(N,x,lAtX,dLdxAtX)
    IMPLICIT NONE
    INTEGER     :: N
    REAL(real64)    :: x
    REAL(real64)    :: lAtX,dLdxAtX
    ! Local
    REAL(real64) :: lNm1,lNm2,dlNm1,dlNm2
    INTEGER  :: i

    IF (N == 0) then

      lAtX = 1.0_real64
      dLdxAtX = 0.0_real64

    ELSEIF (N == 1) then

      lAtX = x
      dLdxAtX = 1.0_real64

    ELSE

      lnM2 = 1.0_real64
      lnM1 = x
      dlnM2 = 0.0_real64
      dlnM1 = 1.0_real64

      DO i = 2,N

        lAtX = ((2.0_real64*REAL(i,real64) - 1.0_real64)*x*lnM1 - &
                (REAL(i,real64) - 1.0_real64)*lnM2)/(REAL(i,real64))

        dldxAtX = dlnM2 + (2.0_real64*REAL(i,real64) - 1.0_real64)*lnM1
        lnM2 = lnM1
        lnM1 = lAtX
        dlnM2 = dlnM1
        dlnM1 = dldxAtX

      END DO

    END IF

  END SUBROUTINE LegendrePolynomial

! =============================================================================================== !
! S/R LegendreQandL
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 65
!   Algorithm 24
! =============================================================================================== !

  SUBROUTINE LegendreQandL(N,x,q,qprime,lN)
    IMPLICIT NONE
    INTEGER    :: N
    REAL(real64) :: x
    REAL(real64) :: lN,q,qprime
    ! Local
    REAL(real64) :: lNm1,lNm2,dlNm1,dlNm2,dlN,lN1,dlN1
    INTEGER    :: i

    lNm2 = 1.0_real64
    lNm1 = x
    dlNm2 = 0.0_real64
    dlNm1 = 1.0_real64

    DO i = 2,N

      lN = (2.0_real64*i - 1.0_real64)/(REAL(i,real64))*x*lNm1 - (REAL(i,real64) - 1.0_real64)/(REAL(i,real64))*lNm2
      dlN = dlNm2 + (2.0_real64*REAL(i,real64) - 1.0_real64)*lNm1
      lNm2 = lNm1
      lNm1 = lN
      dlNm2 = dlNm1
      dlNm1 = dlN

    END DO

    i = N + 1
    lN1 = (2.0_real64*i - 1.0_real64)/(REAL(i,real64))*x*lN - (REAL(i,real64) - 1.0_real64)/(REAL(i,real64))*lNm2
    dlN1 = dlNm2 + (2.0_real64*REAL(i,real64) - 1.0_real64)*lNm1
    q = lN1 - lNm2
    qprime = dlN1 - dlNm2

  END SUBROUTINE LegendreQandL

END MODULE SELF_Quadrature
