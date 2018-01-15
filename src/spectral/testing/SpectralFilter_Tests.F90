! SpectralFilter_Tests.F90
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
!
! This program exercises the 3-D filtering routines using a cutoff filter and a
! rolloff filter. The filters also come equipped with a routine for diagnosing
! modal coefficients. Here, it is used to verify internal consistency of the
! filtering process. The filtering is demonstrated on a test function that is
! the sum of a “smooth” function with low wave number variability and a function
! with high wave number variability.
!
! —> How to verify these procedures ?? Use an identity. Set test function equal
! to the 7th Lagrange interpolating polynomial. Filter that function with a
! modal cutoff filter with the cutoff set to 6. The result should be zero to
! machine precision.
!
! ========================================================================================================= !
PROGRAM SpectralFilter_Tests

USE Quadrature
USE NodalDG_Class
USE NodalDGSolution_3D_Class
USE SpectralFilter_Class

IMPLICIT NONE


  INTEGER, PARAMETER :: polyDegree = 7
 
  TYPE( SpectralFilter )     :: modalFilter
  TYPE( NodalDG )            :: dgStorage
  TYPE( NodalDGSolution_3D ) :: Legendre
  TYPE( NodalDGSolution_3D ) :: filteredFunction
  REAL(prec)                 :: L(0:polyDegree), Ls, transferFunction(1:polyDegree,1:polyDegree)
!  REAL(prec)                 :: modalCoeffs(0:polyDegree,0:polyDegree,0:polyDegree)
  INTEGER                    :: i, j, k, filterID, polyID


    CALL dgStorage % Build( targetPoints  = UniformPoints( -1.0_prec, 1.0_prec, 10 ), &
                            N             = polyDegree, &
                            nTargetPoints = 10, &
                            quadrature    = GAUSS  )

    CALL filteredFunction % Build( N          = polyDegree, &
                                   nEquations = 1, &
                                   nElements  = 1 )

    DO filterID = 1, polyDegree

      CALL modalFilter % Build( dgStorage % interp % interpolationPoints, &
                                dgStorage % quadratureWeights, &
                                polyDegree, filterID, ModalCutoff )
  
      DO polyID = 1, polyDegree
  
        CALL Legendre % Build( N          = polyDegree, &
                               nEquations = 1, &
                               nElements  = 1 )
    
        DO i = 0, polyDegree
          CALL LegendrePolynomial( polyID, dgStorage % interp % interpolationPoints(i), &
                                   L(i), Ls)
        ENDDO
       
        DO k = 0, polyDegree
          DO j = 0, polyDegree
            DO i = 0, polyDegree
    
              Legendre % solution(i,j,k,1,1) = L(i)*L(j)*L(k)
    
            ENDDO
          ENDDO
        ENDDO
    
        CALL modalFilter % Filter3D( f = Legendre % solution,&
                                     filteredF = filteredFunction % solution, &
                                     nVariables = 1, &
                                     nElements = 1 )
     
        transferFunction(polyID,filterID) = MAXVAL( ABS(filteredFunction % solution) )/MAXVAL( ABS(Legendre % solution) )
        CALL Legendre % Trash( )

      ENDDO
  
      CALL modalFilter % Trash( )

    ENDDO

    DO filterID = 1, polyDegree 

      WRITE(*,'(7(E13.6))') transferFunction(1:polyDegree,filterID)

    ENDDO


    CALL filteredFunction % Trash( )
    CALL dgStorage % Trash( )

END PROGRAM SpectralFilter_Tests
