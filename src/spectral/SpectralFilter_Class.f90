! SpectralFilter_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE SpectralFilter_Class

USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE Quadrature

#ifdef HAVE_CUDA
USE cudafor
#endif

IMPLICIT NONE

! ========================================================================================================= !
! SpectralFilter_Class 
!
!   A data-structure for handling Legendre Modal Filtering in 3-D using either a modal cutoff filter or
!   a roll-off filter.
!
!  This module provides a data structure for constructing and storing element-local filter matrices
!  that can be used for polynomial de-aliasing or as part of SGS parameterizations.
!
!   This module was inspired by the paper
! 
!    D. Flad, A. Beck, and C. Munz, (2016) "Simulation of underresolved turbulent flows by adaptive 
!       filtering using the high order discontinuous Galerkin spectral element method", JCP, 313, 
!       pp. 1-12
!
!   It is assumed that data is described by the same polynomial degree in each computational direction
!   and that the filtering is performed using the same cutoff degree in each computational direction.
!
! ========================================================================================================= !


  TYPE SpectralFilter
     INTEGER                 :: N
     REAL(prec), ALLOCATABLE :: filterMat(:,:)
     REAL(prec), ALLOCATABLE :: nodalToModal(:,:)

#ifdef HAVE_CUDA
     INTEGER, ALLOCATABLE, DEVICE    :: N_dev
     REAL(prec), ALLOCATABLE, DEVICE :: filterMat_dev(:,:)
     REAL(prec), ALLOCATABLE, DEVICE :: nodalToModal_dev(:,:)
#endif
      
     CONTAINS
 
     PROCEDURE :: Build => Build_SpectralFilter
     PROCEDURE :: Trash => Trash_SpectralFilter
     
     PROCEDURE :: Filter3D

  END TYPE SpectralFilter


CONTAINS

! ========================================================================================================= !
! 
! Build_SpectralFilter 
! Allocates space for the modal cutoff filter and initializes the attributes of the data structure 
! 
!   Usage :
!     TYPE(SpectralFilter) :: thisFilter
!     INTEGER              :: N, nCutoff, filterType
!     REAL(prec)           :: interpNodes(0:N), interpWeights(0:N)
!     
!       filterType = ModalCutoff ! or set to TanhRollOff
!       CALL thisFilter % Build( s, w, N, nCutoff, filterType )
! 
!   Input/Output :
!   
!    thisFilter (in)
! 
!    interpNodes (in)
!      The interpolation nodes of your data
!
!    interpWeights (in)
!      The quadrature weights for discrete integration associated with the interpolation nodes
!
!    N (in)
!      Polynomial degree associated with the quadrature points and weights
!
!    nCutoff (in)
!      Cutoff polynomial degree for the spectral filter
!
!    filterType (in)
!      An integer flag that can be set to ModalCutoff or TanhRolloff to select filter type.
!   
! ========================================================================================================= !

  SUBROUTINE Build_SpectralFilter( thisFilter, interpNodes, interpWeights, N, nCutoff, filter_a, filter_b, filterType )

    IMPLICIT NONE
    CLASS(SpectralFilter), INTENT(inout) :: thisFilter
    INTEGER, INTENT(in)                  :: N, nCutoff
    REAL(prec), INTENT(in)               :: filter_a, filter_b 
    REAL(prec), INTENT(in)               :: interpNodes(0:N)
    REAL(prec), INTENT(in)               :: interpWeights(0:N)
    INTEGER, INTENT(in)                  :: filterType
    ! Local
    REAL(dp) :: interpNodes_dp(0:N)
    REAL(dp) :: interpWeights_dp(0:N)
    REAL(dp) :: Lnorm, Li, Ls, sc, r
    REAL(dp) :: Pfilt(0:N,0:N), V(0:N,0:N), VInv(0:N,0:N)
    INTEGER  :: row, col 
   
      DO col = 0, N
        interpNodes_dp(col)   = REAL( interpNodes(col), prec )
        interpWeights_dp(col) = REAL( interpWeights(col), prec )
      ENDDO

      thisFilter % N       = N

      ALLOCATE( thisFilter % filterMat(0:N,0:N), &
                thisFilter % nodalToModal(0:N,0:N) )

      thisFilter % filterMat = 0.0_dp

      Pfilt = 0.0_dp
      V     = 0.0_dp
      VInv  = 0.0_dp

      DO row = 0, N 
 
        r = real(row,dp)

        IF( filterType == ModalCutoff )THEN

          IF( row <= nCutoff )THEN
            Pfilt(row,row) = 1.0_dp
          ELSE
            Pfilt(row,row) = 0.0_dp
          ENDIF

        ELSEIF( filterType == TanhRollOff )THEN

          IF( row < N-2 )THEN
            Pfilt(row,row) = 1.0_prec
          ELSE
            Pfilt(row,row) = 0.5_dp*(1.0_dp - tanh( (r- filter_b)*filter_a ) )
          ENDIF

        ELSEIF( filterType == RampFilter )THEN

          IF( row == N-1 )THEN
            Pfilt(row,row) = filter_b 
          ELSEIF( row == N )THEN
            Pfilt(row,row) = filter_a
          ELSE
            Pfilt(row,row) = 1.0_prec
          ENDIF

        ENDIF


        Lnorm = 0.0_dp

        DO col = 0, N

          CALL LegendrePolynomial( row, interpNodes_dp(col), Li, Ls)

          Lnorm = Lnorm + Li*Li*interpWeights_dp(col)
          thisFilter % nodalToModal(row,col) = REAL( Li*interpWeights_dp(col), prec )
          VInv(col,row) = Li
             
        ENDDO

        thisFilter % nodalToModal(row,0:N) = thisFilter % nodalToModal(row,0:N)/REAL( Lnorm, prec )

      ENDDO

      Pfilt = MATMUL( Pfilt, thisFilter % nodalToModal )
      thisFilter % filterMat = REAL( TRANSPOSE( MATMUL( VInv, Pfilt ) ), prec )

#ifdef HAVE_CUDA
      ALLOCATE( thisFilter % N_dev, &
                thisfilter % filterMat_dev(0:N,0:N), &
                thisFilter % nodalToModal_dev(0:N,0:N) )

      thisFilter % filterMat_dev    = thisFilter % filterMat
      thisFilter % N_dev            = N
      thisFilter % nodalToModal_dev = thisFilter % nodalToModal
#endif


  END SUBROUTINE Build_SpectralFilter

! ========================================================================================================= !
!
! Trash_SpectralFilter 
!   Frees memory associated with the modal cutoff filter 
! 
!   Usage :
!     TYPE(SpectralFilter) :: thisFilter
!
!       CALL thisFilter % Trash( )
! 
!   Input/Output
!   
!     thisFilter (in/out)
!       On input, a previously constructed SpectralFilter instance. On output, memory for that instance
!       has been freed.
!
!   
! ========================================================================================================= !

  SUBROUTINE Trash_SpectralFilter( thisFilter )

    IMPLICIT NONE
    CLASS(SpectralFilter), INTENT(inout) :: thisFilter
   
      DEALLOCATE( thisFilter % filterMat, &
                  thisFilter % nodalToModal )

#ifdef HAVE_CUDA
      DEALLOCATE( thisFilter % N_dev, &
                  thisfilter % filterMat_dev, &
                  thisFilter % nodalToModal_dev )
#endif

  END SUBROUTINE Trash_SpectralFilter

! ========================================================================================================= !
!
! Filter3D 
!   Applies the spectral filter to 3-D data for a number of variables over a number of elements.
!   This routine serves as a wrapper that calls either the CPU or CUDA kernel according to the compiler
!   flags 
!
!   Usage :
!     TYPE(SpectralFilter) :: thisFilter
!     INTEGER              :: nVariables, nElements
!     REAL(prec)           :: f(0:N,0:N,0:N,1:nVariables,1:nElements)
!     REAL(prec)           :: filteredF(0:N,0:N,0:N,1:nVariables,1:nElements)
!
!       CALL thisFilter % Filter3D( f, filteredF, nVariables, nElements )
! 
!   * Note that "N" here refers to "thisFilter % N"
!
!   Input/Output
!   
!     thisFilter (in)
!       A Spectral Filter instance
!
!     f (in)
!       3-D data representing a nodal values of a function for a set of variables over a number of elements
!
!     filteredF (out)
!       3-D data after filter has been applied
!
!     nVariables (in)
!       The number of functions/variables being filtered. The fourth dimension of the array
!
!     nElements (in)
!       The number of elements
!   
! ========================================================================================================= !

  SUBROUTINE Filter3D( thisFilter, f, filteredF, nVariables, nElements )
 
    IMPLICIT NONE
    CLASS(SpectralFilter), INTENT(in) :: thisFilter 
   
#ifdef HAVE_CUDA

    INTEGER, DEVICE, INTENT(in)       :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)    :: f(0:thisFilter % N, &
                                           0:thisFilter % N, &
                                           0:thisFilter % N, &
                                           1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out)   :: filteredF(0:thisFilter % N, &
                                                   0:thisFilter % N, &
                                                   0:thisFilter % N, &
                                                   1:nVariables, 1:nElements)
    ! Local
    TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3( 4*(ceiling( REAL(thisFilter % N+1)/4 ) ), &
                     4*(ceiling( REAL(thisFilter % N+1)/4 ) ) , &
                     4*(ceiling( REAL(thisFilter % N+1)/4 ) ) )
      grid = dim3( nVariables, nElements, 1 )
  
      CALL Filter3D_CUDAKernel<<<grid,tBlock>>>( f, filteredF, &
                                                 thisFilter % filterMat_dev, &
                                                 thisFilter % N_dev, nVariables, nElements )
                                                 

#else

    INTEGER, INTENT(in)    :: nVariables, nElements
    REAL(prec), INTENT(in) :: f(0:thisFilter % N, &
                                0:thisFilter % N, &
                                0:thisFilter % N, &
                                1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: filteredF(0:thisFilter % N, &
                                        0:thisFilter % N, &
                                        0:thisFilter % N, &
                                        1:nVariables, 1:nElements)


      filteredF = Filter3D_SpectralFilter( thisFilter, f, nVariables, nElements )

#endif

  END SUBROUTINE Filter3D

! ========================================================================================================= !
!
! CalculateModalCoefficients3D 
!   Applies the nodal-to-modal matrix to 3-D data for a number of variables over a number of elements to 
!   calculate the Legendre modal coefficients.
!
!   This routine serves as a wrapper that calls either the CPU or CUDA kernel according to the compiler
!   flags 
!
!   Usage :
!     TYPE(SpectralFilter) :: thisFilter
!     INTEGER              :: nVariables, nElements
!     REAL(prec)           :: nodalF(0:N,0:N,0:N,1:nVariables,1:nElements)
!     REAL(prec)           :: modalF(0:N,0:N,0:N,1:nVariables,1:nElements)
!
!       CALL thisFilter % Filter3D( nodalF, modalF, nVariables, nElements )
! 
!   * Note that "N" here refers to "thisFilter % N"
!
!   Input/Output
!   
!     thisFilter (in)
!       A Spectral Filter instance
!
!     nodalF (in)
!       3-D data representing a nodal values of a function for a set of variables over a number of elements
!
!     modalF (out)
!       3-D data representing the modal coefficients for each variable and element
!
!     nVariables (in)
!       The number of functions/variables being filtered. The fourth dimension of the array
!
!     nElements (in)
!       The number of elements
!   
! ========================================================================================================= !

  SUBROUTINE CalculateModalCoefficients3D( thisFilter, nodalF, modalF, nVariables, nElements )
 
    IMPLICIT NONE
    CLASS(SpectralFilter), INTENT(in) :: thisFilter 
   
#ifdef HAVE_CUDA

    INTEGER, DEVICE, INTENT(in)       :: nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)    :: nodalF(0:thisFilter % N, &
                                                0:thisFilter % N, &
                                                0:thisFilter % N, &
                                                1:nVariables, 1:nElements)
    REAL(prec), DEVICE, INTENT(out)   :: modalF(0:thisFilter % N, &
                                                0:thisFilter % N, &
                                                0:thisFilter % N, &
                                                1:nVariables, 1:nElements)
    ! Local
    TYPE(dim3) :: grid, tBlock
  
      tBlock = dim3( 4*(ceiling( REAL(thisFilter % N+1)/4 ) ), &
                     4*(ceiling( REAL(thisFilter % N+1)/4 ) ) , &
                     4*(ceiling( REAL(thisFilter % N+1)/4 ) ) )
      grid = dim3( nVariables, nElements, 1 )
  
      CALL Filter3D_CUDAKernel<<<grid,tBlock>>>( nodalF, modalF, &
                                                 thisFilter % nodalToModal_dev, &
                                                 thisFilter % N_dev, nVariables, nElements )
                                                 

#else

    INTEGER, INTENT(in)     :: nVariables, nElements
    REAL(prec), INTENT(in)  :: nodalF(0:thisFilter % N, &
                                      0:thisFilter % N, &
                                      0:thisFilter % N, &
                                      1:nVariables, 1:nElements)
    REAL(prec), INTENT(out) :: modalF(0:thisFilter % N, &
                                      0:thisFilter % N, &
                                      0:thisFilter % N, &
                                      1:nVariables, 1:nElements)


      modalF = CalculateModalCoefficients3D_SpectralFilter( thisFilter, nodalF, nVariables, nElements )

#endif

  END SUBROUTINE CalculateModalCoefficients3D
 
! ================================================================================================ !
! ------------------------------------- PRIVATE ROUTINES ----------------------------------------- !
! ================================================================================================ !

  FUNCTION Filter3D_SpectralFilter( thisFilter, f, nVariables, nElements ) RESULT( filteredF )

    IMPLICIT NONE
    CLASS( SpectralFilter ) :: thisFilter
    INTEGER                 :: nVariables, nElements
    REAL(prec)              :: f(0:thisFilter % N, &
                                 0:thisFilter % N, &
                                 0:thisFilter % N, &
                                 1:nVariables, 1:nElements)
    REAL(prec)              :: filteredF(0:thisFilter % N, &
                                         0:thisFilter % N, &
                                         0:thisFilter % N, &
                                         1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, j, k, iVar, iEl
    INTEGER :: ii, jj, kk
    REAL(prec) :: uij, ui
 
     
      !$OMP PARALLEL
      !$OMP DO PRIVATE( i, j, k, iVar, iEl, ii, jj, kk, uij, ui )
      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, thisFilter % N
            DO j = 0, thisFilter % N
              DO i = 0, thisFilter % N
                  
                filteredF(i,j,k,iVar,iEl) = 0.0_prec

                DO kk = 0, thisFilter % N
                     
                  uij = 0.0_prec

                  DO jj = 0, thisFilter % N
                           
                    ui = 0.0_prec

                    DO ii = 0, thisFilter % N

                      ui = ui + thisFilter % filterMat(ii,i)*f(ii,jj,kk,iVar,iEl)

                    ENDDO
                           
                    uij = uij + thisFilter % filterMat(jj,j)*ui

                  ENDDO
                        
                  filteredF(i,j,k,iVar,iEl) = filteredF(i,j,k,iVar,iEl) + thisFilter % filterMat(kk,k)*uij
                        
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !$OMP ENDDO
      !$OMP END PARALLEL

         
  END FUNCTION Filter3D_SpectralFilter
 
  FUNCTION CalculateModalCoefficients3D_SpectralFilter( thisFilter, nodalF, nVariables, nElements ) RESULT( modalF )
    IMPLICIT NONE
    CLASS( SpectralFilter ) :: thisFilter
    INTEGER                 :: nVariables, nElements
    REAL(prec)              :: nodalf(0:thisFilter % N, &
                                      0:thisFilter % N, &
                                      0:thisFilter % N, &
                                      1:nVariables, 1:nElements)
    REAL(prec)              :: modalF(0:thisFilter % N, &
                                      0:thisFilter % N, &
                                      0:thisFilter % N, &
                                      1:nVariables, 1:nElements)
    ! Local
    INTEGER :: i, j, k, iVar, iEl
    INTEGER :: ii, jj, kk
    REAL(prec) :: uij, ui
                                        

      DO iEl = 1, nElements
        DO iVar = 1, nVariables
          DO k = 0, thisFilter % N
            DO j = 0, thisFilter % N
              DO i = 0, thisFilter % N
                     
                modalF(i,j,k,iVar,iEl) = 0.0_prec

                DO kk = 0, thisFilter % N
                        
                  uij = 0.0_prec

                  DO jj = 0, thisFilter % N
                              
                    ui = 0.0_prec

                    DO ii = 0, thisFilter % N

                      ui = ui + thisFilter % nodalToModal(ii,i)*nodalF(ii,jj,kk,iVar,iEl)

                    ENDDO
                              
                    uij = uij + thisFilter % nodalToModal(jj,j)*ui

                  ENDDO
                           
                  nodalF(i,j,k,iVar,iEl) = modalF(i,j,k,iVar,iEl) + thisFilter % nodalToModal(kk,k)*uij
                           
                ENDDO

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
         
  END FUNCTION CalculateModalCoefficients3D_SpectralFilter

#ifdef HAVE_CUDA
  ATTRIBUTES( Global ) SUBROUTINE Filter3D_CUDAKernel( f, filteredF, filterMatrix, N, nVariables, nElements )
    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nVariables, nElements
    REAL(prec), DEVICE, INTENT(in)  :: f(0:N,0:N,0:N,1:nVariables,1:nElements)
    REAL(prec), DEVICE, INTENT(in)  :: filterMatrix(0:N,0:N)
    REAL(prec), DEVICE, INTENT(out) :: filteredF(0:N,0:N,0:N,1:nVariables,1:nElements)
    ! Local
    INTEGER            :: i, j, k, iVar, iEl
    INTEGER            :: ii, jj, kk
    REAL(prec)         :: uijk, uij, ui
    REAL(prec), SHARED :: fLocal(0:7,0:7,0:7)
   
   
      iVar = blockIdx % x
      iEl  = blockIdx % y
      
      i = threadIdx % x-1
      j = threadIdx % y-1
      k = threadIdx % z-1
      
      IF( i <= N .AND. j <= N .AND. k <= N )THEN

        fLocal(i,j,k) = f(i,j,k,iVar,iEl)

        CALL syncthreads( )
      
        uijk = 0.0_prec

        DO kk = 0, N
        
          uij = 0.0_prec

          DO jj = 0, N
              
            ui = 0.0_prec

            DO ii = 0, N

              ui = ui + filterMatrix(ii,i)*fLocal(ii,jj,kk)

            ENDDO
              
            uij = uij + filterMatrix(jj,j)*ui

          ENDDO
           
          uijk = uijk + filterMatrix(kk,k)*uij
           
        ENDDO
                        
        filteredF(i,j,k,iVar,iEl) = uijk
   
      ENDIF

   
 END SUBROUTINE Filter3D_CUDAKernel
#endif 
 
END MODULE SpectralFilter_Class
