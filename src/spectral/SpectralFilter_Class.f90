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

! SpectralFilter_Class 
! A data-structure for handling Legendre Modal Filtering in 1, 2, and 3
! dimensions using either a modal cutoff filter or a roll-off filter.
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
!  SpectralFilter </H2>
!  Attributes </H3>
!    
!       <tr> <th> N <td> INTEGER  <td> Polynomial degree associated with the filter
!       <tr> <th> nCutoff <td> INTEGER <td> Cutoff polynomial degree indicating which Legendre 
!                                           modal coefficients are made null
!       <tr> <th> nPacked <td> REAL(prec) <td> Upper bound of 2-D array obtained from packing down
!                                              3-D array; used for filtering 3-D data
!       <tr> <th> filterMat(:,:) <td> REAL(prec) <td> The filtering matrix
!    </table>
!
!  Procedures </H3>
!    See \ref SpectralFilter_Class for more information. The first column lists the "call-name" and 
!    the second column lists the name of routine that is aliased onto the call-name.
!    
!       <tr> <th> Build <td> Build_SpectralFilter
!       <tr> <th> Trash <td> Trash_SpectralFilter
!       <tr> <th> Apply1DFilter <td> Apply1DFilter_SpectralFilter
!       <tr> <th> Apply2DFilter <td> Apply2DFilter_SpectralFilter
!       <tr> <th> Apply3DFilter <td> Apply3DFilter_SpectralFilter
!    </table>
!


   TYPE SpectralFilter
      INTEGER                 :: N
      REAL(prec), ALLOCATABLE :: filterMat(:,:)

#ifdef HAVE_CUDA
      INTEGER, ALLOCATABLE, DEVICE    :: N_dev
      REAL(prec), ALLOCATABLE, DEVICE :: filterMat_dev(:,:)
#endif
      
      CONTAINS
  
      PROCEDURE :: Build => Build_SpectralFilter
      PROCEDURE :: Trash => Trash_SpectralFilter
      
      PROCEDURE :: Filter3D

   END TYPE SpectralFilter


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup SpectralFilter_Class 
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_SpectralFilter 
!! Allocates space for the modal cutoff filter and initializes the attributes of the data structure 
!! 
!!  Usage : </H2> 
!! TYPE</B>(SpectralFilter) :: this <BR>
!! INTEGER</B>                 :: N, nCutoff <BR>
!! REAL</B>(prec)              :: s(0:N), w(0:N) <BR>
!!         .... <BR>
!!     CALL</B> this % Build( s, w, N, nCutoff ) <BR>
!! 
!!   Parameters : </H2>
!!   
!!   <tr> <td> out <th> thisFilter <td> SpectralFilter <td> 
!!   <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> The interpolation nodes of your data
!!   <tr> <td> in <th> w(0:N) <td> REAL(prec) <td> The quadrature weights for discrete integration
!!                                                 associated with the interpolation nodes
!!   <tr> <td> in <th> N <td> REAL(prec) <td> Polynomial degree of the data
!!   <tr> <td> in <th> nCutoff <td> REAL(prec) <td> Cutoff polynomial degree
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_SpectralFilter( thisFilter, s, w, N, nCutoff, filterType )

   IMPLICIT NONE
   CLASS(SpectralFilter), INTENT(inout) :: thisFilter
   INTEGER, INTENT(in)                 :: N, nCutoff
   REAL(prec), INTENT(in)              :: s(0:N)
   REAL(prec), INTENT(in)              :: w(0:N)
   INTEGER, INTENT(in)                 :: filterType
   ! Local
   REAL(prec) :: Lnorm, Li, Ls, sc, r
   REAL(prec) :: Pfilt(0:N,0:N), V(0:N,0:N), VInv(0:N,0:N)
   INTEGER    :: row, col 
   
      thisFilter % N       = N

      ALLOCATE( thisFilter % filterMat(0:N,0:N) )

      thisFilter % filterMat = 0.0_prec

      Pfilt = 0.0_prec
      V     = 0.0_prec
      VInv  = 0.0_prec

      DO row = 0, N 
 
         r = real(row,prec)

         IF( filterType == ModalCutoff )THEN
           IF( row <= nCutoff )THEN
              Pfilt(row,row) = 1.0_prec
           ENDIF
         ELSEIF( filterType == TanhRollOff )THEN
           Pfilt(row,row) = 0.5_prec*(1.0_prec - tanh( (r- REAL(nCutoff,prec)) ) )
         ENDIF

         Lnorm = 0.0_prec

         DO col = 0, N

            CALL LegendrePolynomial( row, s(col), Li, Ls)

            Lnorm = Lnorm + Li*Li*w(col)
            V(row,col) = Li*w(col)
            VInv(col,row) = Li
              
         ENDDO
            V(row,0:N) = V(row,0:N)/Lnorm
      ENDDO

      Pfilt = MATMUL( Pfilt, V )
      thisFilter % filterMat = TRANSPOSE( MATMUL( VInv, Pfilt ) )

    
#ifdef HAVE_CUDA
      ALLOCATE( thisFilter % N_dev, &
                thisfilter % filterMat_dev(0:N,0:N) )

      thisFilter % filterMat_dev = thisFilter % filterMat
#endif


 END SUBROUTINE Build_SpectralFilter
!
!> \addtogroup SpectralFilter_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_SpectralFilter 
!! Frees memory associated with the modal cutoff filter 
!! 
!!  Usage : </H2> 
!! TYPE</B>(SpectralFilter) :: this <BR>
!!         .... <BR>
!!     CALL</B> this % Trash(  ) <BR>
!! 
!!   Parameters : </H2>
!!   
!!   <tr> <td> out <th> thisFilter <td> SpectralFilter <td> 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_SpectralFilter( thisFilter )

   IMPLICIT NONE
   CLASS(SpectralFilter), INTENT(inout) :: thisFilter
   
      DEALLOCATE( thisFilter % filterMat )

#ifdef HAVE_CUDA
      DEALLOCATE( thisFilter % N_dev, &
                  thisfilter % filterMat_dev )
#endif

 END SUBROUTINE Trash_SpectralFilter
!
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
                                 ui = ui + thisFilter % filterMat(ii,i)*&
                                           f(ii,jj,kk,iVar,iEl)
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
         
 END FUNCTION Filter3D_SpectralFilter
 
#ifdef HAVE_CUDA
 ATTRIBUTES( Global ) SUBROUTINE Filter3D_CUDAKernel( f, filteredF, filterMatrix, N, nVariables, nElements )
   IMPLICIT NONE
   INTEGER, DEVICE, INTENT(in) :: N, nVariables, nElements
   REAL(prec), DEVICE, INTENT(in) :: f(0:N,0:N,0:N,1:nVariables,1:nElements)
   REAL(prec), DEVICE, INTENT(in) :: filterMatrix(0:N,0:N)
   REAL(prec), DEVICE, INTENT(out) :: filteredF(0:N,0:N,0:N,1:nVariables,1:nElements)
   ! Local
   INTEGER :: i,j,k,iVar,iEl
   INTEGER :: ii,jj,kk
   REAL(prec) :: uijk, uij, ui
   
   
     iVar = blockIdx % x
     iEl  = blockIdx % y
     
     i = threadIdx % x-1
     j = threadIdx % y-1
     k = threadIdx % z-1
     
     IF( i <= N .AND. j <= N .AND. k <= N )THEN
     
     uijk = 0.0_prec
    DO kk = 0, N
    
       uij = 0.0_prec
       DO jj = 0, N
          
          ui = 0.0_prec
          DO ii = 0, N
             ui = ui + filterMatrix(ii,i)*f(ii,jj,kk,iVar,iEl)
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
