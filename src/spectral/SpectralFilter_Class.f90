! SpectralFilter_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE SpectralFilter_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
! src/spectralops/
USE Quadrature


IMPLICIT NONE

!> \addtogroup SpectralFilter_Class 
!! @{

!> \struct SpectralFilter
!! A data-structure for handling Legendre Modal Filtering in 1, 2, and 3 dimensions using a "Roll-Off"
!! Filter
!!
!!  This module provides a data structure for constructing and storing element-local filter matrices
!!  that can be used for polynomial de-aliasing or as part of SGS parameterizations. 1, 2, or 3-D
!!  arrays,representative of the nodal values of an interpolant, can be passed in to the supplied 
!!  procedures "Apply1DFilter", "Apply2DFilter", or "Apply3DFilter"  to return a filtered form of 
!!  the interpolant.  
!!
!!   This module was inspired by the paper
!! 
!!    D. Flad, A. Beck, and C. Munz, (2016) "Simulation of underresolved turbulent flows by adaptive 
!!       filtering using the high order discontinuous Galerkin spectral element method", JCP, 313, 
!!       pp. 1-12
!!
!!
!!   The data-structure provided here can be used with a high-end routine in order to incorporate
!!   "switch-based" filtering of prognostic solution variables. For under-resolved simulations, such 
!!   filtering can provide fine-tuned control of dissipation necessary for stability.
!!
!!   It is assumed that data is described by the same polynomial degree in each computational direction
!!   and that the filtering is performed using the same cutoff degree in each computational direction.
!!
!! <H2> SpectralFilter </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree associated with the filter
!!       <tr> <th> nCutoff <td> INTEGER <td> Cutoff polynomial degree indicating which Legendre 
!!                                           modal coefficients are made null
!!       <tr> <th> nPacked <td> REAL(prec) <td> Upper bound of 2-D array obtained from packing down
!!                                              3-D array; used for filtering 3-D data
!!       <tr> <th> filterMat(:,:) <td> REAL(prec) <td> The filtering matrix
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref SpectralFilter_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_SpectralFilter
!!       <tr> <th> Trash <td> Trash_SpectralFilter
!!       <tr> <th> Apply1DFilter <td> Apply1DFilter_SpectralFilter
!!       <tr> <th> Apply2DFilter <td> Apply2DFilter_SpectralFilter
!!       <tr> <th> Apply3DFilter <td> Apply3DFilter_SpectralFilter
!!    </table>
!!

!>@}

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
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(SpectralFilter) :: this <BR>
!! <B>INTEGER</B>                 :: N, nCutoff <BR>
!! <B>REAL</B>(prec)              :: s(0:N), w(0:N) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( s, w, N, nCutoff ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
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
   REAL(prec)              :: Lnorm, Li, Ls, sc, alpha, r
   REAL(prec), ALLOCATABLE :: Pfilt(:,:), V(:,:), VInv(:,:)
   INTEGER                 :: row, col 
   
      thisFilter % N       = N

      ALLOCATE( thisFilter % filterMat(0:N,0:N) )

      ALLOCATE( Pfilt(0:N,0:N), V(0:N,0:N), VInv(0:N,0:N) )

      thisFilter % filterMat = 0.0_prec

      Pfilt = 0.0_prec
      V     = 0.0_prec
      VInv  = 0.0_prec

      sc = real(nCutoff,prec)
      alpha = log(TWO)/(sc*sc)
      DO row = 0, N 
 
         r = real(row,prec)

         IF( filterType == ModalCutoff )THEN
           IF( row <= nCutoff )THEN
              Pfilt(row,row) = 1.0_prec
           ENDIF
         ELSEIF( filterType == TanhRollOff )THEN
           Pfilt(row,row) = 0.5_prec*(1.0_prec - tanh( (r-nCutoff) ) )
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

    
      DEALLOCATE( V, Pfilt, Vinv )

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
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(SpectralFilter) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
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
!
!==================================================================================================!
!------------------------------------- Type Specific ----------------------------------------------!
!==================================================================================================!
!
! 
END MODULE SpectralFilter_Class
