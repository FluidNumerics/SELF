! RollOffFilter_Class.f90
! 
! Copyright 2017,2018 Joseph Schoonover <schoonover.numerics@gmail.com>
!
! RollOffFilter_Class.f90 is part of the SELF-Fluids.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file RollOffFilter_Class.f90
!! Contains the \ref RollOffFilter_Class module, and <BR>
!! defines the \ref RollOffFilter data-structure.


!> \defgroup RollOffFilter_Class RollOffFilter_Class 
!! This module defines the RollOffFilter data-structure and its associated routines.

MODULE RollOffFilter_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
! src/spectralops/
USE Quadrature


IMPLICIT NONE

!> \addtogroup RollOffFilter_Class 
!! @{

!> \struct RollOffFilter
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
!! <H2> RollOffFilter </H2>
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
!!    See \ref RollOffFilter_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_RollOffFilter
!!       <tr> <th> Trash <td> Trash_RollOffFilter
!!       <tr> <th> Apply1DFilter <td> Apply1DFilter_RollOffFilter
!!       <tr> <th> Apply2DFilter <td> Apply2DFilter_RollOffFilter
!!       <tr> <th> Apply3DFilter <td> Apply3DFilter_RollOffFilter
!!    </table>
!!

!>@}

   TYPE RollOffFilter
      INTEGER                 :: N, nCutoff, nPacked
      REAL(prec), ALLOCATABLE :: filterMat(:,:)
      
      CONTAINS
  
      PROCEDURE :: Build => Build_RollOffFilter
      PROCEDURE :: Trash => Trash_RollOffFilter
      PROCEDURE :: Apply1DFilter => Apply1DFilter_RollOffFilter
      PROCEDURE :: Apply2DFilter => Apply2DFilter_RollOffFilter
      PROCEDURE :: Apply3DFilter => Apply3DFilter_RollOffFilter

   END TYPE RollOffFilter


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup RollOffFilter_Class 
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_RollOffFilter 
!! Allocates space for the modal cutoff filter and initializes the attributes of the data structure 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(RollOffFilter) :: this <BR>
!! <B>INTEGER</B>                 :: N, nCutoff <BR>
!! <B>REAL</B>(prec)              :: s(0:N), w(0:N) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( s, w, N, nCutoff ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> thisFilter <td> RollOffFilter <td> 
!!   <tr> <td> in <th> s(0:N) <td> REAL(prec) <td> The interpolation nodes of your data
!!   <tr> <td> in <th> w(0:N) <td> REAL(prec) <td> The quadrature weights for discrete integration
!!                                                 associated with the interpolation nodes
!!   <tr> <td> in <th> N <td> REAL(prec) <td> Polynomial degree of the data
!!   <tr> <td> in <th> nCutoff <td> REAL(prec) <td> Cutoff polynomial degree
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_RollOffFilter( thisFilter, s, w, N, nCutoff )

   IMPLICIT NONE
   CLASS(RollOffFilter), INTENT(inout) :: thisFilter
   INTEGER, INTENT(in)                 :: N, nCutoff
   REAL(prec), INTENT(in)              :: s(0:N)
   REAL(prec), INTENT(in)              :: w(0:N)
   ! Local
   REAL(prec)              :: Lnorm, Li, Ls, sc, alpha, r
   REAL(prec), ALLOCATABLE :: Pfilt(:,:), V(:,:), VInv(:,:)
   INTEGER                 :: row, col 
   
      thisFilter % N       = N
      thisFilter % nCutoff = nCutoff
      thisFilter % nPacked = (N+1)*(N+1)-1

      ALLOCATE( thisFilter % filterMat(0:N,0:N) )

      ALLOCATE( Pfilt(0:N,0:N), V(0:N,0:N), VInv(0:N,0:N) )

      thisFilter % filterMat = ZERO

      Pfilt = ZERO
      V     = ZERO
      VInv  = ZERO

      sc = real(nCutoff,prec)
      alpha = log(TWO)/(sc*sc)
      DO row = 0, N 
 
         r = real(row,prec)
         !Pfilt(row,row) = exp( -alpha*r**2 )
         Pfilt(row,row) = 0.5_prec*(1.0_prec - tanh( (r-nCutoff) ) )
         Lnorm = ZERO

         DO col = 0, N

            CALL LegendrePolynomial( row, s(col), Li, Ls)

            Lnorm = Lnorm + Li*Li*w(col)
            V(row,col) = Li*w(col)
            ! The inverse of the modal projection matrix is easy to build since the Legendre basis
            ! is an orthogonal basis
            VInv(col,row) = Li
              
         ENDDO
            V(row,0:N) = V(row,0:N)/Lnorm
      ENDDO

      Pfilt = MATMUL( Pfilt, V )
      thisFilter % filterMat = MATMUL( VInv, Pfilt )
    
      DEALLOCATE( V, Pfilt, Vinv )

 END SUBROUTINE Build_RollOffFilter
!
!> \addtogroup RollOffFilter_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_RollOffFilter 
!! Frees memory associated with the modal cutoff filter 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(RollOffFilter) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> thisFilter <td> RollOffFilter <td> 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_RollOffFilter( thisFilter )

   IMPLICIT NONE
   CLASS(RollOffFilter), INTENT(inout) :: thisFilter
   
      DEALLOCATE( thisFilter % filterMat )

 END SUBROUTINE Trash_RollOffFilter
!
!
!==================================================================================================!
!------------------------------------- Type Specific ----------------------------------------------!
!==================================================================================================!
!
! 
!> \addtogroup RollOffFilter_Class 
!! @{ 
! ================================================================================================ !
! S/R Apply3DFilter
! 
!> \fn Apply3DFilter_RollOffFilter 
!! Allocates space for the modal cutoff filter and initializes the attributes of the data structure 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(RollOffFilter) :: this <BR>
!! <B>REAL</B>(prec)              :: sol(0:this % N,0:this % N,0:this % N) <BR>
!! <B>REAL</B>(prec)              :: fsol(0:this % N,0:this % N,0:this % N) <BR>
!!         .... <BR>
!!     fSol = this % Apply3DFilter( sol ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> thisFilter <td> RollOffFilter <td> Previously constructed modal cutoff filter
!!   <tr> <td> in <th> sol(0:this % N,0:this % N,0:this % N) <td> REAL(prec) <td>
!!                     3-D array of data at the interpolation nodes associated with the filter 
!!   <tr> <td> out <th> fsol(0:this % N,0:this % N,0:this % N) <td> REAL(prec) <td> 
!!                     3-D array of filtered data at the interpolation nodes associated with the filter
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION Apply3DFilter_RollOffFilter( thisFilter, sol ) RESULT(fSol)

   IMPLICIT NONE
   CLASS(RollOffFilter) :: thisFilter
   REAL(prec)           :: sol(0:thisFilter % N, 0:thisFilter % N, 0:thisFilter % N)
   REAL(prec)           :: fSol(0:thisFilter % N, 0:thisFilter % N, 0:thisFilter % N)
   ! LOCAL 
   REAL(prec) :: floc(0:thisFilter % N, 0:thisFilter % nPacked)
   REAL(prec) :: fint(0:thisFilter % N, 0:thisFilter % nPacked)
   REAL(prec) :: fmat(0:thisFilter % N, 0:thisFilter % N)
   REAL(prec) :: fUnpacked(0:thisFilter % N, 0:thisFilter % N, 0:thisFilter % N)
   INTEGER    :: N, i, j, k 

      N = thisFilter % N
      fmat = thisFilter % filterMat
      
      ! Pack down the last two dimensions of the solution array
      DO j = 0, N
         DO i = 0, N
            k = i + j*(N+1)
            floc(0:N,k) = sol(0:N,i,j)
         ENDDO
      ENDDO
      
      ! Apply the filter matrix to filter the first computational dimension
      fint = MATMUL( fmat, floc )
      ! Unpack array
      DO j = 0, N
         DO i = 0, N
            k = i + j*(N+1)
            fUnpacked(0:N,i,j) = fint(0:N,k)
         ENDDO
      ENDDO
      ! Repack array along 2nd dimension
      DO j = 0, N
         DO i = 0, N
            k = i + j*(N+1)
            floc(0:N,k) = fUnpacked(i,0:N,j)
         ENDDO
      ENDDO
      
      ! Apply the filter matrix to filter the first computational dimension
      fint = MATMUL( fmat, floc )
      ! Unpack array
      DO j = 0, N
         DO i = 0, N
            k = i + j*(N+1)
            fUnpacked(i,0:N,j) = fint(0:N,k)
         ENDDO
      ENDDO
      ! Repack array along 3rd dimension
      DO j = 0, N
         DO i = 0, N
            k = i + j*(N+1)
            floc(0:N,k) = fUnpacked(i,j,0:N)
         ENDDO
      ENDDO
      
      ! Apply the filter matrix to filter the first computational dimension
      fint = MATMUL( fmat, floc )
      ! Unpack array
      DO j = 0, N
         DO i = 0, N
            k = i + j*(N+1)
            fUnpacked(i,j,0:N) = fint(0:N,k)
         ENDDO
      ENDDO
      
      fSol = fUnpacked
      
 END FUNCTION Apply3DFilter_RollOffFilter
!
!> \addtogroup RollOffFilter_Class 
!! @{ 
! ================================================================================================ !
! S/R Apply2DFilter
! 
!> \fn Apply2DFilter_RollOffFilter 
!! Allocates space for the modal cutoff filter and initializes the attributes of the data structure 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(RollOffFilter) :: this <BR>
!! <B>REAL</B>(prec)          :: sol(0:this % N,0:this % N) <BR>
!! <B>REAL</B>(prec)          :: fsol(0:this % N,0:this % N) <BR>
!!         .... <BR>
!!     fSol = this % Apply2DFilter( sol ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> thisFilter <td> RollOffFilter <td> Previously constructed modal cutoff filter
!!   <tr> <td> in <th> sol(0:this % N,0:this % N) <td> REAL(prec) <td>
!!                     2-D array of data at the interpolation nodes associated with the filter 
!!   <tr> <td> out <th> fsol(0:this % N,0:this % N) <td> REAL(prec) <td> 
!!                     2-D array of filtered data at the interpolation nodes associated with the filter
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION Apply2DFilter_RollOffFilter( thisFilter, sol ) RESULT(fSol)

   IMPLICIT NONE
   CLASS(RollOffFilter) :: thisFilter
   REAL(prec)               :: sol(0:thisFilter % N, 0:thisFilter % N)
   REAL(prec)               :: fSol(0:thisFilter % N, 0:thisFilter % N)
   ! LOCAL 
   REAL(prec) :: floc(0:thisFilter % N, 0:thisFilter % N)
   REAL(prec) :: fint(0:thisFilter % N, 0:thisFilter % N)
   REAL(prec) :: fmat(0:thisFilter % N, 0:thisFilter % N)

      fmat = thisFilter % filterMat

      floc = sol
      ! Apply the filter matrix to filter the first computational dimension
      fint = MATMUL( fmat, floc )
      ! Flip the dimension of the "intermediate" array to prepare for application of filtering
      ! matrix to the second dimension
      floc = TRANSPOSE( fint )
      ! Apply the filter matrix to the second dimension
      fint = MATMUL( fmat, floc )
      ! Transpose again to obtain correct row-column ordering for the filtered solution
      fSol = TRANSPOSE( fint )
      
 END FUNCTION Apply2DFilter_RollOffFilter
!
!> \addtogroup RollOffFilter_Class 
!! @{ 
! ================================================================================================ !
! S/R Apply1DFilter
! 
!> \fn Apply1DFilter_RollOffFilter 
!! Allocates space for the modal cutoff filter and initializes the attributes of the data structure 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(RollOffFilter) :: this <BR>
!! <B>REAL</B>(prec)          :: sol(0:this % N) <BR>
!! <B>REAL</B>(prec)          :: fsol(0:this % N) <BR>
!!         .... <BR>
!!     fSol = this % Apply1DFilter( sol ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> thisFilter <td> RollOffFilter <td> Previously constructed modal cutoff filter
!!   <tr> <td> in <th> sol(0:this % N) <td> REAL(prec) <td>
!!                     1-D array of data at the interpolation nodes associated with the filter 
!!   <tr> <td> out <th> fsol(0:this % N) <td> REAL(prec) <td> 
!!                     1-D array of filtered data at the interpolation nodes associated with the filter
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION Apply1DFilter_RollOffFilter( thisFilter, sol ) RESULT(fSol)

   IMPLICIT NONE
   CLASS(RollOffFilter) :: thisFilter
   REAL(prec)           :: sol(0:thisFilter % N)
   REAL(prec)           :: fSol(0:thisFilter % N)
   ! LOCAL 
   REAL(prec) :: floc(0:thisFilter % N)
   REAL(prec) :: fmat(0:thisFilter % N,0:thisFilter % N)

      fmat = thisFilter % filterMat
      floc = sol
      
      ! Apply the filter matrix to filter the first computational dimension
      fSol = MATMUL( fmat, floc )
      
 END FUNCTION Apply1DFilter_RollOffFilter 
!
END MODULE RollOffFilter_Class
