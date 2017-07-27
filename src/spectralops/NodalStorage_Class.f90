! NodalStorage_Class.f90
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
! NodalStorage_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file NodalStorage_Class.f90
!! Contains the \ref NodalStorage_Class module, and <BR>
!! defines the \ref NodalStorage data-structure.


!> \defgroup NodalStorage_Class NodalStorage_Class 
!! This module defines the NodalStorage data-structure and its associated routines.

MODULE NodalStorage_Class

! src/common/
 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines
! src/interp/
 USE Quadrature
 USE Lagrange_Class

IMPLICIT NONE

!> \addtogroup NodalStorage_Class 
!! @{

!> \struct NodalStorage
!!  The NodalStorage class contains attributes needed for implementing spectral element methods
!!  in 3-D.
!!  
!!  An interpolant is formed that handles mapping between a computational "quadrature" mesh and a
!!  uniform plotting mesh. Quadrature (integration) weights are stored for use with Galerkin type
!!  methods. Galerkin derivative matrices (collocation derivative matrices weighted by appropriate
!!  ratios of the quadrature weights) are stored to facilitate the computation of weak derivatives.
!!  Finally, an interpolation matrix and accompanying subroutine is provided to interpolate 3-D data, 
!!  defined on the quadrature mesh, to the element boundaries.
!!
!! <H2> NodalStorage </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the spectral element method
!!       <tr> <th> nPlot <td> INTEGER <td> Number uniform plotting points
!!       <tr> <th> interp <td> Lagrange <td> Lagrange interpolant
!!       <tr> <th> qWeight(0:N) <td> REAL(prec) <td> Quadrature integration weights
!!       <tr> <th> dMatS(0:N,0:N) <td> REAL(prec) <td> Either the DG or CG derivative matrix
!!       <tr> <th> dMatP(0:N,0:N) <td> REAL(prec) <td> Either the DG or CG derivative matrix transpose
!!       <tr> <th> bMat(0:1,0:N) <td> REAL(prec) <td> Matrix for interpolating data to element boundaries
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref NodalStorage_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_NodalStorage
!!       <tr> <th> Trash <td> Trash_NodalStorage
!!       <tr> <th> CalculateAtBoundaries_1D <td> CalculateAtBoundaries_1D_NodalStorage
!!       <tr> <th> CalculateAtBoundaries_2D <td> CalculateAtBoundaries_2D_NodalStorage
!!       <tr> <th> CalculateAtBoundaries_3D <td> CalculateAtBoundaries_3D_NodalStorage
!!    </table>
!!

!>@}
   TYPE NodalStorage
      INTEGER                 :: N, nPlot
      TYPE(Lagrange)          :: interp      ! Lagrange interpolant
      REAL(prec), ALLOCATABLE :: qWeight(:)  ! Quadrature weights for integration
      REAL(prec), ALLOCATABLE :: dMatS(:,:)  ! Derivative matrix
      REAL(prec), ALLOCATABLE :: dMatP(:,:)  ! Derivative matrix
      REAL(prec), ALLOCATABLE :: bMat(:,:)   ! Matrix for interpolating functions to boundaries of 
                                             ! an element
      CONTAINS

      ! Manual Constructors/Destructors
      PROCEDURE :: Build => Build_NodalStorage
      PROCEDURE :: Trash => Trash_NodalStorage

      ! Type-Specific
      PROCEDURE :: CalculateAtBoundaries_1D => CalculateAtBoundaries_1D_NodalStorage
      PROCEDURE :: CalculateAtBoundaries_2D => CalculateAtBoundaries_2D_NodalStorage
      PROCEDURE :: CalculateAtBoundaries_3D => CalculateAtBoundaries_3D_NodalStorage
      
    END TYPE NodalStorage
    
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup NodalStorage_Class 
!! @{ 
! ================================================================================================ !
! S/R Build 
! 
!> \fn Build_NodalStorage 
!!  Allocates space fills values for the NodalStorage attributes using to the specified 
!!  quadrature and approximation form.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalStorage) :: this <BR>
!! <B>INTEGER</B>               :: N, nPlot
!!         .... <BR>
!!     ! To build a  structure for Continuous Galerkin with Gauss-Lobatto quadrature <BR>
!!     <B>CALL</B> this % Build( N, nPlot, GAUSS_LOBATTO, CG ) <BR>
!!
!!     ! To build a  structure for Discontinuous Galerkin with Gauss quadrature <BR>
!!     <B>CALL</B> this % Build( N, nPlot, GAUSS, DG ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myNodal <td> NodalStorage <td> On output, the attributes of the
!!                                                        NodalStorage data structure are filled
!!                                                        in.
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the method.
!!   <tr> <td> in <th> nPlot <td> INTEGER <td> The number of uniform plotting points in each 
!!                                             computational direction.
!!   <tr> <td> in <th> quadrature <td> INTEGER <td> A flag for specifying the desired type of 
!!                                                  quadrature. Can be set to either GAUSS or
!!                                                  GAUSS_LOBATTO. See \ref ModelFlags.f90 for
!!                                                  flag definitions.
!!   <tr> <td> in <th> approxForm <td> INTEGER <td> A flag for specifying the type of method that 
!!                                                  you are using. Can be set to either CG or DG.
!!                                                  See \ref ModelFlags.f90 for flag definitions.
!! 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_NodalStorage( myNodal, N, nPlot, quadrature, approxForm  )

   IMPLICIT NONE
   CLASS(NodalStorage), INTENT(out) :: myNodal
   INTEGER, INTENT(in)                 :: N
   INTEGER, INTENT(in)                 :: nPlot
   INTEGER, INTENT(in)                 :: quadrature
   INTEGER, INTENT(in)                 :: approxForm
   !LOCAL
   INTEGER                 :: i, j
   REAL(prec), ALLOCATABLE :: tempS(:), tempQ(:), tempUni(:)

      myNodal % N     = N
      myNodal % nPlot = nPlot
      
      ! Allocate space
      ALLOCATE( myNodal % dMatS(0:N,0:N), &
                myNodal % dMatP(0:N,0:N), &
                myNodal % bMat(0:N,0:1), &
                myNodal % qWeight(0:N) )
                
      myNodal % dMatS   = ZERO
      myNodal % dMatP   = ZERO
      myNodal % bMat    = ZERO
      myNodal % qWeight = ZERO

      ALLOCATE( tempS(0:N), tempQ(0:N), tempUni(0:nPlot) )
      tempS   = 0.0_prec
      tempUni = 0.0_prec
      tempQ   = 0.0_prec
      ! Generate the quadrature
      CALL LegendreQuadrature( N, tempS, tempQ, quadrature )
      myNodal % qWeight = tempQ      
       
      ! Build and store the interpolant
      tempUni = UniformPoints( -ONE, ONE, nPlot )

      CALL myNodal % interp % Build( N, nPlot, tempS, tempUni )
   
      ! Calculate and store the interpolants evaluated at the endpoints
      myNodal % bMat(0:N,0) = myNodal % interp % CalculateLagrangePolynomials( -ONE )
      myNodal % bMat(0:N,1) = myNodal % interp % CalculateLagrangePolynomials( ONE )
      
      IF( approxForm == CG )then ! Continuous Galerkin, store the derivative matrices as is

         myNodal % dMatS = myNodal % interp % D
         myNodal % dMatP = myNodal % interp % DTr

      ELSEIF( approxForm == DG )then
      
         ! For Discontinuous Galerkin, the matrix is transposed and multiplied by a ratio of quadrature
         ! weights.
        DO j = 0, N ! loop over the matrix rows
            DO i = 0, N ! loop over the matrix columns

               myNodal % dMatS(i,j) = -myNodal % interp % D(j,i)*&
                                         myNodal % qWeight(j)/&
                                         myNodal % qWeight(i)

            ENDDO
         ENDDO
         
         DO j = 0, N ! loop over the matrix rows
            DO i = 0, N ! loop over the matrix columns

               ! Here, we are purposefully using the transpose of the p-derivative matrix
               ! to conform with a new version of "MappedTimeDerivative"
               myNodal % dMatP(j,i) = myNodal % dMatS(i,j)

            ENDDO
         ENDDO
        
         
      ELSE

         PRINT*,'Module NodalStorage_Class.f90 : S/R BuildNodalStorage : Invalid SEM form. Stopping'
         STOP

      ENDIF
      
      DEALLOCATE( tempS, tempQ, tempUni )

 END SUBROUTINE Build_NodalStorage
!
!> \addtogroup NodalStorage_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_NodalStorage  
!! Frees memory held by the attributes of the NodalStorage class. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalStorage) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myNodal <td> NodalStorage <td>
!!                         On <B>input</B>, a NodalStorage class that has previously been 
!!                         constructed. <BR>
!!                         On <B>output</B>, the memory held by the attributes of this 
!!                         data-structure have been freed.
!!                                                           
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_NodalStorage( myNodal)

   IMPLICIT NONE
   CLASS(NodalStorage), INTENT(inout) :: myNodal

   CALL myNodal % interp % TRASH( )
   DEALLOCATE( myNodal % qWeight, myNodal % dMatS, myNodal % dMatP, myNodal % bMat )


 END SUBROUTINE Trash_NodalStorage
!
!
!==================================================================================================!
!--------------------------- Type Specific Routines -----------------------------------------------!
!==================================================================================================!
!
!
!
!
!> \addtogroup NodalStorage_Class 
!! @{ 
! ================================================================================================ !
! Function CalculateAtBoundaries 
! 
!> \fn CalculateAtBoundaries_NodalStorage
!!  Interpolates a 1-D array of nodal function values, defined on the quadrature mesh,
!!  to the boundaries of the computational element.
!!
!!  Recall that Lagrange-interpolation of data onto a point can be expressed as a matrix-vector 
!!  multiplication (in 1-D).  In 1-D, to interpolate to the boundaries of the computational domain
!!  ( \f$ \xi = -1,1 \f$ ), a \f$ 2 \times N \f$ matrix can be constructed whose columns are the 
!!  Lagrange interpolating polynomials evaluated at \f$ \xi = -1,1 \f$. An array of nodal function
!!  values at the interpolation nodes (identical to the quadrature nodes in a spectral element 
!!  method) can be interpolated to the boundaries through matrix vector multiplication.
!! 
!!  A single matrix-matrix multiplications effectively result in interpolation onto the boundaries.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalStorage) :: this <BR>
!! <B>REAL</B>(prec)         :: f(0:this%N) <BR>
!! <B>REAL</B>(prec)         :: fbound(1:2) <BR>
!!         .... <BR>
!!     fbound = this % CalculateAtBoundaries_1D( f ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myNodal <td> NodalStorage <td> Previously constructed NodalStorage data
!!                                                       structure. 
!!   <tr> <td> in <th> f(0:myNodal%N) <td> REAL(prec) <td> 
!!                     1-D array of nodal function values defined on the quadrature mesh.
!!   <tr> <td> out <th> fBound(1:2) <td> REAL(prec) <td> 
!!                     1-D array of nodal function values derived from interpolation onto the element
!!                     boundaries. The third index cycles over the sides of a 1-D element.
!!                     LEFT=1, RIGHT=2. See \ref ConstantsDictionary.f90 for more details on
!!                     boundary-integer flags.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION CalculateAtBoundaries_1D_NodalStorage( myNodal, f ) RESULT( fBound )

   IMPLICIT NONE
   CLASS(NodalStorage), INTENT(in) :: myNodal
   REAL(prec), INTENT(in)          :: f(0:myNodal % N)
   REAL(prec)                      :: fBound(1:2)
   ! Local
   INTEGER :: i, k   

         fBound(1:2) = 0.0_prec
         DO k = 0, myNodal % N
            ! Left
            fBound(1) = fBound(1) + myNodal % bMat(k,0)*f(k)
            ! Right
            fBound(2) = fBound(2) + myNodal % bMat(k,1)*f(k)
         ENDDO
      
 END FUNCTION CalculateAtBoundaries_1D_NodalStorage
!
!> \addtogroup NodalStorage_Class 
!! @{ 
! ================================================================================================ !
! Function CalculateAtBoundaries 
! 
!> \fn CalculateAtBoundaries_NodalStorage
!!  Interpolates a 2-D array of nodal function values, defined on the quadrature mesh,
!!  to the boundaries of the computational element.
!!
!!  Recall that Lagrange-interpolation of data onto a point can be expressed as a matrix-vector 
!!  multiplication (in 1-D).  In 1-D, to interpolate to the boundaries of the computational domain
!!  ( \f$ \xi = -1,1 \f$ ), a \f$ 2 \times N \f$ matrix can be constructed whose columns are the 
!!  Lagrange interpolating polynomials evaluated at \f$ \xi = -1,1 \f$. An array of nodal function
!!  values at the interpolation nodes (identical to the quadrature nodes in a spectral element 
!!  method) can be interpolated to the boundaries through matrix vector multiplication.
!! 
!!  In 2-D we need to interpolate \f$ 2(N+1) \f$ 1-D arrays to the boundaries. This routine views
!!  the nodal function data as an \f$ N+1 \times (N+1) \f$ matrix. Matrix-matrix multiplcation
!!  between the boundary interpolation matrix and the data-matrix results in boundary data in the
!!  first computational direction. Matrix-matrix multiplcation between the boundary interpolation 
!!  matrix and the data-matrix transpose results in boundary data in the second computational
!!  direction.
!!  Two matrix-matrix multiplications effectively result in interpolation onto the boundaries.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalStorage) :: this <BR>
!! <B>REAL</B>(prec)         :: f(0:this%N, 0:this%N) <BR>
!! <B>REAL</B>(prec)         :: fbound(0:this%N, 1:4) <BR>
!!         .... <BR>
!!     fbound = this % CalculateAtBoundaries_2D( f ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myNodal <td> NodalStorage <td> Previously constructed NodalStorage data
!!                                                       structure. 
!!   <tr> <td> in <th> f(0:myNodal%N, 0:myNodal%N) <td> REAL(prec) <td> 
!!                     2-D array of nodal function values defined on the quadrature mesh.
!!   <tr> <td> out <th> fBound(0:myNodal%N, 1:4) <td> REAL(prec) <td> 
!!                     2-D array of nodal function values derived from interpolation onto the element
!!                     boundaries. The third index cycles over the sides of a 2-D element.
!!                     SOUTH=1, EAST=2, NORTH=3, WEST=4. See \ref ConstantsDictionary.f90
!!                     for more details on boundary-integer flags.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION CalculateAtBoundaries_2D_NodalStorage( myNodal, f ) RESULT( fBound )

   IMPLICIT NONE
   CLASS(NodalStorage) :: myNodal
   REAL(prec)          :: f(0:myNodal % N, 0:myNodal % N)
   REAL(prec)          :: fBound(0:myNodal % N, 1:nQuadEdges)
   ! Local
   INTEGER :: i, k   

      DO i = 0, myNodal % N
         fBound(i,1:4) = 0.0_prec
         DO k = 0, myNodal % N
            
            ! South
            fBound(i,South) = fBound(i,South) + myNodal % bMat(k,0)*f(i,k)
            !West
            fBound(i,West) = fBound(i,West) + myNodal % bMat(k,0)*f(k,i)
            !East
            fBound(i,East) = fBound(i,East) + myNodal % bMat(k,0)*f(k,i)
            !North
            fBound(i,North) = fBound(i,North) + myNodal % bMat(k,0)*f(i,k)
               
         ENDDO
      ENDDO
      
      
 END FUNCTION CalculateAtBoundaries_2D_NodalStorage
!
!> \addtogroup NodalStorage_Class 
!! @{ 
! ================================================================================================ !
! Function CalculateAtBoundaries_3D 
! 
!> \fn CalculateAtBoundaries_NodalStorage
!!  Interpolates a 3-D array of nodal function values, defined on the quadrature mesh,
!!  to the boundaries of the computational element.
!!
!!  Recall that Lagrange-interpolation of data onto a point can be expressed as a matrix-vector 
!!  multiplication (in 1-D).  In 1-D, to interpolate to the boundaries of the computational domain
!!  ( \f$ \xi = -1,1 \f$ ), a \f$ 2 \times N \f$ matrix can be constructed whose columns are the 
!!  Lagrange interpolating polynomials evaluated at \f$ \xi = -1,1 \f$. An array of nodal function
!!  values at the interpolation nodes (identical to the quadrature nodes in a spectral element 
!!  method) can be interpolated to the boundaries through matrix vector multiplication.
!! 
!!  In 3-D we need to interpolate \f$ 3(N+1)^2 \f$ 1-D arrays to the boundaries. This routine forms
!!  3 matrices that are \f$ N+1 \times (N+1)^2 \f$ by collapsing pairs of array indices. The columns
!!  of each matrix correspond to the 1-D nodal function values in a given computational direction.
!!  Three matrix-matrix multiplications effectively result in interpolation onto the boundaries.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalStorage) :: this <BR>
!! <B>REAL</B>(prec)         :: f(0:this%N, 0:this%N, 0:this%N) <BR>
!! <B>REAL</B>(prec)         :: fbound(0:this%N, 0:this%N, 1:6) <BR>
!!         .... <BR>
!!     fbound = this % CalculateAtBoundaries_3D( f ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myNodal <td> NodalStorage <td> Previously constructed NodalStorage data
!!                                                       structure. 
!!   <tr> <td> in <th> f(0:myNodal%N, 0:myNodal%N, 0:myNodal%N) <td> REAL(prec) <td> 
!!                     3-D array of nodal function values defined on the quadrature mesh.
!!   <tr> <td> out <th> fBound(0:myNodal%N, 0:myNodal%N, 1:6) <td> REAL(prec) <td> 
!!                     3-D array of nodal function values derived from interpolation onto the element
!!                     boundaries. The third index cycles over the sides of a 3-D element.
!!                     SOUTH=1, EAST=2, NORTH=3, WEST=4, BOTTOM=5, TOP=6. See \ref ConstantsDictionary.f90
!!                     for more details on boundary-integer flags.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION CalculateAtBoundaries_3D_NodalStorage( myNodal, f ) RESULT( fBound )

   IMPLICIT NONE
   CLASS(NodalStorage) :: myNodal
   REAL(prec)          :: f(0:myNodal % N, 0:myNodal % N, 0:myNodal % N)
   REAL(prec)          :: fBound(0:myNodal % N, 0:myNodal % N, 1:6)
   ! Local
   INTEGER :: i, j, k   

      DO j = 0, myNodal % N
         DO i = 0, myNodal % N
            fBound(i,j,1:6) = 0.0_prec
            DO k = 0, myNodal % N
            
               ! South
               fBound(i,j,South) = fBound(i,j,South) + myNodal % bMat(k,0)*f(k,i,j)
               ! North
               fBound(i,j,North) = fBound(i,j,North) + myNodal % bMat(k,1)*f(k,i,j)
               ! West
               fBound(i,j,West) = fBound(i,j,West) + myNodal % bMat(k,0)*f(i,k,j)
               ! East
               fBound(i,j,East) = fBound(i,j,East) + myNodal % bMat(k,1)*f(i,k,j)
               ! Bottom
               fBound(i,j,Bottom) = fBound(i,j,West) + myNodal % bMat(k,0)*f(i,j,k)
               ! Top
               fBound(i,j,Top) = fBound(i,j,East) + myNodal % bMat(k,1)*f(i,j,k)
            ENDDO
         ENDDO
      ENDDO

 END FUNCTION CalculateAtBoundaries_3D_NodalStorage

END MODULE NodalStorage_Class
