! DGSEMSolutionStorage_2D_Class.f90
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
! DGSEMSolutionStorage_2D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file DGSEMSolutionStorage_2D_Class.f90
!! Contains the \ref DGSEMSolutionStorage_2D_Class module, and <BR>
!! defines the \ref DGSEMSolution_2D data-structure.

!> \defgroup DGSEMSolutionStorage_2D_Class DGSEMSolutionStorage_2D_Class 
!! This module defines the DGSEMSolution_2D data-structure and its associated routines.
 
MODULE DGSEMSolutionStorage_2D_Class
 
! src/common/
USE ConstantsDictionary
USE ModelPrecision
! src/nodal
USE NodalStorage_Class

IMPLICIT NONE

!> \addtogroup DGSEMSolutionStorage_2D_Class 
!! @{

!> \struct DGSEMSolution_2D
!!  The DGSEMSolution_2D class provides attributes for storing a solution and its flux on a
!!  single spectral element.
!!  
!!  When implement a Discontinuous Galerkin Spectral Element method, it is common practice to use
!!  the Gauss quadrature points within each element. To advance the discrete system, fluxes through
!!  the element boundaries must be known. 
!!  This data structure provides attributes for storing the solution at the mesh points 
!!  and at the element boundaries and the fluxes at the element boundaries. Additionally, an array
!!  is provided for storing the solution tendency that can be used in a time-stepping routine. 
!!
!!  A wrapper routine is also provided for interpolating the solution to the element boundaries by
!!  calling a routine from the \ref NodalStorage_2D_Class module.
!!
!! <H2> DGSEMSolution_2D </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the spectral element method
!!       <tr> <th> nEq <td> INTEGER <td> Number of (prognostic) solution variables.
!!       <tr> <th> solution(0:N,0:N,1:nEq) <td> REAL(prec) <td> An array containing the solution variables
!!       <tr> <th> tendency(0:N,0:N,1:nEq) <td> REAL(prec) <td> An array containing the tendency of the 
!!       <tr> <th> tendency(0:N,0:N,1:nEq) <td> REAL(prec) <td> An array containing the tendency of the 
!!                                                          solution variables.
!!       <tr> <th> boundarySolution(0:N,1:4,1:nEq)* <td> REAL(prec) <td>
!!                  An array containing the solution variables at the element boundary
!!       <tr> <th> boundaryFlux(0:N,1:4,1:nEq)* <td> REAL(prec) <td>
!!                  An array containing the flux at the element boundary
!!    </table>
!!
!!  *For the "boundary" arrays, the sides for a quadrilateral element are numbered as SOUTH=1, 
!!  EAST=2, NORTH=3, WEST=4.
!!
!! <H3> Procedures </H3>
!!    See \ref DGSEMSolutionStorage_2D_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_DGSEMSolution_2D
!!       <tr> <th> Trash <td> Trash_DGSEMSolution_2D
!!       <tr> <th> CalculateBoundarySolution <td> CalculateBoundarySolution_DGSEMSolution_2D
!!    </table>
!!

!>@}
    TYPE DGSEMSolution_2D
      INTEGER                 :: nEq, N
      REAL(prec), allocatable :: solution(:,:,:)
      REAL(prec), allocatable :: tendency(:,:,:)
      REAL(prec), allocatable :: boundarySolution(:,:,:) ! Indexed over the boundaries. 
      REAL(prec), allocatable :: boundaryFlux(:,:,:)     ! Indexed over the boundaries

      CONTAINS

      ! CONSTRUCTORS/DESTRUCTORS
      PROCEDURE :: Build => Build_DGSEMSolution_2D
      PROCEDURE :: Trash => Trash_DGSEMSolution_2D

      PROCEDURE :: CalculateBoundarySolution => CalculateBoundarySolution_DGSEMSolution_2D
      
    END TYPE DGSEMSolution_2D

 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup DGSEMSolutionStorage_2D_Class
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_DGSEMSolution_2D_Class  
!! Allocates space for the DGSEMSolution_2D attributes and initializes arrays to zero.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DGSEMSolution_2D) :: this <BR>
!! <B>INTEGER</B>                :: N, nEq <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( N, nEq ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myDGS <td> DGSEMSolution_2D <td> On output, memory has been allocated for
!!                                                       each attribute and each array is initialized
!!                                                       with a value of zero. 
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the DG-method you plan on using
!!   <tr> <td> in <th> nEq <td> INTEGER <td> Number of prognostic variables; number of equations for
!!                                           the system being solved
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_DGSEMSolution_2D( myDGS, N, nEq )

   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: N, nEq

      myDGS % nEq    = nEq
      myDGS % N      = N
      
      ALLOCATE( myDGS % solution(0:N,0:N,1:nEq) )
      ALLOCATE( myDGS % tendency(0:N,0:N,1:nEq) )
      ALLOCATE( myDGS % boundarySolution(0:N,1:nEq,1:4) ) 
      ALLOCATE( myDGS % boundaryFlux(0:N,1:nEq,1:4) ) 

      myDGS % solution = ZERO
      myDGS % tendency = ZERO
      myDGS % boundarySolution = ZERO
      myDGS % boundaryFlux = ZERO


 END SUBROUTINE Build_DGSEMSolution_2D
!> \addtogroup DGSEMSolutionStorage_2D_Class
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_DGSEMSolution_2D
!! Frees memory held by the attributes of the DGSEMSolution_2D data structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DGSEMSolution_2D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myDGS <td> DGSEMSolution_2D <td> On output, memory held by the attributes
!!                                                          of the data structure has been deallocated.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_DGSEMSolution_2D( myDGS )

   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS

      DEALLOCATE( myDGS % Solution )
      DEALLOCATE( myDGS % tendency )
      DEALLOCATE( myDGS % boundarySolution ) 
      DEALLOCATE( myDGS % boundaryFlux ) 

 END SUBROUTINE Trash_DGSEMSolution_2D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup DGSEMSolutionStorage_2D_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateBoundarySolution
! 
!> \fn CalculateBoundarySolution_DGSEMSolution_2D
!! Using the "solution" attribute, this routine calculates the solution at the element boundaries
!! and stores it in the "boundarySolution" data structure.
!! 
!!  This subroutine depends on
!!   Module \ref NodalStorage_2D_Class : Function \ref CalculateBoundarySolution_NodalStorage_2D
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DGSEMSolution_2D) :: this <BR>
!! <B>TYPE</B>(NodalStorage)     :: dgStorage <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateBoundarySolution( dgStorage ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myDGS <td> DGSEMSolution_2D <td> 
!!                         On <B>input</B>, the solution attribute is filled in, <BR>
!!                         On <B>output</B>, the boundarySolution attribute is updated by 
!!                         interpolating the solution to the element boundaries.
!!   <tr> <td> in <th> dgStorage <td> NodalStorage <td>
!!                     Structure that defines the boundary interpolation matrix that corresponds
!!                     to mapping the solution from quadrature points to the element boundaries
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE CalculateBoundarySolution_DGSEMSolution_2D( myDGS, dgStorage )

   IMPLICIT NONE
   CLASS(DGSEMSolution_2D), INTENT(inout) :: myDGS
   TYPE(NodalStorage), INTENT(in)         :: dgStorage
   ! Local
   INTEGER :: iEq
   INTEGER :: i, k   

       
      DO iEq = 1, myDGS % nEq
         DO i = 0, myDGS % N
            myDGS % boundarySolution(i,iEq,1:4) = 0.0_prec
            DO k = 0, myDGS % N
            
               ! South
               myDGS % boundarySolution(i,iEq,South) = myDGS % boundarySolution(i,iEq,South) + &
                                                       dgStorage % bMat(k,0)*myDGS % solution(i,k,iEq)
               !West
               myDGS % boundarySolution(i,iEq,West) = myDGS % boundarySolution(i,iEq,West) + &
                                                      dgStorage % bMat(k,0)*myDGS % solution(k,i,iEq)
               !East
               myDGS % boundarySolution(i,iEq,East) = myDGS % boundarySolution(i,iEq,East) + &
                                                      dgStorage % bMat(k,1)*myDGS % solution(k,i,iEq)
               !North
               myDGS % boundarySolution(i,iEq,North) = myDGS % boundarySolution(i,iEq,North) + &
                                                       dgStorage % bMat(k,1)*myDGS % solution(i,k,iEq)
            ENDDO
         ENDDO
      ENDDO
      
 END SUBROUTINE CalculateBoundarySolution_DGSEMSolution_2D

END MODULE DGSEMSolutionStorage_2D_Class
