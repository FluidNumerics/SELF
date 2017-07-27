! DGSEMSolutionStorage_1D_Class.f90
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
! DGSEMSolutionStorage_1D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file DGSEMSolutionStorage_1D_Class.f90
!! Contains the \ref DGSEMSolutionStorage_1D_Class module, and <BR>
!! defines the \ref DGSEMSolution_1D data-structure.

!> \defgroup DGSEMSolutionStorage_1D_Class DGSEMSolutionStorage_1D_Class 
!! This module defines the DGSEMSolution_1D data-structure and its associated routines.

MODULE DGSEMSolutionStorage_1D_Class

! src/common/
USE CoNtantsDictionary
USE ModelPrecision
! src/nodal/
USE DGSEMSolutionStorage_Class

IMPLICIT NONE

!> \addtogroup DGSEMSolutionStorage_1D_Class 
!! @{

!> \struct DGSEMSolution_1D
!!  The DGSEMSolution_1D class provides attributes for storing a solution and its flux on a
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
!!  calling a routine from the \ref NodalStorage_1D_Class module.
!!
!! <H2> DGSEMSolution_1D </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the spectral element method
!!       <tr> <th> nEq <td> INTEGER <td> Number of (prognostic) solution variables.
!!       <tr> <th> solution(0:N,1:nEq) <td> REAL(prec) <td> An array containing the solution variables
!!       <tr> <th> tendency(0:N,1:nEq) <td> REAL(prec) <td> An array containing the tendency of the 
!!                                                          solution variables.
!!       <tr> <th> boundarySolution(1:nEq,1:2) <td> REAL(prec) <td>
!!                  An array containing the solution variables at the element boundary
!!       <tr> <th> boundaryFlux(1:nEq,1:2) <td> REAL(prec) <td>
!!                  An array containing the flux at the element boundary
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref DGSEMSolutionStorage_1D_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_DGSEMSolution_1D
!!       <tr> <th> Trash <td> Trash_DGSEMSolution_1D
!!       <tr> <th> CalculateBoundarySolution <td> CalculateBoundarySolution_DGSEMSolution_1D
!!    </table>
!!

!>@}
   TYPE DGSEMSolution_1D
      INTEGER                 :: nEq, N
      REAL(prec), ALLOCATABLE :: solution(:,:)
      REAL(prec), ALLOCATABLE :: tendency(:,:)
      REAL(prec), ALLOCATABLE :: boundarySolution(:,:)
      REAL(prec), ALLOCATABLE :: boundaryFlux(:,:)

      CONTAINS

      ! CoNtructors/Destructors
      PROCEDURE :: Build => Build_DGSEMSolution_1D
      PROCEDURE :: Trash => Trash_DGSEMSolution_1D
      ! Type/Specific
      PROCEDURE :: CalculateBoundarySolution => CalculateBoundarySolution_DGSEMSolution_1D
      
   END TYPE DGSEMSolution_1D

 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup DGSEMSolutionStorage_1D_Class
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_DGSEMSolution_1D_Class  
!! Allocates space for the DGSEMSolution_1D attributes and initializes arrays to zero.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DGSEMSolution_1D) :: this <BR>
!! <B>INTEGER</B>                :: N, nEq <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( N, nEq ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myDGS <td> DGSEMSolution_1D <td> On output, memory has been allocated for
!!                                                       each attribute and each array is initialized
!!                                                       with a value of zero. 
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the DG-method you plan on using
!!   <tr> <td> in <th> nEq <td> INTEGER <td> Number of prognostic variables; number of equations for
!!                                           the system being solved
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_DGSEMSolution_1D( myDGS, N, nEq )

   IMPLICIT NONE
   CLASS(DGSEMSolution_1D), INTENT(out) :: myDGS
   INTEGER, INTENT(in)                  :: N, nEq

      myDGS % nEq = nEq
      myDGS % N   = N
      
      ALLOCATE( myDGS % solution(0:N,1:nEq) )
      ALLOCATE( myDGS % tendency(0:N,1:nEq) )
      ALLOCATE( myDGS % boundarySolution(1:nEq,1:2) ) 
      ALLOCATE( myDGS % boundaryFlux(1:nEq,1:2) ) 
      
      myDGS % solution = ZERO
      myDGS % tendency = ZERO
      myDGS % boundarySolution = ZERO
      myDGS % boundaryFlux = ZERO
      
 END SUBROUTINE Build_DGSEMSolution_1D
!
!> \addtogroup DGSEMSolutionStorage_1D_Class
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_DGSEMSolution_1D
!! Frees memory held by the attributes of the DGSEMSolution_1D data structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DGSEMSolution_1D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myDGS <td> DGSEMSolution_1D <td> On output, memory held by the attributes
!!                                                          of the data structure has been deallocated.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_DGSEMSolution_1D( myDGS )

   IMPLICIT NONE
   CLASS(DGSEMSolution_1D), INTENT(inout) :: myDGS

      DEALLOCATE( myDGS % solution )
      DEALLOCATE( myDGS % tendency )
      DEALLOCATE( myDGS % boundarySolution ) 
      DEALLOCATE( myDGS % boundaryFlux ) 

 END SUBROUTINE Trash_DGSEMSolution_1D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup DGSEMSolutionStorage_1D_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateBoundarySolution
! 
!> \fn CalculateBoundarySolution_DGSEMSolution_1D
!! Using the "solution" attribute, this routine calculates the solution at the element boundaries
!! and stores it in the "boundarySolution" data structure.
!! 
!!  This subroutine depends on
!!   Module \ref NodalStorage_1D_Class : Function \ref CalculateSolutionAtBoundaries_NodalStorage_1D
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DGSEMSolution_1D) :: this <BR>
!! <B>TYPE</B>(NodalStorage)     :: dgStorage <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateBoundarySolution( dgStorage ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myDGS <td> DGSEMSolution_1D <td> 
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
 SUBROUTINE CalculateBoundarySolution_DGSEMSolution_1D( myDGS, dgStorage )

   IMPLICIT NONE
   CLASS(DGSEMSolution_1D), INTENT(inout) :: myDGS
   TYPE(NodalStorage), INTENT(in)         :: dgStorage
   ! LOCAL
   INTEGER :: iEq, N, nEq
   REAL(prec) :: sol(0:myDGS % N)
   
      nEq = myDGS % nEq
      N   = myDGS % N

      DO iEq = 1, nEq 
         sol = myDGS % solution(0:N,iEq)
         myDGS % boundarySolution(iEq,1:2) = dgStorage % CalculateSolutionAtBoundaries_1D( sol )
      ENDDO 
         
   
 END SUBROUTINE CalculateBoundarySolution_DGSEMSolution_1D

END MODULE DGSEMSolutionStorage_1D_Class
