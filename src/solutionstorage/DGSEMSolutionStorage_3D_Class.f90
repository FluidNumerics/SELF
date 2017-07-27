! DGSEMSolutionStorage_3D_Class.f90
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
! DGSEMSolutionStorage_3D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file DGSEMSolutionStorage_3D_Class.f90
!! Contains the \ref DGSEMSolutionStorage_3D_Class module, and <BR>
!! defines the \ref DGSEMSolution_3D data-structure.

!> \defgroup DGSEMSolutionStorage_3D_Class DGSEMSolutionStorage_3D_Class 
!! This module defines the DGSEMSolution_3D data-structure and its associated routines.
  
MODULE DGSEMSolutionStorage_3D_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/nodal/
USE NodalStorage_Class

IMPLICIT NONE

!> \addtogroup DGSEMSolutionStorage_3D_Class 
!! @{

!> \struct DGSEMSolution_3D
!!  The DGSEMSolution_3D class provides attributes for storing a solution and its flux on a
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
!!  calling a routine from the \ref NodalStorage_3D_Class module.
!!
!! <H2> DGSEMSolution_3D </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the spectral element method
!!       <tr> <th> nEq <td> INTEGER <td> Number of (prognostic) solution variables.
!!       <tr> <th> solution(0:N,0:N,0:N,1:nEq) <td> REAL(prec) <td> An array containing the solution variables
!!       <tr> <th> tendency(0:N,0:N,0:N,1:nEq) <td> REAL(prec) <td> An array containing the tendency of the 
!!       <tr> <th> tendency(0:N,0:N,0:N,1:nEq) <td> REAL(prec) <td> An array containing the tendency of the 
!!                                                          solution variables.
!!       <tr> <th> boundarySolution(0:N,0:N,1:6,1:nEq)* <td> REAL(prec) <td>
!!                  An array containing the solution variables at the element boundary
!!       <tr> <th> boundaryFlux(0:N,0:N,1:6,1:nEq)* <td> REAL(prec) <td>
!!                  An array containing the flux at the element boundary
!!    </table>
!!
!!  *For the "boundary" arrays, the sides for a quadrilateral element are numbered as SOUTH=1, 
!!  EAST=2, NORTH=3, WEST=4, BOTTOM=5, TOP=6.
!!
!! <H3> Procedures </H3>
!!    See \ref DGSEMSolutionStorage_3D_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_DGSEMSolution_3D
!!       <tr> <th> Trash <td> Trash_DGSEMSolution_3D
!!       <tr> <th> CalculateBoundarySolution <td> CalculateBoundarySolution_DGSEMSolution_3D
!!    </table>
!!

!>@}

    TYPE DGSEMSolution_3D
      INTEGER                 :: nEq, N
      REAL(prec), allocatable :: Solution(:,:,:,:,:)
      REAL(prec), allocatable :: tendency(:,:,:,:,:)
      REAL(prec), allocatable :: boundarySolution(:,:,:,:,:) ! Indexed over the boundaries. 
      REAL(prec), allocatable :: boundaryFlux(:,:,:,:,:)     ! Indexed over the boundaries

      CONTAINS

      ! CONSTRUCTORS/DESTRUCTORS
      PROCEDURE :: Build => Build_DGSEMSolution_3D
      PROCEDURE :: Trash => Trash_DGSEMSolution_3D

!      PROCEDURE :: CalculateBoundarySolution => CalculateBoundarySolution_DGSEMSolution_3D
    END TYPE DGSEMSolution_3D
    
    
    ! Light Storage : No tendency or boundary flux
    TYPE LightDGSEMSolution_3D
      INTEGER                 :: nEq, N
      REAL(prec), allocatable :: solution(:,:,:,:,:)
      REAL(prec), allocatable :: boundarySolution(:,:,:,:,:) ! Indexed over the boundaries. 

      CONTAINS

      ! CONSTRUCTORS/DESTRUCTORS
      PROCEDURE :: Build => Build_LightDGSEMSolution_3D
      PROCEDURE :: Trash => Trash_LightDGSEMSolution_3D

     ! PROCEDURE :: CalculateBoundarySolution => CalculateBoundarySolution_LightDGSEMSolution_3D
    END TYPE LightDGSEMSolution_3D

 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup DGSEMSolutionStorage_3D_Class
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_DGSEMSolution_3D_Class  
!! Allocates space for the DGSEMSolution_3D attributes and initializes arrays to 0.0_prec.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DGSEMSolution_3D) :: this <BR>
!! <B>INTEGER</B>                :: N, nEq <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( N, nEq ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myDGS <td> DGSEMSolution_3D <td> On output, memory has been allocated for
!!                                                       each attribute and each array is initialized
!!                                                       with a value of 0.0_prec. 
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the DG-method you plan on using
!!   <tr> <td> in <th> nEq <td> INTEGER <td> Number of prognostic variables; number of equations for
!!                                           the system being solved
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_DGSEMSolution_3D( myDGS, N, nEq, nElems )

   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                    :: N, nEq, nElems
      
      myDGS % N   = N
      myDGS % nEq = nEq

      ALLOCATE( myDGS % solution(0:N,0:N,0:N,1:nEq,1:nElems) )
      ALLOCATE( myDGS % tendency(0:N,0:N,0:N,1:nEq,1:nElems) )
      ALLOCATE( myDGS % boundarySolution(0:N,0:N,1:nEq,1:nHexFaces,1:nElems) ) 
      ALLOCATE( myDGS % boundaryFlux(0:N,0:N,1:nEq,1:nHexFaces,1:nElems) ) 
      
      myDGS % solution = 0.0_prec
      myDGS % tendency = 0.0_prec
      myDGS % boundarySolution = 0.0_prec
      myDGS % boundaryFlux = 0.0_prec


 END SUBROUTINE Build_DGSEMSolution_3D
!
!> \addtogroup DGSEMSolutionStorage_3D_Class
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_DGSEMSolution_3D
!! Frees memory held by the attributes of the DGSEMSolution_3D data structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DGSEMSolution_3D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myDGS <td> DGSEMSolution_3D <td> On output, memory held by the attributes
!!                                                          of the data structure has been deallocated.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_DGSEMSolution_3D( myDGS )

   IMPLICIT NONE
   CLASS(DGSEMSolution_3D), INTENT(inout) :: myDGS

      DEALLOCATE( myDGS % Solution )
      DEALLOCATE( myDGS % tendency )
      DEALLOCATE( myDGS % boundarySolution ) 
      DEALLOCATE( myDGS % boundaryFlux ) 

 END SUBROUTINE Trash_DGSEMSolution_3D
!
!
!> \addtogroup LightDGSEMSolutionStorage_3D_Class
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_LightDGSEMSolution_3D_Class  
!! Allocates space for the LightDGSEMSolution_3D attributes and initializes arrays to 0.0_prec.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(LightDGSEMSolution_3D) :: this <BR>
!! <B>INTEGER</B>                :: N, nEq <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( N, nEq ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myDGS <td> LightDGSEMSolution_3D <td> On output, memory has been allocated for
!!                                                       each attribute and each array is initialized
!!                                                       with a value of 0.0_prec. 
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the DG-method you plan on using
!!   <tr> <td> in <th> nEq <td> INTEGER <td> Number of prognostic variables; number of equations for
!!                                           the system being solved
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_LightDGSEMSolution_3D( myDGS, N, nEq, nElems )

   IMPLICIT NONE
   CLASS(LightDGSEMSolution_3D), INTENT(inout) :: myDGS
   INTEGER, INTENT(in)                         :: N, nEq, nElems
      
      myDGS % N   = N
      myDGS % nEq = nEq

      ALLOCATE( myDGS % solution(0:N,0:N,0:N,1:nEq,1:nElems) )
      ALLOCATE( myDGS % boundarySolution(0:N,0:N,1:nEq,1:nHexFaces,nElems) ) 
      
      myDGS % solution = 0.0_prec
      myDGS % boundarySolution = 0.0_prec

 END SUBROUTINE Build_LightDGSEMSolution_3D
!
!> \addtogroup LightDGSEMSolutionStorage_3D_Class
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_LightDGSEMSolution_3D
!! Frees memory held by the attributes of the LightDGSEMSolution_3D data structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(LightDGSEMSolution_3D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myDGS <td> LightDGSEMSolution_3D <td> On output, memory held by the attributes
!!                                                          of the data structure has been deallocated.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_LightDGSEMSolution_3D( myDGS )

   IMPLICIT NONE
   CLASS(LightDGSEMSolution_3D), INTENT(inout) :: myDGS

      DEALLOCATE( myDGS % Solution )
      DEALLOCATE( myDGS % boundarySolution ) 

 END SUBROUTINE Trash_LightDGSEMSolution_3D

END MODULE DGSEMSolutionStorage_3D_Class
