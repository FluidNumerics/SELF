! NodalDGSolution_3D_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file NodalDGSolution_3D_Class.f90
!! Contains the \ref NodalDGSolution_3D_Class module, and <BR>
!! defines the \ref NodalDGSolution_3D data-structure.

!> \defgroup NodalDGSolution_3D_Class NodalDGSolution_3D_Class 
!! This module defines the NodalDGSolution_3D data-structure and its associated routines.
  
MODULE NodalDGSolution_3D_Class

USE ModelPrecision

IMPLICIT NONE

!> \addtogroup NodalDGSolution_3D_Class 
!! @{

!> \struct NodalDGSolution_3D
!!  The NodalDGSolution_3D class provides attributes for storing a solution and its flux on a
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
!! <H2> NodalDGSolution_3D </H2>
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
!!    See \ref NodalDGSolution_3D_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_NodalDGSolution_3D
!!       <tr> <th> Trash <td> Trash_NodalDGSolution_3D
!!       <tr> <th> CalculateBoundarySolution <td> CalculateBoundarySolution_NodalDGSolution_3D
!!    </table>
!!

!>@}

  TYPE NodalDGSolution_3D
    INTEGER                 :: N, nEquations
    REAL(prec), ALLOCATABLE :: solution(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: tendency(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: boundarySolution(:,:,:,:,:) ! Indexed over the boundaries. 
    REAL(prec), ALLOCATABLE :: boundaryFlux(:,:,:,:,:)     ! Indexed over the boundaries

    CONTAINS

    PROCEDURE :: Build => Build_NodalDGSolution_3D
    PROCEDURE :: Trash => Trash_NodalDGSolution_3D

  END TYPE NodalDGSolution_3D
    
    
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup NodalDGSolution_3D_Class
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_NodalDGSolution_3D_Class  
!! Allocates space for the NodalDGSolution_3D attributes and initializes arrays to 0.0_prec.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalDGSolution_3D) :: this <BR>
!! <B>INTEGER</B>                :: N, nEq <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( N, nEq ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myDGS <td> NodalDGSolution_3D <td> On output, memory has been allocated for
!!                                                       each attribute and each array is initialized
!!                                                       with a value of 0.0_prec. 
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the DG-method you plan on using
!!   <tr> <td> in <th> nEq <td> INTEGER <td> Number of prognostic variables; number of equations for
!!                                           the system being solved
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_NodalDGSolution_3D( myDGS, N, nEq, nElems )

   IMPLICIT NONE
   CLASS(NodalDGSolution_3D), INTENT(inout) :: myDGS
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


 END SUBROUTINE Build_NodalDGSolution_3D
!
!> \addtogroup NodalDGSolution_3D_Class
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_NodalDGSolution_3D
!! Frees memory held by the attributes of the NodalDGSolution_3D data structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(NodalDGSolution_3D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myDGS <td> NodalDGSolution_3D <td> On output, memory held by the attributes
!!                                                          of the data structure has been deallocated.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_NodalDGSolution_3D( myDGS )

   IMPLICIT NONE
   CLASS(NodalDGSolution_3D), INTENT(inout) :: myDGS

      DEALLOCATE( myDGS % Solution )
      DEALLOCATE( myDGS % tendency )
      DEALLOCATE( myDGS % boundarySolution ) 
      DEALLOCATE( myDGS % boundaryFlux ) 

 END SUBROUTINE Trash_NodalDGSolution_3D


END MODULE NodalDGSolution_3D_Class
