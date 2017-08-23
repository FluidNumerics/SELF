! Surface_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
!> \file Surface_Class.f90
!! Contains the \ref Surface_Class module, and <BR>
!! defines the \ref Surface data-structure.

!> \defgroup Surface_Class Surface_Class 
!! This module defines the Surface data-structure and its associated routines.


MODULE Surface_Class
 
! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/interp/
USE Lagrange_Class

IMPLICIT NONE
!> \addtogroup Surface_Class 
!! @{

!> \struct Surface
!!  The Surface class provides attributes and type-bound procedures for defining and manipulating
!!  surfaces in multiple dimensions.
!!
!!  A surface is a geometric primitive that can be described two free parameters.
!!
!!  As an example, a surface in three-dimensions is represented as
!!  \f[
!!       \vec{x}(\xi^1,\xi^2) = x(\xi^1,\xi^2) \hat{x} + y(\xi^1,\xi^2) \hat{y} + z(\xi^1,\xi^2) \hat{z}
!!  \f]
!!  where \f$ (\xi^1,\xi^2) \f$ are parameters defined on [-1,1]x[-1,1]. 
!!  In the SELF, surfaces in 3-D are primarily used in the generation of mappings between physical 
!!  and computational space for hexahedral elements.
!!
!!  The Surface class permits surfaces that reside in higher dimensions.
!!
!!
!! <H2> Surface </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the interpolant used to describe the surface
!!       <tr> <th> nDim <td> INTEGER  <td> Dimension of the surface
!!       <tr> <th> interp <td> Lagrange <td> Lagrange interpolant that contains the discrete
!!                                              values of the free parameters.
!!       <tr> <th> x(0:interp % N,0:interp % N,1:nDim) <td> REAL(prec) <td> 
!!                  Discrete position vectors of the surface.
!!       <tr> <th> dxds(0:interp % N,0:interp % N,1:nDim) <td> REAL(prec) <td>
!!                 Derivative of the position wrt the first free parameter
!!       <tr> <th> dxdp(0:interp % N,0:interp % N,1:nDim) <td> REAL(prec) <td>
!!                 Derivative of the position wrt the second free parameter
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref surface_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_Surface
!!       <tr> <th> Trash <td> Trash_Surface
!!       <tr> <th> CalculateSlope <td> CalculateSlope_Surface
!!       <tr> <th> Reset <td> Reset_Surface
!!       <tr> <th> Evaluate <td> Evaluate_Surface
!!       <tr> <th> EvaluateSlope <td> EvaluateSlope_Surface
!!    </table>
!!

!>@}
   TYPE Surface 
      INTEGER                 :: N, nDim
      TYPE(Lagrange)          :: interp
      REAL(prec), ALLOCATABLE :: x(:,:,:)
      REAL(prec), ALLOCATABLE :: dxds(:,:,:), dxdp(:,:,:)
      
      CONTAINS

      PROCEDURE :: Build => Build_Surface
      PROCEDURE :: Trash => Trash_Surface

      PROCEDURE :: CalculateSlope => CalculateSlope_Surface
      PROCEDURE :: Reset          => Reset_Surface
      PROCEDURE :: Evaluate       => Evaluate_Surface
      PROCEDURE :: EvaluateSlope  => EvaluateSlope_Surface
      
      
   END TYPE Surface

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Surface_Class 
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_Surface
!! Initializes and assigns the attributes of the Surface class.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref Build_Lagrange <BR> 
!!   Module \ref Surface_Class : S/R \ref CalculateSlope_Surface <BR> 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Surface) :: this <BR>
!! <B>INTEGER</B>       :: N, nDim <BR>
!! <B>REAL</B>(prec)    :: x(0:N,0:N,1:nDim), nodes(0:N) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( x, nodes, N, nDim ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> mySurface <td> Surface <td> On output, the surface interpolant, position vectors,
!!                                              and derivatives are filled in.
!!   <tr> <td> in <th> x(0:N,0:N,1:nDim) <td> REAL(prec) <td> Surface position vectors
!!   <tr> <td> in <th> nodes(0:N)* <td> REAL(prec) <td> Discrete locations of the free parameter.
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the interpolant that describes
!!                                         the surface.
!!   <tr> <td> in <th> nDim <td> INTEGER <td> Number of spatial dimensions that the surface resides in. 
!!  </table>  
!!  * If the surface is being used for Mapped-Geometry element construction, the nodes must be between 
!!    [-1,1]. A tensor product of "nodes" with itself is used to construct the 2-D free parameter
!!    space
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_Surface( mySurface, x, nodes, N, nDim )

   IMPLICIT NONE
   CLASS( Surface ), INTENT(out) :: mySurface
   INTEGER, INTENT(in)           :: N, nDim
   REAL(prec), INTENT(in)        :: x(0:N,0:N,1:nDim), nodes(0:N)
  
      mySurface % N    = N
      mySurface % nDim = nDim
      ALLOCATE( mySurface % x(0:N,0:N,1:nDim), &
                mySurface % dxds(0:N,0:N,1:nDim), &
                mySurface % dxdp(0:N,0:N,1:nDim))

      CALL mySurface % interp % Build( N, N, nodes, nodes )
      mySurface % x = x
      CALL mySurface % CalculateSlope( )
 
 END SUBROUTINE Build_Surface
!
!> \addtogroup Surface_Class
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_Surface
!! Frees memory held by the attributes of the Surface class.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref Trash_Lagrange <BR>  
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Surface) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> mySurface <td> Surface <td> 
!!                         On <B>input</B>, a previously constructed Surface data-structure, <BR>
!!                         On <B>output</B>, the memory held by its attributes is freed. <BR>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_Surface( mySurface )

  IMPLICIT NONE
  CLASS( Surface ), INTENT(inout)     :: mySurface
  
      CALL mySurface % interp % Trash( )
      DEALLOCATE( mySurface % x, mySurface % dxds, mySurface % dxdp )
 
 END SUBROUTINE Trash_Surface
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Surface_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateSlope 
! 
!> \fn CalculateSlope_Surface  
!! Calculates the derivative of the surface position wrt to the free parameter and stores the result
!! within the data structure (attribute "dxds").
!! 
!! The interpolant's derivative matrix is used to quickly compute the surface slope at each of the 
!! points where the surface position is known.
!!
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref ApplyDerivativeMatrix_Lagrange <BR> 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Surface) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateSlope(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> mySurface <td> Surface <td>
!!                         On <B>input</B>, a Surface structure with the position ("x") attribute
!!                         filled in, <BR>
!!                         On <B>output</B>, the surface slope ("dxds") attribute is filled in <BR> 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE CalculateSlope_Surface( mySurface )

   IMPLICIT NONE
   CLASS( Surface ), INTENT(inout) :: mySurface
   ! LOCAL
   REAL(prec) :: dx(0:mySurface % N, 0:mySurface % N,1:2)
   INTEGER    :: nDim, i
   
      nDim = mySurface % nDim
      DO i = 1, nDim
         dx = mySurface % interp % ApplyDerivativeMatrix_2D( mySurface % x(:,:,i) )
         mySurface % dxds(:,:,i) = dx(:,:,1)
         mySurface % dxdp(:,:,i) = dx(:,:,2)
      ENDDO
      
 END SUBROUTINE CalculateSlope_Surface
!
!> \addtogroup Surface_Class 
!! @{ 
! ================================================================================================ !
! S/R Rset
! 
!> \fn Reset_Surface
!! Resets the Surface position data and re-calculates the surface slope data.
!!
!!  If you want to re-use an previously constructed surface, this routine overwrites the surface position
!!  and slope information. Note that the number of position data points must not change.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Surface_Class : S/R \ref CalculateSlope_Surface <BR> 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Surface) :: this <BR>
!! <B>REAL</B>(prec)    :: x(0:this % N,0:this % N,1:this % nDim) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Reset( x ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> inout <th> mySurface <td> Surface <td>
!!                         On <B>input</B>, a previously constructed Surface, <BR>
!!                         On <B>output</B>, new surface positions and surface slopes are filled in <BR> 
!!   <tr> <td> in <th> x(0:mySurface % N,0:mySurface % N, 1: mySurface % nDim) <td> REAL(prec) <td> 
!!                     Surface position vectors 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Reset_Surface( mySurface, x )

   IMPLICIT NONE
   CLASS( Surface ), INTENT(inout) :: mySurface
   REAL(prec), INTENT(in)          :: x(0:mySurface % N, 0:mySurface % N, 1:mySurface % nDim)

      mySurface % x = x
      CALL mySurface % CalculateSlope( )
 
 END SUBROUTINE Reset_Surface
!
!> \addtogroup Surface_Class 
!! @{ 
! ================================================================================================ !
! Function Evaluate 
! 
!> \fn Evaluate_Surface  
!! Estimates the surface position at a given value of the surface parameter. 
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref Interpolate_Lagrange <BR> 
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Surface) :: this <BR>
!! <B>REAL</B>(prec)    :: s(1:2) <BR>
!! <B>REAL</B>(prec)    :: x(1:this % nDim) <BR>
!!         .... <BR>
!!     x = this % Evaluate( s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> mySurface <td> Surface <td> A previously constructed Surface data structure
!!   <tr> <td> in <th> s(1:2) <td> REAL(prec) <td> Value of the surface parameters where the surface position
!!                                            is desired.
!!   <tr> <td> in <th> x(1:mySurface % nDim) <td> REAL(prec) <td> Position of the surface at "s" 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION Evaluate_Surface( mySurface, s ) RESULT( x )

   IMPLICIT NONE
   CLASS( Surface ) :: mySurface
   REAL(prec)       :: s(1:2)
   REAL(prec)       :: x(1:mySurface % nDim)
   ! LOCAL
   INTEGER    :: nDim, i
   
      nDim = mySurface % nDim
      DO i = 1, nDim
         x(i) = mySurface % interp % Interpolate_2D( mySurface % x(:,:,i), s )
      ENDDO
      
 END FUNCTION Evaluate_Surface
!
!> \addtogroup Surface_Class 
!! @{ 
! ================================================================================================ !
! Function EvaluateSlope 
! 
!> \fn EvaluateSlope_Surface  
!! Estimates the surface slope at a given value of the surface parameter. 
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref Differentiate_Lagrange <BR> 
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Surface) :: this <BR>
!! <B>REAL</B>(prec)    :: s(1:2) <BR>
!! <B>REAL</B>(prec)    :: dxds(1:this % nDim,1:2) <BR>
!!         .... <BR>
!!     dxds = this % EvaluateSlope( s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> mySurface <td> Surface <td> A previously constructed Surface data structure
!!   <tr> <td> in <th> s(1:2) <td> REAL(prec) <td> Value of the surface parameters where the surface slope
!!                                            is desired.
!!   <tr> <td> in <th> dxds(1:mySurface % nDim,1:2) <td> REAL(prec) <td> Slope of the surface at "s" 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION EvaluateSlope_Surface( mySurface, s ) RESULT( dxds )

   IMPLICIT NONE
   CLASS( Surface ) :: mySurface
   REAL(prec)       :: s(1:2)
   REAL(prec)       :: dxds(1:mySurface % nDim,1:2)
   ! LOCAL
   INTEGER    :: nDim, i
   
      nDim = mySurface % nDim
      DO i = 1, nDim
         dxds(i,1:2) = mySurface % interp % Differentiate_2D( mySurface % x(:,:,i), s )
      ENDDO

 END FUNCTION EvaluateSlope_Surface
!
END MODULE Surface_Class
