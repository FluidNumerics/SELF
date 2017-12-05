! Curve_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file Curve_Class.f90
!! Contains the \ref Curve_Class module, and <BR>
!! defines the \ref Curve data-structure.

!> \defgroup Curve_Class Curve_Class 
!! This module defines the Curve data-structure and its associated routines.

MODULE Curve_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/interp/
USE Lagrange_Class

IMPLICIT NONE
!> \addtogroup Curve_Class 
!! @{

!> \struct Curve
!!  The Curve class provides attributes and type-bound procedures for defining and manipulating
!!  curves in multiple dimensions.
!!
!!  As an example, a two-dimensional curve is represented as a parametric curve
!!  \f[
!!       \vec{x}(\xi) = x(\xi) \hat{x} + y(\xi) \hat{y}
!!  \f]
!!  where \f$ \xi \f$ is a parameter that varies between [-1,1]. In the SELF, 2-D curves are primarily
!!  used in the generation of mappings between physical and computational space for quadrilateral
!!  elements.
!!
!!  The Curve class permits curves that reside in higher dimensions
!!
!!
!! <H2> Curve </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the interpolant used to describe the curve
!!       <tr> <th> nDim <td> INTEGER  <td> Dimension of the curve
!!       <tr> <th> interp <td> Lagrange <td>  Lagrange interpolant that contains the discrete
!!                                              values of the free parameter.
!!       <tr> <th> x(0:interp % N,1:nDim) <td> REAL(prec) <td> Position vector of the curve
!!       <tr> <th> dxds(0:interp % N,1:nDim) <td> REAL(prec) <td> Derivative of the position wrt to the 
!!                                                         curve parameter
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Curve_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Build <td> Build_Curve
!!       <tr> <th> Trash <td> Trash_Curve
!!       <tr> <th> CalculateSlope <td> CalculateSlope_Curve
!!       <tr> <th> Reset <td> Reset_Curve
!!       <tr> <th> Evaluate <td> Evaluate_Curve
!!       <tr> <th> EvaluateSlope <td> EvaluateSlope_Curve
!!    </table>
!!

!>@}

   TYPE Curve 
      INTEGER                  :: N, nDim
      TYPE(Lagrange)           :: interp
      REAL(prec),  ALLOCATABLE :: x(:,:)
      REAL(prec),  ALLOCATABLE :: dxds(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_Curve
      PROCEDURE :: Trash => Trash_Curve

      PROCEDURE :: CalculateSlope => CalculateSlope_Curve
      PROCEDURE :: Reset          => Reset_Curve
      PROCEDURE :: Evaluate       => Evaluate_Curve
      PROCEDURE :: EvaluateSlope  => EvaluateSlope_Curve
      
      
   END TYPE Curve

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Curve_Class 
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_Curve
!! Initializes and assigns the attributes of the Curve class.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref Build_Lagrange <BR> 
!!   Module \ref Curve_Class : S/R \ref CalculateSlope_Curve <BR> 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Curve) :: this <BR>
!! <B>INTEGER</B>     :: N, nDim <BR>
!! <B>REAL</B>(prec)  :: x(0:N,1:nDim), nodes(0:N) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( x, nodes, N, nDim ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myCurve <td> Curve <td> On output, the curve interpolant, position vectors,
!!                                              and derivatives are filled in.
!!   <tr> <td> in <th> x(0:N,1:nDim) <td> REAL(prec) <td> Curve position vectors
!!   <tr> <td> in <th> nodes(0:N)* <td> REAL(prec) <td> Discrete locations of the free parameter.
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the interpolant that describes
!!                                         the curve.
!!   <tr> <td> in <th> nDim <td> INTEGER <td> Number of spatial dimensions that the curve resides in. 
!!  </table>  
!!  * If the curve is being used for Mapped-Geometry element construction, the nodes must be between 
!!    [-1,1].
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_Curve( myCurve, x, nodes, N, nDim )

   IMPLICIT NONE
   CLASS( Curve ), INTENT(out) :: myCurve
   INTEGER, INTENT(in)         :: N, nDim
   REAL(prec), INTENT(in)      :: x(0:N, 1:nDim), nodes(0:N)

      myCurve % N    = N
      myCurve % nDim = nDim
      ALLOCATE( myCurve % x(0:N,1:nDim), myCurve % dxds(0:N,1:nDim) )

      CALL myCurve % interp % Build( N, N, nodes, nodes )
      myCurve % x = x
      
      CALL myCurve % CalculateSlope( )
 
 END SUBROUTINE Build_Curve
!
!> \addtogroup Curve_Class
!! @{ 
! ================================================================================================ !
! S/R Trash
! 
!> \fn Trash_Curve
!! Frees memory held by the attributes of the Curve class.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref Trash_Lagrange <BR>  
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Curve) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myCurve <td> Curve <td> 
!!                         On <B>input</B>, a previously constructed Curve data-structure, <BR>
!!                         On <B>output</B>, the memory held by its attributes is freed. <BR>
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_Curve( myCurve )

  IMPLICIT NONE
  CLASS( Curve ), INTENT(inout) :: myCurve
  
      CALL myCurve % interp % Trash( )

      DEALLOCATE( myCurve % x, myCurve % dxds )

 END SUBROUTINE Trash_Curve
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup Curve_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateSlope 
! 
!> \fn CalculateSlope_Curve  
!! Calculates the derivative of the curve position wrt to the free parameter and stores the result
!! within the data structure (attribute "dxds").
!! 
!! The interpolant's derivative matrix is used to quickly compute the curve slope at each of the 
!! points where the curve position is known.
!!
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref ApplyDerivativeMatrix_Lagrange <BR> 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Curve) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateSlope(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myCurve <td> Curve <td>
!!                         On <B>input</B>, a Curve structure with the position ("x") attribute
!!                         filled in, <BR>
!!                         On <B>output</B>, the curve slope ("dxds") attribute is filled in <BR> 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE CalculateSlope_Curve( myCurve )

   IMPLICIT NONE
   CLASS( Curve ), INTENT(inout) :: myCurve
   ! Local
   INTEGER :: nDim, i
   
      nDim = myCurve % nDim
      DO i = 1, nDim
         myCurve % dxds(:,i) = myCurve % interp % ApplyDerivativeMatrix_1D( myCurve % x(:,i) )
      ENDDO
      
 END SUBROUTINE CalculateSlope_Curve
!
!> \addtogroup Curve_Class 
!! @{ 
! ================================================================================================ !
! S/R Rset
! 
!> \fn Reset_Curve
!! Resets the Curve position data and re-calculates the curve slope data.
!!
!!  If you want to re-use an previously constructed curve, this routine overwrites the curve position
!!  and slope information. Note that the number of position data points must not change.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Curve_Class : S/R \ref CalculateSlope_Curve <BR> 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Curve) :: this <BR>
!! <B>REAL</B>(prec)  :: x(0:this % N,1:this % nDim) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Reset( x ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> inout <th> myCurve <td> Curve <td>
!!                         On <B>input</B>, a previously constructed Curve, <BR>
!!                         On <B>output</B>, new curve positions and curve slopes are filled in <BR> 
!!   <tr> <td> in <th> x(0:myCurve % N,1:myCurve % nDim) <td> REAL(prec) <td> Curve position vectors 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Reset_Curve( myCurve, x )

   IMPLICIT NONE
   CLASS( Curve ), INTENT(inout) :: myCurve
   REAL(prec), INTENT(in)        :: x(0:myCurve % N, 1:myCurve % nDim)

      myCurve % x = x
      CALL myCurve % CalculateSlope( )
 
 END SUBROUTINE Reset_Curve
!
!> \addtogroup Curve_Class 
!! @{ 
! ================================================================================================ !
! Function Evaluate 
! 
!> \fn Evaluate_Curve  
!! Estimates the curve position at a given value of the curve parameter. 
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref Interpolate_Lagrange <BR> 
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Curve) :: this <BR>
!! <B>REAL</B>(prec)  :: s <BR>
!! <B>REAL</B>(prec)  :: x(1:this % nDim) <BR>
!!         .... <BR>
!!     x = this % Evaluate( s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myCurve <td> Curve <td> A previously constructed Curve data structure
!!   <tr> <td> in <th> s <td> REAL(prec) <td> Value of the curve parameter where the curve position
!!                                            is desired.
!!   <tr> <td> in <th> x(1:myCurve % nDim) <td> REAL(prec) <td> Position of the curve at "s" 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION Evaluate_Curve( myCurve, s ) RESULT( x )

   IMPLICIT NONE
   CLASS( Curve ) :: myCurve
   REAL(prec)     :: s 
   REAL(prec)     :: x(1:myCurve % nDim)
   ! Local
   INTEGER :: nDim, i
   
      nDim = myCurve % nDim
      DO i = 1, nDim
         x(i) = myCurve % interp % Interpolate_1D( myCurve % x(:,i), s )
      ENDDO

 END FUNCTION Evaluate_Curve
!
!> \addtogroup Curve_Class 
!! @{ 
! ================================================================================================ !
! Function EvaluateSlope 
! 
!> \fn EvaluateSlope_Curve  
!! Estimates the curve slope at a given value of the curve parameter. 
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange : S/R \ref Differentiate_1D_Lagrange <BR> 
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(Curve) :: this <BR>
!! <B>REAL</B>(prec)  :: s <BR>
!! <B>REAL</B>(prec)  :: dxds(1:this % nDim) <BR>
!!         .... <BR>
!!     dxds = this % EvaluateSlope( s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myCurve <td> Curve <td> A previously constructed Curve data structure
!!   <tr> <td> in <th> s <td> REAL(prec) <td> Value of the curve parameter where the curve slope
!!                                            is desired.
!!   <tr> <td> in <th> dxds(1:myCurve % nDim) <td> REAL(prec) <td> Slope of the curve at "s" 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION EvaluateSlope_Curve( myCurve, s ) RESULT( dxds )

   IMPLICIT NONE
   CLASS( Curve ) :: myCurve
   REAL(prec)     :: s 
   REAL(prec)     :: dxds(1:myCurve % nDim)
   ! Local
   INTEGER :: nDim, i
   
      nDim = myCurve % nDim
      DO i = 1, nDim
         dxds(i) = myCurve % interp % Differentiate_1D( myCurve % x(:,i), s )
      ENDDO

 END FUNCTION EvaluateSlope_Curve
!

END MODULE Curve_Class
