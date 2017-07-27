! MappedGeometry_2D_Class.f90
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
! MappedGeometry_2D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file MappedGeometry_2D_Class.f90
!! Contains the \ref MappedGeometry_2D_Class module, and <BR>
!! defines the \ref MappedGeometry_2D data-structure.

!> \defgroup MappedGeometry_2D_Class MappedGeometry_2D_Class 
!! This module defines the MappedGeometry_2D data-structure and its associated routines.

MODULE MappedGeometry_2D_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
! src/interp/
USE Lagrange_Class
! src/geom/
USE Curve_Class

IMPLICIT NONE

!> \addtogroup MappedGeometry_2D_Class 
!! @{

!> \struct MappedGeometry_2D
!!  The MappedGeometry_2D class provides attributes and type-bound procedures for handling 
!!  curvilinear mappings between physical space and a reference computational space.
!!
!!  The MappedGeometry_2D class enables the implementation of spectral element methods in complicated
!!  geometry. Given the four bounding curves of an element, the internal geometry (physical node 
!!  positions, covariant basis vectors, and Jacobian) is constructed using transfinite interpolation
!!  with linear blending.
!!
!! <H2> MappedGeometry_2D </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the spectral element method
!!       <tr> <th> nHat(1:2,0:N,1:4) <td> REAL(prec) <td> Outward pointing element boundary normal vectors
!!                                                     along each element edge.
!!       <tr> <th> xBound(0:N,1:4) <td> REAL(prec) <td> x-position of each of the boundary curves
!!       <tr> <th> yBound(0:N,1:4) <td> REAL(prec) <td> y-position of each of the boundary curves
!!       <tr> <th> x(0:N,0:N) <td> REAL(prec) <td> x-position of the interior nodes
!!       <tr> <th> y(0:N,0:N) <td> REAL(prec) <td> y-position of the interior nodes
!!       <tr> <th> J(0:N,0:N) <td> REAL(prec) <td> Jacobian of the mapping at the interior nodes
!!       <tr> <th> dxds(0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial x}{\partial \xi^1} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dxdp(0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial x}{\partial \xi^2} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dyds(0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial y}{\partial \xi^1} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dydp(0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial y}{\partial \xi^2} \f$ at
!!                                                    each of the interior nodes
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref MappedGeometry_2D_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_MappedGeometry_2D
!!       <tr> <th> Build <td> Build_MappedGeometry_2D
!!       <tr> <th> Trash <td> Trash_MappedGeometry_2D
!!       <tr> <th> GenerateMesh <td> GenerateMesh_MappedGeometry_2D
!!       <tr> <th> GenerateMetrics <td> GenerateMetrics_MappedGeometry_2D
!!       <tr> <th> GenerateBoundaryMetrics <td> GenerateBoundaryMetrics_MappedGeometry_2D
!!       <tr> <th> ScaleGeometry <td> ScaleGeometry_MappedGeometry_2D
!!       <tr> <th> CalculateLocation <td> CalculateLocation_MappedGeometry_2D
!!       <tr> <th> CalculateMetrics <td> CalculateMetrics_MappedGeometry_2D
!!       <tr> <th> CalculateComputationalCoordinates <td> CalculateComputationalCoordinates_MappedGeometry_2D
!!       <tr> <th> WriteTecplot <td> WriteTecplot_MappedGeometry_2D
!!    </table>
!!

!>@}
   TYPE MappedGeometry_2D
      INTEGER                    :: N
      REAL(prec), ALLOCATABLE    :: nHat(:,:,:) 
      REAL(prec), ALLOCATABLE    :: xBound(:,:)
      REAL(prec), ALLOCATABLE    :: yBound(:,:) 
      REAL(prec), ALLOCATABLE    :: x(:,:), y(:,:)
      REAL(prec), ALLOCATABLE    :: J(:,:)    
      REAL(prec), ALLOCATABLE    :: dxds(:,:), dxdp(:,:)
      REAL(prec), ALLOCATABLE    :: dyds(:,:), dydp(:,:)

      CONTAINS

      PROCEDURE :: Initialize => Initialize_MappedGeometry_2D
      PROCEDURE :: Build      => Build_MappedGeometry_2D
      PROCEDURE :: Trash      => Trash_MappedGeometry_2D

      PROCEDURE :: GenerateMesh => GenerateMesh_MappedGeometry_2D
      PROCEDURE :: GenerateMetrics => GenerateMetrics_MappedGeometry_2D
      PROCEDURE :: GenerateBoundaryMetrics => GenerateBoundaryMetrics_MappedGeometry_2D
      PROCEDURE :: ScaleGeometry => ScaleGeometry_MappedGeometry_2D
      PROCEDURE :: CalculateLocation => CalculateLocation_MappedGeometry_2D
      PROCEDURE :: CalculateMetrics => CalculateMetrics_MappedGeometry_2D
      PROCEDURE :: CalculateComputationalCoordinates => CalculateComputationalCoordinates_MappedGeometry_2D
      
      PROCEDURE :: WriteTecplot => WriteTecplot_MappedGeometry_2D
   END TYPE MappedGeometry_2D


 PRIVATE :: TransfiniteInterpolation, Unidirectional
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup MappedGeometry_2D_Class
!! @{ 
! ================================================================================================ !
! S/R Initialize
! 
!> \fn Initialize_MappedGeometry_2D  
!! Allocates memory for each of the attributes of the MappedGeometry_2D Class and initializes all
!! arrays to zero.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>INTEGER</B>                 :: N <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize( N ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myGeom <td> MappedGeometry_2D <td> On output, an initialized MappedGeometry_2D
!!                                                         data structure
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the spectral element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_MappedGeometry_2D( myGeom, N )

  IMPLICIT NONE
  CLASS(MappedGeometry_2D), INTENT(out) :: myGeom
  INTEGER, INTENT(in)                   :: N

      myGeom % N = N
       
      ! Allocate space
      ALLOCATE( myGeom % dxds(0:N,0:N), myGeom % dxdp(0:N,0:N) )
      ALLOCATE( myGeom % dyds(0:N,0:N), myGeom % dydp(0:N,0:N) )
      ALLOCATE( myGeom % J(0:N,0:N) )
      ALLOCATE( myGeom % x(0:N,0:N), myGeom % y(0:N,0:N) )
      ALLOCATE( myGeom % xBound(0:N,1:nQuadEdges) )
      ALLOCATE( myGeom % yBound(0:N,1:nQuadEdges) )
      ALLOCATE( myGeom % nHat(1:2,0:N,1:nQuadEdges) )
      
      myGeom % dxds   = ZERO
      myGeom % dxdp   = ZERO
      myGeom % dyds   = ZERO
      myGeom % dydp   = ZERO
      myGeom % J      = ZERO
      myGeom % x      = ZERO
      myGeom % y      = ZERO
      myGeom % xBound = ZERO
      myGeom % yBound = ZERO
      myGeom % nHat   = ZERO

 END SUBROUTINE Initialize_MappedGeometry_2D
!
!> \addtogroup MappedGeometry_2D_Class
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_MappedGeometry_2D  
!! Allocates memory for each of the attributes of the MappedGeometry_2D Class and generates the
!! physical mesh and metric terms.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref MappedGeometry_2D_Class : S/R \ref Initialize_MappedGeometry_2D <BR>
!!   Module \ref MappedGeometry_2D_Class : S/R \ref GenerateMesh_MappedGeometry_2D <BR>
!!   Module \ref MappedGeometry_2D_Class : S/R \ref GenerateMetrics_MappedGeometry_2D <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>TYPE</B>(Curve)             :: boundaries(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( interp, boundaries ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myGeom <td> MappedGeometry_2D <td> On output, an initialized MappedGeometry_2D
!!                                                         data structure with the mesh and metric
!!                                                         terms filled in
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolant that defines the 
!!                                                  reference computational mesh. 
!!   <tr> <td> in <th> myCurves(1:4) <td> Curve <td> 2-D Curve that defines the each of the 
!!                                                      boundaries of the quadrilateral element. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_MappedGeometry_2D( myGeom, interp, myCurves )

   IMPLICIT NONE
   CLASS(MappedGeometry_2D), INTENT(out) :: myGeom
   TYPE(Lagrange), INTENT(in)            :: interp
   TYPE(Curve), INTENT(in)               :: myCurves(1:nQuadEdges)
   !LOCAL
   INTEGER :: N

      N = interp % N
      CALL myGeom % Initialize( N )
      ! Generate the mesh locations, and the mapping metrics
      CALL myGeom % GenerateMesh( interp, myCurves )
      CALL myGeom % GenerateMetrics( interp )
 
 END SUBROUTINE Build_MappedGeometry_2D
!
!> \addtogroup MappedGeometry_2D_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_MappedGeometry_2D_Class  
!! Frees memory held by the attributes of the MappedGeometry_2D data-structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_2D <td>
!!                         On <B>input</B>, a previously constructed MappedGeometry_2D data structure <BR>
!!                         On <B>output</B>, memory held by attributes is freed
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_MappedGeometry_2D( myGeom )

   IMPLICIT NONE
   CLASS(MappedGeometry_2D), INTENT(inout)  :: myGeom

      DEALLOCATE( myGeom % dxds, myGeom % dxdp )
      DEALLOCATE( myGeom % dyds, myGeom % dydp)
      DEALLOCATE( myGeom % J, myGeom % x, myGeom % y )
      DEALLOCATE( myGeom % xBound, myGeom % yBound )
      DEALLOCATE( myGeom % nHat )
 
 END SUBROUTINE Trash_MappedGeometry_2D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup MappedGeometry_2D_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMesh 
! 
!> \fn GenerateMesh_MappedGeometry_2D  
!! Generates the physical interior and boundary positions at each of the computational mesh points.
!! 
!! Given the four boundary curves and an interpolant that stores the computational mesh points,
!! the physical mesh positions are generated using transfinite interpolation with linear blending.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref MappedGeometry_2D_Class : PRIVATE Function TransfiniteInterpolation <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>TYPE</B>(Curve)             :: boundaries(1:4) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateMesh( interp, boundaries ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_2D <td> 
!!                         On <B>input</B>, an initialized MappedGeometry_2D data structure, <BR>
!!                         On <B>output</B>, the mesh physical positions (interior and boundary)
!!                         are filled in.
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolant that defines the 
!!                                                  reference computational mesh. 
!!   <tr> <td> in <th> myCurves(1:4) <td> Curve <td> 2-D Curve that defines the each of the 
!!                                                      boundaries of the quadrilateral element. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GenerateMesh_MappedGeometry_2D( myGeom, interp, theCurves )

   IMPLICIT NONE
   CLASS( MappedGeometry_2D ), INTENT(inout) :: myGeom
   TYPE( Lagrange ), INTENT(in)              :: interp
   TYPE( Curve ), INTENT(in)                 :: theCurves(1:4)
   ! Local
   INTEGER    :: iS, iP, N
   REAL(prec) :: s(0:interp % N), p, x(1:2)
   
      N = interp % N
      s = interp % s
      
      DO iS = 0, N
         DO iP = 0,N
            x = TransfiniteInterpolation( theCurves, s(iS), s(iP) )
            myGeom % x(iS,iP) = x(1)
            myGeom % y(iS,iP) = x(2)
         ENDDO 
      ENDDO 
      
      ! Do the boundary locations
      DO iS = 0, N
      
         p = -ONE  ! south boundary
         x = TransfiniteInterpolation( theCurves, s(iS), p )
         myGeom % xBound(iS,south) = x(1)
         myGeom % yBound(iS,south) = x(2)
         ! west boundary
         x = TransfiniteInterpolation( theCurves, p, s(iS) )
         myGeom % xBound(iS,west) = x(1)
         myGeom % yBound(iS,west) = x(2)
         
         p = ONE  ! north boundary
         x = TransfiniteInterpolation( theCurves, s(iS), p )
         myGeom % xBound(iS,north) = x(1)
         myGeom % yBound(iS,north) = x(2)
         ! east boundary
         x = TransfiniteInterpolation( theCurves, p, s(iS) )
         myGeom % xBound(iS,east) = x(1)
         myGeom % yBound(iS,east) = x(2)
      
      ENDDO

 END SUBROUTINE GenerateMesh_MappedGeometry_2D
!
!> \addtogroup MappedGeometry_2D_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMetrics 
! 
!> \fn GenerateMetrics_MappedGeometry_2D  
!! Generates and stores the covariant basis vector components, Jacobian of the mapping, and the 
!! outward pointing boundary normal vectors.
!! 
!! Once the mesh has been generated, we effectively have
!! \f[ 
!!     \vec{x} = \vec{x}(\vec{\xi})
!! \f]
!! represented as a 2-D Lagrange interpolant. Differentiation of the interpolant allows us to 
!! estimate the covariant basis vectors. The determinant of the covariant metric tensor gives
!! the Jacobian of the mapping at each location in the computational mesh
!! \f[
!!     J = | \frac{\partial \vec{x}}{\partial \vec{\xi}} |
!! \f]
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref ApplyDerivativeMatrix_Lagrange <BR>
!!   Module \ref MappedGeometry_2D_Class : S/R \ref GenerateBoundaryMetrics_MappedGeometry_2D <BR>
!! 
!!  ** To produce meaningful output, GenerateMesh_MappedGeometry_2D must be called prior to calling
!!     this routine.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateMetrics( interp ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_2D <td> 
!!                         On <B>input</B>, an initialized MappedGeometry_2D data structure, <BR>
!!                         On <B>output</B>, the metric terms are filled in, including the
!!                         outward pointing normal vectors on the element boundaries
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolant that defines the 
!!                                                  reference computational mesh 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GenerateMetrics_MappedGeometry_2D( myGeom, interp )

   IMPLICIT NONE
   CLASS( MappedGeometry_2D ), INTENT(inout) :: myGeom
   TYPE( Lagrange ), INTENT(in)              :: interp
   ! Local
   INTEGER    :: iS, iP, N
   REAL(prec) :: cv(1:2,1:2)
   REAL(prec) :: covT(0:interp % N, 0:interp % N, 1:2,1:2)
   
      N = interp % N
      
      covT(0:N,0:N,1,1:2) = interp % ApplyDerivativeMatrix_2D( myGeom % x )
      covT(0:N,0:N,2,1:2) = interp % ApplyDerivativeMatrix_2D( myGeom % y )
      
      DO iS = 0, N
         DO iP = 0, N

            myGeom % dxds(iS,iP) = covT(iS,iP,1,1)
            myGeom % dxdp(iS,iP) = covT(iS,iP,1,2)
            myGeom % dyds(iS,iP) = covT(iS,iP,2,1)
            myGeom % dydp(iS,iP) = covT(iS,iP,2,2)
            
            cv = covT(iS,iP,1:2,1:2)
            myGeom % J(iS,iP) = Determinant( cv, 2 )

         ENDDO
      ENDDO 
      
      CALL myGeom % GenerateBoundaryMetrics( interp )

 END SUBROUTINE GenerateMetrics_MappedGeometry_2D
!
!> \addtogroup MappedGeometry_2D_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMetrics 
! 
!> \fn GenerateMetrics_MappedGeometry_2D  
!! Generates and stores the outward pointing boundary normal vectors.
!!
!!  The outward pointing boundary normal vectors are equivalent to the contravariant basis vectors
!!  evaluated at the element boundaries. These are computed here by differentiating the Lagrange
!!  interpolant of the mesh positions and the computational mesh boundaries.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref Differentiate_2D_Lagrange <BR>
!! 
!!  ** To produce meaningful output, GenerateMesh_MappedGeometry_2D must be called prior to calling
!!     this routine.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateBoundaryMetrics( interp ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_2D <td> 
!!                         On <B>input</B>, an initialized MappedGeometry_2D data structure, <BR>
!!                         On <B>output</B>, the outward pointing normal vectors on the element
!!                         boundaries are filled in
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolant that defines the 
!!                                                  reference computational mesh 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
  SUBROUTINE GenerateBoundaryMetrics_MappedGeometry_2D( myGeom, interp  )

   IMPLICIT NONE
   CLASS( MappedGeometry_2D ), INTENT(inout) :: myGeom
   TYPE( Lagrange ), INTENT(in)              :: interp
   ! Local
   INTEGER    :: iS, N
   REAL(prec) :: s(0:interp % N), p, cv(1:2,1:2)
   REAL(prec) :: J, signJ
   REAL(prec) :: node(1:2)
   
      N = interp % N
      s = interp % s
           ! Do the boundary locations
      DO iS = 0, N
      
         p = -ONE  ! south boundary
         node = (/ s(iS), p /)
         cv(1,1:2) = interp % Differentiate_2D( myGeom % x, node )
         cv(2,1:2) = interp % Differentiate_2D( myGeom % y, node )
         J = Determinant( cv, 2 )
         signJ = abs(J)/J                    
         ! Setting southern boundary normal         -dyds    ,    dxds
         myGeom % nHat(1:2,iS,south)= -signJ*(/ -cv(2,1), cv(1,1) /)
         
         ! west boundary
         node = (/ p, s(iS) /)
         cv(1,1:2) = interp % Differentiate_2D( myGeom % x, node )
         cv(2,1:2) = interp % Differentiate_2D( myGeom % y, node )
         J = Determinant( cv, 2 )
         signJ = abs(J)/J                    
         ! Setting western boundary normal         dydp    ,    -dxdp
         myGeom % nHat(1:2,iS,west) = -signJ*(/ cv(2,2), -cv(1,2) /)
          
         p = ONE  ! north boundary
         node = (/ s(iS), p /)
         cv(1,1:2) = interp % Differentiate_2D( myGeom % x, node )
         cv(2,1:2) = interp % Differentiate_2D( myGeom % y, node )
         J = Determinant( cv, 2 )
         signJ = abs(J)/J                    
         ! Setting southern boundary normal         -dyds    ,    dxds
         myGeom % nHat(1:2,iS,north) = signJ*(/ -cv(2,1), cv(1,1) /)
         
         ! east boundary
         node = (/ p, s(iS) /)
         cv(1,1:2) = interp % Differentiate_2D( myGeom % x, node )
         cv(2,1:2) = interp % Differentiate_2D( myGeom % y, node )
         J = Determinant( cv, 2 )
         signJ = abs(J)/J                    
         ! Setting eastern boundary normal         dydp    ,    -dxdp
         myGeom % nHat(1:2,iS,east) = signJ*(/ cv(2,2), -cv(1,2) /)
         
      ENDDO

 END SUBROUTINE GenerateBoundaryMetrics_MappedGeometry_2D
!
!> \addtogroup MappedGeometry_2D_Class
!! @{ 
! ================================================================================================ !
! Function CalculateLocation 
! 
!> \fn CalculateLocation_MappedGeometry_2D  
!! Given a computational coordinate, the physical coordinate is calculated using Lagrange 
!! interpolation.
!! 
!!  This function depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref Interpolate_Lagrange <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: s(1:2), x(1:2) <BR>
!!         .... <BR>
!!     x = this % CalculateLocation( interp, s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myGeom <td> MappedGeometry_2D <td> An intialized and constructed 
!!                                                        MappedGeometry_2D data structure
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolant data structure that
!!                                                  contains the computational mesh.
!!   <tr> <td> in <th> s(1:2) <td> REAL(prec)  <td> Computational location where the physical 
!!                                                  position is desired.
!!   <tr> <td> out <th> x(1:2) <td> REAL(prec)  <td> Physical location estimated by interpolation
!!                                                  onto the given computational position
!!
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION CalculateLocation_MappedGeometry_2D( myGeom, interp, s ) RESULT( x )

   IMPLICIT NONE
   CLASS( MappedGeometry_2D ) :: myGeom
   TYPE( Lagrange )           :: interp
   REAL(prec)                 :: s(1:2)
   REAL(prec)                 :: x(1:2)
  
      x(1) = interp % Interpolate_2D( myGeom % x, s )
      x(2) = interp % Interpolate_2D( myGeom % y, s  )
  
 END FUNCTION CalculateLocation_MappedGeometry_2D
!
!> \addtogroup MappedGeometry_2D_Class
!! @{ 
! ================================================================================================ !
! Function CalculateMetrics 
! 
!> \fn CalculateMetrics_MappedGeometry_2D  
!! Given a computational coordinate, the covariant metric tensor is estimated by differentiation
!! of the Lagrange interpolant of the mesh positions.
!!
!!  The output of this function is a 2x2 array whose entries are
!!  \f[
!!      covT(1,1) = \frac{\partial x}{\partial \xi^1}
!!  \f]
!!  \f[
!!      covT(1,2) = \frac{\partial x}{\partial \xi^2}
!!  \f]
!!  \f[
!!      covT(2,1) = \frac{\partial y}{\partial \xi^1}
!!  \f]
!!  \f[
!!      covT(2,2) = \frac{\partial y}{\partial \xi^2}
!!  \f]
!! 
!!  This function depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref Differentiate_2D_Lagrange <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>TYPE</B>(Lagrange)       :: interp <BR>
!! <B>REAL</B>(prec)              :: s(1:2), covT(1:2,1:2) <BR>
!!         .... <BR>
!!     x = this % CalculateMetrics( interp, s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myGeom <td> MappedGeometry_2D <td> An intialized and constructed 
!!                                                        MappedGeometry_2D data structure
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolant data structure that
!!                                                  contains the computational mesh.
!!   <tr> <td> in <th> s(1:2) <td> REAL(prec)  <td> Computational position where the covariant metric
!!                                                  tensor is desired.
!!   <tr> <td> out <th> covT(1:2,1:2) <td> REAL(prec)  <td> Covariant metric tensor estimated by 
!!                                                         differentiation of a Lagrange interpolant
!!                                                         of the mesh positions at the given 
!!                                                         computational coordinate
!!
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION CalculateMetrics_MappedGeometry_2D( myGeom, interp, s ) RESULT( covT )

   IMPLICIT NONE
   CLASS( MappedGeometry_2D ) :: myGeom
   TYPE( Lagrange )           :: interp
   REAL(prec)                 :: s(1:2)
   REAL(prec)                 :: covT(1:2,1:2)
 
      covT(1,1:2) = interp % Differentiate_2D( myGeom % x, s )
      covT(2,1:2) = interp % Differentiate_2D( myGeom % y, s )

 END FUNCTION CalculateMetrics_MappedGeometry_2D
!
!> \addtogroup MappedGeometry_2D_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateComputationalCoordinates
! 
!> \fn CalculateComputationalCoordinates_MappedGeometry_2D 
!! Given a physical position, the corresponding computational coordinate is estimated.
!!
!! From the physical coordinates \f$ ( x^*, y^* ) \f$, the mapping
!! \f[
!!      \vec{x}^* = \vec{x}(\vec{\xi})
!! \f] 
!! must be inverted to obtain the computational coordinates. This routine uses Newton's method
!! to approximate the solution. It is assumed that the computational domain is over the square,
!! \f$ [-1,1] \times [-1,1] \f$. If the estimated computational coordinate is outside of this
!! domain, the method is assumed unsuccessful. 
!!
!! Given an initial guess for the computational coordinate,\f$ \vec{\xi}_i \f$, Newton's method
!! proceeds by solving for a correction based on directly solving a linearized form of the mapping,
!! \f[
!!      \vec{\Delta \xi} = C^{-1} \vec{r}_i
!! \f]
!! where \f$ C \f$ is the \f$ 2 \times 2 \f$ covariant metric tensor and 
!! \f[
!!       \vec{r}_i = \vec{x}^* - \vec{x}(\vec{\xi}_i)
!! \f]
!!  is the residual at iterate "i". In this routine, C is inverted exactly. The computational 
!!  coordinate is updated, and the process is repeated until the residual magnitude falls below
!!  a specified tolerance (parameter "newtonTolerance" in \ref ConstantsDictionary.f90 ).
!!
!!  This subroutine depends on <BR>
!!   Module \ref MappedGeometry_2D_Class : Function \ref CalculateLocation_MappedGeometry_2D <BR>
!!   Module \ref MappedGeometry_2D_Class : Function \ref CalculateMetrics_MappedGeometry_2D <BR>
!!   Module \ref CommonRoutines          : Function \ref Invert_2x2 <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: x(1:2), s(1:2) <BR>
!! <B>LOGICAL</B>                 :: successful <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateComputationalCoordinates( interp, x, s, successful ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myGeom <td> MappedGeometry_2D <td> An intialized and constructed 
!!                                                        MappedGeometry_2D data structure
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolant data structure that
!!                                                  contains the computational mesh.
!!   <tr> <td> in <th> x(1:2) <td> REAL(prec)  <td> Physical location where we would like to determine
!!                                                  the computational coordinate
!!   <tr> <td> out <th> s(1:2) <td> REAL(prec)  <td> Computational coordinate corresponding to the
!!                                                  given physical coordinate
!!   <tr> <td> out (optional) <th> success <td> LOGICAL <td> A flag that determines if the Newton's
!!                                                           iteration was successful and returned
!!                                                           computational coordinates within
!!                                                           [-1,1]x[-1,1]
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_2D( myGeom, interp, x, s, success )

   IMPLICIT NONE
   CLASS( MappedGeometry_2D ), INTENT(in) :: myGeom
   TYPE( Lagrange ), INTENT(in)           :: interp
   REAL(prec), INTENT(in)                 :: x(1:2)
   REAL(prec), INTENT(out)                :: s(1:2)
   LOGICAL, INTENT(out), OPTIONAL         :: success
   ! LOCAL
   REAL(prec) :: dr(1:2), ds(1:2), A(1:2,1:2), Ainv(1:2,1:2)
   REAL(prec) :: thisX(1:2), thisS(1:2), resi
   INTEGER    :: i 

      thisS = ZERO ! Initial guess is at the origin
     
      IF( PRESENT(success) )THEN
         success = .FALSE.
      ENDIF
     
      DO i = 1, newtonMax
     
         ! Calculate the physical coordinate associated with the computational coordinate guess
         thisX = myGeom % CalculateLocation( interp, thisS )
     
         ! Calculate the residual
         dr = x - thisX
         resi = SQRT( DOT_PRODUCT( dr, dr ) )
     
         IF( resi < newtonTolerance )THEN
            s = thisS
            IF( PRESENT(success) .AND. &
                ABS(thisS(1))<=ONE .AND. &
                ABS(thisS(2))<=ONE )THEN
               success = .TRUE.
            ENDIF
            RETURN
         ENDIF
        
         A = myGeom % CalculateMetrics( interp, thisS ) ! Calculate the covariant metric tensor
         Ainv = Invert_2x2( A ) ! Invert the covariant metric tensor.
                                ! This matrix is ivertable as long as the Jacobian is non-zero.
         ds = MATMUL( Ainv, dr ) ! calculate the correction in the computational coordinate
         thisS = thisS + ds
     
      ENDDO
     
      ! Calculate the residual
      dr = x - thisX 
      resi = SQRT( DOT_PRODUCT( dr, dr ) )
      s = thisS
      IF( resi < newtonTolerance )THEN
         IF( PRESENT(success) .AND. &
             ABS(thisS(1))<=ONE .AND. &
             ABS(thisS(2))<=ONE )THEN
            success = .TRUE.
         ENDIF
         RETURN
      ELSE
         PRINT*,'Module MappedGeometryClass_2D.f90 : S/R CalculateComputationalCoordinates :'
         PRINT*,'Search for coordinates failed. Final residual norm :', resi
         RETURN
      ENDIF
      
     
 END SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_2D
!
!> \addtogroup MappedGeometry_2D_Class
!! @{ 
! ================================================================================================ !
! S/R ScaleGeometry
! 
!> \fn ScaleGeometry_MappedGeometry_2D 
!! Scales the physical coordinates and metric terms by a given x-scale and y-scale.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: xScale, yScale <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ScaleGeometry( interp, xScale, yScale ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_2D <td> 
!!                         On <B>input</B>, an intialized and constructed MappedGeometry_2D data 
!!                         structure <BR>
!!                         On <B>output</B>, the physical coordinates and metric terms have been
!!                         scaled by the given x and y scales.
!!   <tr> <td> in <th> interp <td> Lagrange <td> 2-D Lagrange interpolant data structure that
!!                                                  contains the computational mesh.
!!   <tr> <td> in <th> xScale <td> REAL(prec)  <td> Factor to multiply the physical x position
!!                                                  and metrics by.
!!   <tr> <td> in <th> yScale <td> REAL(prec)  <td> Factor to multiply the physical y position
!!                                                  and metrics by.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ScaleGeometry_MappedGeometry_2D( myGeom, interp, xScale, yScale )

   IMPLICIT NONE
   CLASS( MappedGeometry_2D ), INTENT(inout) :: myGeom
   TYPE( Lagrange ), INTENT(in)              :: interp
   REAL(prec), INTENT(in)                    :: xScale, yScale
   
         myGeom % x = xScale*( myGeom % x )
         myGeom % y = yScale*( myGeom % y )
         myGeom % xBound = xScale*( myGeom % xBound )
         myGeom % yBound = yScale*( myGeom % yBound )

         myGeom % dxds = xScale*( myGeom % dxds )
         myGeom % dxdp = xScale*( myGeom % dxdp )
         myGeom % dyds = yScale*( myGeom % dyds )
         myGeom % dydp = yScale*( myGeom % dydp )
          
         myGeom % J = xScale*yScale*( myGeom % J )

         ! Update the boundary metrics -- normals and normal lengths
         CALL myGeom % GenerateBoundaryMetrics( interp  )
         
 END SUBROUTINE ScaleGeometry_MappedGeometry_2D
!
!
!==================================================================================================!
!------------------------------------ File I/O Routines -------------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup MappedGeometry_2D_Class 
!! @{ 
! ================================================================================================ !
! S/R WriteTecplot 
! 
!> \fn WriteTecplot_MappedGeometry_2D 
!! Writes a tecplot of the metric terms at the physical coordinates contained in the MappedGeometry_2D
!! data structure.
!! 
!! Given a filename (say filename="foo"), the file written is "foo.tec". If a filename is not given
!! the file is called "LocalGeometry.tec".
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_2D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % WriteTecplot(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myGeom <td> MappedGeometry_2D <td> An intialized and constructed 
!!                                                        MappedGeometry_2D data structure
!!   <tr> <td> in (optional) <th> filename <td> CHARACTER <td> "Base"-name of the tecplot file.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteTecplot_MappedGeometry_2D( myGeom, filename )

  IMPLICIT NONE
  CLASS( MappedGeometry_2D ), INTENT(in) :: myGeom
  CHARACTER(*), INTENT(in), OPTIONAL     :: filename  
  ! Local
  INTEGER :: iX, iY, N, fUnit
  REAL(prec) :: x, y, dxds, dxdp, dyds, dydp, J
  
    N = myGeom % N
    
    IF( PRESENT(filename) )THEN
       OPEN( UNIT   = NEWUNIT(fUnit), &
             FILE   = TRIM(filename)//'.tec', &
             FORM   = 'formatted', & 
             STATUS = 'REPLACE' )
    ELSE
       OPEN( UNIT   = NEWUNIT(fUnit), &
             FILE   = 'LocalGeometry.tec', &
             FORM   = 'formatted', & 
             STATUS = 'REPLACE' )
    ENDIF

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Jacobian", "dxds", "dxdp", "dyds", "dydp" '
    
    WRITE(fUnit,*)  'ZONE T="el0", I=',N+1,', J=', N+1,',F=POINT'
    DO iY = 0, N
       DO iX = 0, N
       
          x = myGeom % x(iX,iY)
          y = myGeom % y(iX,iY)
          J = myGeom % J(iX,iY)
          
          dxds = myGeom % dxds(iX,iY)
          dxdp = myGeom % dxdp(iX,iY)
          dyds = myGeom % dyds(iX,iY)
          dydp = myGeom % dydp(iX,iY)
          
          WRITE(fUnit,*)  x, y, J, dxds, dxdp, dyds, dydp

       ENDDO
    ENDDO

    CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_MappedGeometry_2D
!
!
! ///////////////////////////////////// PRIVATE ////////////////////////////////////////////////// !
 FUNCTION TransfiniteInterpolation( curves, a, b ) RESULT( P )
 ! TransfiniteInterpolation
 !  Takes in the four curves (south, east, north, west) and evaluates the 
 !  bidirectional mapping at xi^1 = a, xi^2 = b. The south and north curves
 !  are assumed to be functions of xi^1, and are located at xi^2 = -1,1 respectively.
 !
 !   An area in 2-D is the assumed geometrical object generated from this evaluation.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Curve ) :: curves(1:4)
   REAL(prec)    :: a, b
   REAL(prec)    :: P(1:2)
   ! LOCAL
   REAL(prec)  :: P1(1:2), P2(1:2)
   REAL(prec)  :: Psouth(1:2), Pnorth(1:2)
   REAL(prec)  :: leftCurve(1:2), rightCurve(1:2)  
   REAL(prec)  :: c1(1:2), c2(1:2), c3(1:2), c4(1:2)
 

     ! Obtain the corner node locations (Tensor product portion of boolean summation)
     c1 = curves(1) % Evaluate( -ONE ) ! southwest
     c2 = curves(1) % Evaluate( ONE ) ! southeast
     c3 = curves(3) % Evaluate( -ONE ) ! northwest
     c4 = curves(3) % Evaluate( ONE ) ! northeast
 
   ! Do the unidirectional interpolation between the east and west curves at xi^2 = b
     leftCurve  = curves(4) % Evaluate( b ) ! west curve
     rightCurve = curves(2) % Evaluate( b ) ! east curve

     P1 = Unidirectional( leftCurve, rightCurve, a )

     leftCurve  = curves(1) % Evaluate( a ) ! south curve
     rightCurve = curves(3) % Evaluate( a ) ! north curve

     P2 = Unidirectional( leftCurve, rightCurve, b )

     ! Use the recursive formula
     Psouth = Unidirectional( c1, c2, a )
     Pnorth = Unidirectional( c3, c4, a )
   
     P = P1 + P2 - Unidirectional( Psouth, Pnorth, b ) 
     

 END FUNCTION TransfiniteInterpolation
!
 FUNCTION Unidirectional( valLeft, valRight, a ) RESULT( P )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: valLeft(1:2), valRight(1:2)
   REAL(prec) :: a
   REAL(prec) :: P(1:2)

       P = HALF*( (ONE - a)*valLeft + (ONE + a)*valRight )
    
 END FUNCTION Unidirectional

END MODULE MappedGeometry_2D_Class
