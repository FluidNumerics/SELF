! MappedGeometry_3D_Class.f90
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
! MappedGeometry_3D_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file MappedGeometry_3D_Class.f90
!! Contains the \ref MappedGeometry_3D_Class module, and <BR>
!! defines the \ref MappedGeometry_3D data-structure.

!> \defgroup MappedGeometry_3D_Class MappedGeometry_3D_Class 
!! This module defines the MappedGeometry_3D data-structure and its associated routines.

MODULE MappedGeometry_3D_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
! src/interp/
USE Lagrange_Class
#ifdef HAVE_CUDA
USE Lagrange_Cuda_Class
#endif
! src/geom/
USE Surface_Class

IMPLICIT NONE

!> \addtogroup MappedGeometry_3D_Class 
!! @{

!> \struct MappedGeometry_3D
!!  The MappedGeometry_3D class provides attributes and type-bound procedures for handling 
!!  curvilinear mappings between physical space and a reference computational space.
!!
!!  The MappedGeometry_3D class enables the implementation of spectral element methods in complicated
!!  geometry. Given the four bounding surfaces of an element, the internal geometry (physical node 
!!  positions, covariant basis vectors, and Jacobian) is constructed using transfinite interpolation
!!  with linear blending.
!!
!! <H2> MappedGeometry_3D </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> N <td> INTEGER  <td> Polynomial degree of the spectral element method
!!       <tr> <th> nHat(1:3,0:N,0:N,1:6) <td> REAL(prec) <td> Outward pointing element boundary normal vectors
!!                                                     along each element edge.
!!       <tr> <th> xBound(0:N,0:N,1:6) <td> REAL(prec) <td> x-position of each of the boundary surfaces
!!       <tr> <th> yBound(0:N,0:N,1:6) <td> REAL(prec) <td> y-position of each of the boundary surfaces
!!       <tr> <th> zBound(0:N,0:N,1:6) <td> REAL(prec) <td> z-position of each of the boundary surfaces
!!       <tr> <th> x(0:N,0:N,0:N) <td> REAL(prec) <td> x-position of the interior nodes
!!       <tr> <th> y(0:N,0:N,0:N) <td> REAL(prec) <td> y-position of the interior nodes
!!       <tr> <th> z(0:N,0:N,0:N) <td> REAL(prec) <td> z-position of the interior nodes
!!       <tr> <th> J(0:N,0:N,0:N) <td> REAL(prec) <td> Jacobian of the mapping at the interior nodes
!!       <tr> <th> dxds(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial x}{\partial \xi^1} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dxdp(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial x}{\partial \xi^2} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dxdq(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial x}{\partial \xi^3} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dyds(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial y}{\partial \xi^1} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dydp(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial y}{\partial \xi^2} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dydq(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial y}{\partial \xi^3} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dzds(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial z}{\partial \xi^1} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dzdp(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial z}{\partial \xi^2} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> dzdq(0:N,0:N,0:N) <td> REAL(prec) <td> \f$ \frac{\partial z}{\partial \xi^3} \f$ at
!!                                                    each of the interior nodes
!!       <tr> <th> Ja(0:N,0:N,0:N,1:3,1:3) <td> REAL(prec) <td> The contravariant metric tensor at
!!                                                         all of the interior nodes.
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref MappedGeometry_3D_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_MappedGeometry_3D
!!       <tr> <th> Build <td> Build_MappedGeometry_3D
!!       <tr> <th> Trash <td> Trash_MappedGeometry_3D
!!       <tr> <th> GenerateMesh <td> GenerateMesh_MappedGeometry_3D
!!       <tr> <th> GenerateMetrics <td> GenerateMetrics_MappedGeometry_3D
!!       <tr> <th> GenerateBoundaryMetrics <td> GenerateBoundaryMetrics_MappedGeometry_3D
!!       <tr> <th> ScaleGeometry <td> ScaleGeometry_MappedGeometry_3D
!!       <tr> <th> CalculateLocation <td> CalculateLocation_MappedGeometry_3D
!!       <tr> <th> CalculateMetrics <td> CalculateMetrics_MappedGeometry_3D
!!       <tr> <th> CalculateComputationalCoordinates <td> CalculateComputationalCoordinates_MappedGeometry_3D
!!       <tr> <th> WriteTecplot <td> WriteTecplot_MappedGeometry_3D
!!    </table>
!!

!>@}
   TYPE MappedGeometry_3D
      INTEGER                    :: N
      REAL(prec), ALLOCATABLE    :: nHat(:,:,:,:) 
      REAL(prec), ALLOCATABLE    :: xBound(:,:,:)
      REAL(prec), ALLOCATABLE    :: yBound(:,:,:) 
      REAL(prec), ALLOCATABLE    :: zBound(:,:,:) 
      REAL(prec), ALLOCATABLE    :: x(:,:,:), y(:,:,:), z(:,:,:)
      REAL(prec), ALLOCATABLE    :: J(:,:,:)    
      REAL(prec), ALLOCATABLE    :: dxds(:,:,:), dxdp(:,:,:), dxdq(:,:,:)
      REAL(prec), ALLOCATABLE    :: dyds(:,:,:), dydp(:,:,:), dydq(:,:,:)
      REAL(prec), ALLOCATABLE    :: dzds(:,:,:), dzdp(:,:,:), dzdq(:,:,:)
      REAL(prec), ALLOCATABLE    :: Ja(:,:,:,:,:)

      CONTAINS

      PROCEDURE :: Initialize => Initialize_MappedGeometry_3D
      PROCEDURE :: Build      => Build_MappedGeometry_3D
      PROCEDURE :: Trash      => Trash_MappedGeometry_3D

      PROCEDURE :: GenerateMesh => GenerateMesh_MappedGeometry_3D
      PROCEDURE :: GenerateMetrics => GenerateMetrics_MappedGeometry_3D
      PROCEDURE :: GenerateBoundaryMetrics => GenerateBoundaryMetrics_MappedGeometry_3D
      PROCEDURE :: ScaleGeometry => ScaleGeometry_MappedGeometry_3D
      PROCEDURE :: CalculateLocation => CalculateLocation_MappedGeometry_3D
      PROCEDURE :: CalculateMetrics => CalculateMetrics_MappedGeometry_3D
      PROCEDURE :: CalculateComputationalCoordinates => CalculateComputationalCoordinates_MappedGeometry_3D
      
      PROCEDURE :: ResetInternalMesh => ResetInternalMesh_MappedGeometry_3D
      
      PROCEDURE :: WriteTecplot => WriteTecplot_MappedGeometry_3D
   END TYPE MappedGeometry_3D


 INTEGER, PRIVATE :: nDims = 3
 PRIVATE :: TransfiniteInterpolation, LinearBlend
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup MappedGeometry_3D_Class
!! @{ 
! ================================================================================================ !
! S/R Initialize
! 
!> \fn Initialize_MappedGeometry_3D  
!! Allocates memory for each of the attributes of the MappedGeometry_3D Class and initializes all
!! arrays to zero.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>INTEGER</B>                 :: N <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize( N ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myGeom <td> MappedGeometry_3D <td> On output, an initialized MappedGeometry_3D
!!                                                         data structure
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the spectral element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_MappedGeometry_3D( myGeom, N )

  IMPLICIT NONE
  CLASS(MappedGeometry_3D), INTENT(out) :: myGeom
  INTEGER, INTENT(in)                   :: N

      myGeom % N = N
       
      ! Allocate space
      ALLOCATE( myGeom % dxds(0:N,0:N,0:N), myGeom % dxdp(0:N,0:N,0:N), myGeom % dxdq(0:N,0:N,0:N) )
      ALLOCATE( myGeom % dyds(0:N,0:N,0:N), myGeom % dydp(0:N,0:N,0:N), myGeom % dydq(0:N,0:N,0:N) )
      ALLOCATE( myGeom % dzds(0:N,0:N,0:N), myGeom % dzdp(0:N,0:N,0:N), myGeom % dzdq(0:N,0:N,0:N) )
      ALLOCATE( myGeom % Ja(0:N,0:N,0:N,1:3,1:3) )
      ALLOCATE( myGeom % J(0:N,0:N,0:N) )
      ALLOCATE( myGeom % x(0:N,0:N,0:N), myGeom % y(0:N,0:N,0:N), myGeom % z(0:N,0:N,0:N) )
      ALLOCATE( myGeom % xBound(0:N,0:N,1:nHexFaces) )
      ALLOCATE( myGeom % yBound(0:N,0:N,1:nHexFaces) )
      ALLOCATE( myGeom % zBound(0:N,0:N,1:nHexFaces) )
      ALLOCATE( myGeom % nHat(1:3,0:N,0:N,1:nHexFaces) )
      
      myGeom % dxds   = ZERO
      myGeom % dxdp   = ZERO
      myGeom % dxdq   = ZERO
      myGeom % dyds   = ZERO
      myGeom % dydp   = ZERO
      myGeom % dydq   = ZERO
      myGeom % dzds   = ZERO
      myGeom % dzdp   = ZERO
      myGeom % dzdq   = ZERO
      myGeom % J      = ZERO
      myGeom % x      = ZERO
      myGeom % y      = ZERO
      myGeom % z      = ZERO
      myGeom % xBound = ZERO
      myGeom % yBound = ZERO
      myGeom % zBound = ZERO

 END SUBROUTINE Initialize_MappedGeometry_3D
!
!> \addtogroup MappedGeometry_3D_Class
!! @{ 
! ================================================================================================ !
! S/R Build
! 
!> \fn Build_MappedGeometry_3D  
!! Allocates memory for each of the attributes of the MappedGeometry_3D Class and generates the
!! physical mesh and metric terms.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref MappedGeometry_3D_Class : S/R \ref Initialize_MappedGeometry_3D <BR>
!!   Module \ref MappedGeometry_3D_Class : S/R \ref GenerateMesh_MappedGeometry_3D <BR>
!!   Module \ref MappedGeometry_3D_Class : S/R \ref GenerateMetrics_MappedGeometry_3D <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>TYPE</B>(Surface)           :: boundaries(1:6) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( interp, boundaries ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myGeom <td> MappedGeometry_3D <td> On output, an initialized MappedGeometry_3D
!!                                                         data structure with the mesh and metric
!!                                                         terms filled in
!!   <tr> <td> in <th> interp <td> Lagrange <td> Lagrange interpolant that defines the 
!!                                                  reference computational mesh. 
!!   <tr> <td> in <th> mySurfaces(1:6)* <td> Surface <td> 3-D Surface that defines the each of the 
!!                                                      boundaries of the quadrilateral element. 
!!  </table>  
!!   
!!  * The boundary surfaces must be specified using the following mapping :
!!    boundaries(1)  -> South Boundary
!!    boundaries(2)  -> East Boundary
!!    boundaries(3)  -> North Boundary
!!    boundaries(4)  -> West Boundary
!!    boundaries(5)  -> Bottom Boundary
!!    boundaries(6)  -> Top Boundary
!!   For convenience, the integers corresponding to the boundary side have been given in the 
!!   file \ref ConstantsDictionary.f90 .
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_MappedGeometry_3D( myGeom, interp, mySurfaces )

   IMPLICIT NONE
   CLASS(MappedGeometry_3D), INTENT(out) :: myGeom
#ifdef HAVE_CUDA
   TYPE(Lagrange_Cuda), INTENT(in)       :: interp
#else
   TYPE(Lagrange), INTENT(in)            :: interp
#endif
   TYPE(Surface), INTENT(in)             :: mySurfaces(1:6)
   !LOCAL
   INTEGER :: N

      N = interp % N
      CALL myGeom % Initialize( N )
      ! Generate the mesh locations, and the mapping metrics
      CALL myGeom % GenerateMesh( interp, mySurfaces )
      CALL myGeom % GenerateMetrics( interp )
 
 END SUBROUTINE Build_MappedGeometry_3D
!
!> \addtogroup MappedGeometry_3D_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_MappedGeometry_3D_Class  
!! Frees memory held by the attributes of the MappedGeometry_3D data-structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_3D <td>
!!                         On <B>input</B>, a previously constructed MappedGeometry_3D data structure <BR>
!!                         On <B>output</B>, memory held by attributes is freed
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_MappedGeometry_3D( myGeom )

   IMPLICIT NONE
   CLASS(MappedGeometry_3D), INTENT(inout)  :: myGeom

      DEALLOCATE( myGeom % dxds, myGeom % dxdp, myGeom % dxdq )
      DEALLOCATE( myGeom % dyds, myGeom % dydp, myGeom % dydq )
      DEALLOCATE( myGeom % dzds, myGeom % dzdp, myGeom % dzdq )
      DEALLOCATE( myGeom % J, myGeom % x, myGeom % y, myGeom % z )
      DEALLOCATE( myGeom % xBound, myGeom % yBound, myGeom % zBound )
      DEALLOCATE( myGeom % nHat )
 
 END SUBROUTINE Trash_MappedGeometry_3D
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup MappedGeometry_3D_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMesh 
! 
!> \fn GenerateMesh_MappedGeometry_3D  
!! Generates the physical interior and boundary positions at each of the computational mesh points.
!! 
!! Given the four boundary surfaces and an interpolant that stores the computational mesh points,
!! the physical mesh positions are generated using transfinite interpolation with linear blending.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref MappedGeometry_3D_Class : PRIVATE Function TransfiniteInterpolation <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>TYPE</B>(Surface)           :: boundaries(1:6) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateMesh( interp, boundaries ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_3D <td> 
!!                         On <B>input</B>, an initialized MappedGeometry_3D data structure, <BR>
!!                         On <B>output</B>, the mesh physical positions (interior and boundary)
!!                         are filled in.
!!   <tr> <td> in <th> interp <td> Lagrange <td> Lagrange interpolant that defines the 
!!                                                  reference computational mesh. 
!!   <tr> <td> in <th> mySurfaces(1:6) <td> Surface <td> 3-D Surface that defines the each of the 
!!                                                      boundaries of the quadrilateral element. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GenerateMesh_MappedGeometry_3D( myGeom, interp, theSurfaces )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in)         :: interp
#else
   TYPE( Lagrange ), INTENT(in)              :: interp
#endif
   TYPE( Surface ), INTENT(in)               :: theSurfaces(1:6)
   ! Local
   INTEGER    :: iS, iP, iQ, N
   REAL(prec) :: s(0:interp % N), p, x(1:3)
   
      N = interp % N
      s = interp % s
      
      DO iQ = 0, N
         DO iP = 0,N
            DO iS = 0,N
               x = TransfiniteInterpolation( theSurfaces, s(iS), s(iP), s(iQ) )
               myGeom % x(iS,iP,iQ) = x(1)
               myGeom % y(iS,iP,iQ) = x(2)
               myGeom % z(iS,iP,iQ) = x(3)
            ENDDO
         ENDDO 
      ENDDO 
!      print*, s
!      print*,'---------- south surface ----------------'
!      print*, theSurfaces(south) % x(1,0,1)
!      print*, theSurfaces(south) % x(1,0,2)
!      print*, theSurfaces(south) % x(1,0,3)
!      print*,'-----------------------------------------'
      ! Do the boundary locations
      DO iP = 0, N
         DO iS = 0, N
            p = -ONE  ! south boundary
            x = TransfiniteInterpolation( theSurfaces, s(iS), p, s(iP) )
            myGeom % xBound(iS,iP,south) = x(1)
            myGeom % yBound(iS,iP,south) = x(2)
            myGeom % zBound(iS,iP,south) = x(3)
            ! west boundary
            x = TransfiniteInterpolation( theSurfaces, p, s(iS), s(iP) )
            myGeom % xBound(iS,iP,west) = x(1)
            myGeom % yBound(iS,iP,west) = x(2)
            myGeom % zBound(iS,iP,west) = x(3)
            ! bottom boundary
            x = TransfiniteInterpolation( theSurfaces, s(iS), s(iP), p )
            myGeom % xBound(iS,iP,bottom) = x(1)
            myGeom % yBound(iS,iP,bottom) = x(2)
            myGeom % zBound(iS,iP,bottom) = x(3)
            
            p = ONE  ! north boundary
            x = TransfiniteInterpolation( theSurfaces, s(iS), p, s(iP) )
            myGeom % xBound(iS,iP,north) = x(1)
            myGeom % yBound(iS,iP,north) = x(2)
            myGeom % zBound(iS,iP,north) = x(3)
            ! east boundary
            x = TransfiniteInterpolation( theSurfaces, p, s(iS), s(iP) )
            myGeom % xBound(iS,iP,east) = x(1)
            myGeom % yBound(iS,iP,east) = x(2)
            myGeom % zBound(iS,iP,east) = x(3)
            ! top boundary
            x = TransfiniteInterpolation( theSurfaces, s(iS), s(iP), p )
            myGeom % xBound(iS,iP,top) = x(1)
            myGeom % yBound(iS,iP,top) = x(2)
            myGeom % zBound(iS,iP,top) = x(3)
         ENDDO
      ENDDO
! PRINT*, myGeom % xBound(:,:,south)
! PRINT*, myGeom % yBound(:,:,south)
! PRINT*, myGeom % zBound(:,:,south)
! STOP
 END SUBROUTINE GenerateMesh_MappedGeometry_3D
!
 SUBROUTINE ResetInternalMesh_MappedGeometry_3D( myGeom, interp )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in)         :: interp
#else
   TYPE( Lagrange ), INTENT(in)              :: interp
#endif
   ! Local
   INTEGER    :: i, j, k
   REAL(prec) :: s(0:interp % N), p, x(1:3)
   
      
      DO k = 0, interp % N
         DO j = 0, interp % N
            DO i = 0, interp % N
               x = TransfiniteInterpolation_Alt( interp, &
                                                 myGeom % xBound, &
                                                 myGeom % yBound, &
                                                 myGeom % zBound, &
                                                 interp % s(i), &
                                                 interp % s(j), &
                                                 interp % s(k) )
               myGeom % x(i,j,k) = x(1)
               myGeom % y(i,j,k) = x(2)
               myGeom % z(i,j,k) = x(3)
            ENDDO
         ENDDO 
      ENDDO 



 END SUBROUTINE ResetInternalMesh_MappedGeometry_3D
!
!> \addtogroup MappedGeometry_3D_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMetrics 
! 
!> \fn GenerateMetrics_MappedGeometry_3D  
!! Generates and stores the covariant and contravariant basis vector components, Jacobian of the 
!! mapping, and the outward pointing boundary normal vectors.
!! 
!! Once the mesh has been generated, we effectively have
!! \f[ 
!!     \vec{x} = \vec{x}(\vec{\xi})
!! \f]
!! represented as a 3-D Lagrange interpolant. Differentiation of the interpolant allows us to 
!! estimate the covariant basis vectors. The determinant of the covariant metric tensor gives
!! the Jacobian of the mapping at each location in the computational mesh
!! \f[
!!     J = | \frac{\partial \vec{x}}{\partial \vec{\xi}} |
!! \f]
!! 
!!  The contravariant basis vectors are computed using the curl-invariant form described in 
!!  D.A. Kopriva, 2006, "Metric Identities and the Discontinuous Spectral Element Method on 
!!  Curvilinear Meshes", J. Sci. Comp., 26-3, DOI: 10.1007/s10915-005-9070-8
!!
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref ApplyDerivativeMatrix_3D_Lagrange <BR>
!!   Module \ref MappedGeometry_3D_Class : S/R \ref GenerateBoundaryMetrics_MappedGeometry_3D <BR>
!! 
!!  ** To produce meaningful output, GenerateMesh_MappedGeometry_3D must be called prior to calling
!!     this routine.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateMetrics( interp ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_3D <td> 
!!                         On <B>input</B>, an initialized MappedGeometry_3D data structure, <BR>
!!                         On <B>output</B>, the metric terms are filled in, including the
!!                         outward pointing normal vectors on the element boundaries
!!   <tr> <td> in <th> interp <td> Lagrange <td>  Lagrange interpolant that defines the 
!!                                                  reference computational mesh 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GenerateMetrics_MappedGeometry_3D( myGeom, interp )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in)              :: interp
#else
   TYPE( Lagrange ), INTENT(in)              :: interp
#endif
   ! Local
   INTEGER    :: iS, iP, iQ, N
   REAL(prec) :: cv(1:3,1:3)
   REAL(prec) :: covT(0:interp % N, 0:interp % N, 0:interp % N, 1:3,1:3)
   REAL(prec) :: v(0:interp % N, 0:interp % N, 0:interp % N,1:3)
   REAL(prec) :: Dv1(0:interp % N, 0:interp % N, 0:interp % N,1:3)
   REAL(prec) :: Dv2(0:interp % N, 0:interp % N, 0:interp % N,1:3)
   REAL(prec) :: Dv3(0:interp % N, 0:interp % N, 0:interp % N,1:3)
   
      N = interp % N
      
      covT(0:N,0:N,0:N,1,1:3) = interp % ApplyDerivativeMatrix_3D( myGeom % x )
      covT(0:N,0:N,0:N,2,1:3) = interp % ApplyDerivativeMatrix_3D( myGeom % y )
      covT(0:N,0:N,0:N,3,1:3) = interp % ApplyDerivativeMatrix_3D( myGeom % z )
      
      DO iQ = 0, N
         DO iP = 0, N
            DO iS = 0, N

               myGeom % dxds(iS,iP,iQ) = covT(iS,iP,iQ,1,1)
               myGeom % dxdp(iS,iP,iQ) = covT(iS,iP,iQ,1,2)
               myGeom % dxdq(iS,iP,iQ) = covT(iS,iP,iQ,1,3)
               myGeom % dyds(iS,iP,iQ) = covT(iS,iP,iQ,2,1)
               myGeom % dydp(iS,iP,iQ) = covT(iS,iP,iQ,2,2)
               myGeom % dydq(iS,iP,iQ) = covT(iS,iP,iQ,2,3)
               myGeom % dzds(iS,iP,iQ) = covT(iS,iP,iQ,3,1)
               myGeom % dzdp(iS,iP,iQ) = covT(iS,iP,iQ,3,2)
               myGeom % dzdq(iS,iP,iQ) = covT(iS,iP,iQ,3,3)
               
               cv = covT(iS,iP,iQ,1:3,1:3)
               myGeom % J(iS,iP,iQ) = Determinant( cv, 3 )
               
            ENDDO
         ENDDO
      ENDDO 
      
      ! Generate the contravariant metric tensor a la Kopriva (2006)
      !Ja_1
      DO iQ = 0, N
         DO iP = 0, N
            DO iS = 0, N
               v(iS,iP,iQ,1:3)  = -myGeom % z(iS,iP,iQ)*covT(iS,iP,iQ,2,1:3) 
            ENDDO
         ENDDO
      ENDDO
      Dv1 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,1) ) ! ( dv1/ds, dv1/dp, dv1/dq )
      Dv2 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,2) ) ! ( dv2/ds, dv2/dp, dv2/dq )
      Dv3 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,3) ) ! ( dv3/ds, dv3/dp, dv3/dq )
      ! Take the curl to obtain the first dimension of each of the contravariant basis vectors
      ! The contravariant metric tensor stores each contravariant basis vector in each column
      ! of the tensor
      myGeom % Ja(0:N,0:N,0:N,1,1) = ( Dv3(0:N,0:N,0:N,2) - Dv2(0:N,0:N,0:N,3) )
      myGeom % Ja(0:N,0:N,0:N,1,2) = -( Dv3(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,3) )
      myGeom % Ja(0:N,0:N,0:N,1,3) = ( Dv2(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,2) )
      !Ja_2
      DO iQ = 0, N
         DO iP = 0, N
            DO iS = 0, N
               v(iS,iP,iQ,1:3)  = -myGeom % x(iS,iP,iQ)*covT(iS,iP,iQ,3,1:3) 
            ENDDO
         ENDDO
      ENDDO
      Dv1 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,1) ) ! ( dv1/ds, dv1/dp, dv1/dq )
      Dv2 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,2) ) ! ( dv2/ds, dv2/dp, dv2/dq )
      Dv3 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,3) ) ! ( dv3/ds, dv3/dp, dv3/dq )
      ! Take the curl to obtain the first dimension of each of the contravariant basis vectors
      ! The contravariant metric tensor stores each contravariant basis vector in each column
      ! of the tensor
      myGeom % Ja(0:N,0:N,0:N,2,1) = ( Dv3(0:N,0:N,0:N,2) - Dv2(0:N,0:N,0:N,3) )
      myGeom % Ja(0:N,0:N,0:N,2,2) = -( Dv3(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,3) )
      myGeom % Ja(0:N,0:N,0:N,2,3) = ( Dv2(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,2) )
      !Ja_3
      DO iQ = 0, N
         DO iP = 0, N
            DO iS = 0, N
               v(iS,iP,iQ,1:3)  = -myGeom % y(iS,iP,iQ)*covT(iS,iP,iQ,1,1:3) 
            ENDDO
         ENDDO
      ENDDO
      Dv1 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,1) ) ! ( dv1/ds, dv1/dp, dv1/dq )
      Dv2 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,2) ) ! ( dv2/ds, dv2/dp, dv2/dq )
      Dv3 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,3) ) ! ( dv3/ds, dv3/dp, dv3/dq )
      ! Take the curl to obtain the first dimension of each of the contravariant basis vectors
      ! The contravariant metric tensor stores each contravariant basis vector in each column
      ! of the tensor
      myGeom % Ja(0:N,0:N,0:N,3,1) = ( Dv3(0:N,0:N,0:N,2) - Dv2(0:N,0:N,0:N,3) )
      myGeom % Ja(0:N,0:N,0:N,3,2) = -( Dv3(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,3) )
      myGeom % Ja(0:N,0:N,0:N,3,3) = ( Dv2(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,2) )
      
      CALL myGeom % GenerateBoundaryMetrics( interp )

 END SUBROUTINE GenerateMetrics_MappedGeometry_3D
!
!> \addtogroup MappedGeometry_3D_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMetrics 
! 
!> \fn GenerateMetrics_MappedGeometry_3D  
!! Generates and stores the outward pointing boundary normal vectors.
!!
!!  The outward pointing boundary normal vectors are equivalent to the contravariant basis vectors
!!  evaluated at the element boundaries. These are computed here by differentiating the Lagrange
!!  interpolant of the mesh positions and the computational mesh boundaries.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref Differentiate_3D_Lagrange <BR>
!! 
!!  ** To produce meaningful output, GenerateMesh_MappedGeometry_3D must be called prior to calling
!!     this routine.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateBoundaryMetrics( interp ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_3D <td> 
!!                         On <B>input</B>, an initialized MappedGeometry_3D data structure, <BR>
!!                         On <B>output</B>, the outward pointing normal vectors on the element
!!                         boundaries are filled in
!!   <tr> <td> in <th> interp <td> Lagrange <td> Lagrange interpolant that defines the 
!!                                                  reference computational mesh 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
  SUBROUTINE GenerateBoundaryMetrics_MappedGeometry_3D( myGeom, interp  )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in)         :: interp
#else
   TYPE( Lagrange ), INTENT(in)              :: interp
#endif
   ! Local
   INTEGER    :: iS, iP, N
   REAL(prec) :: s(0:interp % N), p
   REAL(prec) :: Jain(0:interp % N,0:interp % N,0:interp % N)
   REAL(prec) :: J, signJ, nx, ny, nz
   REAL(prec) :: node(1:3)
   
      N = interp % N
      s = interp % s
           ! Do the boundary locations
      DO iP = 0, N
         DO iS = 0, N
         
            p = -ONE  ! bottom, south, and west boundaries
            
            !bottom boundary
            node = (/s(iS), s(iP), p /)
            J = interp % Interpolate_3D( myGeom % J, node ) !Determinant( cv, nDims )
            signJ = SIGN(ONE,J)                    
            ! Setting bottom boundary normal
            Jain = myGeom % Ja(0:N,0:N,0:N,1,3)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,2,3)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,3,3)
            nz = interp % Interpolate_3D( Jain, node )
            myGeom % nHat(1:nDims,iS,iP,bottom) = -signJ*(/ nx, ny, nz /)
            
            node = (/ s(iS), p, s(iP) /)
            J = interp % Interpolate_3D( myGeom % J, node ) !Determinant( cv, nDims )
            signJ = SIGN(ONE,J)                    
            ! Setting southern boundary normal
            Jain = myGeom % Ja(0:N,0:N,0:N,1,2)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,2,2)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,3,2)
            nz = interp % Interpolate_3D( Jain, node )
            myGeom % nHat(1:nDims,iS,iP,south)= -signJ*(/ nx, ny, nz /)
            
            ! west boundary
            node = (/ p, s(iS), s(iP) /)
            J = interp % Interpolate_3D( myGeom % J, node ) !Determinant( cv, nDims )
            signJ = SIGN(ONE,J)                    
            ! Setting western boundary normal
            Jain = myGeom % Ja(0:N,0:N,0:N,1,1)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,2,1)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,3,1)
            nz = interp % Interpolate_3D( Jain, node )
            myGeom % nHat(1:nDims,iS,iP,west) = -signJ*(/ nx, ny, nz /)
             
            p = ONE  ! top, north, and east boundaries
            
            !top boundary
            node = (/s(iS), s(iP), p /)
            J = interp % Interpolate_3D( myGeom % J, node )!Determinant( cv, nDims )
            signJ = SIGN(ONE,J)      
            ! Setting top boundary normal
            Jain = myGeom % Ja(0:N,0:N,0:N,1,3)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,2,3)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,3,3)
            nz = interp % Interpolate_3D( Jain, node )
            myGeom % nHat(1:nDims,iS,iP,top) = signJ*(/ nx, ny, nz /)
            
            !north boundary
            node = (/ s(iS), p, s(iP) /)
            J = interp % Interpolate_3D( myGeom % J, node ) !Determinant( cv, nDims )
            signJ = SIGN(ONE,J)    
            ! Setting southern boundary normal
            Jain = myGeom % Ja(0:N,0:N,0:N,1,2)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,2,2)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,3,2)
            nz = interp % Interpolate_3D( Jain, node )
            myGeom % nHat(1:nDims,iS,iP,north) = signJ*(/ nx, ny, nz /)
            
            ! east boundary
            node = (/ p, s(iS), s(iP) /)
            J = interp % Interpolate_3D( myGeom % J, node ) !Determinant( cv, nDims )
            signJ = SIGN(ONE,J)                    
            ! Setting eastern boundary normal
            Jain = myGeom % Ja(0:N,0:N,0:N,1,1)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,2,1)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myGeom % Ja(0:N,0:N,0:N,3,1)
            nz = interp % Interpolate_3D( Jain, node )
            myGeom % nHat(1:nDims,iS,iP,east) = signJ*(/ nx, ny, nz /)
            
         ENDDO
      ENDDO

 END SUBROUTINE GenerateBoundaryMetrics_MappedGeometry_3D
!
!> \addtogroup MappedGeometry_3D_Class
!! @{ 
! ================================================================================================ !
! Function CalculateLocation 
! 
!> \fn CalculateLocation_MappedGeometry_3D  
!! Given a computational coordinate, the physical coordinate is calculated using Lagrange 
!! interpolation.
!! 
!!  This function depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref Interpolate_3D_Lagrange <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: s(1:2), x(1:2) <BR>
!!         .... <BR>
!!     x = this % CalculateLocation( interp, s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myGeom <td> MappedGeometry_3D <td> An intialized and constructed 
!!                                                        MappedGeometry_3D data structure
!!   <tr> <td> in <th> interp <td> Lagrange <td> Lagrange interpolant data structure that
!!                                                  contains the computational mesh.
!!   <tr> <td> in <th> s(1:3) <td> REAL(prec)  <td> Computational location where the physical 
!!                                                  position is desired.
!!   <tr> <td> out <th> x(1:3) <td> REAL(prec)  <td> Physical location estimated by interpolation
!!                                                  onto the given computational position
!!
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION CalculateLocation_MappedGeometry_3D( myGeom, interp, s ) RESULT( x )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ) :: myGeom
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda )      :: interp
#else
   TYPE( Lagrange )           :: interp
#endif
   REAL(prec)                 :: s(1:3)
   REAL(prec)                 :: x(1:3)
  
      x(1) = interp % Interpolate_3D( myGeom % x, s )
      x(2) = interp % Interpolate_3D( myGeom % y, s )
      x(3) = interp % Interpolate_3D( myGeom % z, s )
  
 END FUNCTION CalculateLocation_MappedGeometry_3D
!
!> \addtogroup MappedGeometry_3D_Class
!! @{ 
! ================================================================================================ !
! Function CalculateMetrics 
! 
!> \fn CalculateMetrics_MappedGeometry_3D  
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
!!      covT(1,3) = \frac{\partial x}{\partial \xi^3}
!!  \f]
!!  \f[
!!      covT(2,1) = \frac{\partial y}{\partial \xi^1}
!!  \f]
!!  \f[
!!      covT(2,2) = \frac{\partial y}{\partial \xi^2}
!!  \f]
!!  \f[
!!      covT(2,3) = \frac{\partial y}{\partial \xi^3}
!!  \f]
!!  \f[
!!      covT(3,1) = \frac{\partial z}{\partial \xi^1}
!!  \f]
!!  \f[
!!      covT(3,2) = \frac{\partial z}{\partial \xi^2}
!!  \f]
!!  \f[
!!      covT(3,3) = \frac{\partial z}{\partial \xi^3}
!!  \f]
!! 
!!  This function depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref Differentiate_3D_Lagrange <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: s(1:3), covT(1:3,1:3) <BR>
!!         .... <BR>
!!     x = this % CalculateMetrics( interp, s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myGeom <td> MappedGeometry_3D <td> An intialized and constructed 
!!                                                        MappedGeometry_3D data structure
!!   <tr> <td> in <th> interp <td> Lagrange <td> Lagrange interpolant data structure that
!!                                                  contains the computational mesh.
!!   <tr> <td> in <th> s(1:3) <td> REAL(prec)  <td> Computational position where the covariant metric
!!                                                  tensor is desired.
!!   <tr> <td> out <th> covT(1:3,1:3) <td> REAL(prec)  <td> Covariant metric tensor estimated by 
!!                                                         differentiation of a Lagrange interpolant
!!                                                         of the mesh positions at the given 
!!                                                         computational coordinate
!!
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION CalculateMetrics_MappedGeometry_3D( myGeom, interp, s ) RESULT( covT )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ) :: myGeom
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda )      :: interp
#else
   TYPE( Lagrange )           :: interp
#endif
   REAL(prec)                 :: s(1:3)
   REAL(prec)                 :: covT(1:3,1:3)
 
      covT(1,1:3) = interp % Differentiate_3D( myGeom % x, s )
      covT(2,1:3) = interp % Differentiate_3D( myGeom % y, s )
      covT(3,1:3) = interp % Differentiate_3D( myGeom % z, s )

 END FUNCTION CalculateMetrics_MappedGeometry_3D
!
!> \addtogroup MappedGeometry_3D_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateComputationalCoordinates
! 
!> \fn CalculateComputationalCoordinates_MappedGeometry_3D 
!! Given a physical position, the corresponding computational coordinate is estimated.
!!
!! From the physical coordinates \f$ ( x^*, y^*, z^* ) \f$, the mapping
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
!! where \f$ C \f$ is the \f$ 3 \times 3 \f$ covariant metric tensor and 
!! \f[
!!       \vec{r}_i = \vec{x}^* - \vec{x}(\vec{\xi}_i)
!! \f]
!!  is the residual at iterate "i". In this routine, C is inverted exactly. The computational 
!!  coordinate is updated, and the process is repeated until the residual magnitude falls below
!!  a specified tolerance (parameter "newtonTolerance" in \ref ConstantsDictionary.f90 ).
!!
!!  This subroutine depends on <BR>
!!   Module \ref MappedGeometry_3D_Class : Function \ref CalculateLocation_MappedGeometry_3D <BR>
!!   Module \ref MappedGeometry_3D_Class : Function \ref CalculateMetrics_MappedGeometry_3D <BR>
!!   Module \ref CommonRoutines          : Function \ref Invert_2x2 <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: x(1:3), s(1:3) <BR>
!! <B>LOGICAL</B>                 :: successful <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateComputationalCoordinates( interp, x, s, successful ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myGeom <td> MappedGeometry_3D <td> An intialized and constructed 
!!                                                        MappedGeometry_3D data structure
!!   <tr> <td> in <th> interp <td> Lagrange <td> Lagrange interpolant data structure that
!!                                                  contains the computational mesh.
!!   <tr> <td> in <th> x(1:3) <td> REAL(prec)  <td> Physical location where we would like to determine
!!                                                  the computational coordinate
!!   <tr> <td> out <th> s(1:3) <td> REAL(prec)  <td> Computational coordinate corresponding to the
!!                                                  given physical coordinate
!!   <tr> <td> out (optional) <th> success <td> LOGICAL <td> A flag that determines if the Newton's
!!                                                           iteration was successful and returned
!!                                                           computational coordinates within
!!                                                           [-1,1]x[-1,1]
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_3D( myGeom, interp, x, s, success )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ), INTENT(in) :: myGeom
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda )                  :: interp
#else
   TYPE( Lagrange ), INTENT(in)           :: interp
#endif
   REAL(prec), INTENT(in)                 :: x(1:3)
   REAL(prec), INTENT(out)                :: s(1:3)
   LOGICAL, INTENT(out), OPTIONAL         :: success
   ! LOCAL
   REAL(prec) :: dr(1:3), ds(1:3), A(1:3,1:3), Ainv(1:3,1:3)
   REAL(prec) :: thisX(1:3), thisS(1:3), resi
   INTEGER    :: i 

      thisS = 0.0_prec ! Initial guess is at the origin
     
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
                ABS(thisS(1))<=ONE+newtonTolerance .AND. &
                ABS(thisS(2))<=ONE+newtonTolerance .AND. &
                ABS(thisS(3))<=ONE+newtonTolerance )THEN
               success = .TRUE.
            ENDIF
            RETURN
         ENDIF
        
         A = myGeom % CalculateMetrics( interp, thisS ) ! Calculate the covariant metric tensor
       !  print*, A
         Ainv = Invert_3x3( A ) ! Invert the covariant metric tensor.
                                ! This matrix is ivertable as long as the Jacobian is non-zero.
       !  print*, Ainv
       !  STOP                       
         ds = MATMUL( Ainv, dr ) ! calculate the correction in the computational coordinate
         thisS = thisS + ds
     
      ENDDO
     
      ! Calculate the residual
      dr = x - thisX 
      resi = SQRT( DOT_PRODUCT( dr, dr ) )
      s = thisS
      IF( resi < newtonTolerance )THEN
         IF( PRESENT(success) .AND. &
             ABS(thisS(1))<=ONE+newtonTolerance .AND. &
             ABS(thisS(2))<=ONE+newtonTolerance .AND. &
             ABS(thisS(3))<=ONE+newtonTolerance )THEN
            success = .TRUE.
         ENDIF
         RETURN
      ELSE
        ! PRINT*,'Module MappedGeometryClass_3D.f90 : S/R CalculateComputationalCoordinates :'
        ! PRINT*,'Search for coordinates failed. Final residual norm :', resi
         RETURN
      ENDIF
      
     
 END SUBROUTINE CalculateComputationalCoordinates_MappedGeometry_3D
!
!> \addtogroup MappedGeometry_3D_Class
!! @{ 
! ================================================================================================ !
! S/R ScaleGeometry
! 
!> \fn ScaleGeometry_MappedGeometry_3D 
!! Scales the physical coordinates and metric terms by a given x-scale and y-scale.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: xScale, yScale <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ScaleGeometry( interp, xScale, yScale ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myGeom <td> MappedGeometry_3D <td> 
!!                         On <B>input</B>, an intialized and constructed MappedGeometry_3D data 
!!                         structure <BR>
!!                         On <B>output</B>, the physical coordinates and metric terms have been
!!                         scaled by the given x and y scales.
!!   <tr> <td> in <th> interp <td> Lagrange <td> 3-D Lagrange interpolant data structure that
!!                                                  contains the computational mesh.
!!   <tr> <td> in <th> xScale <td> REAL(prec)  <td> Factor to multiply the physical x position
!!                                                  and metrics by.
!!   <tr> <td> in <th> yScale <td> REAL(prec)  <td> Factor to multiply the physical y position
!!                                                  and metrics by.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ScaleGeometry_MappedGeometry_3D( myGeom, interp, xScale, yScale, zScale )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ), INTENT(inout) :: myGeom
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ), INTENT(in)         :: interp
#else
   TYPE( Lagrange ), INTENT(in)              :: interp
#endif
   REAL(prec), INTENT(in)                    :: xScale, yScale, zScale
   
         myGeom % x = xScale*( myGeom % x )
         myGeom % y = yScale*( myGeom % y )
         myGeom % z = zScale*( myGeom % z )
         myGeom % xBound = xScale*( myGeom % xBound )
         myGeom % yBound = yScale*( myGeom % yBound )
         myGeom % zBound = zScale*( myGeom % zBound )

!         myGeom % dxds = xScale*( myGeom % dxds )
!         myGeom % dxdp = xScale*( myGeom % dxdp )
!         myGeom % dxdq = xScale*( myGeom % dxdq )
!         myGeom % dyds = yScale*( myGeom % dyds )
!         myGeom % dydp = yScale*( myGeom % dydp )
!         myGeom % dydq = yScale*( myGeom % dydq )
!         myGeom % dzds = zScale*( myGeom % dzds )
!         myGeom % dzdp = zScale*( myGeom % dzdp )
!         myGeom % dzdq = zScale*( myGeom % dzdq )
          
!         myGeom % J = xScale*yScale*zScale*( myGeom % J )

         ! Update the boundary metrics -- normals and normal lengths
         CALL myGeom % GenerateMetrics( interp )
         CALL myGeom % GenerateBoundaryMetrics( interp  )
         
 END SUBROUTINE ScaleGeometry_MappedGeometry_3D
!
!
!==================================================================================================!
!------------------------------------ File I/O Routines -------------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup MappedGeometry_3D_Class 
!! @{ 
! ================================================================================================ !
! S/R WriteTecplot 
! 
!> \fn WriteTecplot_MappedGeometry_3D 
!! Writes a tecplot of the metric terms at the physical coordinates contained in the MappedGeometry_3D
!! data structure.
!! 
!! Given a filename (say filename="foo"), the file written is "foo.tec". If a filename is not given
!! the file is called "LocalGeometry.tec".
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MappedGeometry_3D) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % WriteTecplot(  ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myGeom <td> MappedGeometry_3D <td> An intialized and constructed 
!!                                                        MappedGeometry_3D data structure
!!   <tr> <td> in (optional) <th> filename <td> CHARACTER <td> "Base"-name of the tecplot file.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WriteTecplot_MappedGeometry_3D( myGeom, filename )

   IMPLICIT NONE
   CLASS( MappedGeometry_3D ), INTENT(in) :: myGeom
   CHARACTER(*), INTENT(in), OPTIONAL     :: filename  
   ! Local
   INTEGER :: iX, iY, iZ, N, fUnit
   REAL(prec) :: x, y, z, dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq, J
  
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

      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "Jacobian", "dxds", "dxdp", "dxdq",'//&
                                 '"dyds", "dydp", "dydq", "dzds", "dzdp", "dzdq"'
      WRITE(fUnit,*)  'ZONE T="el0", I=',N+1,', J=', N+1,', K=', N+1,',F=POINT'
      
      DO iZ = 0, N
         DO iY = 0, N
            DO iX = 0, N
          
               x = myGeom % x(iX,iY,iZ)
               y = myGeom % y(iX,iY,iZ)
               z = myGeom % z(iX,iY,iZ)
               J = myGeom % J(iX,iY,iZ)
             
               dxds = myGeom % dxds(iX,iY,iZ)
               dxdp = myGeom % dxdp(iX,iY,iZ)
               dxdq = myGeom % dxdq(iX,iY,iZ)
               dyds = myGeom % dyds(iX,iY,iZ)
               dydp = myGeom % dydp(iX,iY,iZ)
               dydq = myGeom % dydq(iX,iY,iZ)
               dzds = myGeom % dzds(iX,iY,iZ)
               dzdp = myGeom % dzdp(iX,iY,iZ)
               dzdq = myGeom % dzdq(iX,iY,iZ)
             
               WRITE(fUnit,*)  x, y, z, J, dxds, dxdp, dxdq, dyds, dydp, dydq, dzds, dzdp, dzdq
   
            ENDDO
         ENDDO
      ENDDO

      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = 'SouthFace.tec', &
            FORM   = 'formatted', & 
            STATUS = 'REPLACE' )
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "nx", "ny", "nz" '
      WRITE(fUnit,*)  'ZONE T="el0", I=',N+1,', J=', N+1,',F=POINT'
      DO iY = 0, N
         DO iX = 0, N
            WRITE(fUnit,*)  myGeom % xBound(iX,iY,South), &
                            myGeom % yBound(iX,iY,South), &
                            myGeom % zBound(iX,iY,South), &
                            myGeom % nHat(1:3,iX,iY,South)
         ENDDO
      ENDDO
      CLOSE(UNIT=fUnit)
      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = 'NorthFace.tec', &
            FORM   = 'formatted', & 
            STATUS = 'REPLACE' )
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "nx", "ny", "nz" '
      WRITE(fUnit,*)  'ZONE T="el0", I=',N+1,', J=', N+1,',F=POINT'
      DO iY = 0, N
         DO iX = 0, N

            WRITE(fUnit,*)  myGeom % xBound(iX,iY,North), &
                            myGeom % yBound(iX,iY,North), &
                            myGeom % zBound(iX,iY,North), &
                            myGeom % nHat(1:3,iX,iY,North)
         ENDDO
      ENDDO
      CLOSE(UNIT=fUnit)
      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = 'WestFace.tec', &
            FORM   = 'formatted', & 
            STATUS = 'REPLACE' )
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "nx", "ny", "nz" '
      WRITE(fUnit,*)  'ZONE T="el0", I=',N+1,', J=', N+1,',F=POINT'
      DO iY = 0, N
         DO iX = 0, N
            WRITE(fUnit,*)  myGeom % xBound(iX,iY,West), &
                            myGeom % yBound(iX,iY,West), &
                            myGeom % zBound(iX,iY,West), &
                            myGeom % nHat(1:3,iX,iY,West)
         ENDDO
      ENDDO
      CLOSE(UNIT=fUnit)
      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = 'EastFace.tec', &
            FORM   = 'formatted', & 
            STATUS = 'REPLACE' )
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "nx", "ny", "nz" '
      WRITE(fUnit,*)  'ZONE T="el0", I=',N+1,', J=', N+1,',F=POINT'
      DO iY = 0, N
         DO iX = 0, N
            WRITE(fUnit,*)  myGeom % xBound(iX,iY,East), &
                            myGeom % yBound(iX,iY,East), &
                            myGeom % zBound(iX,iY,East), &
                            myGeom % nHat(1:3,iX,iY,East)
         ENDDO
      ENDDO
      CLOSE(UNIT=fUnit)
      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = 'BottomFace.tec', &
            FORM   = 'formatted', & 
            STATUS = 'REPLACE' )
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "nx", "ny", "nz" '
      WRITE(fUnit,*)  'ZONE T="el0", I=',N+1,', J=', N+1,',F=POINT'
      DO iY = 0, N
         DO iX = 0, N
            WRITE(fUnit,*)  myGeom % xBound(iX,iY,Bottom), &
                            myGeom % yBound(iX,iY,Bottom), &
                            myGeom % zBound(iX,iY,Bottom), &
                            myGeom % nHat(1:3,iX,iY,Bottom)
         ENDDO
      ENDDO
      CLOSE(UNIT=fUnit)
      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = 'TopFace.tec', &
            FORM   = 'formatted', & 
            STATUS = 'REPLACE' )
      WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "nx", "ny", "nz" '
      WRITE(fUnit,*)  'ZONE T="el0", I=',N+1,', J=', N+1,',F=POINT'
      DO iY = 0, N
         DO iX = 0, N
            WRITE(fUnit,*)  myGeom % xBound(iX,iY,Top), &
                            myGeom % yBound(iX,iY,Top), &
                            myGeom % zBound(iX,iY,Top), &
                            myGeom % nHat(1:3,iX,iY,Top)
         ENDDO
      ENDDO
      CLOSE(UNIT=fUnit)

    RETURN

 END SUBROUTINE WriteTecplot_MappedGeometry_3D
!
!
! ///////////////////////////////////// PRIVATE ////////////////////////////////////////////////// !
 FUNCTION TransfiniteInterpolation( surfaces, a, b, c ) RESULT( P )
 ! TransfiniteInterpolation
 !  Takes in the six surfaces (south, east, north, west, bottom, top) and evaluates the 
 !  bidirectional mapping at xi^1 = a, xi^2 = b, xi^3 = c. The boundary of the computational 
 !  coordinate system is assumed to be at +/- 1 in each direction.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( Surface )  :: surfaces(1:6)
   REAL(prec)       :: a, b, c
   REAL(prec)       :: P(1:nDims)
   ! LOCAL
   REAL(prec)  :: P1(1:nDims), P2(1:nDims), P3(1:nDims)
   REAL(prec)  :: sSurf(1:nDims), nSurf(1:nDims), eSurf(1:nDims), wSurf(1:nDims), bSurf(1:nDims), tSurf(1:nDims)
   REAL(prec)  :: l1(1:2), l2(1:2), l3(1:2)
   REAL(prec)  :: ref(1:2)
   INTEGER     :: i, j
   
      ref = (/ -ONE, ONE /)

      ! Transfinite interpolation with linear blending uses linear lagrange interpolating polynomials
      ! to blend the bounding surfaces.
      ! The linear blending weights in the first computational direction are calculated.
      l1 = LinearBlend( a )
      l2 = LinearBlend( b )
      l3 = LinearBlend( c )
!      print*,a,b,c
      ! The bounding surfaces need to be evaluated at the provided computational coordinates
      wSurf = surfaces(west) % Evaluate( (/b, c/) )   ! west
      eSurf = surfaces(east) % Evaluate( (/b, c/) )   ! east
      sSurf = surfaces(south) % Evaluate( (/a, c/) )  ! south
      nSurf = surfaces(north) % Evaluate( (/a, c/) )  ! north
      bSurf = surfaces(bottom) % Evaluate( (/a, b/) ) ! bottom
      tSurf = surfaces(top) % Evaluate( (/a, b/) )    ! top
!      print*, 'xwest : ', wSurf
!      print*, 'xeast : ',  eSurf
!      print*, 'xsouth : ', sSurf
!      print*, 'xnorth : ', nSurf
!      print*, 'xbottom : ', bSurf
!      print*, 'xtop : ', tSurf
!      print*,'-----------------------'
      ! P1 contains the interpolation in the first computational coordinate
      ! The first computational coordinate is assumed to vary between the (computational) east and
      ! west boundaries.
      P1 = l1(1)*wSurf + l1(2)*eSurf
      ! P2 contains the interpolation in the second computational coordinate
      ! The second computational coordinate is assumed to vary between the (computational) south and
      ! north boundaries.
      P2 = l2(1)*sSurf + l2(2)*nSurf
      ! P3 contains the interpolation in the first computational coordinate
      ! The first computational coordinate is assumed to vary between the (computational) bottom and
      ! top boundaries.
      P3 = l3(1)*bSurf + l3(2)*tSurf
!      print*, 'P1 : ', P1
!      print*, 'P2 : ', P2
!      print*, 'P3 : ', P3
!      STOP
      DO i = 1, 2
         ! Now we need to compute the tensor product of the first and second computational direction 
         ! interpolants and subtract from P1.
         wSurf = surfaces(west) % Evaluate( (/ref(i), c/) )
         eSurf = surfaces(east) % Evaluate( (/ref(i), c/) )
         !CALL surfaces(west) % Evaluate( ref(i), c, wSurf(1), wSurf(2), wSurf(3) )
         !CALL surfaces(east) % Evaluate( ref(i), c, eSurf(1), eSurf(2), eSurf(3) )
         P1 = P1 - l2(i)*( wSurf*l1(1) + eSurf*l1(2) )
         
         ! Now we need to compute the tensor product of the first and third computational direction 
         ! interpolants and subtract from P1.
         wSurf = surfaces(west) % Evaluate( (/b, ref(i)/) )
         eSurf = surfaces(east) % Evaluate( (/b, ref(i)/) )
        ! CALL surfaces(west) % Evaluate( b, ref(i), wSurf(1), wSurf(2), wSurf(3) )
        ! CALL surfaces(east) % Evaluate( b, ref(i), eSurf(1), eSurf(2), eSurf(3) )

         P1 = P1 - l3(i)*( wSurf*l1(1) + eSurf*l1(2) )
      
         ! Now we need to compute the tensor product of the second and third computational direction 
         ! interpolants and subtract from P2.
         sSurf = surfaces(south) % Evaluate( (/a, ref(i)/) )
         nSurf = surfaces(north) % Evaluate( (/a, ref(i)/) )
        ! CALL surfaces(south) % Evaluate( a, ref(i), sSurf(1), sSurf(2), sSurf(3) )
        ! CALL surfaces(north) % Evaluate( a, ref(i), nSurf(1), nSurf(2), nSurf(3) )

         P2 = P2 - l3(i)*( sSurf*l2(1) + nSurf*l2(2) )
      
      ENDDO

      ! Next, the compounded tensor product is computed and added to P3.
      DO j = 1,2
         DO i = 1,2
         
            wSurf = surfaces(west) % Evaluate( (/ref(i), ref(j)/) )
            eSurf = surfaces(east) % Evaluate( (/ref(i), ref(j)/) )
            !CALL surfaces(west) % Evaluate( ref(i), ref(j), wSurf(1), wSurf(2), wSurf(3) )
            !CALL surfaces(east) % Evaluate( ref(i), ref(j), eSurf(1), eSurf(2), eSurf(3) )
            P3 = P3 + l2(i)*l3(j)*( wSurf*l1(1) + eSurf*l1(2) )
         
         ENDDO
      ENDDO
      
      !Finally, the sum the interpolants is computed to yield the computational coordinate
      P = P1 + P2 + P3
      
 END FUNCTION TransfiniteInterpolation
!
 FUNCTION TransfiniteInterpolation_Alt( interp, x, y, z, a, b, c ) RESULT( P )
 ! TransfiniteInterpolation
 !  Takes in the six surfaces (south, east, north, west, bottom, top) and evaluates the 
 !  bidirectional mapping at xi^1 = a, xi^2 = b, xi^3 = c. The boundary of the computational 
 !  coordinate system is assumed to be at +/- 1 in each direction.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
#ifdef HAVE_CUDA
   TYPE( Lagrange_Cuda ) :: interp
#else
   TYPE( Lagrange ) :: interp
#endif
   REAL(prec)       :: x(0:interp % N,0:interp % N,1:6), y(0:interp % N,0:interp % N,1:6), z(0:interp % N,0:interp % N,1:6)
   REAL(prec)       :: a, b, c
   REAL(prec)       :: P(1:nDims)
   ! LOCAL
   REAL(prec)  :: P1(1:nDims), P2(1:nDims), P3(1:nDims)
   REAL(prec)  :: sSurf(1:nDims), nSurf(1:nDims), eSurf(1:nDims), wSurf(1:nDims), bSurf(1:nDims), tSurf(1:nDims)
   REAL(prec)  :: l1(1:2), l2(1:2), l3(1:2)
   REAL(prec)  :: ref(1:2)
   INTEGER     :: i, j
   
      ref = (/ -ONE, ONE /)

      ! Transfinite interpolation with linear blending uses linear lagrange interpolating polynomials
      ! to blend the bounding surfaces.
      ! The linear blending weights in the first computational direction are calculated.
      l1 = LinearBlend( a )
      l2 = LinearBlend( b )
      l3 = LinearBlend( c )

      ! The bounding surfaces need to be evaluated at the provided computational coordinates
      wSurf(1) = interp % Interpolate_2D( x(:,:,west), (/b, c/) ) ! west
      wSurf(2) = interp % Interpolate_2D( y(:,:,west), (/b, c/) )
      wSurf(3) = interp % Interpolate_2D( z(:,:,west), (/b, c/) )
      
      eSurf(1) = interp % Interpolate_2D( x(:,:,east), (/b, c/) ) ! east
      eSurf(2) = interp % Interpolate_2D( y(:,:,east), (/b, c/) )
      eSurf(3) = interp % Interpolate_2D( z(:,:,east), (/b, c/) )

      sSurf(1) = interp % Interpolate_2D( x(:,:,south), (/a, c/) ) ! south
      sSurf(2) = interp % Interpolate_2D( y(:,:,south), (/a, c/) )
      sSurf(3) = interp % Interpolate_2D( z(:,:,south), (/a, c/) )

      nSurf(1) = interp % Interpolate_2D( x(:,:,north), (/a, c/) ) ! north
      nSurf(2) = interp % Interpolate_2D( y(:,:,north), (/a, c/) )
      nSurf(3) = interp % Interpolate_2D( z(:,:,north), (/a, c/) )
      
      bSurf(1) = interp % Interpolate_2D( x(:,:,bottom), (/a, b/) ) ! bottom
      bSurf(2) = interp % Interpolate_2D( y(:,:,bottom), (/a, b/) )
      bSurf(3) = interp % Interpolate_2D( z(:,:,bottom), (/a, b/) )
   
      tSurf(1) = interp % Interpolate_2D( x(:,:,top), (/a, b/) ) ! top
      tSurf(2) = interp % Interpolate_2D( y(:,:,top), (/a, b/) )
      tSurf(3) = interp % Interpolate_2D( z(:,:,top), (/a, b/) )
      
      ! P1 contains the interpolation in the first computational coordinate
      ! The first computational coordinate is assumed to vary between the (computational) east and
      ! west boundaries.
      P1 = l1(1)*wSurf + l1(2)*eSurf
      ! P2 contains the interpolation in the second computational coordinate
      ! The second computational coordinate is assumed to vary between the (computational) south and
      ! north boundaries.
      P2 = l2(1)*sSurf + l2(2)*nSurf
      ! P3 contains the interpolation in the first computational coordinate
      ! The first computational coordinate is assumed to vary between the (computational) bottom and
      ! top boundaries.
      P3 = l3(1)*bSurf + l3(2)*tSurf

      DO i = 1, 2
         ! Now we need to compute the tensor product of the first and second computational direction 
         ! interpolants and subtract from P1.
         wSurf(1) = interp % Interpolate_2D( x(:,:,west), (/ref(i), c/) ) ! west
         wSurf(2) = interp % Interpolate_2D( y(:,:,west), (/ref(i), c/) )
         wSurf(3) = interp % Interpolate_2D( z(:,:,west), (/ref(i), c/) )
         
         eSurf(1) = interp % Interpolate_2D( x(:,:,east), (/ref(i), c/) ) ! east
         eSurf(2) = interp % Interpolate_2D( y(:,:,east), (/ref(i), c/) )
         eSurf(3) = interp % Interpolate_2D( z(:,:,east), (/ref(i), c/) )

         P1 = P1 - l2(i)*( wSurf*l1(1) + eSurf*l1(2) )
         
         ! Now we need to compute the tensor product of the first and third computational direction 
         ! interpolants and subtract from P1.
         wSurf(1) = interp % Interpolate_2D( x(:,:,west), (/b, ref(i)/) ) ! west
         wSurf(2) = interp % Interpolate_2D( y(:,:,west), (/b, ref(i)/) )
         wSurf(3) = interp % Interpolate_2D( z(:,:,west), (/b, ref(i)/) )
         
         eSurf(1) = interp % Interpolate_2D( x(:,:,east), (/b, ref(i)/) ) ! east
         eSurf(2) = interp % Interpolate_2D( y(:,:,east), (/b, ref(i)/) )
         eSurf(3) = interp % Interpolate_2D( z(:,:,east), (/b, ref(i)/) )

         P1 = P1 - l3(i)*( wSurf*l1(1) + eSurf*l1(2) )
      
         ! Now we need to compute the tensor product of the second and third computational direction 
         ! interpolants and subtract from P2.
         sSurf(1) = interp % Interpolate_2D( x(:,:,south), (/a, ref(i)/) ) ! south
         sSurf(2) = interp % Interpolate_2D( y(:,:,south), (/a, ref(i)/) )
         sSurf(3) = interp % Interpolate_2D( z(:,:,south), (/a, ref(i)/) )
         
         nSurf(1) = interp % Interpolate_2D( x(:,:,north), (/a, ref(i)/) ) ! north
         nSurf(2) = interp % Interpolate_2D( y(:,:,north), (/a, ref(i)/) )
         nSurf(3) = interp % Interpolate_2D( z(:,:,north), (/a, ref(i)/) )

         P2 = P2 - l3(i)*( sSurf*l2(1) + nSurf*l2(2) )
      
      ENDDO

      ! Next, the compounded tensor product is computed and added to P3.
      DO j = 1,2
         DO i = 1,2
         
     !       wSurf = surfaces(west) % Evaluate( (/ref(i), ref(j)/) )
     !       eSurf = surfaces(east) % Evaluate( (/ref(i), ref(j)/) )
            
            wSurf(1) = interp % Interpolate_2D( x(:,:,west), (/ref(i), ref(j)/) ) ! west
            wSurf(2) = interp % Interpolate_2D( y(:,:,west), (/ref(i), ref(j)/) )
            wSurf(3) = interp % Interpolate_2D( z(:,:,west), (/ref(i), ref(j)/) )
         
            eSurf(1) = interp % Interpolate_2D( x(:,:,east), (/ref(i), ref(j)/) ) ! east
            eSurf(2) = interp % Interpolate_2D( y(:,:,east), (/ref(i), ref(j)/) )
            eSurf(3) = interp % Interpolate_2D( z(:,:,east), (/ref(i), ref(j)/) )
            P3 = P3 + l2(i)*l3(j)*( wSurf*l1(1) + eSurf*l1(2) )
         
         ENDDO
      ENDDO
      
      !Finally, the sum the interpolants is computed to yield the computational coordinate
      P = P1 + P2 + P3
      
 END FUNCTION TransfiniteInterpolation_Alt
!
 FUNCTION Unidirectional( valLeft, valRight, a ) RESULT( P )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   REAL(prec) :: valLeft(1:3), valRight(1:3)
   REAL(prec) :: a
   REAL(prec) :: P(1:3)

       P = HALF*( (ONE - a)*valLeft + (ONE + a)*valRight )
    
 END FUNCTION Unidirectional
 FUNCTION LinearBlend( a ) RESULT( weights )

   IMPLICIT NONE
   REAL(prec) :: a
   REAL(prec) :: weights(1:2)

       weights(1) = HALF*(ONE - a)
       weights(2) = HALF*(ONE + a)
    
 END FUNCTION LinearBlend
!
END MODULE MappedGeometry_3D_Class
