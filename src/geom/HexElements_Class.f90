! HexElement_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE HexElement_Class
 
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
USE Lagrange_Class
USE Surface_Class


IMPLICIT NONE


!  The HexElement data structure defines attributes needed to describe a quadrilateral element.
!
!  HexElements, elements, edges, and faces form the foundation of describing an unstructured mesh. 
!  The relationship betweens nodes and elements and nodes and edges (or faces in 3-D) define the 
!  connectivity in an unstructured mesh. In this data structure a HexElement is defined through an
!  integer ID, its four corner nodes, and its neighbors. 
!
!   --- COMBINE DOCUMENTATION --- 
!
!  The HexElements class provides attributes and type-bound procedures for handling 
!  curvilinear mappings between physical space and a reference computational space.
!
!  The HexElements class enables the implementation of spectral element methods in complicated
!  geometry. Given the four bounding surfaces of an element, the internal geometry (physical node 
!  positions, covariant basis vectors, and Jacobian) is constructed using transfinite interpolation
!  with linear blending.
!  
!

  TYPE HexElements
    INTEGER                 :: N, nElements
    INTEGER, ALLOCATABLE    :: elementID(:)
    INTEGER, ALLOCATABLE    :: nodeIDs(:,:)
    INTEGER, ALLOCATABLE    :: neighbors(:,:)
    REAL(prec), ALLOCATABLE :: nHat(:,:,:,:,:) 
    REAL(prec), ALLOCATABLE :: xBound(:,:,:,:)
    REAL(prec), ALLOCATABLE :: yBound(:,:,:,:) 
    REAL(prec), ALLOCATABLE :: zBound(:,:,:,:) 
    REAL(prec), ALLOCATABLE :: x(:,:,:,:), y(:,:,:,:), z(:,:,:,:)
    REAL(prec), ALLOCATABLE :: J(:,:,:,:)    
    REAL(prec), ALLOCATABLE :: dxds(:,:,:,:), dxdp(:,:,:,:), dxdq(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dyds(:,:,:,:), dydp(:,:,:,:), dydq(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dzds(:,:,:,:), dzdp(:,:,:,:), dzdq(:,:,:,:)
    REAL(prec), ALLOCATABLE :: Ja(:,:,:,:,:)
  
#ifdef HAVE_CUDA
    INTEGER, ALLOCATABLE    :: N_dev, nElements_dev
    INTEGER, ALLOCATABLE    :: elementID_dev(:)
    INTEGER, ALLOCATABLE    :: nodeIDs_dev(:,:)
    INTEGER, ALLOCATABLE    :: neighbors_dev(:,:)
    REAL(prec), ALLOCATABLE :: nHat_dev(:,:,:,:,:) 
    REAL(prec), ALLOCATABLE :: xBound_dev(:,:,:,:)
    REAL(prec), ALLOCATABLE :: yBound_dev(:,:,:,:) 
    REAL(prec), ALLOCATABLE :: zBound_dev(:,:,:,:) 
    REAL(prec), ALLOCATABLE :: x_dev(:,:,:,:), y_dev(:,:,:,:), z_dev(:,:,:,:)
    REAL(prec), ALLOCATABLE :: J_dev(:,:,:,:)    
    REAL(prec), ALLOCATABLE :: dxds_dev(:,:,:,:), dxdp_dev(:,:,:,:), dxdq_dev(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dyds_dev(:,:,:,:), dydp_dev(:,:,:,:), dydq_dev(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dzds_dev(:,:,:,:), dzdp_dev(:,:,:,:), dzdq_dev(:,:,:,:)
    REAL(prec), ALLOCATABLE :: Ja_dev(:,:,:,:,:)
#endif

    CONTAINS
      PROCEDURE :: Build      => Build_HexElements
      PROCEDURE :: Trash      => Trash_HexElements

      PROCEDURE :: GenerateMesh => GenerateMesh_HexElements
      PROCEDURE :: GenerateMetrics => GenerateMetrics_HexElements
      PROCEDURE :: GenerateBoundaryMetrics => GenerateBoundaryMetrics_HexElements
      PROCEDURE :: ScaleGeometry => ScaleGeometry_HexElements
      PROCEDURE :: CalculateLocation => CalculateLocation_HexElements
      PROCEDURE :: CalculateMetrics => CalculateMetrics_HexElements
      PROCEDURE :: CalculateComputationalCoordinates => CalculateComputationalCoordinates_HexElements
      
      PROCEDURE :: ResetInternalMesh => ResetInternalMesh_HexElements
      
      PROCEDURE :: WriteTecplot => WriteTecplot_HexElements
  
  END TYPE HexElements


 PRIVATE :: TransfiniteInterpolation, LinearBlend
 
 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup HexElements_Class
!! @{ 
! ================================================================================================ !
! S/R Initialize
! 
!> \fn Initialize_HexElements  
!! Allocates memory for each of the attributes of the HexElements Class and initializes all
!! arrays to zero.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElements) :: this <BR>
!! <B>INTEGER</B>                 :: N <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize( N ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myHex <td> HexElements <td> On output, an initialized HexElements
!!                                                         data structure
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the spectral element 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_HexElements( myHex, N, nElements )

  IMPLICIT NONE
  CLASS(HexElements), INTENT(out) :: myHex
  INTEGER, INTENT(in)             :: N, nElements

      myHex % N = N
      myHex % nElements = nElements
       
      ! Allocate space
      ALLOCATE( myHex % elementID(1:nElements), myHex % nodeIDs(1:8,1:nElements), myHex % neighbors(1:6,1:nElements) )
      ALLOCATE( myHex % dxds(0:N,0:N,0:N,1:nElements), myHex % dxdp(0:N,0:N,0:N,1:nElements), myHex % dxdq(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % dyds(0:N,0:N,0:N,1:nElements), myHex % dydp(0:N,0:N,0:N,1:nElements), myHex % dydq(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % dzds(0:N,0:N,0:N,1:nElements), myHex % dzdp(0:N,0:N,0:N,1:nElements), myHex % dzdq(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % Ja(0:N,0:N,0:N,1:3,1:3,1:nElements) )
      ALLOCATE( myHex % J(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % x(0:N,0:N,0:N,1:nElements), myHex % y(0:N,0:N,0:N,1:nElements), myHex % z(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % xBound(0:N,0:N,1:nHexFaces,1:nElements) )
      ALLOCATE( myHex % yBound(0:N,0:N,1:nHexFaces,1:nElements) )
      ALLOCATE( myHex % zBound(0:N,0:N,1:nHexFaces,1:nElements) )
      ALLOCATE( myHex % nHat(1:3,0:N,0:N,1:nHexFaces,1:nElements) )
      
      myHex % dxds   = 0.0_prec
      myHex % dxdp   = 0.0_prec
      myHex % dxdq   = 0.0_prec
      myHex % dyds   = 0.0_prec
      myHex % dydp   = 0.0_prec
      myHex % dydq   = 0.0_prec
      myHex % dzds   = 0.0_prec
      myHex % dzdp   = 0.0_prec
      myHex % dzdq   = 0.0_prec
      myHex % J      = 0.0_prec
      myHex % x      = 0.0_prec
      myHex % y      = 0.0_prec
      myHex % z      = 0.0_prec
      myHex % xBound = 0.0_prec
      myHex % yBound = 0.0_prec
      myHex % zBound = 0.0_prec

#ifdef HAVE_CUDA

      ALLOCATE( myHex % N_dev, myHex % nElements_dev )
      ALLOCATE( myHex % elementID_dev(1:nElements), myHex % nodeIDs_dev(1:8,1:nElements), myHex % neighbors_dev(1:6,1:nElements) )
      ALLOCATE( myHex % dxds_dev(0:N,0:N,0:N,1:nElements), myHex % dxdp_dev(0:N,0:N,0:N,1:nElements), myHex % dxdq_dev(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % dyds_dev(0:N,0:N,0:N,1:nElements), myHex % dydp_dev(0:N,0:N,0:N,1:nElements), myHex % dydq_dev(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % dzds_dev(0:N,0:N,0:N,1:nElements), myHex % dzdp_dev(0:N,0:N,0:N,1:nElements), myHex % dzdq_dev(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % Ja_dev(0:N,0:N,0:N,1:3,1:3,1:nElements) )
      ALLOCATE( myHex % J_dev(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % x_dev(0:N,0:N,0:N,1:nElements), myHex % y_dev(0:N,0:N,0:N,1:nElements), myHex % z_dev(0:N,0:N,0:N,1:nElements) )
      ALLOCATE( myHex % xBound_dev(0:N,0:N,1:nHexFaces,1:nElements) )
      ALLOCATE( myHex % yBound_dev(0:N,0:N,1:nHexFaces,1:nElements) )
      ALLOCATE( myHex % zBound_dev(0:N,0:N,1:nHexFaces,1:nElements) )
      ALLOCATE( myHex % nHat_dev(1:3,0:N,0:N,1:nHexFaces,1:nElements) )

      myHex % N_dev =  N
      myHex % nElements_dev = nElements

#endif

 END SUBROUTINE Build_HexElements
!
!> \addtogroup HexElements_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_HexElements_Class  
!! Frees memory held by the attributes of the HexElements data-structure.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElements) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myHex <td> HexElements <td>
!!                         On <B>input</B>, a previously constructed HexElements data structure <BR>
!!                         On <B>output</B>, memory held by attributes is freed
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_HexElements( myHex )

   IMPLICIT NONE
   CLASS(HexElements), INTENT(inout)  :: myHex

      DEALLOCATE( myHex % elementID, myHex % nodeIDs, myHex % neighbors )
      DEALLOCATE( myHex % dxds, myHex % dxdp, myHex % dxdq )
      DEALLOCATE( myHex % dyds, myHex % dydp, myHex % dydq )
      DEALLOCATE( myHex % dzds, myHex % dzdp, myHex % dzdq )
      DEALLOCATE( myHex % J, myHex % Ja, myHex % x, myHex % y, myHex % z )
      DEALLOCATE( myHex % xBound, myHex % yBound, myHex % zBound )
      DEALLOCATE( myHex % nHat )

#ifdef HAVE_CUDA
      
      DEALLOCATE( myHex % elementID_dev, myHex % nodeIDs_dev, myHex % neighbors_dev )
      DEALLOCATE( myHex % dxds_dev, myHex % dxdp_dev, myHex % dxdq_dev )
      DEALLOCATE( myHex % dyds_dev, myHex % dydp_dev, myHex % dydq_dev )
      DEALLOCATE( myHex % dzds_dev, myHex % dzdp_dev, myHex % dzdq_dev )
      DEALLOCATE( myHex % J_dev, myHex % Ja_dev, myHex % x_dev, myHex % y_dev, myHex % z_dev )
      DEALLOCATE( myHex % xBound_dev, myHex % yBound_dev, myHex % zBound_dev )
      DEALLOCATE( myHex % nHat_dev )

#endif
 
 END SUBROUTINE Trash_HexElements


SUBROUTINE UpdateDevice_HexElements

END SUBROUTINE UpdateDevice_HexElements

SUBROUTINE UpdateHost_HexElements

END SUBROUTINE UpdateHost_HexElements

!
!> \addtogroup HexElements_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMesh 
! 
!> \fn GenerateMesh_HexElements  
!! Generates the physical interior and boundary positions at each of the computational mesh points.
!! 
!! Given the four boundary surfaces and an interpolant that stores the computational mesh points,
!! the physical mesh positions are generated using transfinite interpolation with linear blending.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref HexElements_Class : PRIVATE Function TransfiniteInterpolation <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElements) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>TYPE</B>(Surface)           :: boundaries(1:6) <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateMesh( interp, boundaries ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myHex <td> HexElements <td> 
!!                         On <B>input</B>, an initialized HexElements data structure, <BR>
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
 SUBROUTINE GenerateMesh_HexElements( myHex, interp, theSurfaces )

   IMPLICIT NONE
   CLASS( HexElements ), INTENT(inout) :: myHex
   TYPE( Lagrange ), INTENT(in)              :: interp
   TYPE( Surface ), INTENT(in)               :: theSurfaces(1:6)
   ! Local
   INTEGER    :: i, j, k, N
   REAL(prec) :: s(0:interp % N), p, x(1:3)
   
      N = interp % N
      s = interp % s
      
      DO k = 0, N
         DO j = 0,N
            DO i = 0,N
               x = TransfiniteInterpolation( theSurfaces, s(i), s(j), s(k) )
               myHex % x(i,j,k) = x(1)
               myHex % y(i,j,k) = x(2)
               myHex % z(i,j,k) = x(3)
            ENDDO
         ENDDO 
      ENDDO 

      ! Do the boundary locations
      DO j = 0, N
         DO i = 0, N
            p = -ONE  ! south boundary
            x = TransfiniteInterpolation( theSurfaces, s(i), p, s(j) )
            myHex % xBound(i,j,south) = x(1)
            myHex % yBound(i,j,south) = x(2)
            myHex % zBound(i,j,south) = x(3)
            ! west boundary
            x = TransfiniteInterpolation( theSurfaces, p, s(i), s(j) )
            myHex % xBound(i,j,west) = x(1)
            myHex % yBound(i,j,west) = x(2)
            myHex % zBound(i,j,west) = x(3)
            ! bottom boundary
            x = TransfiniteInterpolation( theSurfaces, s(i), s(j), p )
            myHex % xBound(i,j,bottom) = x(1)
            myHex % yBound(i,j,bottom) = x(2)
            myHex % zBound(i,j,bottom) = x(3)
            
            p = ONE  ! north boundary
            x = TransfiniteInterpolation( theSurfaces, s(i), p, s(j) )
            myHex % xBound(i,j,north) = x(1)
            myHex % yBound(i,j,north) = x(2)
            myHex % zBound(i,j,north) = x(3)
            ! east boundary
            x = TransfiniteInterpolation( theSurfaces, p, s(i), s(j) )
            myHex % xBound(i,j,east) = x(1)
            myHex % yBound(i,j,east) = x(2)
            myHex % zBound(i,j,east) = x(3)
            ! top boundary
            x = TransfiniteInterpolation( theSurfaces, s(i), s(j), p )
            myHex % xBound(i,j,top) = x(1)
            myHex % yBound(i,j,top) = x(2)
            myHex % zBound(i,j,top) = x(3)
         ENDDO
      ENDDO

 END SUBROUTINE GenerateMesh_HexElements
!
 SUBROUTINE ResetInternalMesh_HexElements( myHex, interp )

   IMPLICIT NONE
   CLASS( HexElements ), INTENT(inout) :: myHex
   TYPE( Lagrange ), INTENT(in)              :: interp
   ! Local
   INTEGER    :: i, j, k
   REAL(prec) :: s(0:interp % N), p, x(1:3)
   
      
      DO k = 0, interp % N
         DO j = 0, interp % N
            DO i = 0, interp % N
               x = TransfiniteInterpolation_Alt( interp, &
                                                 myHex % xBound, &
                                                 myHex % yBound, &
                                                 myHex % zBound, &
                                                 interp % s(i), &
                                                 interp % s(j), &
                                                 interp % s(k) )
               myHex % x(i,j,k) = x(1)
               myHex % y(i,j,k) = x(2)
               myHex % z(i,j,k) = x(3)
            ENDDO
         ENDDO 
      ENDDO 



 END SUBROUTINE ResetInternalMesh_HexElements
!
!> \addtogroup HexElements_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMetrics 
! 
!> \fn GenerateMetrics_HexElements  
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
!!   Module \ref HexElements_Class : S/R \ref GenerateBoundaryMetrics_HexElements <BR>
!! 
!!  ** To produce meaningful output, GenerateMesh_HexElements must be called prior to calling
!!     this routine.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElements) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateMetrics( interp ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myHex <td> HexElements <td> 
!!                         On <B>input</B>, an initialized HexElements data structure, <BR>
!!                         On <B>output</B>, the metric terms are filled in, including the
!!                         outward pointing normal vectors on the element boundaries
!!   <tr> <td> in <th> interp <td> Lagrange <td>  Lagrange interpolant that defines the 
!!                                                  reference computational mesh 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GenerateMetrics_HexElements( myHex, interp )

   IMPLICIT NONE
   CLASS( HexElements ), INTENT(inout) :: myHex
   TYPE( Lagrange ), INTENT(in)              :: interp
   ! Local
   INTEGER    :: i, j, k, N
   REAL(prec) :: cv(1:3,1:3)
   REAL(prec) :: covT(0:interp % N, 0:interp % N, 0:interp % N, 1:3,1:3)
   REAL(prec) :: v(0:interp % N, 0:interp % N, 0:interp % N,1:3)
   REAL(prec) :: Dv1(0:interp % N, 0:interp % N, 0:interp % N,1:3)
   REAL(prec) :: Dv2(0:interp % N, 0:interp % N, 0:interp % N,1:3)
   REAL(prec) :: Dv3(0:interp % N, 0:interp % N, 0:interp % N,1:3)
   
      N = interp % N
      
      covT(0:N,0:N,0:N,1,1:3) = interp % ApplyDerivativeMatrix_3D( myHex % x )
      covT(0:N,0:N,0:N,2,1:3) = interp % ApplyDerivativeMatrix_3D( myHex % y )
      covT(0:N,0:N,0:N,3,1:3) = interp % ApplyDerivativeMatrix_3D( myHex % z )
      
      DO k = 0, N
         DO j = 0, N
            DO i = 0, N

               myHex % dxds(i,j,k) = covT(i,j,k,1,1)
               myHex % dxdp(i,j,k) = covT(i,j,k,1,2)
               myHex % dxdq(i,j,k) = covT(i,j,k,1,3)
               myHex % dyds(i,j,k) = covT(i,j,k,2,1)
               myHex % dydp(i,j,k) = covT(i,j,k,2,2)
               myHex % dydq(i,j,k) = covT(i,j,k,2,3)
               myHex % dzds(i,j,k) = covT(i,j,k,3,1)
               myHex % dzdp(i,j,k) = covT(i,j,k,3,2)
               myHex % dzdq(i,j,k) = covT(i,j,k,3,3)
               
               cv = covT(i,j,k,1:3,1:3)
               myHex % J(i,j,k) = Determinant( cv, 3 )
               
            ENDDO
         ENDDO
      ENDDO 
      
      ! Generate the contravariant metric tensor a la Kopriva (2006)
      !Ja_1
      DO k = 0, N
         DO j = 0, N
            DO i = 0, N
               v(i,j,k,1:3)  = -myHex % z(i,j,k)*covT(i,j,k,2,1:3) 
            ENDDO
         ENDDO
      ENDDO
      Dv1 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,1) ) ! ( dv1/ds, dv1/dp, dv1/dq )
      Dv2 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,2) ) ! ( dv2/ds, dv2/dp, dv2/dq )
      Dv3 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,3) ) ! ( dv3/ds, dv3/dp, dv3/dq )
      ! Take the curl to obtain the first dimension of each of the contravariant basis vectors
      ! The contravariant metric tensor stores each contravariant basis vector in each column
      ! of the tensor
      myHex % Ja(0:N,0:N,0:N,1,1) = ( Dv3(0:N,0:N,0:N,2) - Dv2(0:N,0:N,0:N,3) )
      myHex % Ja(0:N,0:N,0:N,1,2) = -( Dv3(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,3) )
      myHex % Ja(0:N,0:N,0:N,1,3) = ( Dv2(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,2) )
      !Ja_2
      DO k = 0, N
         DO j = 0, N
            DO i = 0, N
               v(i,j,k,1:3)  = -myHex % x(i,j,k)*covT(i,j,k,3,1:3) 
            ENDDO
         ENDDO
      ENDDO
      Dv1 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,1) ) ! ( dv1/ds, dv1/dp, dv1/dq )
      Dv2 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,2) ) ! ( dv2/ds, dv2/dp, dv2/dq )
      Dv3 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,3) ) ! ( dv3/ds, dv3/dp, dv3/dq )
      ! Take the curl to obtain the first dimension of each of the contravariant basis vectors
      ! The contravariant metric tensor stores each contravariant basis vector in each column
      ! of the tensor
      myHex % Ja(0:N,0:N,0:N,2,1) = ( Dv3(0:N,0:N,0:N,2) - Dv2(0:N,0:N,0:N,3) )
      myHex % Ja(0:N,0:N,0:N,2,2) = -( Dv3(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,3) )
      myHex % Ja(0:N,0:N,0:N,2,3) = ( Dv2(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,2) )
      !Ja_3
      DO k = 0, N
         DO j = 0, N
            DO i = 0, N
               v(i,j,k,1:3)  = -myHex % y(i,j,k)*covT(i,j,k,1,1:3) 
            ENDDO
         ENDDO
      ENDDO
      Dv1 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,1) ) ! ( dv1/ds, dv1/dp, dv1/dq )
      Dv2 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,2) ) ! ( dv2/ds, dv2/dp, dv2/dq )
      Dv3 = interp % ApplyDerivativeMatrix_3D( v(0:N,0:N,0:N,3) ) ! ( dv3/ds, dv3/dp, dv3/dq )
      ! Take the curl to obtain the first dimension of each of the contravariant basis vectors
      ! The contravariant metric tensor stores each contravariant basis vector in each column
      ! of the tensor
      myHex % Ja(0:N,0:N,0:N,3,1) = ( Dv3(0:N,0:N,0:N,2) - Dv2(0:N,0:N,0:N,3) )
      myHex % Ja(0:N,0:N,0:N,3,2) = -( Dv3(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,3) )
      myHex % Ja(0:N,0:N,0:N,3,3) = ( Dv2(0:N,0:N,0:N,1) - Dv1(0:N,0:N,0:N,2) )
      
      CALL myHex % GenerateBoundaryMetrics( interp )

 END SUBROUTINE GenerateMetrics_HexElements
!
!> \addtogroup HexElements_Class 
!! @{ 
! ================================================================================================ !
! S/R GenerateMetrics 
! 
!> \fn GenerateMetrics_HexElements  
!! Generates and stores the outward pointing boundary normal vectors.
!!
!!  The outward pointing boundary normal vectors are equivalent to the contravariant basis vectors
!!  evaluated at the element boundaries. These are computed here by differentiating the Lagrange
!!  interpolant of the mesh positions and the computational mesh boundaries.
!! 
!!  This subroutine depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref Differentiate_3D_Lagrange <BR>
!! 
!!  ** To produce meaningful output, GenerateMesh_HexElements must be called prior to calling
!!     this routine.
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElements) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!!         .... <BR>
!!     <B>CALL</B> this % GenerateBoundaryMetrics( interp ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myHex <td> HexElements <td> 
!!                         On <B>input</B>, an initialized HexElements data structure, <BR>
!!                         On <B>output</B>, the outward pointing normal vectors on the element
!!                         boundaries are filled in
!!   <tr> <td> in <th> interp <td> Lagrange <td> Lagrange interpolant that defines the 
!!                                                  reference computational mesh 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
  SUBROUTINE GenerateBoundaryMetrics_HexElements( myHex, interp  )

   IMPLICIT NONE
   CLASS( HexElements ), INTENT(inout) :: myHex
   TYPE( Lagrange ), INTENT(in)              :: interp
   ! Local
   INTEGER    :: i, j, N
   REAL(prec) :: s(0:interp % N), p
   REAL(prec) :: Jain(0:interp % N,0:interp % N,0:interp % N)
   REAL(prec) :: J, signJ, nx, ny, nz
   REAL(prec) :: node(1:3)
   
      N = interp % N
      s = interp % s
           ! Do the boundary locations
      DO j = 0, N
         DO i = 0, N
         
            p = -ONE  ! bottom, south, and west boundaries
            
            !bottom boundary
            node = (/s(i), s(j), p /)
            J = interp % Interpolate_3D( myHex % J, node ) !Determinant( cv, 3 )
            signJ = SIGN(ONE,J)                    
            ! Setting bottom boundary normal
            Jain = myHex % Ja(0:N,0:N,0:N,1,3)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,2,3)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,3,3)
            nz = interp % Interpolate_3D( Jain, node )
            myHex % nHat(1:3,i,j,bottom) = -signJ*(/ nx, ny, nz /)
            
            node = (/ s(i), p, s(j) /)
            J = interp % Interpolate_3D( myHex % J, node ) !Determinant( cv, 3 )
            signJ = SIGN(ONE,J)                    
            ! Setting southern boundary normal
            Jain = myHex % Ja(0:N,0:N,0:N,1,2)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,2,2)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,3,2)
            nz = interp % Interpolate_3D( Jain, node )
            myHex % nHat(1:3,i,j,south)= -signJ*(/ nx, ny, nz /)
            
            ! west boundary
            node = (/ p, s(i), s(j) /)
            J = interp % Interpolate_3D( myHex % J, node ) !Determinant( cv, 3 )
            signJ = SIGN(ONE,J)                    
            ! Setting western boundary normal
            Jain = myHex % Ja(0:N,0:N,0:N,1,1)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,2,1)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,3,1)
            nz = interp % Interpolate_3D( Jain, node )
            myHex % nHat(1:3,i,j,west) = -signJ*(/ nx, ny, nz /)
             
            p = ONE  ! top, north, and east boundaries
            
            !top boundary
            node = (/s(i), s(j), p /)
            J = interp % Interpolate_3D( myHex % J, node )!Determinant( cv, 3 )
            signJ = SIGN(ONE,J)      
            ! Setting top boundary normal
            Jain = myHex % Ja(0:N,0:N,0:N,1,3)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,2,3)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,3,3)
            nz = interp % Interpolate_3D( Jain, node )
            myHex % nHat(1:3,i,j,top) = signJ*(/ nx, ny, nz /)
            
            !north boundary
            node = (/ s(i), p, s(j) /)
            J = interp % Interpolate_3D( myHex % J, node ) !Determinant( cv, 3 )
            signJ = SIGN(ONE,J)    
            ! Setting southern boundary normal
            Jain = myHex % Ja(0:N,0:N,0:N,1,2)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,2,2)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,3,2)
            nz = interp % Interpolate_3D( Jain, node )
            myHex % nHat(1:3,i,j,north) = signJ*(/ nx, ny, nz /)
            
            ! east boundary
            node = (/ p, s(i), s(j) /)
            J = interp % Interpolate_3D( myHex % J, node ) !Determinant( cv, 3 )
            signJ = SIGN(ONE,J)                    
            ! Setting eastern boundary normal
            Jain = myHex % Ja(0:N,0:N,0:N,1,1)
            nx = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,2,1)
            ny = interp % Interpolate_3D( Jain, node )
            Jain = myHex % Ja(0:N,0:N,0:N,3,1)
            nz = interp % Interpolate_3D( Jain, node )
            myHex % nHat(1:3,i,j,east) = signJ*(/ nx, ny, nz /)
            
         ENDDO
      ENDDO

 END SUBROUTINE GenerateBoundaryMetrics_HexElements
!
!> \addtogroup HexElements_Class
!! @{ 
! ================================================================================================ !
! Function CalculateLocation 
! 
!> \fn CalculateLocation_HexElements  
!! Given a computational coordinate, the physical coordinate is calculated using Lagrange 
!! interpolation.
!! 
!!  This function depends on <BR>
!!   Module \ref Lagrange_Class : Function \ref Interpolate_3D_Lagrange <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElements) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: s(1:2), x(1:2) <BR>
!!         .... <BR>
!!     x = this % CalculateLocation( interp, s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myHex <td> HexElements <td> An intialized and constructed 
!!                                                        HexElements data structure
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
 FUNCTION CalculateLocation_HexElements( myHex, interp, s ) RESULT( x )

   IMPLICIT NONE
   CLASS( HexElements ) :: myHex
   TYPE( Lagrange )           :: interp
   REAL(prec)                 :: s(1:3)
   REAL(prec)                 :: x(1:3)
  
      x(1) = interp % Interpolate_3D( myHex % x, s )
      x(2) = interp % Interpolate_3D( myHex % y, s )
      x(3) = interp % Interpolate_3D( myHex % z, s )
  
 END FUNCTION CalculateLocation_HexElements
!
!> \addtogroup HexElements_Class
!! @{ 
! ================================================================================================ !
! Function CalculateMetrics 
! 
!> \fn CalculateMetrics_HexElements  
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
!! <B>TYPE</B>(HexElements) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: s(1:3), covT(1:3,1:3) <BR>
!!         .... <BR>
!!     x = this % CalculateMetrics( interp, s ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myHex <td> HexElements <td> An intialized and constructed 
!!                                                        HexElements data structure
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
 FUNCTION CalculateMetrics_HexElements( myHex, interp, s ) RESULT( covT )

   IMPLICIT NONE
   CLASS( HexElements ) :: myHex
   TYPE( Lagrange )           :: interp
   REAL(prec)                 :: s(1:3)
   REAL(prec)                 :: covT(1:3,1:3)
 
      covT(1,1:3) = interp % Differentiate_3D( myHex % x, s )
      covT(2,1:3) = interp % Differentiate_3D( myHex % y, s )
      covT(3,1:3) = interp % Differentiate_3D( myHex % z, s )

 END FUNCTION CalculateMetrics_HexElements
!
!> \addtogroup HexElements_Class 
!! @{ 
! ================================================================================================ !
! S/R CalculateComputationalCoordinates
! 
!> \fn CalculateComputationalCoordinates_HexElements 
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
!!   Module \ref HexElements_Class : Function \ref CalculateLocation_HexElements <BR>
!!   Module \ref HexElements_Class : Function \ref CalculateMetrics_HexElements <BR>
!!   Module \ref CommonRoutines          : Function \ref Invert_2x2 <BR>
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElements) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: x(1:3), s(1:3) <BR>
!! <B>LOGICAL</B>                 :: successful <BR>
!!         .... <BR>
!!     <B>CALL</B> this % CalculateComputationalCoordinates( interp, x, s, successful ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myHex <td> HexElements <td> An intialized and constructed 
!!                                                        HexElements data structure
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
 SUBROUTINE CalculateComputationalCoordinates_HexElements( myHex, interp, x, s, success )

   IMPLICIT NONE
   CLASS( HexElements ), INTENT(in) :: myHex
   TYPE( Lagrange ), INTENT(in)           :: interp
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
         thisX = myHex % CalculateLocation( interp, thisS )
     
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
        
         A = myHex % CalculateMetrics( interp, thisS ) ! Calculate the covariant metric tensor
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
      
     
 END SUBROUTINE CalculateComputationalCoordinates_HexElements
!
!> \addtogroup HexElements_Class
!! @{ 
! ================================================================================================ !
! S/R ScaleGeometry
! 
!> \fn ScaleGeometry_HexElements 
!! Scales the physical coordinates and metric terms by a given x-scale and y-scale.
!! 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HexElements) :: this <BR>
!! <B>TYPE</B>(Lagrange)          :: interp <BR>
!! <B>REAL</B>(prec)              :: xScale, yScale <BR>
!!         .... <BR>
!!     <B>CALL</B> this % ScaleGeometry( interp, xScale, yScale ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myHex <td> HexElements <td> 
!!                         On <B>input</B>, an intialized and constructed HexElements data 
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
 SUBROUTINE ScaleGeometry_HexElements( myHex, interp, xScale, yScale, zScale )

   IMPLICIT NONE
   CLASS( HexElements ), INTENT(inout) :: myHex
   TYPE( Lagrange ), INTENT(in)              :: interp
   REAL(prec), INTENT(in)                    :: xScale, yScale, zScale
   
         myHex % x = xScale*( myHex % x )
         myHex % y = yScale*( myHex % y )
         myHex % z = zScale*( myHex % z )
         myHex % xBound = xScale*( myHex % xBound )
         myHex % yBound = yScale*( myHex % yBound )
         myHex % zBound = zScale*( myHex % zBound )

         ! Update the boundary metrics -- normals and normal lengths
         CALL myHex % GenerateMetrics( interp )
         CALL myHex % GenerateBoundaryMetrics( interp  )
         
 END SUBROUTINE ScaleGeometry_HexElements
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
   REAL(prec)       :: P(1:3)
   ! LOCAL
   REAL(prec)  :: P1(1:3), P2(1:3), P3(1:3)
   REAL(prec)  :: sSurf(1:3), nSurf(1:3), eSurf(1:3), wSurf(1:3), bSurf(1:3), tSurf(1:3)
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
   REAL(prec)       :: P(1:3)
   ! LOCAL
   REAL(prec)  :: P1(1:3), P2(1:3), P3(1:3)
   REAL(prec)  :: sSurf(1:3), nSurf(1:3), eSurf(1:3), wSurf(1:3), bSurf(1:3), tSurf(1:3)
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
END MODULE HexElements_Class
