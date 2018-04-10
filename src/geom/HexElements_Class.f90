! HexElements_CLASS.f90
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE HexElements_CLASS

  USE ModelPrecision
  USE ConstantsDictionary
  USE CommonRoutines
  USE Lagrange_CLASS
  USE Surfaces_CLASS


  IMPLICIT NONE


!  The HexElement DATA structure defines attributes needed to describe a quadrilateral element.
!
!  HexElements, elements, edges, and faces form the foundation of describing an unstructured mesh.
!  The relationship betweens nodes and elements and nodes and edges (or faces in 3-D) define the
!  connectivity in an unstructured mesh. In this DATA structure a HexElement is defined through an
!  INTEGER ID, its four corner nodes, and its neighbors.
!
!   --- COMBINE DOCUMENTATION ---
!
!  The HexElements CLASS provides attributes and TYPE-bound procedures for handling
!  curvilinear mappings between physical space and a reference computational space.
!
!  The HexElements CLASS enables the implementation of spectral element methods in complicated
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
    REAL(prec), ALLOCATABLE :: xBound(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: x(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: J(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dxds(:,:,:,:), dxdp(:,:,:,:), dxdq(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dyds(:,:,:,:), dydp(:,:,:,:), dydq(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dzds(:,:,:,:), dzdp(:,:,:,:), dzdq(:,:,:,:)
    REAL(prec), ALLOCATABLE :: Ja(:,:,:,:,:,:)

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE    :: N_dev, nElements_dev
    INTEGER, DEVICE, ALLOCATABLE    :: elementID_dev(:)
    INTEGER, DEVICE, ALLOCATABLE    :: nodeIDs_dev(:,:)
    INTEGER, DEVICE, ALLOCATABLE    :: neighbors_dev(:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: nHat_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: xBound_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: x_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: J_dev(:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: dxds_dev(:,:,:,:), dxdp_dev(:,:,:,:), dxdq_dev(:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: dyds_dev(:,:,:,:), dydp_dev(:,:,:,:), dydq_dev(:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: dzds_dev(:,:,:,:), dzdp_dev(:,:,:,:), dzdq_dev(:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: Ja_dev(:,:,:,:,:,:)
#endif

  CONTAINS

    PROCEDURE :: Build      => Build_HexElements
    PROCEDURE :: Trash      => Trash_HexElements

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateDevice => UpdateDevice_HexElements
    PROCEDURE :: UpdateHost   => UpdateHost_HexElements
#endif

    PROCEDURE :: GenerateMesh => GenerateMesh_HexElements
    PROCEDURE :: GenerateMetrics => GenerateMetrics_HexElements
    PROCEDURE :: GenerateBoundaryMetrics => GenerateBoundaryMetrics_HexElements
    PROCEDURE :: ScaleGeometry => ScaleGeometry_HexElements
    PROCEDURE :: ResetInternalMesh => ResetInternalMesh_HexElements

    PROCEDURE :: CalculateComputationalCoordinates

  END TYPE HexElements


  PRIVATE :: TransfiniteInterpolation, LinearBlend

CONTAINS

!
!  Build_HexElements
! Allocates memory for each of the attributes of the HexElements CLASS and initializes all
! arrays to zero.
!
! <H2> Usage : </H2>
! <B>TYPE</B>(HexElements) :: this <BR>
! <B>INTEGER</B>                 :: N <BR>
!         .... <BR>
!     <B>CALL</B> this % Build( N ) <BR>
!
!  <H2> Parameters : </H2>
!  <table>
!   <tr> <td> out <th> myElements <td> HexElements <td> On output, an initialized HexElements
!                                                         DATA structure
!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree of the spectral element
!  </table>
!
! ================================================================================================ !

  SUBROUTINE Build_HexElements( myElements, N, nElements )

    IMPLICIT NONE
    CLASS(HexElements), INTENT(out) :: myElements
    INTEGER, INTENT(in)             :: N, nElements

    myElements % N = N
    myElements % nElements = nElements

    ALLOCATE( myElements % elementID(1:nElements), &
      myElements % nodeIDs(1:8,1:nElements), &
      myElements % neighbors(1:6,1:nElements), &
      myElements % dxds(0:N,0:N,0:N,1:nElements), &
      myElements % dxdp(0:N,0:N,0:N,1:nElements), &
      myElements % dxdq(0:N,0:N,0:N,1:nElements), &
      myElements % dyds(0:N,0:N,0:N,1:nElements), &
      myElements % dydp(0:N,0:N,0:N,1:nElements), &
      myElements % dydq(0:N,0:N,0:N,1:nElements), &
      myElements % dzds(0:N,0:N,0:N,1:nElements), &
      myElements % dzdp(0:N,0:N,0:N,1:nElements), &
      myElements % dzdq(0:N,0:N,0:N,1:nElements), &
      myElements % Ja(0:N,0:N,0:N,1:3,1:3,1:nElements), &
      myElements % J(0:N,0:N,0:N,1:nElements), &
      myElements % x(0:N,0:N,0:N,1:3,1:nElements), &
      myElements % xBound(0:N,0:N,1:3,1:6,1:nElements), &
      myElements % nHat(1:3,0:N,0:N,1:6,1:nElements) )

    myElements % dxds   = 0.0_prec
    myElements % dxdp   = 0.0_prec
    myElements % dxdq   = 0.0_prec
    myElements % dyds   = 0.0_prec
    myElements % dydp   = 0.0_prec
    myElements % dydq   = 0.0_prec
    myElements % dzds   = 0.0_prec
    myElements % dzdp   = 0.0_prec
    myElements % dzdq   = 0.0_prec
    myElements % J      = 0.0_prec
    myElements % Ja     = 0.0_prec
    myElements % x      = 0.0_prec
    myElements % xBound = 0.0_prec

#ifdef HAVE_CUDA

    ALLOCATE( myElements % N_dev, myElements % nElements_dev )

    ALLOCATE( myElements % elementID_dev(1:nElements), &
      myElements % nodeIDs_dev(1:8,1:nElements), &
      myElements % neighbors_dev(1:6,1:nElements), &
      myElements % dxds_dev(0:N,0:N,0:N,1:nElements), &
      myElements % dxdp_dev(0:N,0:N,0:N,1:nElements), &
      myElements % dxdq_dev(0:N,0:N,0:N,1:nElements), &
      myElements % dyds_dev(0:N,0:N,0:N,1:nElements), &
      myElements % dydp_dev(0:N,0:N,0:N,1:nElements), &
      myElements % dydq_dev(0:N,0:N,0:N,1:nElements), &
      myElements % dzds_dev(0:N,0:N,0:N,1:nElements), &
      myElements % dzdp_dev(0:N,0:N,0:N,1:nElements), &
      myElements % dzdq_dev(0:N,0:N,0:N,1:nElements), &
      myElements % Ja_dev(0:N,0:N,0:N,1:3,1:3,1:nElements), &
      myElements % J_dev(0:N,0:N,0:N,1:nElements), &
      myElements % x_dev(0:N,0:N,0:N,1:3,1:nElements), &
      myElements % xBound_dev(0:N,0:N,1:3,1:6,1:nElements), &
      myElements % nHat_dev(1:3,0:N,0:N,1:6,1:nElements) )

    myElements % N_dev =  N
    myElements % nElements_dev = nElements

    myElements % dxds_dev   = 0.0_prec
    myElements % dxdp_dev   = 0.0_prec
    myElements % dxdq_dev   = 0.0_prec
    myElements % dyds_dev   = 0.0_prec
    myElements % dydp_dev   = 0.0_prec
    myElements % dydq_dev   = 0.0_prec
    myElements % dzds_dev   = 0.0_prec
    myElements % dzdp_dev   = 0.0_prec
    myElements % dzdq_dev   = 0.0_prec
    myElements % J_dev      = 0.0_prec
    myElements % Ja_dev     = 0.0_prec
    myElements % x_dev      = 0.0_prec
    myElements % xBound_dev = 0.0_prec

#endif

  END SUBROUTINE Build_HexElements
!
!> \addtogroup HexElements_CLASS
!! @{
! ================================================================================================ !
! S/R Trash
!
!> \fn Trash_HexElements_CLASS
!! Frees memory held by the attributes of the HexElements DATA-structure.
!!
!! <H2> Usage : </H2>
!! <B>TYPE</B>(HexElements) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!!
!!  <H2> Parameters : </H2>
!!  <table>
!!   <tr> <td> in/out <th> myElements <td> HexElements <td>
!!                         On <B>input</B>, a previously constructed HexElements DATA structure <BR>
!!                         On <B>output</B>, memory held by attributes is freed
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE Trash_HexElements( myElements )

    IMPLICIT NONE
    CLASS(HexElements), INTENT(inout)  :: myElements

    DEALLOCATE( myElements % elementID, &
      myElements % nodeIDs, &
      myElements % neighbors, &
      myElements % dxds, &
      myElements % dxdp, &
      myElements % dxdq, &
      myElements % dyds, &
      myElements % dydp, &
      myElements % dydq, &
      myElements % dzds, &
      myElements % dzdp, &
      myElements % dzdq, &
      myElements % Ja, &
      myElements % J, &
      myElements % x, &
      myElements % xBound, &
      myElements % nHat )

#ifdef HAVE_CUDA

    ALLOCATE( myElements % N_dev, myElements % nElements_dev )

    DEALLOCATE( myElements % elementID_dev, &
      myElements % nodeIDs_dev, &
      myElements % neighbors_dev, &
      myElements % dxds_dev, &
      myElements % dxdp_dev, &
      myElements % dxdq_dev, &
      myElements % dyds_dev, &
      myElements % dydp_dev, &
      myElements % dydq_dev, &
      myElements % dzds_dev, &
      myElements % dzdp_dev, &
      myElements % dzdq_dev, &
      myElements % Ja_dev, &
      myElements % J_dev, &
      myElements % x_dev, &
      myElements % xBound_dev, &
      myElements % nHat_dev )

#endif



  END SUBROUTINE Trash_HexElements

#ifdef HAVE_CUDA
  SUBROUTINE UpdateDevice_HexElements( myElements )
    IMPLICIT NONE
    CLASS( HexElements ), INTENT(inout) :: myElements

    myElements % elementID_dev = myElements % elementID
    myElements % nodeIDs_dev   = myElements % nodeIDs
    myElements % neighbors_dev = myElements % neighbors
    myElements % dxds_dev = myElements % dxds
    myElements % dxdp_dev = myElements % dxdp
    myElements % dxdq_dev = myElements % dxdq
    myElements % dyds_dev = myElements % dyds
    myElements % dydp_dev = myElements % dydp
    myElements % dydq_dev = myElements % dydq
    myElements % dzds_dev = myElements % dzds
    myElements % dzdp_dev = myElements % dzdp
    myElements % dzdq_dev = myElements % dzdq
    myElements % Ja_dev = myElements % Ja
    myElements % J_dev  = myElements % J
    myElements % x_dev = myElements % x
    myElements % xBound_dev = myElements % xBound
    myElements % nHat_dev   = myElements % nHat


  END SUBROUTINE UpdateDevice_HexElements

  SUBROUTINE UpdateHost_HexElements( myElements )
    IMPLICIT NONE
    CLASS( HexElements ), INTENT(inout) :: myElements

    myElements % elementID = myElements % elementID_dev
    myElements % nodeIDs   = myElements % nodeIDs_dev
    myElements % neighbors = myElements % neighbors_dev
    myElements % dxds = myElements % dxds_dev
    myElements % dxdp = myElements % dxdp_dev
    myElements % dxdq = myElements % dxdq_dev
    myElements % dyds = myElements % dyds_dev
    myElements % dydp = myElements % dydp_dev
    myElements % dydq = myElements % dydq_dev
    myElements % dzds = myElements % dzds_dev
    myElements % dzdp = myElements % dzdp_dev
    myElements % dzdq = myElements % dzdq_dev
    myElements % Ja = myElements % Ja_dev
    myElements % J  = myElements % J_dev
    myElements % x = myElements % x_dev
    myElements % xBound = myElements % xBound_dev
    myElements % nHat   = myElements % nHat_dev


  END SUBROUTINE UpdateHost_HexElements
#endif

!
!> \addtogroup HexElements_CLASS
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
!!  This SUBROUTINE depends on <BR>
!!   Module \ref HexElements_CLASS : PRIVATE FUNCTION TransfiniteInterpolation <BR>
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
!!   <tr> <td> in/out <th> myElements <td> HexElements <td>
!!                         On <B>input</B>, an initialized HexElements DATA structure, <BR>
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
  SUBROUTINE GenerateMesh_HexElements( myElements, interp, theSurfaces )

    IMPLICIT NONE
    CLASS( HexElements ), INTENT(inout) :: myElements
    TYPE( Lagrange ), INTENT(in)        :: interp
    TYPE( Surfaces ), INTENT(in)        :: theSurfaces
    ! Local
    INTEGER    :: i, j, k, iEl
    REAL(prec) :: s(0:interp % N), p, x(1:3)

    s = interp % interpolationPoints

    DO iEl = 1, myElements % nElements

      DO k = 0, interp % N
        DO j = 0,interp % N
          DO i = 0,interp % N
            x = TransfiniteInterpolation( theSurfaces, iEl, s(i), s(j), s(k) )
            myElements % x(i,j,k,1:3,iEl) = x
          ENDDO
        ENDDO
      ENDDO

      ! DO the boundary locations
      DO j = 0, interp % N
        DO i = 0, interp % N

          p = -1.0_prec ! south boundary
          myElements % xBound(i,j,1:3,south,iEl) = TransfiniteInterpolation( theSurfaces, iEl, s(i), p, s(j) )

          ! west boundary
          myElements % xBound(i,j,1:3,west,iEl) = TransfiniteInterpolation( theSurfaces, iEl, p, s(i), s(j) )

          ! bottom boundary
          myElements % xBound(i,j,1:3,bottom,iEl) = TransfiniteInterpolation( theSurfaces, iEl, s(i), s(j), p )

          p = 1.0_prec  ! north boundary
          myElements % xBound(i,j,1:3,north,iEl) = TransfiniteInterpolation( theSurfaces, iEl, s(i), p, s(j) )

          ! east boundary
          myElements % xBound(i,j,1:3,east,iEl) = TransfiniteInterpolation( theSurfaces, iEl, p, s(i), s(j) )

          ! top boundary
          myElements % xBound(i,j,1:3,top,iEl) = TransfiniteInterpolation( theSurfaces, iEl, s(i), s(j), p )

        ENDDO
      ENDDO

    ENDDO

  END SUBROUTINE GenerateMesh_HexElements
!
  SUBROUTINE ResetInternalMesh_HexElements( myElements, interp )

    IMPLICIT NONE
    CLASS( HexElements ), INTENT(inout) :: myElements
    TYPE( Lagrange ), INTENT(in)        :: interp
    ! Local
    INTEGER    :: i, j, k, iEl


    DO iEl = 1, myElements % nElements
      DO k = 0, interp % N
        DO j = 0, interp % N
          DO i = 0, interp % N
            myElements % x(i,j,k,1:3,iEl) = TransfiniteInterpolation_Alt( interp, &
              myElements % xBound(0:interp % N, 0:interp % N, 1, 1:6,iEl), &
              myElements % xBound(0:interp % N, 0:interp % N, 2, 1:6,iEl), &
              myElements % xBound(0:interp % N, 0:interp % N, 3, 1:6,iEl), &
              interp % interpolationPoints(i), &
              interp % interpolationPoints(j), &
              interp % interpolationPoints(k) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO



  END SUBROUTINE ResetInternalMesh_HexElements
!
!> \addtogroup HexElements_CLASS
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
!! represented as a 3-D Lagrange interpolant. DIFferentiation of the interpolant allows us to
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
!!  This SUBROUTINE depends on <BR>
!!   Module \ref Lagrange_CLASS : FUNCTION \ref ApplyDerivativeMatrix_3D_Lagrange <BR>
!!   Module \ref HexElements_CLASS : S/R \ref GenerateBoundaryMetrics_HexElements <BR>
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
!!   <tr> <td> in/out <th> myElements <td> HexElements <td>
!!                         On <B>input</B>, an initialized HexElements DATA structure, <BR>
!!                         On <B>output</B>, the metric terms are filled in, including the
!!                         outward pointing normal vectors on the element boundaries
!!   <tr> <td> in <th> interp <td> Lagrange <td>  Lagrange interpolant that defines the
!!                                                  reference computational mesh
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE GenerateMetrics_HexElements( myElements, interp )

    IMPLICIT NONE
    CLASS( HexElements ), INTENT(inout) :: myElements
    TYPE( Lagrange ), INTENT(in)        :: interp
    ! Local
    INTEGER    :: i, j, k, iEl, N
    REAL(prec) :: cv(1:3,1:3)
    REAL(prec) :: xGradient(1:3,0:interp % N,0:interp % N,0:interp % N,1:3,1:myElements % nElements)
    REAL(prec) :: v(0:interp % N, 0:interp % N, 0:interp % N,1:3,1:myElements % nElements)
    REAL(prec) :: Dv(1:3,0:interp % N, 0:interp % N, 0:interp % N,1:3,1:myElements % nElements)

    N = interp % N

    ! Forced call to the CPU Kernel for calculating gradient
    xGradient = CalculateGradient_3D_Lagrange( interp, &
                                               myElements % x,&
                                               3, myElements % nElements )

    DO iEl = 1, myElements % nElements

      DO k = 0, N
        DO j = 0, N
          DO i = 0, N

            myElements % dxds(i,j,k,iEl) = xGradient(1,i,j,k,1,iEl) 
            myElements % dxdp(i,j,k,iEl) = xGradient(2,i,j,k,1,iEl) 
            myElements % dxdq(i,j,k,iEl) = xGradient(3,i,j,k,1,iEl) 
            myElements % dyds(i,j,k,iEl) = xGradient(1,i,j,k,2,iEl) 
            myElements % dydp(i,j,k,iEl) = xGradient(2,i,j,k,2,iEl) 
            myElements % dydq(i,j,k,iEl) = xGradient(3,i,j,k,2,iEl) 
            myElements % dzds(i,j,k,iEl) = xGradient(1,i,j,k,3,iEl) 
            myElements % dzdp(i,j,k,iEl) = xGradient(2,i,j,k,3,iEl) 
            myElements % dzdq(i,j,k,iEl) = xGradient(3,i,j,k,3,iEl) 

            cv(1,1:3) = xGradient(1:3,i,j,k,1,iEl)
            cv(2,1:3) = xGradient(1:3,i,j,k,2,iEl)
            cv(3,1:3) = xGradient(1:3,i,j,k,3,iEl)
            myElements % J(i,j,k,iEl) = Determinant( cv, 3 )

          ENDDO
        ENDDO
      ENDDO

    ENDDO

#ifdef SKEW_METRICS

    !Ja_1
    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            myElements % Ja(i,j,k,1,1,iEl) = myElements % dydp(i,j,k,iEl)*myElements % dzdq(i,j,k,iEl) - &
                                             myElements % dzdp(i,j,k,iEl)*myElements % dydq(i,j,k,iEl)

            myElements % Ja(i,j,k,2,1,iEl) = myElements % dzdp(i,j,k,iEl)*myElements % dxdq(i,j,k,iEl) - &
                                             myElements % dxdp(i,j,k,iEl)*myElements % dzdq(i,j,k,iEl)

            myElements % Ja(i,j,k,3,1,iEl) = myElements % dxdp(i,j,k,iEl)*myElements % dydq(i,j,k,iEl) - &
                                             myElements % dydp(i,j,k,iEl)*myElements % dxdq(i,j,k,iEl)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    !Ja_2
    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            myElements % Ja(i,j,k,1,2,iEl) = myElements % dydq(i,j,k,iEl)*myElements % dzds(i,j,k,iEl) - &
                                             myElements % dzdq(i,j,k,iEl)*myElements % dyds(i,j,k,iEl)

            myElements % Ja(i,j,k,2,2,iEl) = myElements % dzdq(i,j,k,iEl)*myElements % dxds(i,j,k,iEl) - &
                                             myElements % dxdq(i,j,k,iEl)*myElements % dzds(i,j,k,iEl)

            myElements % Ja(i,j,k,3,2,iEl) = myElements % dxdq(i,j,k,iEl)*myElements % dyds(i,j,k,iEl) - &
                                             myElements % dydq(i,j,k,iEl)*myElements % dxds(i,j,k,iEl)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    !Ja_3
    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            myElements % Ja(i,j,k,1,3,iEl) = myElements % dyds(i,j,k,iEl)*myElements % dzdp(i,j,k,iEl) - &
                                             myElements % dzds(i,j,k,iEl)*myElements % dydp(i,j,k,iEl)

            myElements % Ja(i,j,k,2,3,iEl) = myElements % dzds(i,j,k,iEl)*myElements % dxdp(i,j,k,iEl) - &
                                             myElements % dxds(i,j,k,iEl)*myElements % dzdp(i,j,k,iEl)

            myElements % Ja(i,j,k,3,3,iEl) = myElements % dxds(i,j,k,iEl)*myElements % dydp(i,j,k,iEl) - &
                                             myElements % dyds(i,j,k,iEl)*myElements % dxdp(i,j,k,iEl)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
#else
    ! Generate the contravariant basis vectors using the curl form ( Kopriva, 2006 )
    !Ja_1
    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            v(i,j,k,1,iEl)  = -myElements % x(i,j,k,3,iEl)*myElements % dyds(i,j,k,iEl)
            v(i,j,k,2,iEl)  = -myElements % x(i,j,k,3,iEl)*myElements % dydp(i,j,k,iEl)
            v(i,j,k,3,iEl)  = -myElements % x(i,j,k,3,iEl)*myElements % dydq(i,j,k,iEl)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Forced call to the CPU Kernel for calculating gradient
    Dv = CalculateGradient_3D_Lagrange( interp, v, 3, myElements % nElements )
    ! Ja(i,j,k,a,b,iEl) -- a-th direction of the b-th contravariant basis vector
    ! Take the curl to obtain the first dimension of each of the contravariant basis vectors
    ! The contravariant metric tensor stores each contravariant basis vector in each column
    ! of the tensor
    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            myElements % Ja(i,j,k,1,1,iEl) = ( Dv(2,i,j,k,3,iEl) - Dv(3,i,j,k,2,iEl) )
            myElements % Ja(i,j,k,1,2,iEl) = -( Dv(1,i,j,k,3,iEl) - Dv(3,i,i,k,1,iEl) )
            myElements % Ja(i,j,k,1,3,iEl) = ( Dv(1,i,j,k,2,iEl) - Dv(2,i,j,k,1,iEl) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    !Ja_2
    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            v(i,j,k,1,iEl)  = -myElements % x(i,j,k,1,iEl)*myElements % dzds(i,j,k,iEl)
            v(i,j,k,2,iEl)  = -myElements % x(i,j,k,1,iEl)*myElements % dzdp(i,j,k,iEl)
            v(i,j,k,3,iEl)  = -myElements % x(i,j,k,1,iEl)*myElements % dzdq(i,j,k,iEl)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    Dv = CalculateGradient_3D_Lagrange( interp, v, 3, myElements % nElements )

    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            myElements % Ja(i,j,k,2,1,iEl) = ( Dv(2,i,j,k,3,iEl) - Dv(3,i,j,k,2,iEl) )
            myElements % Ja(i,j,k,2,2,iEl) = -( Dv(1,i,j,k,3,iEl) - Dv(3,i,i,k,1,iEl) )
            myElements % Ja(i,j,k,2,3,iEl) = ( Dv(1,i,j,k,2,iEl) - Dv(2,i,j,k,1,iEl) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    !Ja_3
    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            v(i,j,k,1,iEl)  = -myElements % x(i,j,k,2,iEl)*myElements % dxds(i,j,k,iEl)
            v(i,j,k,2,iEl)  = -myElements % x(i,j,k,2,iEl)*myElements % dxdp(i,j,k,iEl)
            v(i,j,k,3,iEl)  = -myElements % x(i,j,k,2,iEl)*myElements % dxdq(i,j,k,iEl)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    Dv = CalculateGradient_3D_Lagrange( interp, v, 3, myElements % nElements )

    DO iEl = 1, myElements % nElements
      DO k = 0, N
        DO j = 0, N
          DO i = 0, N
            myElements % Ja(i,j,k,3,1,iEl) = ( Dv(2,i,j,k,3,iEl) - Dv(3,i,j,k,2,iEl) )
            myElements % Ja(i,j,k,3,2,iEl) = -( Dv(1,i,j,k,3,iEl) - Dv(3,i,i,k,1,iEl) )
            myElements % Ja(i,j,k,3,3,iEl) = ( Dv(1,i,j,k,2,iEl) - Dv(2,i,j,k,1,iEl) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO
#endif

    CALL myElements % GenerateBoundaryMetrics( interp )

  END SUBROUTINE GenerateMetrics_HexElements
!
!> \addtogroup HexElements_CLASS
!! @{
! ================================================================================================ !
! S/R GenerateMetrics
!
!> \fn GenerateMetrics_HexElements
!! Generates and stores the outward pointing boundary normal vectors.
!!
!!  The outward pointing boundary normal vectors are equivalent to the contravariant basis vectors
!!  evaluated at the element boundaries. These are computed here by dIFferentiating the Lagrange
!!  interpolant of the mesh positions and the computational mesh boundaries.
!!
!!  This SUBROUTINE depends on <BR>
!!   Module \ref Lagrange_CLASS : FUNCTION \ref DIFferentiate_3D_Lagrange <BR>
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
!!   <tr> <td> in/out <th> myElements <td> HexElements <td>
!!                         On <B>input</B>, an initialized HexElements DATA structure, <BR>
!!                         On <B>output</B>, the outward pointing normal vectors on the element
!!                         boundaries are filled in
!!   <tr> <td> in <th> interp <td> Lagrange <td> Lagrange interpolant that defines the
!!                                                  reference computational mesh
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE GenerateBoundaryMetrics_HexElements( myElements, interp  )

    IMPLICIT NONE
    CLASS( HexElements ), INTENT(inout) :: myElements
    TYPE( Lagrange ), INTENT(in)              :: interp
    ! Local
    INTEGER    :: i, j, N, iEl
    REAL(prec) :: s(0:interp % N), p
    REAL(prec) :: Jain(0:interp % N,0:interp % N,0:interp % N)
    REAL(prec) :: Jac, signJ, nx, ny, nz
    REAL(prec) :: node(1:3)

    N = interp % N
    s = interp % interpolationPoints
    ! DO the boundary locations
    DO iEl = 1, myElements % nElements

      DO j = 0, N
        DO i = 0, N

          p = -1.0_prec  ! bottom, south, and west boundaries

          !bottom boundary
          node = (/s(i), s(j), p /)
          Jac = interp % Interpolate_3D( myElements % J(0:N,0:N,0:N,iEl), node ) !Determinant( cv, 3 )
          signJ = SIGN(1.0_prec,Jac)
          ! Setting bottom boundary normal
          Jain = myElements % Ja(0:N,0:N,0:N,1,3,iEl)
          nx = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,2,3,iEl)
          ny = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,3,3,iEl)
          nz = interp % Interpolate_3D( Jain, node )
          myElements % nHat(1:3,i,j,bottom,iEl) = -signJ*(/ nx, ny, nz /)

          node = (/ s(i), p, s(j) /)
          Jac = interp % Interpolate_3D( myElements % J(0:N,0:N,0:N,iEl), node ) !Determinant( cv, 3 )
          signJ = SIGN(1.0_prec,Jac)
          ! Setting southern boundary normal
          Jain = myElements % Ja(0:N,0:N,0:N,1,2,iEl)
          nx = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,2,2,iEl)
          ny = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,3,2,iEl)
          nz = interp % Interpolate_3D( Jain, node )
          myElements % nHat(1:3,i,j,south,iEl)= -signJ*(/ nx, ny, nz /)

          ! west boundary
          node = (/ p, s(i), s(j) /)
          Jac = interp % Interpolate_3D( myElements % J(0:N,0:N,0:N,iEl), node ) !Determinant( cv, 3 )
          signJ = SIGN(1.0_prec,Jac)
          ! Setting western boundary normal
          Jain = myElements % Ja(0:N,0:N,0:N,1,1,iEl)
          nx = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,2,1,iEl)
          ny = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,3,1,iEl)
          nz = interp % Interpolate_3D( Jain, node )
          myElements % nHat(1:3,i,j,west,iEl) = -signJ*(/ nx, ny, nz /)

          p = 1.0_prec  ! top, north, and east boundaries

          !top boundary
          node = (/s(i), s(j), p /)
          Jac = interp % Interpolate_3D( myElements % J(0:N,0:N,0:N,iEl), node )!Determinant( cv, 3 )
          signJ = SIGN(1.0_prec,Jac)
          ! Setting top boundary normal
          Jain = myElements % Ja(0:N,0:N,0:N,1,3,iEl)
          nx = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,2,3,iEl)
          ny = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,3,3,iEl)
          nz = interp % Interpolate_3D( Jain, node )
          myElements % nHat(1:3,i,j,top,iEl) = signJ*(/ nx, ny, nz /)

          !north boundary
          node = (/ s(i), p, s(j) /)
          Jac = interp % Interpolate_3D( myElements % J(0:N,0:N,0:N,iEl), node ) !Determinant( cv, 3 )
          signJ = SIGN(1.0_prec,Jac)
          ! Setting southern boundary normal
          Jain = myElements % Ja(0:N,0:N,0:N,1,2,iEl)
          nx = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,2,2,iEl)
          ny = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,3,2,iEl)
          nz = interp % Interpolate_3D( Jain, node )
          myElements % nHat(1:3,i,j,north,iEl) = signJ*(/ nx, ny, nz /)

          ! east boundary
          node = (/ p, s(i), s(j) /)
          Jac = interp % Interpolate_3D( myElements % J(0:N,0:N,0:N,iEl), node ) !Determinant( cv, 3 )
          signJ = SIGN(1.0_prec,Jac)
          ! Setting eastern boundary normal
          Jain = myElements % Ja(0:N,0:N,0:N,1,1,iEl)
          nx = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,2,1,iEl)
          ny = interp % Interpolate_3D( Jain, node )
          Jain = myElements % Ja(0:N,0:N,0:N,3,1,iEl)
          nz = interp % Interpolate_3D( Jain, node )
          myElements % nHat(1:3,i,j,east,iEl) = signJ*(/ nx, ny, nz /)

        ENDDO
      ENDDO

    ENDDO

  END SUBROUTINE GenerateBoundaryMetrics_HexElements
!
!> \addtogroup HexElements_CLASS
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
!!   <tr> <td> in/out <th> myElements <td> HexElements <td>
!!                         On <B>input</B>, an intialized and constructed HexElements DATA
!!                         structure <BR>
!!                         On <B>output</B>, the physical coordinates and metric terms have been
!!                         scaled by the given x and y scales.
!!   <tr> <td> in <th> interp <td> Lagrange <td> 3-D Lagrange interpolant DATA structure that
!!                                                  CONTAINS the computational mesh.
!!   <tr> <td> in <th> xScale <td> REAL(prec)  <td> Factor to multiply the physical x position
!!                                                  and metrics by.
!!   <tr> <td> in <th> yScale <td> REAL(prec)  <td> Factor to multiply the physical y position
!!                                                  and metrics by.
!!  </table>
!!
! ================================================================================================ !
!>@}
  SUBROUTINE ScaleGeometry_HexElements( myElements, interp, xScale, yScale, zScale )

    IMPLICIT NONE
    CLASS( HexElements ), INTENT(inout) :: myElements
    TYPE( Lagrange ), INTENT(in)              :: interp
    REAL(prec), INTENT(in)                    :: xScale, yScale, zScale
#ifdef HAVE_CUDA
    INTEGER :: istat
#endif

    myElements % x(:,:,:,1,:) = xScale*( myElements % x(:,:,:,1,:) )
    myElements % x(:,:,:,2,:) = yScale*( myElements % x(:,:,:,2,:) )
    myElements % x(:,:,:,3,:) = zScale*( myElements % x(:,:,:,3,:) )
    myElements % xBound(:,:,1,:,:) = xScale*( myElements % xBound(:,:,1,:,:) )
    myElements % xBound(:,:,2,:,:) = yScale*( myElements % xBound(:,:,2,:,:) )
    myElements % xBound(:,:,3,:,:) = zScale*( myElements % xBound(:,:,3,:,:) )

    ! Update the boundary metrics -- normals and normal lengths
    CALL myElements % GenerateMetrics( interp )
    CALL myElements % GenerateBoundaryMetrics( interp  )

#ifdef HAVE_CUDA
    CALL myElements % UpdateDevice( )
    istat = cudaDeviceSynchronize( )
#endif

  END SUBROUTINE ScaleGeometry_HexElements
!
  SUBROUTINE CalculateComputationalCoordinates( myElements, interp, x, s, elements, nCoordinates )
    IMPLICIT NONE
    CLASS( HexElements ), INTENT(in) :: myElements
    TYPE( Lagrange ), INTENT(in)     :: interp
    INTEGER, INTENT(in)              :: nCoordinates
    REAL(prec), INTENT(in)           :: x(1:3,1:nCoordinates)
    REAL(prec), INTENT(out)          :: s(1:3,1:nCoordinates)
    INTEGER, INTENT(out)             :: elements(1:nCoordinates)
    ! Local
    INTEGER    :: i, iEl, kIt, row, col
    REAL(prec) :: r(1:3), si(1:3), xi(1:3), delta(1:3), A(1:3,1:3), Ainv(1:3,1:3), dMag 
    LOGICAL    :: converged

      elements = 0
      s = 0.0_prec

      DO i = 1, nCoordinates

        DO iEl = 1, myElements % nElements

          si(1:3) = 0.0_prec

          converged = .FALSE.

          DO kIt = 1, newtonMax

            xi(1) = interp % Interpolate_3D( myElements % x(1,:,:,:,iEl), si )
            xi(2) = interp % Interpolate_3D( myElements % x(2,:,:,:,iEl), si )
            xi(3) = interp % Interpolate_3D( myElements % x(3,:,:,:,iEl), si )
              
            A(1,1) = interp % Interpolate_3D( myElements % dxds(:,:,:,iEl), si )
            A(1,2) = interp % Interpolate_3D( myElements % dxdp(:,:,:,iEl), si )
            A(1,3) = interp % Interpolate_3D( myElements % dxdq(:,:,:,iEl), si )

            A(2,1) = interp % Interpolate_3D( myElements % dyds(:,:,:,iEl), si )
            A(2,2) = interp % Interpolate_3D( myElements % dydp(:,:,:,iEl), si )
            A(2,3) = interp % Interpolate_3D( myElements % dydq(:,:,:,iEl), si )

            A(3,1) = interp % Interpolate_3D( myElements % dzds(:,:,:,iEl), si )
            A(3,2) = interp % Interpolate_3D( myElements % dzdp(:,:,:,iEl), si )
            A(3,3) = interp % Interpolate_3D( myElements % dzdq(:,:,:,iEl), si )

            Ainv = Invert_3x3( A )
            
            r(1:3) = x(1:3,i) - xi(1:3) 

            dMag = 0.0_prec
            DO row = 1, 3
              delta(row) = 0.0_prec
              DO col = 1, 3
                delta(row) = delta(row) + Ainv(row,col)*r(col)
              ENDDO
              dMag = dMag + delta(row)**2
            ENDDO

            si(1:3) = si(1:3) + delta(1:3)

            IF( SQRT(dMag) <= TOL*MAXVAL(si) ) THEN
              converged = .TRUE.
              EXIT
            ENDIF

          ENDDO

          IF( converged .AND. MAXVAL( ABS(si) ) <= 1.0_prec  )THEN

            s(1:3,i)    = si(1:3)
            elements(i) = iEl
            EXIT
            
          ENDIF
      
        ENDDO

      ENDDO

  END SUBROUTINE CalculateComputationalCoordinates
!
  FUNCTION TransfiniteInterpolation( boundingSurfaces, iEl, a, b, c ) RESULT( P )
    ! TransfiniteInterpolation
    !  Takes in the six surfaces (south, east, north, west, bottom, top) and evaluates the
    !  bidirectional mapping at xi^1 = a, xi^2 = b, xi^3 = c. The boundary of the computational
    !  coordinate system is assumed to be at +/- 1 in each direction.
    !
    ! =============================================================================================== !
    ! DECLARATIONS
    IMPLICIT NONE
    TYPE( Surfaces )  :: boundingSurfaces
    INTEGER           :: iEl
    REAL(prec)       :: a, b, c
    REAL(prec)       :: P(1:3)
    ! LOCAL
    REAL(prec)  :: P1(1:3), P2(1:3), P3(1:3)
    REAL(prec)  :: sSurf(1:3), nSurf(1:3), eSurf(1:3), wSurf(1:3), bSurf(1:3), tSurf(1:3)
    REAL(prec)  :: l1(1:2), l2(1:2), l3(1:2)
    REAL(prec)  :: ref(1:2)
    INTEGER     :: i, j, iSurf

    ref = (/ -1.0_prec, 1.0_prec /)

    ! Transfinite interpolation with linear blending USEs linear lagrange interpolating polynomials
    ! to blend the bounding surfaces.
    ! The linear blending weights in the first computational direction are calculated.

    l1 = LinearBlend( a )
    l2 = LinearBlend( b )
    l3 = LinearBlend( c )

    ! The bounding surfaces need to be evaluated at the provided computational coordinates

    wSurf = boundingSurfaces % Evaluate( (/b, c/), west + (iEl-1)*6 )   ! west
    eSurf = boundingSurfaces % Evaluate( (/b, c/), east + (iEl-1)*6 )   ! east
    sSurf = boundingSurfaces % Evaluate( (/a, c/), south + (iEl-1)*6 )  ! south
    nSurf = boundingSurfaces % Evaluate( (/a, c/), north + (iEl-1)*6 )  ! north
    bSurf = boundingSurfaces % Evaluate( (/a, b/), bottom + (iEl-1)*6 ) ! bottom
    tSurf = boundingSurfaces % Evaluate( (/a, b/), top + (iEl-1)*6 )    ! top

    ! P1 CONTAINS the interpolation in the first computational coordinate
    ! The first computational coordinate is assumed to vary between the (computational) east and
    ! west boundaries.

    P1 = l1(1)*wSurf + l1(2)*eSurf

    ! P2 CONTAINS the interpolation in the second computational coordinate
    ! The second computational coordinate is assumed to vary between the (computational) south and
    ! north boundaries.

    P2 = l2(1)*sSurf + l2(2)*nSurf

    ! P3 CONTAINS the interpolation in the first computational coordinate
    ! The first computational coordinate is assumed to vary between the (computational) bottom and
    ! top boundaries.

    P3 = l3(1)*bSurf + l3(2)*tSurf

    DO i = 1, 2

      ! Now we need to compute the tensor product of the first and second computational direction
      ! interpolants and subtract from P1.

      wSurf = boundingsurfaces % Evaluate( (/ref(i), c/), west + (iEl-1)*6 )
      eSurf = boundingsurfaces % Evaluate( (/ref(i), c/), east + (iEl-1)*6 )
      P1 = P1 - l2(i)*( wSurf*l1(1) + eSurf*l1(2) )

      ! Now we need to compute the tensor product of the first and third computational direction
      ! interpolants and subtract from P1.

      wSurf = boundingsurfaces % Evaluate( (/b, ref(i)/), west + (iEl-1)*6 )
      eSurf = boundingsurfaces % Evaluate( (/b, ref(i)/), east + (iEl-1)*6 )

      P1 = P1 - l3(i)*( wSurf*l1(1) + eSurf*l1(2) )

      ! Now we need to compute the tensor product of the second and third computational direction
      ! interpolants and subtract from P2.

      sSurf = boundingsurfaces % Evaluate( (/a, ref(i)/), south + (iEl-1)*6 )
      nSurf = boundingsurfaces % Evaluate( (/a, ref(i)/), north + (iEl-1)*6 )

      P2 = P2 - l3(i)*( sSurf*l2(1) + nSurf*l2(2) )

    ENDDO

    ! Next, the compounded tensor product is computed and added to P3.
    DO j = 1,2
      DO i = 1,2

        wSurf = boundingsurfaces % Evaluate( (/ref(i), ref(j)/), west + (iEl-1)*6 )
        eSurf = boundingsurfaces % Evaluate( (/ref(i), ref(j)/), east + (iEl-1)*6 )
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
    TYPE( Lagrange ) :: interp
    REAL(prec)       :: x(0:interp % N,0:interp % N,1:6), y(0:interp % N,0:interp % N,1:6), z(0:interp % N,0:interp % N,1:6)
    REAL(prec)       :: a, b, c
    REAL(prec)       :: P(1:3)
    ! LOCAL
    REAL(prec)  :: P1(1:3), P2(1:3), P3(1:3)
    REAL(prec)  :: sSurf(1:3), nSurf(1:3), eSurf(1:3), wSurf(1:3), bSurf(1:3), tSurf(1:3)
    REAL(prec)  :: l1(1:2), l2(1:2), l3(1:2)
    REAL(prec)  :: ref(1:2)
    INTEGER     :: i, j

    ref = (/ -1.0_prec, 1.0_prec /)

    ! Transfinite interpolation with linear blending USEs linear lagrange interpolating polynomials
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

    ! P1 CONTAINS the interpolation in the first computational coordinate
    ! The first computational coordinate is assumed to vary between the (computational) east and
    ! west boundaries.

    P1 = l1(1)*wSurf + l1(2)*eSurf

    ! P2 CONTAINS the interpolation in the second computational coordinate
    ! The second computational coordinate is assumed to vary between the (computational) south and
    ! north boundaries.

    P2 = l2(1)*sSurf + l2(2)*nSurf

    ! P3 CONTAINS the interpolation in the first computational coordinate
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

    P = 0.5_prec*( (1.0_prec - a)*valLeft + (1.0_prec + a)*valRight )

  END FUNCTION Unidirectional
  FUNCTION LinearBlend( a ) RESULT( weights )

    IMPLICIT NONE
    REAL(prec) :: a
    REAL(prec) :: weights(1:2)

    weights(1) = 0.5_prec*(1.0_prec - a)
    weights(2) = 0.5_prec*(1.0_prec + a)

  END FUNCTION LinearBlend
!
END MODULE HexElements_CLASS
