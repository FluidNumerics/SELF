! NodalDGSolution_3D_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE NodalDGSolution_3D_Class

USE ModelPrecision
USE NodalDG_Class
USE HexMesh_Class
USE HDF5

IMPLICIT NONE

! ================================================================================================ !
!
!  The NodalDGSolution_3D addtogroup provides attributes for storing a solution and the flux and source
!  terms associated with a conservation law. 
!
!  Wrapper routines are provided that call routines in the
!  NodalDG_Class to interpolate the solution array to the element boundaries and calculate the
!  flux divergence in both the strong and weak form.
!  
!  When building a Spectral Element solver, this addtogroup can be used to cleanly calculate the tendency in a 
!  PDE.
!
!  >> *For the "boundary" arrays, the faces of a hexahedral element are numbered as SOUTH=1, 
!      EAST=2, NORTH=3, WEST=4, BOTTOM=5, TOP=6.
!
! ================================================================================================ !

  TYPE NodalDGSolution_3D

    INTEGER                 :: N, nEquations, nElements, nBoundaryFaces

    REAL(prec), ALLOCATABLE :: solution(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: boundarySolution(:,:,:,:,:) 

    REAL(prec), ALLOCATABLE :: prescribedState(:,:,:,:)
    REAL(prec), ALLOCATABLE :: externalState(:,:,:,:)

    REAL(prec), ALLOCATABLE :: solutionGradient(:,:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: boundaryGradientFlux(:,:,:,:,:,:) 
   
    REAL(prec), ALLOCATABLE :: flux(:,:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: boundaryFlux(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: fluxDivergence(:,:,:,:,:)

    REAL(prec), ALLOCATABLE :: source(:,:,:,:,:)

    REAL(prec), ALLOCATABLE :: tendency(:,:,:,:,:)


#ifdef HAVE_CUDA

    INTEGER, DEVICE, ALLOCATABLE    :: N_dev, nEquations_dev, nElements_dev, nBoundaryFaces_dev

    REAL(prec), DEVICE, ALLOCATABLE :: solution_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: boundarySolution_dev(:,:,:,:,:) 

    REAL(prec), DEVICE, ALLOCATABLE :: prescribedState_dev(:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: externalState_dev(:,:,:,:)

    REAL(prec), DEVICE, ALLOCATABLE :: solutionGradient_dev(:,:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: boundaryGradientFlux_dev(:,:,:,:,:,:)  

    REAL(prec), DEVICE, ALLOCATABLE :: flux_dev(:,:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: boundaryFlux_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: fluxDivergence_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: source_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: tendency_dev(:,:,:,:,:)

#endif

    CONTAINS

      PROCEDURE :: Build => Build_NodalDGSolution_3D
      PROCEDURE :: Trash => Trash_NodalDGSolution_3D
  
      PROCEDURE :: Calculate_Solution_At_Boundaries
      PROCEDURE :: Calculate_Weak_Flux_Divergence
      PROCEDURE :: Calculate_Strong_Flux_Divergence

      PROCEDURE :: Mapped_BassiRebay_Gradient
      PROCEDURE :: Mapped_Strong_Gradient

#ifdef HAVE_CUDA

      PROCEDURE :: UpdateDevice => UpdateDevice_NodalDGSolution_3D
      PROCEDURE :: UpdateHost   => UpdateHost_NodalDGSolution_3D

#endif

  END TYPE NodalDGSolution_3D

CONTAINS

! ================================================================================================ !
! Build_NodalDGSolution_3D  
!   Allocates space for the NodalDGSolution_3D attributes and initializes arrays to 0.0_prec.
! 
!   Usage :
! 
!     TYPE(NodalDGSolution_3D) :: nodalSolution 
!     INTEGER                  :: N, nEquations, nElements
!
!       CALL nodalSolution % Build( N, nEquations, nElements )
!
!   Input/Output :
!
!     nodalSolution (out)
!       NodalDGSolution structure. All array attributes are allocated and initialized to 0.0_prec
!       and the polynomial degree, number of equations, and the number of elements are set upon exit.
!
!     N (in)
!       Polynomial degree of the spectral element approximation
!
!     nEquations (in)
!       The number of equations/variables to be held underneath this data structure
!
!     nElements (in)
!       The number of elements in the spectral element mesh.
!
! ================================================================================================ ! 

  SUBROUTINE Build_NodalDGSolution_3D( myDGS, N, nEquations, nElements, nBoundaryFaces )
    IMPLICIT NONE
    CLASS(NodalDGSolution_3D), INTENT(inout) :: myDGS
    INTEGER, INTENT(in)                      :: N, nEquations, nElements, nBoundaryFaces
    ! Local
    INTEGER :: nUse

     nUse = N
#ifdef HAVE_CUDA

     IF( N > 7 )THEN

       PRINT*, 'Module NodalDGSolution_3D_Class : S/R Build_NodalDGSolution_3D :'
       PRINT*, '  WARNING : Polynomial degree > 7 with CUDA not permitted. '
       PRINT*, '            Forcing polynomial degree = 7.'
       nUse = 7

     ENDIF
       
#endif
      
      myDGS % N          = Nuse
      myDGS % nEquations = nEquations
      myDGS % nElements  = nElements
      myDGS % nBoundaryFaces = nBoundaryFaces

      ALLOCATE( myDGS % solution(0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % boundarySolution(0:nUse,0:nUse,1:nEquations,1:6,1:nElements), &
                myDGS % prescribedState(0:nUse,0:nUse,1:nEquations,1:nBoundaryFaces), &
                myDGS % externalState(0:nUse,0:nUse,1:nEquations,1:nBoundaryFaces), &
                myDGS % solutionGradient(1:3,0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % boundaryGradientFlux(1:3,0:nUse,0:nUse,1:nEquations,1:6,1:nElements), &
                myDGS % flux(1:3,0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % boundaryFlux(0:nUse,0:nUse,1:nEquations,1:6,1:nElements), &
                myDGS % fluxDivergence(0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % source(0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % tendency(0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements) )
      
      myDGS % solution             = 0.0_prec
      myDGS % boundarySolution     = 0.0_prec
      myDGS % prescribedState      = 0.0_prec
      myDGS % externalState        = 0.0_prec
      myDGS % solutionGradient     = 0.0_prec
      myDGS % boundaryGradientFlux = 0.0_prec
      myDGS % flux                 = 0.0_prec
      myDGS % boundaryFlux         = 0.0_prec
      myDGS % fluxDivergence       = 0.0_prec
      myDGS % source               = 0.0_prec
      myDGS % tendency             = 0.0_prec

#ifdef HAVE_CUDA
      ALLOCATE( myDGS % N_dev, myDGS % nEquations_dev, myDGS % nElements_dev, myDGS % nBoundaryFaces_dev )
      myDGS % N_dev          = nUse
      myDGS % nEquations_dev = nEquations
      myDGS % nElements_dev  = nElements
      myDGS % nBoundaryFaces_dev  = nBoundaryFaces
      
      ALLOCATE( myDGS % solution_dev(0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % boundarySolution_dev(0:nUse,0:nUse,1:nEquations,1:6,1:nElements), &
                myDGS % prescribedState_dev(0:nUse,0:nUse,1:nEquations,1:nBoundaryFaces), &
                myDGS % externalState_dev(0:nUse,0:nUse,1:nEquations,1:nBoundaryFaces), &
                myDGS % solutionGradient_dev(1:3,0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % boundaryGradientFlux_dev(1:3,0:nUse,0:nUse,1:nEquations,1:6,1:nElements), &
                myDGS % flux_dev(1:3,0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % boundaryFlux_dev(0:nUse,0:nUse,1:nEquations,1:6,1:nElements), &
                myDGS % fluxDivergence_dev(0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % source_dev(0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements), &
                myDGS % tendency_dev(0:nUse,0:nUse,0:nUse,1:nEquations,1:nElements) )

      myDGS % solution_dev             = 0.0_prec
      myDGS % boundarySolution_dev     = 0.0_prec
      myDGS % prescribedState_dev      = 0.0_prec
      myDGS % externalState_dev        = 0.0_prec
      myDGS % solutionGradient_dev     = 0.0_prec
      myDGS % boundaryGradientFlux_dev = 0.0_prec
      myDGS % flux_dev                 = 0.0_prec
      myDGS % boundaryFlux_dev         = 0.0_prec
      myDGS % fluxDivergence_dev       = 0.0_prec
      myDGS % source_dev               = 0.0_prec
      myDGS % tendency_dev             = 0.0_prec

#endif

  END SUBROUTINE Build_NodalDGSolution_3D

! ================================================================================================ !
! Trash_NodalDGSolution_3D
!   Frees memory held by the attributes of the NodalDGSolution_3D data structure.
! 
!   Usage :
! 
!     TYPE(NodalDGSolution_3D) :: nodalSolution
!
!     CALL nodalSolution % Trash( )
! 
!   Input/Output :
!
!     nodalSolution (in/out) 
!       On output, the memory held by the array attributes has been deallocated. If CUDA is enabled
!       all DEVICE memory used by this structure is deallocated.
!   
! ================================================================================================ ! 

  SUBROUTINE Trash_NodalDGSolution_3D( myDGS )
    IMPLICIT NONE
    CLASS(NodalDGSolution_3D), INTENT(inout) :: myDGS

      DEALLOCATE( myDGS % solution, &
                  myDGS % boundarySolution, &
                  myDGS % prescribedState, &
                  myDGS % externalState, &
                  myDGS % solutionGradient, &
                  myDGS % boundaryGradientFlux, &
                  myDGS % flux, &
                  myDGS % boundaryFlux, &
                  myDGS % fluxDivergence, &
                  myDGS % source, &
                  myDGS % tendency )


#ifdef HAVE_CUDA
      DEALLOCATE( myDGS % N_dev, myDGS % nEquations_dev, myDGS % nElements_dev, myDGS % nBoundaryFaces_dev )
      DEALLOCATE( myDGS % solution_dev, &
                  myDGS % boundarySolution_dev, &
                  myDGS % prescribedState_dev, &
                  myDGS % externalState_dev, &
                  myDGS % solutionGradient_dev, &
                  myDGS % boundaryGradientFlux_dev, &
                  myDGS % flux_dev, &
                  myDGS % boundaryFlux_dev, &
                  myDGS % fluxDivergence_dev, &
                  myDGS % source_dev, &
                  myDGS % tendency_dev )

#endif

  END SUBROUTINE Trash_NodalDGSolution_3D

! ================================================================================================ !
! Calculate_Solution_At_Boundaries
!   Passes the solution attribute to the routine \ref CalculateFunctionsAtBoundaries_3D to calculate
!   the boundarySolution attribute. This wrapper routine passes the device attributes when CUDA is 
!   enabled. 
! 
!  Usage : </H2> 
!    TYPE</B>(NodalDGSolution_3D) :: nodalSolution 
!    TYPE</B>(NodalDG)            :: dgStorage 
!
!      CALL</B> nodalSolution % Calculate_Solution_At_Boundaries( dgStorage ) 
! 
!  Input/Output : </H2>
!
!    nodalSolution  (in/out) 
!      Instance of the \ref NodalDGSolution_3D structure that has been constructed. On output, the 
!      boundarySolution attribute is updated by interpolating the solution attribute to element
!      faces.
!
!    dgStorage  (in) 
!      Instance of the \ref NodalDG structure.
!
!   
! ================================================================================================ ! 

  SUBROUTINE Calculate_Solution_At_Boundaries( myDGS, dgStorage )
    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: myDGS
    TYPE( NodalDG ), INTENT(in)                :: dgStorage

#ifdef HAVE_CUDA

      CALL CalculateFunctionsAtBoundaries_3D( dgStorage, &
                                              myDGS % solution_dev, &
                                              myDGS % boundarySolution_dev, &
                                              myDGS % nEquations_dev, & 
                                              myDGS % nElements_dev )
#else
      CALL CalculateFunctionsAtBoundaries_3D( dgStorage, &
                                              myDGS % solution, &
                                              myDGS % boundarySolution, &
                                              myDGS % nEquations, & 
                                              myDGS % nElements )

#endif

  END SUBROUTINE Calculate_Solution_At_Boundaries

! ================================================================================================ !
! Calculate_Weak_Gradient
!   Passes the solution and boundarySolution attributes to the routine \ref DG_Gradient_3D to calculate
!   the solutionGradient attribute. This wrapper routine passes the device attributes when CUDA is 
!   enabled. This routine assumes no mapping to physical space.
!
!   Usage : </H2> 
!
!     TYPE</B>(NodalDGSolution_3D) :: nodalSolution 
!     TYPE</B>(NodalDG)            :: dgStorage 
!
!       CALL</B> nodalSolution % Calculate_Weak_Gradient( dgStorage ) 
! 
!   Input/Output : </H2>
!
!     nodalSolution  (in/out) 
!       Instance of the \ref NodalDGSolution_3D structure that has been constructed. On output, the 
!       solutionGradient attribute is updated.
!
!     dgStorage  (in) 
!       Instance of the \ref NodalDG structure.
!
! ================================================================================================ ! 

  SUBROUTINE Calculate_Weak_Gradient( myDGS, dgStorage )
    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: myDGS
    TYPE( NodalDG ), INTENT(in)                :: dgStorage

#ifdef HAVE_CUDA
      CALL DG_Gradient_3D( dgStorage, &
                           myDGS % solution_dev, &
                           myDGS % boundarySolution_dev, &
                           myDGS % solutionGradient_dev, &
                           myDGS % nEquations_dev, & 
                           myDGS % nElements_dev )
#else
      CALL DG_Gradient_3D( dgStorage, &
                           myDGS % solution, &
                           myDGS % boundarySolution, &
                           myDGS % solutionGradient, &
                           myDGS % nEquations, & 
                           myDGS % nElements )

#endif

  END SUBROUTINE Calculate_Weak_Gradient
  
! ================================================================================================ !
! Calculate_Weak_Flux_Divergence
!   Passes the flux and boundaryFlux attributes to the routine \ref DG_Divergence_3D to calculate
!   the fluxDivergence attribute. This wrapper routine passes the device attributes when CUDA is 
!   enabled. This routine calculates the flux using the weak formulation of the divergence operator
!   in computational coordinates
!
!   Usage : </H2> 
!
!     TYPE</B>(NodalDGSolution_3D) :: nodalSolution 
!     TYPE</B>(NodalDG)            :: dgStorage 
!
!       CALL</B> nodalSolution % Calculate_Weak_Flux_Divergence( dgStorage ) 
! 
!   Input/Output : </H2>
!
!     nodalSolution  (in/out) 
!       Instance of the \ref NodalDGSolution_3D structure that has been constructed. On output, the 
!       boundarySolution attribute is updated by interpolating the solution attribute to element
!       faces.
!
!     dgStorage  (in) 
!       Instance of the \ref NodalDG structure.
!
! ================================================================================================ ! 

  SUBROUTINE Calculate_Weak_Flux_Divergence( myDGS, dgStorage )
    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: myDGS
    TYPE( NodalDG ), INTENT(in)                :: dgStorage

#ifdef HAVE_CUDA
      CALL DG_Divergence_3D( dgStorage, &
                             myDGS % flux_dev, &
                             myDGS % boundaryFlux_dev, &
                             myDGS % fluxDivergence_dev, &
                             myDGS % nEquations_dev, & 
                             myDGS % nElements_dev )
#else
      CALL DG_Divergence_3D( dgStorage, &
                             myDGS % flux, &
                             myDGS % boundaryFlux, &
                             myDGS % fluxDivergence, &
                             myDGS % nEquations, & 
                             myDGS % nElements )

#endif


  END SUBROUTINE Calculate_Weak_Flux_Divergence

! ================================================================================================ !
! Calculate_Strong_Flux_Divergence
!   Passes the flux attribute to the routine \ref CalculateDivergence to calculate
!   the fluxDivergence attribute. This wrapper routine passes the device attributes when CUDA is 
!   enabled. This routine calculates the flux using the strong formulation of the divergence operator
!   in computational coordinates
!
!   Usage : </H2> 
!
!     TYPE</B>(NodalDGSolution_3D) :: nodalSolution 
!     TYPE</B>(NodalDG)            :: dgStorage 
!
!       CALL</B> nodalSolution % Calculate_Strong_Flux_Divergence( dgStorage ) 
! 
!   Input/Output : </H2>
!
!     nodalSolution  (in/out) 
!       Instance of the \ref NodalDGSolution_3D structure that has been constructed. On output, the 
!       boundarySolution attribute is updated by interpolating the solution attribute to element
!       faces.
!
!     dgStorage  (in) 
!       Instance of the \ref NodalDG structure.
!
! ================================================================================================ ! 

  SUBROUTINE Calculate_Strong_Flux_Divergence( myDGS, dgStorage )
    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: myDGS
    TYPE( NodalDG ), INTENT(in)                :: dgStorage

#ifdef HAVE_CUDA
      CALL CalculateDivergence_3D( dgStorage % interp, &
                                   myDGS % flux_dev, &
                                   myDGS % fluxDivergence_dev, &
                                   myDGS % nEquations_dev, & 
                                   myDGS % nElements_dev )
#else
      CALL CalculateDivergence_3D( dgStorage % interp, &
                                   myDGS % flux, &
                                   myDGS % fluxDivergence, &
                                   myDGS % nEquations, & 
                                   myDGS % nElements )

#endif

  END SUBROUTINE Calculate_Strong_Flux_Divergence

! Algorithms for computing mapped derivative operations !

  SUBROUTINE Mapped_Strong_FluxDivergence( nodalSolution, dgStorage, mesh )
   
    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: nodalSolution 
    TYPE( NodalDG ), INTENT(in)                :: dgStorage
    TYPE( HexMesh ), INTENT(in)                :: mesh
    ! Local
    INTEGER    :: ii, i, j, k, iEq, iEl, idir
    REAL(prec) :: f(1:3,0:dgStorage % N,0:dgStorage % N, 0:dgStorage % N), df

      DO iEl = 1, nodalSolution % nElements
        DO iEq = 1, nodalSolution % nEquations

          f = 0.0_prec
          DO idir = 1, 3
          
            DO k = 0, dgStorage % N
              DO j = 0, dgStorage % N
                DO i = 0, dgStorage % N   
 
                  f(1,i,j,k) = f(1,i,j,k) + nodalSolution % flux(idir,i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,1,iEl)
                  f(2,i,j,k) = f(2,i,j,k) + nodalSolution % flux(idir,i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,2,iEl)
                  f(3,i,j,k) = f(3,i,j,k) + nodalSolution % flux(idir,i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,3,iEl)

                ENDDO
              ENDDO
            ENDDO

          ENDDO

            DO k = 0, dgStorage % N
              DO j = 0, dgStorage % N
                DO i = 0, dgStorage % N   

                  df = 0.0_prec
                  DO ii = 0, dgStorage % N
                    df = df + dgStorage % interp % derivativeMatrixTranspose(ii,i)*f(1,ii,j,k) + &
                              dgStorage % interp % derivativeMatrixTranspose(ii,j)*f(2,i,ii,k) + &
                              dgStorage % interp % derivativeMatrixTranspose(ii,k)*f(3,i,j,ii)
                  ENDDO
                   
                  nodalSolution % fluxDivergence(i,j,k,iEq,iEl) = df/mesh % elements % J(i,j,k,iEl)

                ENDDO  
              ENDDO
            ENDDO
          ENDDO

      ENDDO
                            

  END SUBROUTINE Mapped_Strong_FluxDivergence

  SUBROUTINE Mapped_Strong_Gradient( nodalSolution, dgStorage, mesh ) 

    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: nodalSolution 
    TYPE( NodalDG ), INTENT(in)                :: dgStorage
    TYPE( HexMesh ), INTENT(in)                :: mesh
    !
    INTEGER    :: ii, i, j, k, iEq, iEl, idir
    REAL(prec) :: f(1:3,0:dgStorage % N,0:dgStorage % N, 0:dgStorage % N), df

      DO iEl = 1, nodalSolution % nElements
        DO iEq = 1, nodalSolution % nEquations
          DO idir = 1, 3
          
            DO k = 0, dgStorage % N
              DO j = 0, dgStorage % N
                DO i = 0, dgStorage % N   
 
                  f(1,i,j,k) = nodalSolution % solution(i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,1,iEl)
                  f(2,i,j,k) = nodalSolution % solution(i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,2,iEl)
                  f(3,i,j,k) = nodalSolution % solution(i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,3,iEl)

                ENDDO
              ENDDO
            ENDDO

            DO k = 0, dgStorage % N
              DO j = 0, dgStorage % N
                DO i = 0, dgStorage % N   

                  df = 0.0_prec
                  DO ii = 0, dgStorage % N
                    df = df + dgStorage % interp % derivativeMatrixTranspose(ii,i)*f(1,ii,j,k) + &
                              dgStorage % interp % derivativeMatrixTranspose(ii,j)*f(2,i,ii,k) + &
                              dgStorage % interp % derivativeMatrixTranspose(ii,k)*f(3,i,j,ii)
                  ENDDO

                  nodalSolution % solutionGradient(idir,i,j,k,iEq,iEl) = df/mesh % elements % J(i,j,k,iEl)

                ENDDO
              ENDDO
            ENDDO

          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE Mapped_Strong_Gradient

  SUBROUTINE Mapped_BassiRebay_Gradient( nodalSolution, dgStorage, mesh )

    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: nodalSolution 
    TYPE( NodalDG ), INTENT(in)                :: dgStorage
    TYPE( HexMesh ), INTENT(in)                :: mesh


      CALL BassiRebay_Gradient_InternalFaceFlux( nodalSolution, mesh )

      CALL BassiRebay_Gradient_ExternalFaceFlux( nodalSolution, mesh )

      CALL MappedGradientFluxDivergence( nodalSolution, dgStorage, mesh )


  END SUBROUTINE Mapped_BassiRebay_Gradient


  SUBROUTINE MappedGradientFluxDivergence( nodalSolution, dgStorage, mesh )
   
    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: nodalSolution 
    TYPE( NodalDG ), INTENT(in)                :: dgStorage
    TYPE( HexMesh ), INTENT(in)                :: mesh
    ! Local
    INTEGER    :: ii, i, j, k, iEq, iEl, idir
    REAL(prec) :: f(1:3,0:dgStorage % N,0:dgStorage % N, 0:dgStorage % N), df

      DO iEl = 1, nodalSolution % nElements
        DO iEq = 1, nodalSolution % nEquations
          DO idir = 1, 3
          
            DO k = 0, dgStorage % N
              DO j = 0, dgStorage % N
                DO i = 0, dgStorage % N   
 
                  f(1,i,j,k) = nodalSolution % solution(i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,1,iEl)
                  f(2,i,j,k) = nodalSolution % solution(i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,2,iEl)
                  f(3,i,j,k) = nodalSolution % solution(i,j,k,iEq,iEl)*mesh % elements % Ja(i,j,k,idir,3,iEl)

                ENDDO
              ENDDO
            ENDDO

            DO k = 0, dgStorage % N
              DO j = 0, dgStorage % N
                DO i = 0, dgStorage % N   

                  df = 0.0_prec
                  DO ii = 0, dgStorage % N
                    df = df + dgStorage % dgDerivativeMatrixTranspose(ii,i)*f(1,ii,j,k) + &
                              dgStorage % dgDerivativeMatrixTranspose(ii,j)*f(2,i,ii,k) + &
                              dgStorage % dgDerivativeMatrixTranspose(ii,k)*f(3,i,j,ii)
                  ENDDO
                   
                  nodalSolution % solutionGradient(idir,i,j,k,iEq,iEl) =  ( df+ ( nodalSolution % boundaryGradientFlux(idir,i,k,iEq,1,iEl)*dgStorage % boundaryInterpolationMatrix(j,0) + &
                                                                                nodalSolution % boundaryGradientFlux(idir,i,k,iEq,3,iEl)*dgStorage % boundaryInterpolationMatrix(j,1) )/&
                                                                              dgStorage % quadratureWeights(j) + &
                                                                              ( nodalSolution % boundaryGradientFlux(idir,j,k,iEq,4,iEl)*dgStorage % boundaryInterpolationMatrix(i,0) + &
                                                                                nodalSolution % boundaryGradientFlux(idir,j,k,iEq,2,iEl)*dgStorage % boundaryInterpolationMatrix(i,1) )/&
                                                                              dgStorage % quadratureWeights(i) + &
                                                                              ( nodalSolution % boundaryGradientFlux(idir,i,j,iEq,5,iEl)*dgStorage % boundaryInterpolationMatrix(k,0) + &
                                                                                nodalSolution % boundaryGradientFlux(idir,i,j,iEq,6,iEl)*dgStorage % boundaryInterpolationMatrix(k,1) )/&
                                                                              dgStorage % quadratureWeights(k) )/mesh % elements % J(i,j,k,iEl)

                ENDDO  
              ENDDO
            ENDDO
          ENDDO

        ENDDO
      ENDDO
                            

  END SUBROUTINE MappedGradientFluxDivergence

! ================================================================================================ !
! BassiRebay_Gradient_InternalFaceFlux
!
!  Uses the boundarySolution and the mesh boundary normal vectors to build the boundaryGradientFlux
!  attribute on all internal faces of the mesh. At the element faces, the gradient-flux is taken 
!  as the average of the two boundary states that are share the face.
!
! ================================================================================================ !

  SUBROUTINE BassiRebay_Gradient_InternalFaceFlux( nodalSolution, mesh )

    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: nodalSolution
    TYPE( HexMesh ), INTENT(in)                :: mesh
    ! Local
    INTEGER :: iFace, e1, e2, s1, s2, iEq, i, j, ii, jj, idir, jEq
   

    DO iFace = 1, mesh % faces % nFaces

      e1 = mesh % faces % elementIDs(1,iFace)
      e2 = mesh % faces % elementIDs(2,iFace)
      s1 = mesh % faces % elementSides(1,iFace)
      s2 = ABS(mesh % faces % elementSides(2,iFace))

      IF( e2 > 0 )THEN

        DO iEq = 1, nodalSolution % nEquations
          DO j = 0, nodalSolution % N
            DO i = 0, nodalSolution % N
          
              ii = mesh % faces % iMap(i,j,iFace)
              jj = mesh % faces % jMap(i,j,iFace)

              DO idir = 1, 3        
                nodalSolution % boundaryGradientFlux(idir,i,j,iEq,s1,e1) = 0.5_prec*( nodalSolution % boundarySolution(i,j,iEq,s1,e1) +&
                                                                                 nodalSolution % boundarySolution(ii,jj,iEq,s2,e2) )*&
                                                                                 mesh % elements % nHat(idir,i,j,s1,e1)
               
                nodalSolution % boundaryGradientFlux(idir,ii,jj,iEq,s2,e2) = -nodalSolution % boundaryGradientFlux(idir,i,j,iEq,s1,e1)
              ENDDO

            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ENDDO


  END SUBROUTINE BassiRebay_Gradient_InternalFaceFlux

! ================================================================================================ !
! BassiRebay_Gradient_ExternalFaceFlux
!
!  Uses the boundarySolution, externalState, and the mesh boundary normal vectors to build the 
!  boundaryGradientFlux attribute on all internal faces of the mesh. At the element faces, the 
!  gradient-flux is taken as the average of the two boundary states that are share the face.
!
! ================================================================================================ !
  SUBROUTINE BassiRebay_Gradient_ExternalFaceFlux( nodalSolution, mesh )

    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: nodalSolution
    TYPE( HexMesh ), INTENT(in)                :: mesh
    ! Local
    INTEGER :: iFace, e1, e2, s1, iEq, i, j, ii, jj, idir, bID


    DO iFace = 1, mesh % faces % nFaces

      e1 = mesh % faces % elementIDs(1,iFace)
      e2 = mesh % faces % elementIDs(2,iFace)
      s1 = mesh % faces % elementSides(1,iFace)
      bID = ABS(mesh % faces % boundaryID(iFace))

      IF( e2 < 0 )THEN

        IF( bID == 0 )THEN
          DO iEq = 1, nodalSolution % nEquations
            DO j = 0, nodalSolution % N
              DO i = 0, nodalSolution % N
            
                ii = mesh % faces % iMap(i,j,iFace)
                jj = mesh % faces % jMap(i,j,iFace)
  
                DO idir = 1, 3        
  
                  nodalSolution % boundaryGradientFlux(idir,i,j,iEq,s1,e1) =  nodalSolution % boundarySolution(i,j,iEq,s1,e1)*&
                                                                                   mesh % elements % nHat(idir,i,j,s1,e1)
  
                ENDDO
                 
  
              ENDDO
            ENDDO
          ENDDO

        ELSE

          DO iEq = 1, nodalSolution % nEquations
            DO j = 0, nodalSolution % N
              DO i = 0, nodalSolution % N
            
                ii = mesh % faces % iMap(i,j,iFace)
                jj = mesh % faces % jMap(i,j,iFace)
  
                DO idir = 1, 3        
                  nodalSolution % boundaryGradientFlux(idir,i,j,iEq,s1,e1) = 0.5_prec*( nodalSolution % boundarySolution(i,j,iEq,s1,e1) +&
                                                                                   nodalSolution % externalState(ii,jj,iEq,bID) )*&
                                                                                   mesh % elements % nHat(idir,i,j,s1,e1)
                 
                ENDDO
                 
  
              ENDDO
            ENDDO
          ENDDO

        ENDIF

      ENDIF

    ENDDO


  END SUBROUTINE BassiRebay_Gradient_ExternalFaceFlux

#ifdef HAVE_CUDA

  ATTRIBUTES( Global ) SUBROUTINE Mapped_Strong_FluxDivergence_CUDAKernel( flux, Ja, Jac, derivativeMatrixTranspose, fluxDivergence, N, nEq, nEl )

    IMPLICIT NONE
    INTEGER, DEVICE, INTENT(in)     :: N, nEq, nEl
    REAL(prec), DEVICE, INTENT(in)  :: flux(1:3,0:N,0:N,0:N,1:nEq,1:nEl) 
    REAL(prec), DEVICE, INTENT(in)  :: Ja(0:N,0:N,0:N,1:3,1:3,1:nEl) 
    REAL(prec), DEVICE, INTENT(in)  :: Jac(0:N,0:N,0:N,1:nEl) 
    REAL(prec), DEVICE, INTENT(in)  :: derivativeMatrixTranspose(0:N,0:N) 
    REAL(prec), DEVICE, INTENT(out) :: fluxDivergence(0:N,0:N,0:N,1:nEq,1:nEl) 
    ! Local
    INTEGER            :: iEl, iEq, i, j, k, idir, ii
    REAL(prec)         :: df
    REAL(prec), SHARED :: f(1:3,0:7,0:7,0:7)

      iEl = blockIDx % x
      iEq = blockIDx % y

      i = threadIDx % x-1
      j = threadIDx % y-1
      k = threadIDx % z-1

      f(1,i,j,k) = 0.0_prec
      f(2,i,j,k) = 0.0_prec
      f(3,i,j,k) = 0.0_prec

      DO idir = 1, 3

        f(1,i,j,k) = f(1,i,j,k) + flux(idir,i,j,k,iEq,iEl)*Ja(i,j,k,idir,1,iEl)
        f(2,i,j,k) = f(2,i,j,k) + flux(idir,i,j,k,iEq,iEl)*Ja(i,j,k,idir,2,iEl)
        f(3,i,j,k) = f(3,i,j,k) + flux(idir,i,j,k,iEq,iEl)*Ja(i,j,k,idir,3,iEl)

      ENDDO

      CALL syncthreads( )

      df = 0.0_prec
      DO ii = 0, N
        df = df + derivativeMatrixTranspose(ii,i)*f(1,ii,j,k) + &
                  derivativeMatrixTranspose(ii,j)*f(2,i,ii,k) + &
                  derivativeMatrixTranspose(ii,k)*f(3,i,j,ii)
      ENDDO
                   
      fluxDivergence(i,j,k,iEq,iEl) = df/Jac(i,j,k,iEl)


  END SUBROUTINE Mapped_Strong_FluxDivergence_CUDAKernel

! ================================================================================================ !
! UpdateDevice_NodalDGSolution_3D 
!
!   Copies the host-side (CPU memory) attributes to the device-side (GPU memory) attributes. This
!   routine is only available if CUDA is enabled. The memory copy is done on the default stream
!   and is blocking on the GPU.
!
!   Usage : </H2> 
!
!     TYPE</B>(NodalDGSolution_3D) :: nodalSolution 
!
!       CALL</B> nodalSolution % UpdateDevice( ) 
! 
!   Input/Output : </H2>
!
!     nodalSolution  (in/out) 
!       Instance of the \ref NodalDGSolution_3D structure that has been constructed.  
!       On output, the device-side attributes are updated wth the host-side attributes
!
! ================================================================================================ ! 

  SUBROUTINE UpdateDevice_NodalDGSolution_3D( myDGS )
    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: myDGS

      myDGS % solution_dev             = myDGS % solution
      myDGS % boundarySolution_dev     = myDGS % boundarySolution
      myDGS % prescribedState_dev      = myDGS % prescribedState
      myDGS % externalState_dev        = myDGS % externalState
      myDGS % solutionGradient_dev     = myDGS % solutionGradient
      myDGS % boundaryGradientFlux_dev = myDGS % boundaryGradientFlux
      myDGS % flux_dev                 = myDGS % flux
      myDGS % boundaryFlux_dev         = myDGS % boundaryFlux
      myDGS % fluxDivergence_dev       = myDGS % fluxDivergence
      myDGS % source_dev               = myDGS % source
      myDGS % tendency_dev             = myDGS % tendency

  END SUBROUTINE UpdateDevice_NodalDGSolution_3D

! ================================================================================================ !
! UpdateHost_NodalDGSolution_3D 
!
!   Copies the device-side (GPU memory) attributes to the host-side (CPU memory) attributes. This
!   routine is only available if CUDA is enabled. The memory copy is done on the default stream
!   and is blocking on the GPU.
!
!   Usage : </H2> 
!
!     TYPE</B>(NodalDGSolution_3D) :: nodalSolution 
!
!       CALL</B> nodalSolution % UpdateDevice( ) 
! 
!   Input/Output : </H2>
!
!     nodalSolution  (in/out) 
!       Instance of the \ref NodalDGSolution_3D structure that has been constructed.  
!       On output, the host-side attributes are updated wth the device-side attributes
!
! ================================================================================================ ! 

  SUBROUTINE UpdateHost_NodalDGSolution_3D( myDGS )
    IMPLICIT NONE
    CLASS( NodalDGSolution_3D ), INTENT(inout) :: myDGS

      myDGS % solution             = myDGS % solution_dev
      myDGS % boundarySolution     = myDGS % boundarySolution_dev
      myDGS % prescribedState      = myDGS % prescribedState_dev
      myDGS % externalState        = myDGS % externalState_dev
      myDGS % solutionGradient     = myDGS % solutionGradient_dev
      myDGS % boundaryGradientFlux = myDGS % boundaryGradientFlux_dev
      myDGS % flux                 = myDGS % flux_dev
      myDGS % boundaryFlux         = myDGS % boundaryFlux_dev
      myDGS % fluxDivergence       = myDGS % fluxDivergence_dev
      myDGS % source               = myDGS % source_dev
      myDGS % tendency             = myDGS % tendency_dev

  END SUBROUTINE UpdateHost_NodalDGSolution_3D

#endif

END MODULE NodalDGSolution_3D_Class

