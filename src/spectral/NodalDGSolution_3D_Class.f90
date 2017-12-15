! NodalDGSolution_3D_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE NodalDGSolution_3D_Class

USE ModelPrecision
USE NodalDG_Class

IMPLICIT NONE

! ================================================================================================ !
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
!  Attributes
!    
!     INTEGER  ::  N  <td> Polynomial degree of the spectral element method
!
!     INTEGER  ::  nEquations <td> Number of solution variables.
!
!     INTEGER  ::  nElements <td> Number of elements that make up the spatial domain.
!
!     REAL(prec)  ::  solution(0:N,0:N,0:N,1:nEquations,1:nElements)
!
!     REAL(prec)  ::  flux(1:3,0:N,0:N,0:N,1:nEquations,1:nElements)
!
!     REAL(prec)  ::  boundarySolution(0:N,0:N,1:6,1:nEquations)*
!
!     REAL(prec)  ::  boundaryFlux(0:N,0:N,1:6,1:nEquations)*
!
!     REAL(prec)  ::  fluxDivergence(0:N,0:N,0:N,1:nEquations,1:nElements)
!
!     REAL(prec)  ::  source(0:N,0:N,0:N,1:nEquations,1:nElements)
!
!     REAL(prec)  ::  tendency(0:N,0:N,0:N,1:nEquations)
!
!
!     >> *For the "boundary" arrays, the faces of a hexahedral element are numbered as SOUTH=1, 
!         EAST=2, NORTH=3, WEST=4, BOTTOM=5, TOP=6.
!
!  Procedures
!
!      nodalSolution % Build( N, nEquations, nElements ) 
!      nodalSolution % Trash( )
!      nodalSolution % Calculate_Solution_at_Boundaries( nodalDGStorage )
!      nodalSolution % Calculate_Weak_Flux_Divergence( nodalDGStorage )
!      nodalSolution % Calculate_Strong_Flux_Divergence( nodalDGStorage ) 
!
!      >> nodalDGStorage is an instance of the NodalDG class.
!
!      >> If CUDA is enabled, additional routines are provided for updating the host-side arrays
!         (CPU) and the device-side arrays (GPU)
!
!      nodalSolution % UpdateDevice( )
!      nodalSolution % UpdateHost( )
!
! ================================================================================================ !

  TYPE NodalDGSolution_3D

    INTEGER                 :: N, nEquations, nElements
    REAL(prec), ALLOCATABLE :: solution(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: flux(:,:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: boundarySolution(:,:,:,:,:) 
    REAL(prec), ALLOCATABLE :: boundaryFlux(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: fluxDivergence(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: source(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: tendency(:,:,:,:,:)

#ifdef HAVE_CUDA

    INTEGER, DEVICE, ALLOCATABLE    :: N_dev, nEquations_dev, nElements_dev
    REAL(prec), DEVICE, ALLOCATABLE :: solution_dev(:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: flux_dev(:,:,:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: boundarySolution_dev(:,:,:,:,:) 
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

  SUBROUTINE Build_NodalDGSolution_3D( myDGS, N, nEquations, nElements )
    IMPLICIT NONE
    CLASS(NodalDGSolution_3D), INTENT(inout) :: myDGS
    INTEGER, INTENT(in)                      :: N, nEquations, nElements
      
      myDGS % N          = N
      myDGS % nEquations = nEquations
      myDGS % nElements  = nElements

      ALLOCATE( myDGS % solution(0:N,0:N,0:N,1:nEquations,1:nElements), &
                myDGS % flux(1:3,0:N,0:N,0:N,1:nEquations,1:nElements), &
                myDGS % boundarySolution(0:N,0:N,1:nEquations,1:6,1:nElements), &
                myDGS % boundaryFlux(0:N,0:N,1:nEquations,1:6,1:nElements), &
                myDGS % fluxDivergence(0:N,0:N,0:N,1:nEquations,1:nElements), &
                myDGS % source(0:N,0:N,0:N,1:nEquations,1:nElements), &
                myDGS % tendency(0:N,0:N,0:N,1:nEquations,1:nElements) )
      
      myDGS % solution         = 0.0_prec
      myDGS % flux             = 0.0_prec
      myDGS % boundarySolution = 0.0_prec
      myDGS % boundaryFlux     = 0.0_prec
      myDGS % fluxDivergence   = 0.0_prec
      myDGS % source           = 0.0_prec
      myDGS % tendency         = 0.0_prec

#ifdef HAVE_CUDA
      ALLOCATE( myDGS % N_dev, myDGS % nEquations_dev, myDGS % nElements_dev )
      myDGS % N_dev          = N
      myDGS % nEquations_dev = nEquations
      myDGS % nElements_dev  = nElements
      
      ALLOCATE( myDGS % solution_dev(0:N,0:N,0:N,1:nEquations,1:nElements), &
                myDGS % flux_dev(1:3,0:N,0:N,0:N,1:nEquations,1:nElements), &
                myDGS % boundarySolution_dev(0:N,0:N,1:nEquations,1:6,1:nElements), &
                myDGS % boundaryFlux_dev(0:N,0:N,1:nEquations,1:6,1:nElements), &
                myDGS % fluxDivergence_dev(0:N,0:N,0:N,1:nEquations,1:nElements), &
                myDGS % source_dev(0:N,0:N,0:N,1:nEquations,1:nElements), &
                myDGS % tendency_dev(0:N,0:N,0:N,1:nEquations,1:nElements) )

      myDGS % solution_dev         = 0.0_prec
      myDGS % flux_dev             = 0.0_prec
      myDGS % boundarySolution_dev = 0.0_prec
      myDGS % boundaryFlux_dev     = 0.0_prec
      myDGS % fluxDivergence_dev   = 0.0_prec
      myDGS % source_dev           = 0.0_prec
      myDGS % tendency_dev         = 0.0_prec

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
                  myDGS % flux, &
                  myDGS % boundarySolution, &
                  myDGS % boundaryFlux, &
                  myDGS % fluxDivergence, &
                  myDGS % source, &
                  myDGS % tendency )


#ifdef HAVE_CUDA
      DEALLOCATE( myDGS % N_dev, myDGS % nEquations_dev, myDGS % nElements_dev )
      DEALLOCATE( myDGS % solution_dev, &
                  myDGS % flux_dev, &
                  myDGS % boundarySolution_dev, &
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

#ifdef HAVE_CUDA

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

      myDGS % solution_dev         = myDGS % solution
      myDGS % flux_dev             = myDGS % flux
      myDGS % boundarySolution_dev = myDGS % boundarySolution
      myDGS % boundaryFlux_dev     = myDGS % boundaryFlux
      myDGS % fluxDivergence_dev   = myDGS % fluxDivergence
      myDGS % source_dev           = myDGS % source
      myDGS % tendency_dev         = myDGS % tendency

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

      myDGS % solution         = myDGS % solution_dev
      myDGS % flux             = myDGS % flux_dev
      myDGS % boundarySolution = myDGS % boundarySolution_dev
      myDGS % boundaryFlux     = myDGS % boundaryFlux_dev
      myDGS % fluxDivergence   = myDGS % fluxDivergence_dev
      myDGS % source           = myDGS % source_dev
      myDGS % tendency         = myDGS % tendency_dev

  END SUBROUTINE UpdateHost_NodalDGSolution_3D

#endif

END MODULE NodalDGSolution_3D_Class
