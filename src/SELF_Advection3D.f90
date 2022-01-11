MODULE SELF_Advection3D
!
! Copyright 2021-2022 Fluid Numerics LLC
!
! ============================================================
!
! A class for solving the advection-diffusion equation in 3-D
!
! Environment variables that influence behavior
!
!  SELF_INSTALL_DIR - the path to the SELF installation
!
! ============================================================
USE SELF_Constants
USE SELF_SupportRoutines
USE SELF_Mesh
USE SELF_DG
USE FEQParse
USE FLAP
USE SELF_CLI

! Needed for Fortran-C interoperability
! Helps expose HIP kernels to Fortran
USE ISO_C_BINDING


  TYPE,EXTENDS(DG3D), PUBLIC :: Advection3D

    TYPE(MappedVector3D),PUBLIC :: velocity
    TYPE(Vector3D),PUBLIC :: plotVelocity
    TYPE(Vector3D),PUBLIC :: plotX

    TYPE(EquationParser), ALLOCATABLE :: boundaryConditionEqn(:)
    TYPE(EquationParser), ALLOCATABLE :: solutionEqn(:)
    TYPE(EquationParser), ALLOCATABLE :: sourceEqn(:)

    REAL(prec) :: simulationTime

    ! Model Settings !
    REAL(prec) :: Lx, Ly, Lz ! Domain lengths
    REAL(prec) :: dt ! Default time step size
    REAL(prec) :: outputInterval
    REAL(prec) :: initialTime
    REAL(prec) :: endTime
    REAL(prec) :: diffusivity
    LOGICAL :: diffusiveFlux
    INTEGER :: controlDegree
    INTEGER :: targetDegree
    INTEGER :: controlQuadrature ! ENUMS in SELF_Constants.f90
    INTEGER :: targetQuadrature ! ENUMS in SELF_Constants.f90
    CHARACTER(LEN=self_FileNameLength) :: icFile
    CHARACTER(LEN=self_FileNameLength) :: meshFile
    INTEGER :: nxElements
    INTEGER :: nyElements
    INTEGER :: nzElements
    INTEGER :: integrator ! ENUMS needed in SELF_Constants.f90 !! TO DO !!
    CHARACTER(LEN=self_EquationLength) :: velEqnX ! Velocity Equation (x-direction)
    CHARACTER(LEN=self_EquationLength) :: velEqnY ! Velocity Equation (y-direction)
    CHARACTER(LEN=self_EquationLength) :: velEqnZ ! Velocity Equation (z-direction)
    CHARACTER(LEN=self_EquationLength) :: icEqn ! Initial condition Equation
    CHARACTER(LEN=self_EquationLength) :: bcEqn ! Boundary condition Equation
    LOGICAL :: enableMPI
    LOGICAL :: gpuAccel
    

    CONTAINS

      PROCEDURE,PUBLIC :: Init => Init_Advection3D

      PROCEDURE,PUBLIC :: InitFromCLI => InitFromCLI_Advection3D

      PROCEDURE,PUBLIC :: Free => Free_Advection3D


      GENERIC, PUBLIC :: SetSolution => SetSolutionFromEquation_Advection3D
      PROCEDURE, PRIVATE :: SetSolutionFromEquation_Advection3D

      GENERIC, PUBLIC :: SetSource => SetSourceFromEquation_Advection3D
      PROCEDURE, PRIVATE :: SetSourceFromEquation_Advection3D

      GENERIC, PUBLIC :: SetVelocity => SetVelocityFromEquation_Advection3D
      PROCEDURE, PRIVATE :: SetVelocityFromEquation_Advection3D

      GENERIC, PUBLIC :: SetBoundaryCondition => SetBoundaryConditionFromEquation_Advection3D
      PROCEDURE, PRIVATE :: SetBoundaryConditionFromEquation_Advection3D

      PROCEDURE, PUBLIC :: WriteTecplot => WriteTecplot_Advection3D
      PROCEDURE, PUBLIC :: WritePickup => WritePickup_Advection3D

      PROCEDURE, PUBLIC :: ForwardStep => ForwardStep_Advection3D
      PROCEDURE, PUBLIC :: TimeStepRK3 => TimeStepRK3_Advection3D

      PROCEDURE, PUBLIC :: Tendency => Tendency_Advection3D
      PROCEDURE, PUBLIC :: InternalFlux => InternalFlux_Advection3D
      PROCEDURE, PUBLIC :: SideFlux => SideFlux_Advection3D

  END TYPE Advection3D

!  PRIVATE :: GetCLIParameters

  ! Interfaces to GPU kernels !
  INTERFACE
    SUBROUTINE InternalFlux_Advection3D_gpu_wrapper(flux, solution, velocity, dsdx, N, nVar, nEl) &
      BIND(c,name="InternalFlux_Advection3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, solution, velocity, dsdx
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE InternalFlux_Advection3D_gpu_wrapper 
  END INTERFACE

  INTERFACE
    SUBROUTINE InternalDiffusiveFlux_Advection3D_gpu_wrapper(flux, solutionGradient, dsdx, diffusivity, N, nVar, nEl) &
      BIND(c,name="InternalDiffusiveFlux_Advection3D_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, solutionGradient, dsdx
      REAL(c_prec), VALUE :: diffusivity
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE InternalDiffusiveFlux_Advection3D_gpu_wrapper 
  END INTERFACE

  INTERFACE
    SUBROUTINE SideFlux_Advection3D_gpu_wrapper(flux, boundarySol, extSol, velocity, nHat, nScale, N, nVar, nEl) &
      BIND(c,name="SideFlux_Advection3D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, boundarySol, extSol, velocity, nHat, nScale
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SideFlux_Advection3D_gpu_wrapper 
  END INTERFACE

  INTERFACE
    SUBROUTINE SideDiffusiveFlux_Advection3D_gpu_wrapper(flux, boundarySolGradient, extSolGradient, &
                                                         nHat, nScale, diffusivity, N, nVar, nEl) &
      BIND(c,name="SideDiffusiveFlux_Advection3D_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, boundarySolGradient, extSolGradient, nHat, nScale
      REAL(c_prec), VALUE :: diffusivity
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SideDiffusiveFlux_Advection3D_gpu_wrapper 
  END INTERFACE

  INTERFACE
    SUBROUTINE UpdateGRK3_Advection3D_gpu_wrapper(gRK3, solution, dSdt, rk3A, rk3G, dt, N, nVar, nEl) &
      BIND(c,name="UpdateGRK3_Advection3D_gpu_wrapper")
      USE ISO_C_BINDING
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: gRK3, solution, dSdt
      REAL(c_prec),VALUE :: rk3A, rk3G, dt
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE UpdateGRK3_Advection3D_gpu_wrapper 
  END INTERFACE

CONTAINS

  SUBROUTINE Init_Advection3D(this,cqType,tqType,cqDegree,tqDegree,nvar,enableMPI,spec)
    IMPLICIT NONE
    CLASS(Advection3D),INTENT(out) :: this
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nvar
    LOGICAL,INTENT(in) :: enableMPI
    TYPE(MeshSpec),INTENT(in) :: spec

    CALL this % decomp % Init(enableMPI)

    ! Load Mesh
    IF (enableMPI)THEN
      CALL this % mesh % Load(spec,this % decomp)
    ELSE
      CALL this % mesh % Load(spec)
      CALL this % decomp % setElemToRank(this % mesh % nGlobalElem)
    ENDIF

    CALL this % decomp % SetMaxMsg(this % mesh % nUniqueSides)


    ! Create geometry from mesh
    CALL this % geometry % GenerateFromMesh(&
            this % mesh,cqType,tqType,cqDegree,tqDegree)

    CALL this % plotSolution % Init(&
            tqDegree,tqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % dSdt % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % solution % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % solutionGradient % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % flux % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % velocity % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % plotVelocity % Init(&
            tqDegree,tqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % plotX % Init(&
            tqDegree,tqType,tqDegree,tqType,1,&
            this % mesh % nElem)

    CALL this % source % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
                this % mesh % nElem)

    CALL this % fluxDivergence % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % workScalar % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % workVector % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % workTensor % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % compFlux % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    ALLOCATE (this % solutionMetaData(1:nvar))
    ALLOCATE (this % boundaryConditionEqn(1:nvar))
    ALLOCATE (this % solutionEqn(1:nvar))
    ALLOCATE (this % sourceEqn(1:nvar))


  END SUBROUTINE Init_Advection3D

  SUBROUTINE InitFromCLI_Advection3D(this)
    IMPLICIT NONE
    CLASS(Advection3D),INTENT(inout) :: this
    ! Local
    TYPE(CLI) :: cliConf
    TYPE(MeshSpec) :: spec
    CHARACTER(self_QuadratureTypeCharLength) :: cqTypeChar
    CHARACTER(self_QuadratureTypeCharLength) :: tqTypeChar
    CHARACTER(self_IntegratorTypeCharLength) :: integratorChar
    REAL(prec) :: Lx, Ly, Lz ! Domain lengths
    REAL(prec) :: dt ! Default time step size
    REAL(prec) :: diffusivity
    INTEGER :: controlDegree
    INTEGER :: targetDegree
    INTEGER :: controlQuadrature ! ENUMS in SELF_Constants.f90
    INTEGER :: targetQuadrature ! ENUMS in SELF_Constants.f90
    CHARACTER(LEN=self_FileNameLength) :: meshFile
    INTEGER :: nxElements
    INTEGER :: nyElements
    INTEGER :: nzElements
    INTEGER :: integrator ! ENUMS needed in SELF_Constants.f90 !! TO DO !!
    CHARACTER(LEN=self_EquationLength) :: velEqnX ! Velocity Equation (x-direction)
    CHARACTER(LEN=self_EquationLength) :: velEqnY ! Velocity Equation (y-direction)
    CHARACTER(LEN=self_EquationLength) :: velEqnZ ! Velocity Equation (z-direction)
    CHARACTER(LEN=self_EquationLength) :: icEqn ! Initial condition Equation
    CHARACTER(LEN=self_EquationLength) :: bcEqn ! Boundary condition Equation
    CHARACTER(LEN=self_EquationLength) :: sourceEqn ! Boundary condition Equation
    CHARACTER(LEN=self_FileNameLength) :: selfInstallDir
    LOGICAL :: enableMPI
    LOGICAL :: enableGPU
    LOGICAL :: diffusiveFlux
    REAL(prec) :: outputInterval
    REAL(prec) :: initialTime
    REAL(prec) :: endTime
    TYPE(EquationParser) :: eqn(1)
    TYPE(EquationParser) :: velEqn(1:3)
    LOGICAL :: fileExists

 
    CALL GET_ENVIRONMENT_VARIABLE('SELF_INSTALL_DIR', selfInstallDir)
    IF( TRIM(selfInstallDir) == "" )THEN
      PRINT*, "SELF_INSTALL_DIR environment variable unset."
      PRINT*, "Trying /opt/view/"
      selfInstallDir = "/opt/view"
    ENDIF

    INQUIRE(FILE=TRIM(selfInstallDir)//'/etc/sadv3d.json', EXIST=fileExists)
    IF( .NOT. fileExists )THEN
      PRINT*, TRIM(selfInstallDir)//'/etc/sadv3d.json'//' configuration file not found!'
      STOP "ERROR"
    ENDIF

    CALL cliConf % Init(TRIM(selfInstallDir)//'/etc/sadv3d.json')
    CALL cliConf % LoadFromCLI() 

    ! TO DO : Swap out for cliConf generic "get" based on
    !  cliConf % opts json object
    !
    ! Set the CLI parameters !
    CALL cliConf % cliObj % get(val=enableMPI,switch='--mpi')
    CALL cliConf % cliObj % get(val=enableGPU,switch='--gpu')
    CALL cliConf % cliObj % get(val=meshfile,switch='--mesh')
    CALL cliConf % cliObj % get(val=dt,switch="--time-step")
    CALL cliConf % cliObj % get(val=outputInterval,switch="--output-interval")
    CALL cliConf % cliObj % get(val=initialTime,switch="--initial-time")
    CALL cliConf % cliObj % get(val=endTime,switch="--end-time")
    CALL cliConf % cliObj % get(val=controlDegree,switch="--control-degree")
    CALL cliConf % cliObj % get(val=targetDegree,switch="--target-degree")
    CALL cliConf % cliObj % get(val=cqTypeChar,switch="--control-quadrature")
    CALL cliConf % cliObj % get(val=tqTypeChar,switch="--target-quadrature")
    CALL cliConf % cliObj % get(val=meshFile,switch="--mesh")
    CALL cliConf % cliObj % get(val=nxElements,switch="--nxelements")
    CALL cliConf % cliObj % get(val=nyElements,switch="--nyelements")
    CALL cliConf % cliObj % get(val=nzElements,switch="--nzelements")
    CALL cliConf % cliObj % get(val=Lx, switch="--xlength")
    CALL cliConf % cliObj % get(val=Ly, switch="--ylength")
    CALL cliConf % cliObj % get(val=Lz, switch="--zlength")
    CALL cliConf % cliObj % get(val=velEqnX,switch="--velocity-x")
    CALL cliConf % cliObj % get(val=velEqnY,switch="--velocity-y")
    CALL cliConf % cliObj % get(val=velEqnZ,switch="--velocity-z")
    CALL cliConf % cliObj % get(val=icEqn,switch="--initial-condition")
    CALL cliConf % cliObj % get(val=bcEqn,switch="--boundary-condition")
    CALL cliConf % cliObj % get(val=sourceEqn,switch="--source")
    CALL cliConf % cliObj % get(val=integratorChar,switch="--integrator")
    CALL cliConf % cliObj % get(val=diffusivity,switch="--diffusivity")

    diffusiveFlux = .TRUE.
    IF( diffusivity == 0.0_prec ) THEN
      diffusiveFlux = .FALSE.
    ELSEIF( diffusivity < 0.0_prec ) THEN
      IF( dt > 0.0_prec )THEN
        PRINT*, 'Negative diffusivity provably unstable for forward stepping'
        PRINT*, 'Invalid diffusivity value. Stopping'
        STOP
      ENDIF
    ENDIF

    IF (TRIM(UpperCase(cqTypeChar)) == 'GAUSS') THEN
      controlQuadrature = GAUSS
    ELSEIF (TRIM(UpperCase(cqTypeChar)) == 'GAUSS-LOBATTO') THEN
      controlQuadrature = GAUSS_LOBATTO
    ELSEIF (TRIM(UpperCase(cqTypeChar)) == 'CHEBYSHEV-GAUSS') THEN
      controlQuadrature = CHEBYSHEV_GAUSS
    ELSEIF (TRIM(UpperCase(cqTypeChar)) == 'CHEBYSHEV-GAUSS-LOBATTO') THEN
      controlQuadrature = CHEBYSHEV_GAUSS_LOBATTO
    ELSE
      PRINT *, 'Invalid Control Quadrature'
      STOP - 1
    END IF

    IF (TRIM(UpperCase(tqTypeChar)) == 'UNIFORM') THEN
      targetQuadrature = UNIFORM
    ELSEIF (TRIM(UpperCase(tqTypeChar)) == 'GAUSS') THEN
      targetQuadrature = GAUSS
    ELSEIF (TRIM(UpperCase(tqTypeChar)) == 'GAUSS-LOBATTO') THEN
      targetQuadrature = GAUSS_LOBATTO
    ELSEIF (TRIM(UpperCase(tqTypeChar)) == 'CHEBYSHEV-GAUSS') THEN
      targetQuadrature = CHEBYSHEV_GAUSS
    ELSEIF (TRIM(UpperCase(tqTypeChar)) == 'CHEBYSHEV-GAUSS-LOBATTO') THEN
      targetQuadrature = CHEBYSHEV_GAUSS_LOBATTO
    ELSE
      PRINT *, 'Invalid Target Quadrature'
      STOP - 1
    END IF

    IF (TRIM(UpperCase(integratorChar)) == 'EULER') THEN
      integrator = EULER
    ELSEIF (TRIM(UpperCase(integratorChar)) == 'RK3') THEN
      integrator = RK3
    ELSE
      PRINT *, 'Invalid time integrator'
      STOP - 1
    END IF

    IF( TRIM(meshfile) == "" )THEN
      spec % blockMesh = .TRUE.
    ELSE
      spec % blockMesh = .FALSE.     
    ENDIF
    spec % filename = meshfile
    spec % filetype = SELF_MESH_ISM_V2_3D

    spec % blockMesh_nGeo = 1
    spec % blockMesh_x0 = 0.0_prec
    spec % blockMesh_x1 = Lx
    spec % blockMesh_y0 = 0.0_prec
    spec % blockMesh_y1 = Ly
    spec % blockMesh_z0 = 0.0_prec
    spec % blockMesh_z1 = Lz
    spec % blockMesh_nElemX = nxElements
    spec % blockMesh_nElemY = nyElements
    spec % blockMesh_nElemZ = nzElements

    CALL this % Init(controlQuadrature, &
                     targetQuadrature, &
                     controlDegree, &
                     targetDegree, &
                     1,enableMPI, &
                     spec)

    this % simulationTime = 0.0_prec
    this % Lx = Lx
    this % Ly = Ly ! Domain lengths
    this % Lz = Lz ! Domain lengths
    this % dt = dt ! Default time step size
    this % initialTime = initialTime
    this % simulationTime = initialTime
    this % endTime = endTime
    this % outputInterval = outputInterval
    this % controlDegree = controlDegree
    this % targetDegree = targetDegree
    this % controlQuadrature = controlQuadrature ! ENUMS in SELF_Constants.f90
    this % targetQuadrature = targetQuadrature ! ENUMS in SELF_Constants.f90
    this % meshFile = meshFile
    this % nxElements = nxElements
    this % nyElements = nyElements
    this % nzElements = nzElements
    this % integrator = integrator ! ENUMS needed in SELF_Constants.f90 !! TO DO !!
    this % velEqnX = velEqnX ! Velocity Equation (x-direction)
    this % velEqnY = velEqnY ! Velocity Equation (y-direction)
    this % velEqnZ = velEqnZ ! Velocity Equation (y-direction)
    this % icEqn = icEqn ! Initial condition Equation
    this % bcEqn = bcEqn ! Boundary condition Equation
    this % enableMPI = enableMPI
    this % gpuAccel = enableGPU
    this % diffusivity = diffusivity
    this % diffusiveFlux = diffusiveFlux        

    eqn(1) = EquationParser( icEqn, (/'x','y','z','t'/))
    CALL this % SetSolution( eqn )
    this % solutionEqn(1) = EquationParser( icEqn, (/'x','y','z','t'/))

    velEqn(1) = EquationParser(velEqnX, (/'x','y','z'/))
    velEqn(2) = EquationParser(velEqnY, (/'x','y','z'/))
    velEqn(3) = EquationParser(velEqnZ, (/'x','y','z'/))
    CALL this % SetVelocity( velEqn )

    this % boundaryConditionEqn(1) = EquationParser( bcEqn, (/'x','y','z','t'/))
    CALL this % SetBoundaryCondition( this % boundaryConditionEqn )

    this % sourceEqn(1) = EquationParser( sourceEqn, (/'x','y','z','t'/))

  END SUBROUTINE InitFromCLI_Advection3D

!  SUBROUTINE GetCLIParameters( cli )
!    TYPE(COMMAND_LINE_INTERFACE), INTENT(inout) :: cli
!
!    CALL cli % init(progname="sadv3d", &
!                    version="v0.0.0", &
!                    description="SELF Advection in 3-D", &
!                    license="ANTI-CAPITALIST SOFTWARE LICENSE (v 1.4)", &
!                    authors="Joseph Schoonover (Fluid Numerics LLC)")
!
!    CALL cli % add(switch="--mpi", &
!                   help="Enable MPI", &
!                   act="store_true", &
!                   def="false", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--gpu", &
!                   help="Enable GPU acceleration", &
!                   act="store_true", &
!                   def="false", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--time-step", &
!                   switch_ab="-dt", &
!                   help="The time step size for the time integrator", &
!                   def="0.001", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--initial-time", &
!                   switch_ab="-t0", &
!                   help="The initial time level", &
!                   def="0.0", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--output-interval", &
!                   switch_ab="-oi", &
!                   help="The time between file output", &
!                   def="0.5", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--end-time", &
!                   switch_ab="-tn", &
!                   help="The final time level", &
!                   def="1.0", &
!                   required=.FALSE.)
!
!    ! Get the control degree
!    CALL cli % add(switch="--control-degree", &
!                   switch_ab="-c", &
!                   help="The polynomial degree of the control points."//NEW_LINE("A"), &
!                   def="7", &
!                   required=.FALSE.)
!
!    ! Get the target degree (assumed for plotting)
!    CALL cli % add(switch="--target-degree", &
!                   switch_ab="-t", &
!                   help="The polynomial degree for the"//&
!                  &" target points for interpolation."//&
!                  &" Typically used for plotting"//NEW_LINE("A"), &
!                   def="14", &
!                   required=.FALSE.)
!
!    ! Get the control quadrature
!    ! Everyone know Legendre-Gauss Quadrature is the best...
!    CALL cli % add(switch="--control-quadrature", &
!                   switch_ab="-cq", &
!                   def="gauss", &
!                   help="The quadrature type for control points."//NEW_LINE("A"), &
!                   choices="gauss,gauss-lobatto,chebyshev-gauss,chebyshev-gauss-lobatto", &
!                   required=.FALSE.)
!
!
!    ! Set the target grid quadrature
!    ! Default to uniform (assumed for plotting)
!    CALL cli % add(switch="--target-quadrature", &
!                   switch_ab="-tq", &
!                   def="uniform", &
!                   help="The quadrature type for target points."//NEW_LINE("A"), &
!                   choices="gauss,gauss-lobatto,uniform", &
!                   required=.FALSE.)
!
!    ! (Optional) Provide a file for a mesh
!    ! Assumed in HOPR or ISM-v2 format
!    CALL cli % add(switch="--mesh", &
!                   switch_ab="-m", &
!                   help="Path to a mesh file for control mesh."//NEW_LINE("A"), &
!                   def="", &
!                   required=.FALSE.)
!
!    ! (Optional) If a mesh is not provided, you
!    ! can request a structured grid to be generated
!    ! just set the nxelement, nyelements..
!    CALL cli % add(switch="--nxelements", &
!                   switch_ab="-nx", &
!                   help="The number of elements in the x-direction for structured mesh generation.", &
!                   def="5", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--nyelements", &
!                   switch_ab="-ny", &
!                   help="The number of elements in the y-direction for structured mesh generation.", &
!                   def="5", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--nzelements", &
!                   switch_ab="-nz", &
!                   help="The number of elements in the z-direction for structured mesh generation.", &
!                   def="5", &
!                   required=.FALSE.)
!
!    ! Alright... now tell me some physical mesh dimensions
!    CALL cli % add(switch="--xlength", &
!                   switch_ab="-lx", &
!                   help="The physical x-scale for structured mesh generation."//&
!                   " Ignored if a mesh file is provided", &
!                   def="1.0", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--ylength", &
!                   switch_ab="-ly", &
!                   help="The physical y-scale for structured mesh generation."//&
!                   " Ignored if a mesh file is provided", &
!                   def="1.0", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--zlength", &
!                   switch_ab="-lz", &
!                   help="The physical z-scale for structured mesh generation."//&
!                   " Ignored if a mesh file is provided", &
!                   def="1.0", &
!                   required=.FALSE.)
!
!    ! Set the velocity field
!    CALL cli % add(switch="--velocity-x", &
!                   switch_ab="-vx", &
!                   help="Equation for the x-component of the velocity field",&
!                   def="vx=1.0", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--velocity-y", &
!                   switch_ab="-vy", &
!                   help="Equation for the y-component of the velocity field",&
!                   def="vy=1.0", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--velocity-z", &
!                   switch_ab="-vz", &
!                   help="Equation for the z-component of the velocity field",&
!                   def="vz=1.0", &
!                   required=.FALSE.)
!
!    ! Tracer diffusivity
!    CALL cli % add(switch="--diffusivity", &
!                   switch_ab="-nu", &
!                   help="Tracer diffusivity (applied to all tracers)", &
!                   def="0.0", &
!                   required=.FALSE.)
!
!    ! Set the initial conditions
!    ! .. TO DO .. 
!    !  > How to handle multiple tracer fields ??
!    CALL cli % add(switch="--initial-condition", &
!                   switch_ab="-ic", &
!                   help="Equation for the initial tracer distributions",&
!                   def="f = exp( -( ((x-t)-0.5)^2 + ((y-t)-0.5)^2 + ((z-t)-0.5)^2)/0.01 )", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--boundary-condition", &
!                   switch_ab="-bc", &
!                   help="Equation for the boundary tracer distributions (can be time dependent!)", &
!                   def="f = exp( -( ((x-t)-0.5)^2 + ((y-t)-0.5)^2 + ((z-t)-0.5)^2)/0.01 )", &
!                   required=.FALSE.)
!
!    CALL cli % add(switch="--source", &
!                   switch_ab="-s", &
!                   help="Equation for the source term (can be time dependent!)", &
!                   def="s = 0.0", &
!                   required=.FALSE.)
!
!    ! Give me a time integrator
!    CALL cli % add(switch="--integrator", &
!                   switch_ab="-int", &
!                   help="Sets the time integration method. Only 'euler' or 'williamson_rk3'", &
!                   def="williamson_rk3", &
!                   required=.FALSE.)
!
!  END SUBROUTINE GetCLIParameters

  SUBROUTINE Free_Advection3D(this)
    IMPLICIT NONE
    CLASS(Advection3D),INTENT(inout) :: this

    CALL this % mesh % Free()
    CALL this % decomp % Free()
    CALL this % geometry % Free()
    CALL this % solution % Free()
    CALL this % dSdt % Free()
    CALL this % plotSolution % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()
    CALL this % workScalar % Free()
    CALL this % workVector % Free()
    CALL this % workTensor % Free()
    CALL this % compFlux % Free()
    CALL this % velocity % Free()
    CALL this % plotVelocity % Free()
    CALL this % plotX % Free()
    DEALLOCATE (this % solutionMetaData)
    DEALLOCATE (this % boundaryConditionEqn)
    DEALLOCATE (this % solutionEqn)
    DEALLOCATE (this % sourceEqn)
    CALL this % decomp % Finalize()

  END SUBROUTINE Free_Advection3D

  SUBROUTINE SetSolutionFromEquation_Advection3D( this, eqn )
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    TYPE(EquationParser), INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: i, j, k, iEl, iVar
    REAL(prec) :: x
    REAL(prec) :: y
    REAL(prec) :: z
    REAL(prec) :: t


    DO iEl = 1,this % solution % nElem
      DO iVar = 1, this % solution % nVar
        DO k = 0, this % solution % N
          DO j = 0, this % solution % N
            DO i = 0, this % solution % N

               ! Get the mesh positions
               x = this % geometry % x % interior % hostData(1,i,j,k,1,iEl)
               y = this % geometry % x % interior % hostData(2,i,j,k,1,iEl)
               z = this % geometry % x % interior % hostData(3,i,j,k,1,iEl)
               t = this % simulationTime

               this % solution % interior % hostData(i,j,k,iVar,iEl) = &
                 eqn(iVar) % Evaluate((/x, y, z, t/))


            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF( this % gpuAccel )THEN
      CALL this % solution % interior % UpdateDevice()
    ENDIF

  END SUBROUTINE SetSolutionFromEquation_Advection3D

  SUBROUTINE SetSourceFromEquation_Advection3D( this, eqn )
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    TYPE(EquationParser), INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: i, j, k, iEl, iVar
    REAL(prec) :: x
    REAL(prec) :: y
    REAL(prec) :: z
    REAL(prec) :: t


    DO iEl = 1,this % source % nElem
      DO iVar = 1, this % source % nVar
        DO k = 0, this % source % N
          DO j = 0, this % source % N
            DO i = 0, this % source % N

               ! Get the mesh positions
               x = this % geometry % x % interior % hostData(1,i,j,k,1,iEl)
               y = this % geometry % x % interior % hostData(2,i,j,k,1,iEl)
               z = this % geometry % x % interior % hostData(3,i,j,k,1,iEl)
               t = this % simulationTime

               this % source % interior % hostData(i,j,k,iVar,iEl) = &
                 eqn(iVar) % Evaluate((/x, y, z, t/))


            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF( this % gpuAccel )THEN
      CALL this % source % interior % UpdateDevice()
    ENDIF

  END SUBROUTINE SetSourceFromEquation_Advection3D

  SUBROUTINE SetVelocityFromEquation_Advection3D( this, eqn )
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    TYPE(EquationParser), INTENT(in) :: eqn(1:3)
    ! Local
    INTEGER :: i, j, k, iEl, iVar, iSide
    REAL(prec) :: x
    REAL(prec) :: y
    REAL(prec) :: z


    DO iEl = 1,this % solution % nElem

      ! Set the velocity at the element interiors
      DO k = 0, this % solution % N
        DO j = 0, this % solution % N
          DO i = 0, this % solution % N

             ! Get the mesh positions
             x = this % geometry % x % interior % hostData(1,i,j,k,1,iEl)
             y = this % geometry % x % interior % hostData(2,i,j,k,1,iEl)
             z = this % geometry % x % interior % hostData(3,i,j,k,1,iEl)

             ! Set the velocity in the x-direction
             this % velocity % interior % hostData(1,i,j,k,1,iEl) = &
               eqn(1) % Evaluate((/x, y, z/))

             ! Set the velocity in the y-direction
             this % velocity % interior % hostData(2,i,j,k,1,iEl) = &
               eqn(2) % Evaluate((/x, y, z/))

             ! Set the velocity in the z-direction
             this % velocity % interior % hostData(3,i,j,k,1,iEl) = &
               eqn(3) % Evaluate((/x, y, z/))

          ENDDO
        ENDDO
      ENDDO

      ! Set the velocity at element faces
      DO iSide = 1, 6
        DO j = 0, this % solution % N
          DO i = 0, this % solution % N

             ! Get the mesh positions
             x = this % geometry % x % boundary % hostData(1,i,j,1,iSide,iEl)
             y = this % geometry % x % boundary % hostData(2,i,j,1,iSide,iEl)
             z = this % geometry % x % boundary % hostData(3,i,j,1,iSide,iEl)

             ! Set the velocity in the x-direction
             this % velocity % boundary % hostData(1,i,j,1,iSide,iEl) = &
               eqn(1) % Evaluate((/x, y, z/))

             ! Set the velocity in the y-direction
             this % velocity % boundary % hostData(2,i,j,1,iSide,iEl) = &
               eqn(2) % Evaluate((/x, y, z/))

             ! Set the velocity in the z-direction
             this % velocity % boundary % hostData(3,i,j,1,iSide,iEl) = &
               eqn(3) % Evaluate((/x, y, z/))


          ENDDO
        ENDDO
      ENDDO

    ENDDO

    IF( this % gpuAccel )THEN
      CALL this % velocity % interior % UpdateDevice()
      CALL this % velocity % boundary % UpdateDevice()
    ENDIF

  END SUBROUTINE SetVelocityFromEquation_Advection3D

  SUBROUTINE SetBoundaryConditionFromEquation_Advection3D( this, eqn )
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    TYPE(EquationParser), INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: i, j, iEl, iVar, iSide
    REAL(prec) :: x
    REAL(prec) :: y
    REAL(prec) :: z


    DO iEl = 1,this % solution % nElem
      DO iSide = 1, 6
        DO iVar = 1, this % solution % nvar
          DO j = 0, this % solution % N
            DO i = 0, this % solution % N

               ! If this element's side has no neighbor assigned
               ! it is assumed to be a physical boundary.
               ! In this case, we want to assign the external boundary
               ! condition.
               IF( this % mesh % self_sideInfo % hostData(3,iSide,iEl) == 0 )THEN
                 ! Get the mesh positions
                 x = this % geometry % x % boundary % hostData(1,i,j,1,iSide,iEl)
                 y = this % geometry % x % boundary % hostData(2,i,j,1,iSide,iEl)
                 z = this % geometry % x % boundary % hostData(3,i,j,1,iSide,iEl)

                 ! Set the external boundary condition
                 this % solution % extBoundary % hostData(i,j,iVar,iSide,iEl) = &
                   eqn(iVar) % Evaluate((/x, y, z, this % simulationTime/))
               ENDIF


            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF( this % gpuAccel )THEN
      ! Copy data to the GPU
      CALL this % solution % extBoundary % UpdateDevice()
    ENDIF

  END SUBROUTINE SetBoundaryConditionFromEquation_Advection3D

  SUBROUTINE ForwardStep_Advection3D( this, endTime )
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    REAL(prec), INTENT(in) :: endTime
    ! Local
    INTEGER :: nSteps
    REAL(prec) :: dt
    REAL(prec) :: t1, t2

    IF( this % integrator == RK3 )THEN
    
      ! Step forward
      dt = this % dt
      nSteps = INT(( endTime - this % simulationTime )/dt)
      CALL CPU_TIME(t1)
      CALL this % TimeStepRK3( nSteps )
      CALL CPU_TIME(t2)
      PRINT*, nSteps, ' steps took ', (t2-t1), ' seconds'

      ! Take any additional steps to reach desired endTime
      this % dt = endTime - this % simulationTime
      IF( this % dt > 0 )THEN
        nSteps = 1
        CALL this % TimeStepRK3( nSteps )
      ENDIF

      ! Reset the time step
      this % dt = dt

    ENDIF

  END SUBROUTINE ForwardStep_Advection3D

  SUBROUTINE TimeStepRK3_Advection3D( this, nSteps )
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    INTEGER, INTENT(in) :: nSteps
    ! Local
    INTEGER :: m, iStep
    INTEGER :: iEl
    INTEGER :: iVar
    INTEGER :: i, j, k
    TYPE(hfReal_r5) :: gRK3
    REAL(prec) :: t0
    REAL(prec) :: dt
    REAL(prec) :: rk3A
    REAL(prec) :: rk3G
   
      CALL gRK3 % Alloc(loBound=(/0,0,0,1,1/), &
                        upBound=(/this % solution % N,&
                                  this % solution % N,&
                                  this % solution % N,&
                                  this % solution % nVar, &
                                  this % solution % nElem/))

      dt = this % dt

      DO iStep = 1, nSteps

        t0 = this % simulationTime

        gRK3 % hostData = 0.0_prec
        DO m = 1, 3 ! Loop over RK3 steps

          CALL this % Tendency( )

          IF( this % gpuAccel )THEN

            rk3A = rk3_a(m)
            rk3G = rk3_g(m)

            CALL UpdateGRK3_Advection3D_gpu_wrapper( gRK3 % deviceData, &
                             this % solution % interior % deviceData, &
                             this % dSdt % interior % deviceData, &
                             rk3A, rk3G, dt, &
                             this % solution % N, &
                             this % solution % nVar, &
                             this % solution % nElem )
          ELSE


            DO iEl = 1, this % solution % nElem
              DO iVar = 1, this % solution % nVar
                DO k = 0, this % solution % N
                  DO j = 0, this % solution % N
                    DO i = 0, this % solution % N

                      gRK3 % hostData(i,j,k,iVar,iEl) = rk3_a(m)*gRK3 % hostData(i,j,k,iVar,iEl) + &
                              this % dSdt % interior % hostData(i,j,k,iVar,iEl)


                      this % solution % interior % hostData(i,j,k,iVar,iEl) = &
                              this % solution % interior % hostData(i,j,k,iVar,iEl) + &
                              rk3_g(m)*dt*gRK3 % hostData(i,j,k,iVar,iEl)

                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDIF

          this % simulationTime = this % simulationTime + rk3_b(m)*dt

        ENDDO

        this % simulationTime = t0 + dt

      ENDDO

      CALL gRK3 % Free()

  END SUBROUTINE TimeStepRK3_Advection3D

  SUBROUTINE Tendency_Advection3D( this ) 
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this

      CALL this % solution % BoundaryInterp( this % gpuAccel )
      
      IF (this % diffusiveFlux) THEN
        CALL this % CalculateSolutionGradient( this % gpuAccel )
      ENDIF

      ! Internal Flux calculates both the advective and diffusive flux -- need diffusivity 
      CALL this % InternalFlux( )

      ! Exchange side information between neighboring cells
      CALL this % solution % SideExchange( this % mesh, &
                                           this % decomp, &
                                           this % gpuAccel )

      IF (this % diffusiveFlux) THEN
        CALL this % solutionGradient % SideExchange( this % mesh, &
                                                     this % decomp, &
                                                     this % gpuAccel )
   
      ENDIF

      CALL this % SideFlux( )

      CALL this % CalculateFluxDivergence( this % gpuAccel )

      CALL this % CalculatedSdt( this % gpuAccel )

  END SUBROUTINE Tendency_Advection3D

  SUBROUTINE SideFlux_Advection3D( this )
    !! Calculates the Advective Flux on element sides using a Lax-Friedrich's upwind Riemann Solver
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iSide,iVar,iEl
    REAL(prec) :: nhat(1:3)
    REAL(prec) :: nmag
    REAL(prec) :: un
    REAL(prec) :: extState
    REAL(prec) :: intState


      IF( this % gpuAccel )THEN

        CALL SideFlux_Advection3D_gpu_wrapper( this % flux % boundaryNormal % deviceData, &
                                               this % solution % boundary % deviceData, &
                                               this % solution % extBoundary % deviceData, &
                                               this % velocity % boundary % deviceData, &
                                               this % geometry % nHat % boundary % deviceData, &
                                               this % geometry % nScale % boundary % deviceData, &
                                               this % solution % N, &
                                               this % solution % nVar, &
                                               this % solution % nElem )
        IF (this % diffusiveFlux) THEN
          CALL SideDiffusiveFlux_Advection3D_gpu_wrapper( this % flux % boundaryNormal % deviceData, &
                                                 this % solutionGradient % boundary % deviceData, &
                                                 this % solutionGradient % extBoundary % deviceData, &
                                                 this % geometry % nHat % boundary % deviceData, &
                                                 this % geometry % nScale % boundary % deviceData, &
                                                 this % diffusivity, &
                                                 this % solution % N, &
                                                 this % solution % nVar, &
                                                 this % solution % nElem )
        ENDIF

      ELSE

        DO iEl = 1, this % solution % nElem
          DO iSide = 1, 6
            DO iVar = 1, this % solution % nVar
              DO j = 0, this % solution % N
                DO i = 0, this % solution % N

                   ! Get the boundary normals on cell edges from the mesh geometry
                   nhat(1:3) = this % geometry % nHat % boundary % hostData(1:3,i,j,1,iSide,iEl)

                   ! Calculate the normal velocity at the cell edges
                   un = this % velocity % boundary % hostData(1,i,j,1,iSide,iEl)*nHat(1)+&
                        this % velocity % boundary % hostData(2,i,j,1,iSide,iEl)*nHat(2)+&
                        this % velocity % boundary % hostData(3,i,j,1,iSide,iEl)*nHat(3)

                   ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
                   extState = this % solution % extBoundary % hostData(i,j,iVar,iSide,iEl)
                   intState = this % solution % boundary % hostData(i,j,iVar,iSide,iEl)
                   nmag = this % geometry % nScale % boundary % hostData(i,j,1,iSide,iEl)

                   ! Calculate the flux
                   this % flux % boundaryNormal % hostData(i,j,iVar,iSide,iEl) = 0.5_prec*&
                       ( un*(intState + extState) - abs(un)*(extState - intState) )*nmag

                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        IF (this % diffusiveFlux) THEN
          DO iEl = 1, this % solution % nElem
            DO iSide = 1, 6
              DO iVar = 1, this % solution % nVar
                DO j = 0, this % solution % N
                  DO i = 0, this % solution % N

                    nhat(1:3) = this % geometry % nHat % boundary % hostData(1:3,i,j,1,iSide,iEl)
                    nmag = this % geometry % nScale % boundary % hostData(i,j,1,iSide,iEl)

                    !  Calculate \nabla{f} \cdot \hat{n} on the cell sides
                    extState = this % solutionGradient % extBoundary % hostData(1,i,j,iVar,iSide,iEl)*nHat(1)+&
                               this % solutionGradient % extBoundary % hostData(2,i,j,iVar,iSide,iEl)*nHat(2)+&
                               this % solutionGradient % extBoundary % hostData(3,i,j,iVar,iSide,iEl)*nHat(3)

                    intState = this % solutionGradient % boundary % hostData(1,i,j,iVar,iSide,iEl)*nHat(1)+&
                               this % solutionGradient % boundary % hostData(2,i,j,iVar,iSide,iEl)*nHat(2)+&
                               this % solutionGradient % boundary % hostData(3,i,j,iVar,iSide,iEl)*nHat(3)

                    ! Bassi-Rebay flux is the average of the internal and external diffusive flux vectors.
                    this % flux % boundaryNormal % hostData(i,j,iVar,iSide,iEl) = &
                      this % flux % boundaryNormal % hostData(i,j,iVar,iSide,iEl) -&
                      0.5_prec*this % diffusivity*(extState + intState)*nmag

                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! Diffusivity

      ENDIF ! GPU Acceleration

  END SUBROUTINE SideFlux_Advection3D

  SUBROUTINE InternalFlux_Advection3D( this )
    !! Calculates the advective flux using the provided velocity
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,k,iVar,iEl
    REAL(prec) :: Fx, Fy, Fz

    IF( this % gpuAccel )THEN

      ! When GPU acceleration is enabled (requested by the user)
      ! we call the gpu wrapper interface, which will call the
      ! HIP kernel "under the hood"
      ! 
      ! TO DO : Pass the contravariant basis vector to GPU kernel
      CALL InternalFlux_Advection3D_gpu_wrapper(this % flux % interior % deviceData,&
                                                this % solution % interior % deviceData, &
                                                this % velocity % interior % deviceData, &
                                                this % geometry % dsdx % interior % deviceData, &
                                                this % solution % N, & 
                                                this % solution % nVar, &
                                                this % solution % nElem )

      IF (this % diffusiveFlux) THEN
        CALL InternalDiffusiveFlux_Advection3D_gpu_wrapper(this % flux % interior % deviceData,&
                                                  this % solutionGradient % interior % deviceData, &
                                                  this % geometry % dsdx % interior % deviceData, &
                                                  this % diffusivity, &
                                                  this % solution % N, & 
                                                  this % solution % nVar, &
                                                  this % solution % nElem )
      ENDIF

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO k = 0, this % solution % N
            DO j = 0, this % solution % N
              DO i = 0, this % solution % N

                Fx = this % velocity % interior % hostData(1,i,j,k,1,iEl)*&
                     this % solution % interior % hostData(i,j,k,iVar,iEl)

                Fy = this % velocity % interior % hostData(2,i,j,k,1,iEl)*&
                     this % solution % interior % hostData(i,j,k,iVar,iEl)

                Fz = this % velocity % interior % hostData(3,i,j,k,1,iEl)*&
                     this % solution % interior % hostData(i,j,k,iVar,iEl)

                this % flux % interior % hostData(1,i,j,k,iVar,iEl) = &
                  this % geometry % dsdx % interior % hostData(1,1,i,j,k,1,iel)*Fx + &
                  this % geometry % dsdx % interior % hostData(2,1,i,j,k,1,iel)*Fy + &
                  this % geometry % dsdx % interior % hostData(3,1,i,j,k,1,iel)*Fz

                this % flux % interior % hostData(2,i,j,k,iVar,iEl) = &
                  this % geometry % dsdx % interior % hostData(1,2,i,j,k,1,iel)*Fx + &
                  this % geometry % dsdx % interior % hostData(2,2,i,j,k,1,iel)*Fy + &
                  this % geometry % dsdx % interior % hostData(3,2,i,j,k,1,iel)*Fz

                this % flux % interior % hostData(3,i,j,k,iVar,iEl) = &
                  this % geometry % dsdx % interior % hostData(1,3,i,j,k,1,iel)*Fx + &
                  this % geometry % dsdx % interior % hostData(2,3,i,j,k,1,iel)*Fy + &
                  this % geometry % dsdx % interior % hostData(3,3,i,j,k,1,iel)*Fz


              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! When diffusivity == 0, then we don't bother calculating the diffusive flux
      IF (this % diffusiveFlux) THEN
        ! Otherwise, we add the diffusive flux to to the flux vector

        DO iEl = 1,this % solution % nElem
          DO iVar = 1, this % solution % nVar
            DO k = 0, this % solution % N
              DO j = 0, this % solution % N
                DO i = 0, this % solution % N

                  ! Diffusive flux is diffusivity coefficient mulitplied by 
                  ! solution gradient
                  Fx = this % solutionGradient % interior % hostData(1,i,j,k,iVar,iEl)*&
                       this % diffusivity

                  Fy = this % solutionGradient % interior % hostData(2,i,j,k,iVar,iEl)*&
                       this % diffusivity

                  Fz = this % solutionGradient % interior % hostData(3,i,j,k,iVar,iEl)*&
                       this % diffusivity

                  ! Project the diffusive flux vector onto computational coordinates
                  this % flux % interior % hostData(1,i,j,k,iVar,iEl) = &
                    this % flux % interior % hostData(1,i,j,k,iVar,iEl) - &
                    ( this % geometry % dsdx % interior % hostData(1,1,i,j,k,1,iel)*Fx + &
                      this % geometry % dsdx % interior % hostData(2,1,i,j,k,1,iel)*Fy + & 
                      this % geometry % dsdx % interior % hostData(3,1,i,j,k,1,iel)*Fz) 

                  this % flux % interior % hostData(2,i,j,k,iVar,iEl) = &
                    this % flux % interior % hostData(2,i,j,k,iVar,iEl) - &
                    ( this % geometry % dsdx % interior % hostData(1,2,i,j,k,1,iel)*Fx + &
                      this % geometry % dsdx % interior % hostData(2,2,i,j,k,1,iel)*Fy + & 
                      this % geometry % dsdx % interior % hostData(3,2,i,j,k,1,iel)*Fz)

                  this % flux % interior % hostData(3,i,j,k,iVar,iEl) = &
                    this % flux % interior % hostData(3,i,j,k,iVar,iEl) - &
                    ( this % geometry % dsdx % interior % hostData(1,3,i,j,k,1,iel)*Fx + &
                      this % geometry % dsdx % interior % hostData(2,3,i,j,k,1,iel)*Fy + & 
                      this % geometry % dsdx % interior % hostData(3,3,i,j,k,1,iel)*Fz )
  
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF   ! DiffusiveFlux

    ENDIF ! GPU Acceleration

  END SUBROUTINE InternalFlux_Advection3D

  SUBROUTINE WriteTecplot_Advection3D(self, filename)
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: self
    CHARACTER(*), INTENT(in), OPTIONAL :: filename
    ! Local
    CHARACTER(8) :: zoneID
    INTEGER :: fUnit
    INTEGER :: iEl, i, j, k 
    CHARACTER(LEN=self_FileNameLength) :: tecFile
    CHARACTER(13) :: timeStampString
    CHARACTER(5) :: rankString

    IF( PRESENT(filename) )THEN
      tecFile = filename
    ELSE
      timeStampString = TimeStamp(self % simulationTime, 's')

      IF( self % decomp % mpiEnabled )THEN
        WRITE(rankString,'(I5.5)') self % decomp % rankId 
        tecFile = 'solution.'//rankString//'.'//timeStampString//'.tec'
      ELSE
        tecFile = 'solution.'//timeStampString//'.tec'
      ENDIF

    ENDIF
                      
    IF( self % gpuAccel )THEN
      ! Copy data to the CPU
      CALL self % solution % interior % UpdateHost()
    ENDIF

    ! Map the mesh positions to the target grid
    CALL self % geometry % x % GridInterp(self % plotX, gpuAccel=.FALSE.)

    ! Map the solution to the target grid
    CALL self % solution % GridInterp(self % plotSolution,gpuAccel=.FALSE.)

    ! Map the velocity to the target grid 
    CALL self % velocity % GridInterp(self % plotVelocity,gpuAccel=.FALSE.)
   
    ! Let's write some tecplot!! 
     OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= TRIM(tecFile), &
      FORM='formatted', &
      STATUS='replace')

    ! TO DO :: Adjust for multiple tracer fields
    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "tracer","u","v","w"'

    DO iEl = 1, self % solution % nElem

      ! TO DO :: Get the global element ID 
      WRITE(zoneID,'(I8.8)') iEl
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',self % solution % M+1,&
                                                 ', J=',self % solution % M+1,&
                                                 ', K=',self % solution % M+1,',F=POINT'

      DO k = 0, self % solution % M
        DO j = 0, self % solution % M
          DO i = 0, self % solution % M

            WRITE(fUnit,'(7(E15.7,1x))') self % plotX % interior % hostData(1,i,j,k,1,iEl), &
                                         self % plotX % interior % hostData(2,i,j,k,1,iEl), &
                                         self % plotX % interior % hostData(3,i,j,k,1,iEl), &
                                         self % plotSolution % interior % hostData(i,j,k,1,iEl),&
                                         self % plotVelocity % interior % hostData(1,i,j,k,1,iEl),&
                                         self % plotVelocity % interior % hostData(2,i,j,k,1,iEl),&
                                         self % plotVelocity % interior % hostData(3,i,j,k,1,iEl)

          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)

  END SUBROUTINE WriteTecplot_Advection3D

  SUBROUTINE WritePickup_Advection3D(self, filename)
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: self
    CHARACTER(*), INTENT(in), OPTIONAL :: filename
    ! Local
    CHARACTER(LEN=self_FileNameLength) :: pickupFile
    CHARACTER(13) :: timeStampString

    IF( PRESENT(filename) )THEN
      pickupFile = filename
    ELSE
      timeStampString = TimeStamp(self % simulationTime, 's')
      pickupFile = 'solution.'//timeStampString//'.h5'

    ENDIF

    IF( self % gpuAccel )THEN
      CALL self % solution % UpdateHost()
    ENDIF


    CALL self % Write(pickupFile)

  END SUBROUTINE WritePickup_Advection3D

END MODULE SELF_Advection3D
