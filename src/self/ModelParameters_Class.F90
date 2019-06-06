! ModelParameters_CLASS.f90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE ModelParameters_CLASS

! src/COMMON/
  USE ModelPrecision
  USE CommonRoutines
  USE ConstantsDictionary


  IMPLICIT NONE


  TYPE ModelParameters
    ! TimeManagement
    REAL(prec)    :: dt
    REAL(prec)    :: startTime
    REAL(prec)    :: endTime
    REAL(prec)    :: outputFrequency
    INTEGER       :: nDumps
    INTEGER       :: nStepsPerDump
    REAL(prec)    :: jacobianStepSize
    ! SpaceManagement
    CHARACTER(50) :: SpecMeshFile
    CHARACTER(50) :: SELFMeshFile
    CHARACTER(50) :: UCDMeshFile
    INTEGER       :: MeshTYPE
    INTEGER       :: TopographicShape
    INTEGER       :: QuadTYPE
    INTEGER       :: polyDeg
    INTEGER       :: nXElem
    INTEGER       :: nYElem
    INTEGER       :: nZElem
    INTEGER       :: nProc
    INTEGER       :: nProcX
    INTEGER       :: nProcY
    INTEGER       :: nProcZ
    INTEGER       :: nPlot
    REAL(prec)    :: dxPlot
    REAL(prec)    :: xScale
    REAL(prec)    :: yScale
    REAL(prec)    :: zScale
    ! SubsgridScale
    INTEGER       :: SubGridModel
    LOGICAL       :: SpectralFilter
    INTEGER       :: filterTYPE
    REAL(prec)    :: viscosity
    REAL(prec)    :: viscLengthScale
    INTEGER       :: nCutoff
    REAL(prec)    :: filter_a, filter_b
    ! Physical
    REAL(prec)    :: fRotX  ! coriolis parameter (x-component)
    REAL(prec)    :: fRotY  ! "                " (y-component)
    REAL(prec)    :: fRotZ  ! "                " (z-component)
    REAL(prec)    :: Cd
    REAL(prec)    :: dragscale
    REAL(prec)    :: g  ! gravitational acceleration
    REAL(prec)    :: Cv ! Heat Capacity at constant volume
    REAL(prec)    :: R  ! "Ideal gas constant"
    REAL(prec)    :: T0 ! Reference Temperature
    REAL(prec)    :: dTdz ! Linear temperature stratIFication
    REAL(prec)    :: rho0 ! Reference density
    REAL(prec)    :: P0   ! reference pressure
    REAL(prec)    :: v0
    REAL(prec)    :: hCapRatio
    REAL(prec)    :: rC


  CONTAINS

    PROCEDURE :: Build => Build_ModelParameters

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateDevice => UpdateDevice_ModelParameters
#endif

  END TYPE ModelParameters

  !=================================================================!
  ! ------ SubGridScale Model (Fluid) ---------- !
  !=================================================================!
  INTEGER, PARAMETER :: Laplacian         = 100
  INTEGER, PARAMETER :: SpectralEKE       = 101
  INTEGER, PARAMETER :: SpectralFiltering = 102
  !=================================================================!

#ifdef HAVE_CUDA
    INTEGER, CONSTANT       :: polyDeg_dev
    REAL(prec), CONSTANT    :: dt_dev
    REAL(prec), CONSTANT    :: viscosity_dev
    REAL(prec), CONSTANT    :: viscLengthScale_dev
    ! Physical
    REAL(prec), CONSTANT    :: fRotX_dev  ! coriolis parameter (x-component)
    REAL(prec), CONSTANT    :: fRotY_dev   ! "                " (y-component)
    REAL(prec), CONSTANT    :: fRotZ_dev   ! "                " (z-component)
    REAL(prec), CONSTANT    :: Cd_dev
    REAL(prec), CONSTANT    :: dragscale_dev
    REAL(prec), CONSTANT    :: g_dev   ! gravitational acceleration
    REAL(prec), CONSTANT    :: Cv_dev  ! Heat Capacity at constant volume
    REAL(prec), CONSTANT    :: R_dev   ! "Ideal gas constant"
    REAL(prec), CONSTANT    :: T0_dev  ! Reference Temperature
    REAL(prec), CONSTANT    :: dTdz_dev  ! Linear temperature stratIFication
    REAL(prec), CONSTANT    :: rho0_dev  ! Reference density
    REAL(prec), CONSTANT    :: P0_dev    ! reference pressure
    REAL(prec), CONSTANT    :: v0_dev
    REAL(prec), CONSTANT    :: hCapRatio_dev
    REAL(prec), CONSTANT    :: rC_dev

#endif
CONTAINS


  SUBROUTINE Build_ModelParameters( params, paramFile, readSuccess )

    CLASS( ModelParameters ), INTENT(out) :: params
    CHARACTER(*), INTENT(in)              :: paramFile
    LOGICAL, INTENT(out)                  :: readSuccess
    ! LOCAL
    INTEGER :: nUnit
    ! TimeManagement
    CHARACTER(1)  :: units
    REAL(prec)    :: dt
    REAL(prec)    :: startTime
    REAL(prec)    :: endTime
    REAL(prec)    :: outputFrequency
    INTEGER       :: nStepsPerDump
    REAL(prec)    :: jacobianStepSize
    ! SpaceManagement
    CHARACTER(50) :: SpecMeshFile
    CHARACTER(50) :: SELFMeshFile
    CHARACTER(50) :: UCDMeshFile
    CHARACTER(20) :: MeshTYPE
    CHARACTER(20) :: TopographicShape
    INTEGER       :: QuadTYPE
    INTEGER       :: polyDeg
    INTEGER       :: nXElem
    INTEGER       :: nYElem
    INTEGER       :: nZElem
    INTEGER       :: nProc
    INTEGER       :: nProcX
    INTEGER       :: nProcY
    INTEGER       :: nProcZ
    INTEGER       :: nPlot
    REAL(prec)    :: xScale
    REAL(prec)    :: yScale
    REAL(prec)    :: zScale
    ! SubgridScale
    LOGICAL       :: SpectralFilter
    CHARACTER(20) :: SubGridModel
    CHARACTER(20) :: filterTYPE
    REAL(prec)    :: viscosity
    REAL(prec)    :: viscLengthScale
    INTEGER       :: nCutoff
    REAL(prec)    :: filter_a, filter_b
    ! Physical
    REAL(prec)    :: fRotX
    REAL(prec)    :: fRotY
    REAL(prec)    :: fRotZ
    REAL(prec)    :: Cd
    REAL(prec)    :: dragscale
    REAL(prec)    :: g  ! gravitational acceleration
    REAL(prec)    :: Cv ! Heat Capacity at constant volume
    REAL(prec)    :: R  ! "Ideal gas constant"
    REAL(prec)    :: T0 ! Reference Temperature
    REAL(prec)    :: dTdz
    REAL(prec)    :: rho0 ! Reference Density
    REAL(prec)    :: P0   ! reference pressure
    REAL(prec)    :: v0

    LOGICAL :: fileExists


    NAMELIST / TimeManagement / units, dt, startTime, endTime, outputFrequency, jacobianStepSize
    NAMELIST / SpaceManagement / SpecMeshFile, SELFMeshFile, UCDMeshFile, MeshTYPE, topographicShape, QuadTYPE, polyDeg, &
      nXElem, nYElem, nZElem, nProc, nProcX, nProcY, nProcZ, &
      nPlot, xScale, yScale, zScale
    NAMELIST / SubgridScale / SubGridModel, SpectralFilter, filterTYPE, viscosity, viscLengthScale, nCutoff, filter_a, filter_b
    NAMELIST / PhysicalConstants / fRotX, fRotY, fRotZ, Cd, dragscale, g, Cv, R, T0, dTdz, rho0, P0, v0

    readSuccess = .FALSE.

    ! TimeManagement
    units           = 's'
    dt              = 0.1_prec
    startTime       = 0.0_prec
    endTime         = 1.0_prec
    outputFrequency = 0.5_prec
    jacobianStepSize = 1.0_prec*10.0_prec**(-6)
    ! SpaceManagement
    SpecMeshFile  = nada
    SELFMeshFile = nada
    UCDMeshFile   = nada
    MeshTYPE      = 'Default'
    topographicShape = 'Default'
    QuadTYPE      = GAUSS
    polyDeg       = 5
    nXElem        = 5
    nYElem        = 5
    nZElem        = 5
    nProc         = 1
    nProcX        = 0
    nProcY        = 0
    nProcZ        = 0
    nPlot         = 10
    xScale        = 1.0_prec
    yScale        = 1.0_prec
    zScale        = 1.0_prec
    ! SubgridScale
    SpectralFilter  = .FALSE.
    SubGridModel    = 'Laplacian'
    filterTYPE      = 'TanhRolloff'
    viscosity       = 0.0_prec  ! (m^2/s)
    viscLengthScale = 1.0_prec  ! (m)
    nCutoff = 5
    filter_a = 8.0_prec
    filter_b = 1.0_prec
    ! PhysicalConstants
    fRotX          = 0.0_prec ! 1/s
    fRotY          = 0.0_prec ! 1/s
    fRotZ          = 0.0_prec ! 1/s
    Cd = 0.0_prec
    dragscale = 1.0_prec
    g  = 9.81_prec   ! m/s^2     (Earth)
    Cv = 717.5_prec  ! J/(kg*K)  (Dry air)
    R  = 287.14_prec ! J/(kg*K)  (Dry air)
    T0 = 273.0_prec  ! K         (Room Temperature)
    dTdz = 0.0_prec  ! K/m
    rho0 = 2.0_prec  ! kg/m^3    (Surface density)
    P0   = rho0*R*T0
    v0   = 10.0_prec


    INQUIRE( FILE = TRIM(paramFile), EXIST = fileExists )

    IF( fileExists )THEN
      OPEN( UNIT = NEWUNIT(nUnit), FILE = TRIM(paramFile))
      READ( UNIT = nUnit, NML = TimeManagement )
      READ( UNIT = nUnit, NML = SpaceManagement )
      READ( UNIT = nUnit, NML = SubgridScale )
      READ( UNIT = nUnit, NML = PhysicalConstants )
      CLOSE( UNIT = nUnit )

      ! TimeManagement
      ! Internally, units of seconds are always USEd
      IF( units(1:1) == 'h' )THEN

        params % dt        = dt/3600.0_prec
        params % startTime = startTime/3600.0_prec
        params % endTime   = endTime/3600.0_prec
        params % outputFrequency = outputFrequency/3600.0_prec

      ELSEIF( units(1:1) == 'm' )THEN

        params % dt        = dt/60.0_prec
        params % startTime = startTime/60.0_prec
        params % endTime   = endTime/60.0_prec
        params % outputFrequency = outputFrequency/60.0_prec

      ELSEIF( units(1:1) == 's' )THEN

        params % dt        = dt
        params % startTime = startTime
        params % endTime   = endTime
        params % outputFrequency = outputFrequency

      ELSE

        PRINT(MsgFMT), '  Unknown units specified for time in '//TRIM(paramFile)
        RETURN

      ENDIF

      params % nStepsPerDump = INT( params % outputFrequency/params % dt )
      params % nDumps        = INT( (params % endTime-params % startTime)/params % outputFrequency )


      params % jacobianStepSize = jacobianStepSize
      ! SpaceManagement
      params % SpecMeshFile = SpecMeshFile
      params % SELFMeshFile = SELFMeshFile
      params % UCDMeshFile = UCDMeshFile
      IF( TRIM( UpperCASE( MeshTYPE ) )=='DEFAULT' )THEN
        params % MeshTYPE = DefaultMesh
      ELSEIF( TRIM( UpperCASE( MeshTYPE ) )=='DOUBLYPERIODIC' )THEN
        params % MeshTYPE = DoublyPeriodic
      ELSE
        PRINT*, '   Invalid mesh type : '//UpperCASE(MeshTYPE)
        PRINT*, '   Valid options are "Default" or "DoublyPeriodic"'
      ENDIF

      params % topographicShape = DefaultMesh

      params % QuadTYPE      = QuadTYPE
      params % polyDeg = polyDeg
      params % nXElem = nXElem
      params % nYElem = nYElem
      params % nZElem = nZElem

      IF( nProcX == 0 .AND. nProcY == 0 .AND. nProcZ == 0)THEN
        params % nProc  = nProc
        params % nProcX = nProcX
        params % nProcY = nProcY
        params % nProcZ = nProcZ
      ELSE
        params % nProc  = nProcX*nProcY*nProcZ
        params % nProcX = nProcX
        params % nProcY = nProcY
        params % nProcZ = nProcZ
      ENDIF

      params % nPlot = nPlot
      params % dxPlot = 2.0_prec/REAL(nPlot,prec)
      params % xScale = xScale
      params % yScale = yScale
      params % zScale = zScale

      ! SubgridScale
      IF( TRIM( UpperCASE( SubGridModel) )=='LAPLACIAN' )THEN
        params % SubGridModel = Laplacian ! convert to the INTEGER flag
      ELSEIF( TRIM(UpperCASE( SubGridModel) )=='SPECTRALEKE' )THEN
        params % SubGridModel = SpectralEKE
      ELSEIF( TRIM(UpperCASE( SubGridModel) )=='SPECTRALFILTERING' )THEN
        params % SubGridModel = SpectralFiltering
      ELSE
        PRINT*, '   Invalid subgrid-scale model : '//UpperCASE(SubGridModel)
        PRINT*, '   Valid options are "Laplacian", "SpectralEKE", or "SpectralFiltering"'
        RETURN
      ENDIF

      IF( TRIM( UpperCASE( FilterTYPE ) )=='TANHROLLOFF' )THEN
        params % filterTYPE = tanhRollOff
      ELSEIF( TRIM( UpperCASE( FilterTYPE ) )=='MODALCUTOFF' )THEN
        params % filterTYPE = modalCutoff
      ELSEIF( TRIM( UpperCASE( FilterTYPE ) )=='RAMPFILTER' )THEN
        params % filterTYPE = rampFilter
      ENDIF

      !params % subGridModel = subGridModel
      params % SpectralFilter  = SpectralFilter
      params % viscosity       = viscosity
      params % viscLengthScale = viscLengthScale
      params % nCutoff         = nCutoff
      params % filter_a        = filter_a
      params % filter_b        = filter_b
      ! PhysicalConstants
      params % fRotX  = fRotX
      params % fRotY  = fRotY
      params % fRotZ  = fRotZ
      params % Cd     = Cd
      params % dragScale = dragScale
      params % g  = g
      params % Cv = Cv
      params % R  = R
      params % T0 = T0
      params % dTdz = dTdz
      params % rho0 = rho0
      params % P0   = rho0*R*T0
      params % v0   = v0

      readSuccess = .TRUE.

    ELSE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = TRIM(paramFile), ACTION = 'WRITE' )
      WRITE( UNIT = nUnit, NML = TimeManagement )
      WRITE( UNIT = nUnit, NML = SpaceManagement )
      WRITE( UNIT = nUnit, NML = SubgridScale )
      WRITE( UNIT = nUnit, NML = PhysicalConstants )
      CLOSE( UNIT = nUnit )

      PRINT(MsgFMT), TRIM(paramFile)//'not found.'
      PRINT(MsgFMT), 'A sample namelist file has been generated   '
      PRINT(MsgFMT), 'for you in your current directory.'
      RETURN 

    ENDIF

    params % hCapRatio = ( params % R + params % Cv ) / params % Cv
    params % rC        =   params % R / ( params % R + params % Cv )

#ifdef HAVE_CUDA
    CALL params % UpdateDevice( )
#endif


  END SUBROUTINE Build_ModelParameters


#ifdef HAVE_CUDA
  SUBROUTINE UpdateDevice_ModelParameters( params )
    IMPLICIT NONE
    CLASS( ModelParameters ), INTENT(inout) :: params

    polyDeg_dev         = params % polyDeg
    dt_dev              = params % dt
    viscosity_dev       = params % viscosity
    viscLengthScale_dev = params % viscLengthScale

    fRotX_dev     = params % fRotX
    fRotY_dev     = params % fRotY
    fRotZ_dev     = params % fRotZ
    Cd_dev        = params % Cd
    dragScale_dev = params % dragScale
    g_dev         = params % g
    Cv_dev        = params % Cv
    R_dev         = params % R
    T0_dev        = params % T0
    dTdz_dev      = params % dTdz
    rho0_dev      = params % rho0
    P0_dev        = params % P0
    v0_dev        = params % v0

    hCapRatio_dev = params % hCapRatio
    rC_dev        = params % rC

  END SUBROUTINE UpdateDevice_ModelParameters
#endif

END MODULE ModelParameters_CLASS
