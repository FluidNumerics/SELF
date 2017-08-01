! FluidParams_Class.f90
! 
! Copyright 2017 Joseph Schoonover <schoonover.numerics@gmail.com>
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
MODULE FluidParams_Class

! src/common/
USE ModelPrecision
USE CommonRoutines
USE ConstantsDictionary
 

 IMPLICIT NONE


    TYPE FluidParams
      ! TimeManagement
      REAL(prec)    :: dt
      REAL(prec)    :: Time
      INTEGER       :: iterInit
      INTEGER       :: nTimeSteps
      INTEGER       :: dumpFreq
      ! SpaceManagement
      CHARACTER(50) :: SpecMeshFile
      CHARACTER(50) :: PeaceMeshFile
      CHARACTER(50) :: UCDMeshFile
      INTEGER       :: MeshType
      INTEGER       :: TopographicShape
      INTEGER       :: QuadType
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
      REAL(prec)    :: viscosity
      REAL(prec)    :: viscLengthScale
      INTEGER       :: nCutoff
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
      REAL(prec)    :: dTdz ! Linear temperature stratification
      REAL(prec)    :: dTdx ! Linear temperature lateral gradient (x-component, for initial condition, NOT background)
      REAL(prec)    :: dTdy ! Linear temperature lateral gradient (y-component, for initial condition, NOT background)
      REAL(prec)    :: rho0 ! Reference density
      REAL(prec)    :: P0   ! reference pressure
      REAL(prec)    :: v0

      CONTAINS

      PROCEDURE :: Build => Build_FluidParams

    END TYPE FluidParams 

  !=================================================================!
  ! ------ SubGridScale Model (Fluid) ---------- !
  !=================================================================!
   INTEGER, PARAMETER :: Laplacian         = 100
   INTEGER, PARAMETER :: SpectralEKE       = 101
   INTEGER, PARAMETER :: SpectralFiltering = 102
  !=================================================================!
  
 CONTAINS


 SUBROUTINE Build_FluidParams( thisParam )

   CLASS( FluidParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
      ! TimeManagement
      REAL(prec)    :: dt
      REAL(prec)    :: Time
      INTEGER       :: iterInit
      INTEGER       :: nTimeSteps
      INTEGER       :: dumpFreq
      ! SpaceManagement
      CHARACTER(50) :: SpecMeshFile
      CHARACTER(50) :: PeaceMeshFile
      CHARACTER(50) :: UCDMeshFile
      CHARACTER(20) :: MeshType
      CHARACTER(20) :: TopographicShape
      INTEGER       :: QuadType
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
      CHARACTER(20) :: SubGridModel
      REAL(prec)    :: viscosity
      REAL(prec)    :: viscLengthScale
      INTEGER       :: nCutoff
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
      REAL(prec)    :: dTdx ! Linear temperature lateral gradient (x-component, for initial condition, NOT background)
      REAL(prec)    :: dTdy ! Linear temperature lateral gradient (y-component, for initial condition, NOT background)
      REAL(prec)    :: rho0 ! Reference Density
      REAL(prec)    :: P0   ! reference pressure
      REAL(prec)    :: v0
       
       
      NAMELIST / TimeManagement / dt, Time, iterInit, nTimeSteps, dumpFreq
      NAMELIST / SpaceManagement / SpecMeshFile, PeaceMeshFile, UCDMeshFile, MeshType, topographicShape, QuadType, polyDeg, &
                                    nXElem, nYElem, nZElem, nProc, nProcX, nProcY, nProcZ, &
                                   nPlot, xScale, yScale, zScale
      NAMELIST / SubgridScale / SubGridModel, viscosity, viscLengthScale, nCutoff
      NAMELIST / PhysicalConstants / fRotX, fRotY, fRotZ, Cd, dragscale, g, Cv, R, T0, dTdz, dTdx, dTdy, rho0, P0, v0
      
      ! TimeManagement
      dt = 1.0_prec
      Time = -1.0_prec
      iterInit = 0
      nTimeSteps = 0
      dumpFreq = 1
      ! SpaceManagement
      SpecMeshFile  = nada
      PeaceMeshFile = nada
      UCDMeshFile   = nada
      MeshType      = 'Default'
      topographicShape = 'Default'
      QuadType      = GAUSS
      polyDeg       = 5
      nXElem        = 5
      nYElem        = 5
      nZElem        = 5
      nProc         = 1
      nProcX        = 0
      nProcY        = 0
      nProcZ        = 0
      nPlot         = 10
      xScale        = ONE
      yScale        = ONE 
      zScale        = ONE
      ! SubgridScale
      SubGridModel    = 'Laplacian'
      viscosity       = 0.0_prec  ! (m^2/s)
      viscLengthScale = 1.0_prec  ! (m)
      nCutoff = 5
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
      dTdx = 0.0_prec  ! K/m
      dTdy = 0.0_prec  ! K/m
      rho0 = 2.0_prec  ! kg/m^3    (Surface density)
      P0   = rho0*R*T0
      v0   = 10.0_prec
      
       
      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = TimeManagement )
         READ( UNIT = nUnit, NML = SpaceManagement )
         READ( UNIT = nUnit, NML = SubgridScale )
         READ( UNIT = nUnit, NML = PhysicalConstants )
      CLOSE( UNIT = nUnit ) 

      ! TimeManagement
      thisParam % dt = dt
      thisParam % iterInit = iterInit
      IF( Time <= 0.0_prec )THEN
         thisParam % nTimeSteps = nTimeSteps
      ELSE
         nTimeSteps             = Time/dt
         thisParam % nTimeSteps = Time/dt
      ENDIF
      thisParam % dumpFreq = dumpFreq
      ! SpaceManagement 
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % PeaceMeshFile = PeaceMeshFile
      thisParam % UCDMeshFile = UCDMeshFile
      IF( TRIM( UpperCase( MeshType ) )=='DEFAULT' )THEN
         thisParam % MeshType = DefaultMesh
      ELSEIF( TRIM( UpperCase( MeshType ) )=='DOUBLYPERIODIC' )THEN
         thisParam % MeshType = DoublyPeriodic
      ELSE
         PRINT*, 'Module FluidParams_Class.f90 : S/R Build '
         PRINT*, '   Invalid MeshType : '//UpperCase(MeshType)
         PRINT*, '   Valid options are "Default" or "DoublyPeriodic"'
         STOP 'STOPPING!'
      ENDIF
      
      IF( TRIM( UpperCase( topographicShape ) ) =='DEFAULT' )THEN
         thisParam % topographicShape = DefaultMesh
      ELSEIF( TRIM( UpperCase( topographicShape ) ) =='GAUSSIANHILL' )THEN
         thisParam % topographicShape = Gaussian
      ENDIF
      
      thisParam % QuadType      = QuadType
      thisParam % polyDeg = polyDeg
      thisParam % nXElem = nXElem
      thisParam % nYElem = nYElem 
      thisParam % nZElem = nZElem 
      IF( nProcX == 0 .AND. nProcY == 0 .AND. nProcZ == 0)THEN
         thisParam % nProc  = nProc
         thisParam % nProcX = nProcX
         thisParam % nProcY = nProcY
         thisParam % nProcZ = nProcZ
      ELSE
         thisParam % nProc  = nProcX*nProcY*nProcZ
         thisParam % nProcX = nProcX
         thisParam % nProcY = nProcY
         thisParam % nProcZ = nProcZ
      ENDIF
      thisParam % nPlot = nPlot
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale = xScale
      thisParam % yScale = yScale
      thisParam % zScale = zScale
      ! SubgridScale
      IF( TRIM( UpperCase( SubGridModel) )=='LAPLACIAN' )THEN
      thisParam % SubGridModel = Laplacian ! convert to the integer flag
      ELSEIF( TRIM(UpperCase( SubGridModel) )=='SPECTRALEKE' )THEN
         thisParam % SubGridModel = SpectralEKE
      ELSEIF( TRIM(UpperCase( SubGridModel) )=='SPECTRALFILTERING' )THEN
         thisParam % SubGridModel = SpectralFiltering
      ELSE
         PRINT*, 'Module FluidParams_Class.f90 : S/R Build '
         PRINT*, '   Invalid SubGridModel : '//UpperCase(SubGridModel)
         PRINT*, '   Valid options are "Laplacian", "SpectralEKE", or "SpectralFiltering"'
         STOP 'STOPPING!'
      ENDIF
      !thisParam % subGridModel = subGridModel
      thisParam % viscosity       = viscosity
      thisParam % viscLengthScale = viscLengthScale
      thisParam % nCutoff         = nCutoff
      ! PhysicalConstants
      thisParam % fRotX  = fRotX
      thisParam % fRotY  = fRotY
      thisParam % fRotZ  = fRotZ
      thisParam % Cd     = Cd
      thisParam % dragScale = dragScale
      thisParam % g  = g
      thisParam % Cv = Cv
      thisParam % R  = R
      thisParam % T0 = T0
      thisParam % dTdz = dTdz
      thisParam % dTdx = dTdx
      thisParam % dTdy = dTdy
      thisParam % rho0 = rho0
      thisParam % P0   = P0
      thisParam % v0   = v0
      
      
            ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = TimeManagement )
      WRITE( UNIT = *, NML = SpaceManagement )
      WRITE( UNIT = *, NML = SubgridScale )
      WRITE( UNIT = *, NML = PhysicalConstants )
      
 END SUBROUTINE Build_FluidParams

END MODULE FluidParams_Class
