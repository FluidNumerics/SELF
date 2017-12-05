! FluidParams_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
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
      REAL(prec)    :: startTime
      REAL(prec)    :: endTime
      REAL(prec)    :: outputFrequency
      INTEGER       :: nDumps
      INTEGER       :: nStepsPerDump
      REAL(prec)    :: jacobianStepSize
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


 SUBROUTINE Build_FluidParams( thisParam, readSuccess )

   CLASS( FluidParams ), intent(out) :: thisParam
   LOGICAL, INTENT(out)              :: readSuccess
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
      
      LOGICAL :: fileExists
       
       
      NAMELIST / TimeManagement / units, dt, startTime, endTime, outputFrequency, jacobianStepSize
      NAMELIST / SpaceManagement / SpecMeshFile, PeaceMeshFile, UCDMeshFile, MeshType, topographicShape, QuadType, polyDeg, &
                                    nXElem, nYElem, nZElem, nProc, nProcX, nProcY, nProcZ, &
                                   nPlot, xScale, yScale, zScale
      NAMELIST / SubgridScale / SubGridModel, viscosity, viscLengthScale, nCutoff
      NAMELIST / PhysicalConstants / fRotX, fRotY, fRotZ, Cd, dragscale, g, Cv, R, T0, dTdz, dTdx, dTdy, rho0, P0, v0
      
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
      
       
      INQUIRE( FILE = 'runtime.params', EXIST = fileExists )

      IF( fileExists )THEN 
         OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
            READ( UNIT = nUnit, NML = TimeManagement )
            READ( UNIT = nUnit, NML = SpaceManagement )
            READ( UNIT = nUnit, NML = SubgridScale )
            READ( UNIT = nUnit, NML = PhysicalConstants )
         CLOSE( UNIT = nUnit ) 

         ! TimeManagement
         ! Internally, units of seconds are always used
         IF( units(1:1) == 'h' )THEN

            thisParam % dt        = dt/3600.0_prec
            thisParam % startTime = startTime/3600.0_prec
            thisParam % endTime   = endTime/3600.0_prec
            thisParam % outputFrequency = outputFrequency/3600.0_prec

         ELSEIF( units(1:1) == 'm' )THEN

            thisParam % dt        = dt/60.0_prec
            thisParam % startTime = startTime/60.0_prec
            thisParam % endTime   = endTime/60.0_prec
            thisParam % outputFrequency = outputFrequency/60.0_prec

         ELSEIF( units(1:1) == 's' )THEN
 
            thisParam % dt        = dt
            thisParam % startTime = startTime
            thisParam % endTime   = endTime
            thisParam % outputFrequency = outputFrequency
 
         ELSE

            PRINT(MsgFMT), 'S/R Build_Params : Unknown units specified.'
            RETURN

         ENDIF

         thisParam % nStepsPerDump = INT( thisParam % outputFrequency/thisParam % dt )
         thisParam % nDumps        = INT( (thisParam % endTime-thisParam % startTime)/thisParam % outputFrequency )
         
         PRINT(MsgFMT), 'S/R Build_Params : Estimated Number of Time Steps :'
         PRINT('(4x,I10)'), thisParam % nStepsPerDump*thisParam % nDumps 
        
         thisParam % jacobianStepSize = jacobianStepSize
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

         readSuccess = .TRUE.
      ELSE   

         OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params', ACTION = 'WRITE' )
            WRITE( UNIT = nUnit, NML = TimeManagement )
            WRITE( UNIT = nUnit, NML = SpaceManagement )
            WRITE( UNIT = nUnit, NML = SubgridScale )
            WRITE( UNIT = nUnit, NML = PhysicalConstants )
         CLOSE( UNIT = nUnit ) 

         PRINT(MsgFMT), 'S/R Build_FluidParams : runtime.params not found.' 
         PRINT(MsgFMT), 'S/R Build_FluidParams : A sample runtime.params namelist file has been'
         PRINT(MsgFMT), 'generated for you in your current directory.'
         readSuccess = .FALSE.

      ENDIF
      
      
 END SUBROUTINE Build_FluidParams

END MODULE FluidParams_Class
