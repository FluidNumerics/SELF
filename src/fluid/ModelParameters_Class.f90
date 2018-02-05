! ModelParameters_Class.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
MODULE ModelParameters_Class

! src/common/
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
      INTEGER       :: filterType
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
      REAL(prec)    :: rho0 ! Reference density
      REAL(prec)    :: P0   ! reference pressure
      REAL(prec)    :: v0

#ifdef HAVE_CUDA
      ! Physical
      REAL(prec), DEVICE, ALLOCATABLE    :: fRotX_dev  ! coriolis parameter (x-component)
      REAL(prec), DEVICE, ALLOCATABLE    :: fRotY_dev   ! "                " (y-component)
      REAL(prec), DEVICE, ALLOCATABLE    :: fRotZ_dev   ! "                " (z-component)
      REAL(prec), DEVICE, ALLOCATABLE    :: Cd_dev 
      REAL(prec), DEVICE, ALLOCATABLE    :: dragscale_dev 
      REAL(prec), DEVICE, ALLOCATABLE    :: g_dev   ! gravitational acceleration
      REAL(prec), DEVICE, ALLOCATABLE    :: Cv_dev  ! Heat Capacity at constant volume
      REAL(prec), DEVICE, ALLOCATABLE    :: R_dev   ! "Ideal gas constant"
      REAL(prec), DEVICE, ALLOCATABLE    :: T0_dev  ! Reference Temperature
      REAL(prec), DEVICE, ALLOCATABLE    :: dTdz_dev  ! Linear temperature stratification
      REAL(prec), DEVICE, ALLOCATABLE    :: rho0_dev  ! Reference density
      REAL(prec), DEVICE, ALLOCATABLE    :: P0_dev    ! reference pressure
      REAL(prec), DEVICE, ALLOCATABLE    :: v0_dev 

#endif

      CONTAINS

      PROCEDURE :: Build => Build_ModelParameters
      PROCEDURE :: Trash => Trash_ModelParameters

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
  
 CONTAINS


 SUBROUTINE Build_ModelParameters( params, readSuccess )

   CLASS( ModelParameters ), intent(out) :: params
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
      CHARACTER(50) :: SELFMeshFile
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
      CHARACTER(20) :: filterType
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
      REAL(prec)    :: rho0 ! Reference Density
      REAL(prec)    :: P0   ! reference pressure
      REAL(prec)    :: v0
      
      LOGICAL :: fileExists
       
       
      NAMELIST / TimeManagement / units, dt, startTime, endTime, outputFrequency, jacobianStepSize
      NAMELIST / SpaceManagement / SpecMeshFile, SELFMeshFile, UCDMeshFile, MeshType, topographicShape, QuadType, polyDeg, &
                                    nXElem, nYElem, nZElem, nProc, nProcX, nProcY, nProcZ, &
                                   nPlot, xScale, yScale, zScale
      NAMELIST / SubgridScale / SubGridModel, filterType, viscosity, viscLengthScale, nCutoff
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
      xScale        = 1.0_prec
      yScale        = 1.0_prec
      zScale        = 1.0_prec
      ! SubgridScale
      SubGridModel    = 'Laplacian'
      filterType      = 'TanhRollOff'
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

            PRINT(MsgFMT), 'S/R Build_Params : Unknown units specified.'
            RETURN

         ENDIF

         params % nStepsPerDump = INT( params % outputFrequency/params % dt )
         params % nDumps        = INT( (params % endTime-params % startTime)/params % outputFrequency )
         
         PRINT(MsgFMT), 'S/R Build_Params : Estimated Number of Time Steps :'
         PRINT('(4x,I10)'), params % nStepsPerDump*params % nDumps 
        
         params % jacobianStepSize = jacobianStepSize
         ! SpaceManagement 
         params % SpecMeshFile = SpecMeshFile
         params % SELFMeshFile = SELFMeshFile
         params % UCDMeshFile = UCDMeshFile
         IF( TRIM( UpperCase( MeshType ) )=='DEFAULT' )THEN
            params % MeshType = DefaultMesh
         ELSEIF( TRIM( UpperCase( MeshType ) )=='DOUBLYPERIODIC' )THEN
            params % MeshType = DoublyPeriodic
         ELSE
            PRINT*, 'Module ModelParameters_Class.f90 : S/R Build '
            PRINT*, '   Invalid MeshType : '//UpperCase(MeshType)
            PRINT*, '   Valid options are "Default" or "DoublyPeriodic"'
            STOP 'STOPPING!'
         ENDIF
         
         IF( TRIM( UpperCase( topographicShape ) ) =='DEFAULT' )THEN
            params % topographicShape = DefaultMesh
         ELSEIF( TRIM( UpperCase( topographicShape ) ) =='GAUSSIANHILL' )THEN
            params % topographicShape = Gaussian
         ENDIF
         
         params % QuadType      = QuadType
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
         IF( TRIM( UpperCase( SubGridModel) )=='LAPLACIAN' )THEN
         params % SubGridModel = Laplacian ! convert to the integer flag
         ELSEIF( TRIM(UpperCase( SubGridModel) )=='SPECTRALEKE' )THEN
            params % SubGridModel = SpectralEKE
         ELSEIF( TRIM(UpperCase( SubGridModel) )=='SPECTRALFILTERING' )THEN
            params % SubGridModel = SpectralFiltering
         ELSE
            PRINT*, 'Module ModelParameters_Class.f90 : S/R Build '
            PRINT*, '   Invalid SubGridModel : '//UpperCase(SubGridModel)
            PRINT*, '   Valid options are "Laplacian", "SpectralEKE", or "SpectralFiltering"'
            STOP 'STOPPING!'
         ENDIF
         IF( TRIM( UpperCase( FilterType ) )=='TANHROLLOFF' )THEN
           params % filterType = tanhRollOff
         ELSEIF( TRIM( UpperCase( FilterType ) )=='MODALCUTOFF' )THEN
           params % filterType = modalCutoff
         ENDIF
         !params % subGridModel = subGridModel
         params % viscosity       = viscosity
         params % viscLengthScale = viscLengthScale
         params % nCutoff         = nCutoff
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
         params % P0   = P0
         params % v0   = v0

         readSuccess = .TRUE.
      ELSE   

         OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params', ACTION = 'WRITE' )
            WRITE( UNIT = nUnit, NML = TimeManagement )
            WRITE( UNIT = nUnit, NML = SpaceManagement )
            WRITE( UNIT = nUnit, NML = SubgridScale )
            WRITE( UNIT = nUnit, NML = PhysicalConstants )
         CLOSE( UNIT = nUnit ) 

         PRINT(MsgFMT), 'S/R Build_ModelParameters : runtime.params not found.' 
         PRINT(MsgFMT), 'S/R Build_ModelParameters : A sample runtime.params namelist file has been'
         PRINT(MsgFMT), 'generated for you in your current directory.'
         readSuccess = .FALSE.

      ENDIF

#ifdef HAVE_CUDA

      ALLOCATE( params % fRotX_dev, & 
                params % fRotY_dev, &
                params % fRotZ_dev, &
                params % Cd_dev, &
                params % dragscale_dev, & 
                params % g_dev, &
                params % Cv_dev, &
                params % R_dev, &
                params % T0_dev, &
                params % dTdz_dev, &
                params % rho0_dev, &
                params % P0_dev, & 
                params % v0_dev ) 

     CALL params % UpdateDevice( )
#endif
      
      
 END SUBROUTINE Build_ModelParameters

 SUBROUTINE Trash_ModelParameters( params )
 IMPLICIT NONE
 CLASS( ModelParameters ), INTENT(inout) :: params

#ifdef HAVE_CUDA

      DEALLOCATE( params % fRotX_dev, & 
                  params % fRotY_dev, &
                  params % fRotZ_dev, &
                  params % Cd_dev, &
                  params % dragscale_dev, & 
                  params % g_dev, &
                  params % Cv_dev, &
                  params % R_dev, &
                  params % T0_dev, &
                  params % dTdz_dev, &
                  params % rho0_dev, &
                  params % P0_dev, & 
                  params % v0_dev ) 

#endif

 END SUBROUTINE Trash_ModelParameters
 

#ifdef HAVE_CUDA
 SUBROUTINE UpdateDevice_ModelParameters( params )
 IMPLICIT NONE
 CLASS( ModelParameters ), INTENT(inout) :: params

 params % fRotX_dev     = params % fRotX
 params % fRotY_dev     = params % fRotY
 params % fRotZ_dev     = params % fRotZ
 params % Cd_dev        = params % Cd
 params % dragScale_dev = params % dragScale
 params % g_dev         = params % g
 params % Cv_dev        = params % Cv
 params % R_dev         = params % R
 params % T0_dev        = params % T0
 params % dTdz_dev      = params % dTdz
 params % rho0_dev      = params % rho0
 params % P0_dev        = params % P0
 params % v0_dev        = params % v0

 END SUBROUTINE UpdateDevice_ModelParameters
#endif

END MODULE ModelParameters_Class
