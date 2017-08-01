PROGRAM SampleAlongRays

 ! src/common/
 USE ModelPrecision
 ! src/fluid/
 USE Fluid_Class
 ! src/viz/ 
 USE Ray_Class
 
    INTEGER                 :: nSamples
    TYPE( Ray )             :: thisRay
    TYPE( Fluid )           :: myCNS
    REAL(prec)              :: x, y, z
    REAL(prec), ALLOCATABLE :: elementMarkers(:)
    REAL(prec), ALLOCATABLE :: rho(:,:,:,:)
    REAL(prec), TARGET      :: vantagePoint(1:3)
    INTEGER                 :: i, j, k, iel
    INTEGER                 :: iter0, nT, dFreq
    CHARACTER(10)           :: iterChar
    LOGICAL                 :: setupFail


      setupFail = .FALSE.
      CALL InitializeFromCommandLine( )

      IF( setupFail )THEN
        PRINT*, '  sars setup failed!'
        STOP
      ENDIF
    
      CALL myCNS % Build( 0, 1 )
      
      
      ALLOCATE( elementMarkers(1:myCNS % mesh % nElems), &
                rho(0:myCNS % params % polyDeg, &
                  0:myCNS % params % polyDeg, &
                  0:myCNS % params % polyDeg, &
                  1:myCNS % mesh % nElems) )

      elementMarkers = 0.0_prec
      rho = 0.0_prec
                                                
      CALL thisRay % Build( nSamples )
      
      ! Set the vantage point 
      ! Hard Coded to be centered in (x,y) and below the bottom in (z)
      vantagePoint(1) = 0.5_prec*myCNS % params % xScale ! x
      vantagePoint(2) = 0.5_prec*myCNS % params % yScale ! y
      vantagePoint(3) = -10.0_prec
      
      thisRay % sourcePoint => vantagePoint
      ! Set the ray direction so that it points upwards
      thisRay % direction(1) = 0.0_prec
      thisRay % direction(2) = 0.0_prec
      thisRay % direction(3) = 1.0_prec
      
      CALL thisRay % FindMeshEndPoints( myCNS % mesh, &
                                        myCNS % extcomm, &
                                        myCNS % dgStorage % interp )
      ! Locate the element ID's and the computational coordinates for all of the ray points
      CALL thisRay % FindComputationalCoordinates( myCNS % mesh, &
                                                   myCNS % dgStorage % interp )
      
      DO i = 1, thisRay % nSamples
         IF( thisRay % elementIDs(i) > 0 )THEN
            elementMarkers( thisRay % elementIDs(i) ) = 1.0_prec
         ENDIF
      ENDDO
      
      ! Set the vantage point
      vantagePoint(1) = 0.5_prec ! x
      vantagePoint(2) = -1.5_prec ! y
      vantagePoint(3) = 0.5_prec
      
      thisRay % sourcePoint => vantagePoint
      ! Set the ray direction
      thisRay % direction(1) = 0.0_prec
      thisRay % direction(2) = 1.5_prec
      thisRay % direction(3) = 0.0_prec
      
      CALL thisRay % FindMeshEndPoints( myCNS % mesh, &
                                        myCNS % extcomm, &
                                        myCNS % dgStorage % interp )
      ! Locate the element ID's and the computational coordinates for all of the ray points
      CALL thisRay % FindComputationalCoordinates( myCNS % mesh, &
                                                   myCNS % dgStorage % interp )
     
     
      DO i = 1, thisRay % nSamples
         IF( thisRay % elementIDs(i) > 0 )THEN
            elementMarkers( thisRay % elementIDs(i) ) = 2.0_prec
         ENDIF
      ENDDO
      
      ! This mesh is helpful for diagnosing possible issues with the process            
      CALL myCNS % mesh % WriteMaterialTecplot( elementMarkers )
      iter0  = myCNS % params % iterInit
      nT     = myCNS % params % nTimeSteps
      dFreq  = myCNS % params % dumpFreq

      DO iT = iter0, iter0+nT-1, dFreq ! Loop over time-steps

         CALL myCNS % ReadPickup( iT, 0 )
 
         ! Extract the density from the Fluid's data structure
         DO iel = 1, myCNS % mesh % nelems
            DO k = 0, myCNS % params % polyDeg
               DO j = 0, myCNS % params % polyDeg
                  DO i = 0, myCNS % params % polyDeg
                    rho(i,j,k,iel) = myCNS % state % solution(i,j,k,4,iEl)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         CALL thisRay % SampleScalar( rho, &
                                      myCNS % dgStorage % interp, &
                                      myCNS %  mesh % nElems, 0.0_prec, 1.0_prec )
         WRITE(iterChar,'(I10.10)')iT         
         OPEN( UNIT=2, FILE = 'rho_sampled.'//iterChar//'.curve' )
         WRITE(2, *) '#sampled'
         DO i = 1, thisRay % nSamples
            WRITE(2, *) i, thisRay % samples(i)
         ENDDO  
         CLOSE(UNIT=2)    

         OPEN( UNIT = 2, &
               FILE = 'rho_sampled.'//iterChar//'.bin', &
               FORM = 'UNFORMATTED', &
               ACCESS = 'DIRECT', & 
               STATUS = 'REPLACE', &
               RECL = prec*thisRay % nSamples )
         WRITE( 2, REC=1 ) thisRay % samples
         CLOSE(UNIT=2) 

      ENDDO
      ! Clean up !
      
      CALL myCNS % Trash( )
      CALL thisRay % Trash( )
      
      DEALLOCATE( elementMarkers )
CONTAINS

 SUBROUTINE InitializeFromCommandLine( )
   IMPLICIT NONE

   INTEGER       :: nArg, argID
   CHARACTER(50) :: argname
   LOGICAL       :: nGiven


     nGiven   = .FALSE.
     nSamples = 0
     nArg = command_argument_count( )

     IF( nArg > 0 )THEN
 
        DO argID = 1, nArg
  
          CALL get_command_argument( argID, argName )

          SELECT CASE( TRIM(argName) )

             CASE("--nSamples")
                nGiven = .TRUE.

             CASE DEFAULT

               IF( nGiven )THEN
               
                  READ( argName, '(I7)') nSamples
                  nGiven = .FALSE.
 
               ENDIF

          END SELECT 
        ENDDO

        

        IF( nSamples == 0 )THEN
          PRINT*, ' Number of ray sample points needs to be specified.'
          setupFail = .TRUE.
        ENDIF
        
     ELSE

        ! List possible options
        PRINT*, '  sars : Sample Along Rays'
        PRINT*, '    A tool for sampling SELF-CNS data along a ray.'
        PRINT*, '--------------------------------------------------------------'
        PRINT*, '  Usage : sars --nSamples < number of samples along the ray > '
        PRINT*, '--------------------------------------------------------------'
        PRINT*, ' This executable samples SELF-CNS density along a ray centered'
        PRINT*, ' on the (x,y) center and from z=0 to the top of the domain.   '
        PRINT*, '--------------------------------------------------------------'
        STOP
        
     ENDIF   


 END SUBROUTINE InitializeFromCommandLine
      
END PROGRAM SampleAlongRays
