PROGRAM TestRay_Class



 ! src/spectralOps/
 USE Lagrange_Class
 USE NodalStorage_Class
 ! src/templates
 USE Params_Class
 ! src/geom/
 USE HexMesh_Class
 USE BoundaryCommunicator_Class
 ! src/viz/ 
 USE Ray_Class
 
    INTEGER, PARAMETER   :: nSamples = 64
    TYPE( Ray )          :: thisRay
    TYPE( NodalStorage ) :: dgStorage
    TYPE( RunParams )    :: params
    TYPE( HexMesh )      :: mesh
    TYPE( BoundaryCommunicator ) :: extComm
    REAL(prec) :: x, y, z
    REAL(prec), ALLOCATABLE :: elementMarkers(:), f(:,:,:,:)
    REAL(prec), TARGET :: vantagePoint(1:3)
    INTEGER :: i,j,k,iel
    
      CALL params % Build( )
      
      CALL dGStorage % Build( params % polyDeg, params % nPlot, GAUSS, DG )
      
      PRINT*, 'Reading mesh from '//trim(params % PeaceMeshFile)//'.'
      CALL mesh % ReadPeaceMeshFile( params % PeaceMeshFile )
      
      ALLOCATE( elementMarkers(1:mesh % nElems), &
                f(0:params % polyDeg, 0:params % polyDeg, 0:params % polyDeg, 1:mesh % nElems) )
      elementMarkers = 0.0_prec
      f = 0.0_prec
                                                
      CALL extComm % ReadPickup( 'ExtComm' )
      
      CALL thisRay % Build( nSamples )
      
      ! Assume a domain on [0,1]^3
      ! Set the vantage point
      vantagePoint(1) = 0.1_prec ! x
      vantagePoint(2) = -1.5_prec ! y
      vantagePoint(3) = 1.5_prec
      
      thisRay % sourcePoint => vantagePoint
      ! Set the ray direction
      thisRay % direction(1) = 0.4_prec
      thisRay % direction(2) = 1.5_prec
      thisRay % direction(3) = -1.0_prec
      
      CALL thisRay % FindMeshEndPoints( mesh, extcomm, dgStorage % interp )
      ! Locate the element ID's and the computational coordinates for all of the ray points
      CALL thisRay % FindComputationalCoordinates( mesh, dgStorage % interp )
      
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
      
      CALL thisRay % FindMeshEndPoints( mesh, extcomm, dgStorage % interp )
      ! Locate the element ID's and the computational coordinates for all of the ray points
      CALL thisRay % FindComputationalCoordinates( mesh, dgStorage % interp )
     
     
      DO i = 1, thisRay % nSamples
         IF( thisRay % elementIDs(i) > 0 )THEN
            elementMarkers( thisRay % elementIDs(i) ) = 2.0_prec
         ENDIF
      ENDDO
      
            
      CALL mesh % WriteMaterialTecplot( elementMarkers )
      
      ! Fill in function value... use a radially symmetric gaussian
      DO iel = 1, mesh % nelems
         DO k = 0, params % polyDeg
            DO j = 0, params % polyDeg
               DO i = 0, params % polyDeg
                 
                 x = mesh % geom(iel) % x(i,j,k)
                 y = mesh % geom(iel) % y(i,j,k)
                 z = mesh % geom(iel) % z(i,j,k)
                 f(i,j,k,iel) = exp( -( (x-0.5_prec)**2 + (y-0.5_prec)**2 + (z-0.5_prec)**2 )/&
                                       (2.0_prec*0.01) )
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL thisRay % SampleScalar( f, dgStorage % interp, mesh % nElems, 0.0_prec, 1.0_prec )
      
      OPEN( UNIT=2, FILE = 'sampled.curve' )
      WRITE(2, *) '#sampled'
      DO i = 1, thisRay % nSamples
         WRITE(2, *) i, thisRay % samples(i)
      ENDDO  
      CLOSE(UNIT=2)    
      ! Clean up !
      
      CALL dgStorage % Trash( )
      CALL mesh % Trash( )
      CALL extComm % Trash( )
      CALL thisRay % Trash( )
      
      DEALLOCATE( elementMarkers )
      

END PROGRAM TestRay_Class
