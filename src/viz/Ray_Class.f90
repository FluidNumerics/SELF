
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.





MODULE Ray_Class


 ! src/common
 USE ModelPrecision
 USE CommonRoutines
 ! src/spectralops/
 USE Lagrange_Class
 ! src/geom/
 USE BoundaryCommunicator_Class
 USE MappedGeometry_3D_Class 
 USE HexMesh_Class


 IMPLICIT NONE
 
    ! A ray is defined by its starting point and direction. For 3-D graphics rendering, typically
    !
    TYPE Ray
       INTEGER                 :: nSamples
       REAL(prec), POINTER     :: sourcePoint(:)
       REAL(prec)              :: entryPoint(1:3), exitPoint(1:3)
       REAL(prec)              :: direction(1:3)
       REAL(prec)              :: emissionWeight
       REAL(prec), ALLOCATABLE :: physicalEmissionPoints(:,:) ! x,y,z by nSamples
       REAL(prec), ALLOCATABLE :: computationalEmissionPoints(:,:)  ! (s,p,q) by nSamples
       INTEGER, ALLOCATABLE    :: elementIDs(:) ! nSamples
       REAL(prec), ALLOCATABLE :: samples(:) ! sampled scalar field on physicalEmissionPoints
       REAL(SP), ALLOCATABLE   :: rgba(:,:) !(R,G,B,alpha) by nSamples (single precision
    
       CONTAINS
       PROCEDURE :: Build => Build_Ray
       PROCEDURE :: Trash => Trash_Ray
    
       PROCEDURE :: FindMeshEndPoints => FindMeshEndPoints_Ray
       PROCEDURE :: FindComputationalCoordinates => FindComputationalCoordinates_Ray
       PROCEDURE :: SampleScalar => SampleScalar_Ray
       
    END TYPE Ray
 
 
 CONTAINS
 
 SUBROUTINE Build_Ray( myRay, nSamples )
   ! Allocates space for single ray and initializes attributes to default values
   IMPLICIT NONE
   CLASS( Ray ), INTENT(inout) :: myRay
   INTEGER, INTENT(in)         :: nSamples
   
      myRay % nSamples = nSamples
      
      myRay % sourcePoint => NULL( )
      myRay % entryPoint    = 0.0_prec
      myRay % exitPoint     = 0.0_prec
      myRay % direction     = 0.0_prec
      myRay % emissionWeight = 0.0_prec
      
      ALLOCATE( myRay % physicalEmissionPoints(1:3,1:nSamples), &
                myRay % computationalEmissionPoints(1:3,1:nSamples), &
                myRay % elementIDs(1:nSamples), &
                myRay % samples(1:nSamples), &
                myRay % rgba(1:4,1:nSamples) )
                
      myRay % physicalEmissionPoints      = 0.0_prec
      myRay % computationalEmissionPoints = 0.0_prec
      myRay % elementIDs                  = 0
      myRay % samples                     = 0.0_prec
      myRay % rgba                        = 0.0_prec
      
 END SUBROUTINE Build_Ray
!
 SUBROUTINE Trash_Ray( myRay )

   IMPLICIT NONE
   CLASS( Ray ), INTENT(inout) :: myRay
   
      DEALLOCATE( myRay % physicalEmissionPoints, &
                  myRay % computationalEmissionPoints, &
                  myRay % elementIDs, &
                  myRay % samples, &
                  myRay % rgba )
                
 END SUBROUTINE Trash_Ray
!
!
!
 SUBROUTINE FindMeshEndPoints_Ray( myRay, mesh, extComm, interp )
  ! This subroutine uses the array starting point and the array direction to cast a ray towards
  ! the volume of space described by the unstructured HexMesh and the boundary communicator.
  ! The boundary communicator is used to only perform searches through the boundary faces of the 
  ! mesh.
  ! 
  ! The routine cycles through the boundary faces and determines if an intersection is possible 
  ! through the boundary face. We assume that there are at most 2 faces where this is possible, 
  ! one entry point and one exit point.
  !
  ! The ray is parameterized x(t) = x0 + d*t, where x0 is the "sourcePoint" and d is the "direction"
  ! Each bounding surface is a plane with computational coordinates between [-1,1]x[-1,1].
  ! what we want to test for each plane is if there is a (t,p,q) such that
  !   x(t) = x0 + d*t = xs(p,q)
  !
  ! This system is solved using Newton's method for each bounding surface. This method is terminated 
  ! early if (p,q) venture outside of [-1,1]x[-1,1] --> in this case, it is likely that the ray 
  ! does not interesect the bounded surface. If a solution is found, the (t,p,q) are recorded in
  ! addition to the bounding face ID. This is used to initialize the ray for ray-casting.
  !
  ! The relative size of the "t" parameters found determine the entry and exit points through
  ! the volume.
  !        
  ! This routine assumes one entry and on exit point (no holes in the mesh).
  !
   CLASS( Ray ), INTENT(inout) :: myRay
   TYPE( HexMesh ), INTENT(in) :: mesh
   TYPE( BoundaryCommunicator ), INTENT(in) :: extComm
   TYPE( Lagrange ), INTENT(in)             :: interp
   ! Local
   INTEGER :: bFace, iFace, i, j, k, e1, s1, nsurf
   INTEGER :: bFaces(1:2)
   REAL(prec) :: s(1:3), xs(1:3), xr(1:3), ds(1:3), t, resi
   REAL(prec) :: xsurf(0:interp % N, 0:interp % N,1:3)
   REAL(prec) :: dFdp(1:3,1:3), Jinv(1:3,1:3)
   REAL(prec) :: coords(1:2)
   
   
      dFdp = 0.0_prec
      dFdp(:,1) = myRay % direction
      
      nsurf = 0
      coords = 0.0_prec
      bFaces = 0
   
      DO bFace = 1, extComm % nBoundaries
      
         IF( nsurf < 2 )THEN
            ! Initial gues
            s(1) = 1.0_prec
            s(2:3) = 0.0_prec
        
            iFace = extComm % boundaryIDs(bFace) ! Get the global face ID for this boundary face  
            e1    = mesh % faces(iFace) % elementIDs(1) ! Obtain the element ID associated with this face  
            s1    = mesh % faces(iFace) % elementSides(1) ! Obtain the local side ID for this face
            
            ! For now, copy the surface positions for this face
            xsurf(:,:,1) = mesh % geom(e1) % xBound(:,:,s1)
            xsurf(:,:,2) = mesh % geom(e1) % yBound(:,:,s1)
            xsurf(:,:,3) = mesh % geom(e1) % zBound(:,:,s1)
            
            DO i = 1, newtonMax
        
               ! Calculate the ray position
               xr = myRay % sourcePoint + myRay % direction*s(1) 
               
               ! Evaluate the surface position and gradients
               DO j = 1,3
                 xs(j) = interp % Interpolate_2D( xsurf(:,:,j), s(2:3) )
                 dFdp(j,2:3) = -interp % Differentiate_2D( xsurf(:,:,j), s(2:3) )
               ENDDO
        
               ! Calculate the inverse of the matrix system
               Jinv = Invert_3x3( dFdp )
               
               ! Update the solution
               DO j = 1, 3
                  ds(j) = 0.0_prec
                  DO k = 1,3
                     ds(j) = ds(j) - Jinv(j,k)*( xr(k) - xs(k) )
                  ENDDO
                  s(j) = s(j) + ds(j)
               ENDDO
               
               ! If the surface parameters venture outside of [-1,1], we terminate this loop
               IF( ABS(s(2)) > 1.01_prec .OR. ABS(s(3)) > 1.01_prec )THEN
                  EXIT
               ENDIF
               resi = SQRT( DOT_PRODUCT( ds, ds ) )
        
               IF( resi < newtonTolerance )THEN
                  nsurf = nsurf + 1
                !  PRINT*, 'Found new surface!'
                 
                  bfaces(nsurf) = bface
                  coords(nsurf) = s(1)
                  !PRINT*, s
                  !PRINT*, xs
                  EXIT
               ENDIF
        
            ENDDO
     
         ENDIF
         
      ENDDO
      
     ! PRINT*, 'Ray_Class : S/R FindMeshEndPoints : N intersecting mesh surfaces :',nSurf
      IF( nSurf < 2 )THEN
        PRINT*, 'Ray_Class : S/R FindMeshEndPoints : N intersecting mesh surfaces :',nSurf
        PRINT*, 'Missing Either entry or exit face'
      ENDIF
      
      myRay % entryPoint(1:3) = myRay % sourcePoint + myRay % direction*MINVAL( coords )
      myRay % exitPoint(1:3)  = myRay % sourcePoint + myRay % direction*MAXVAL( coords )
      
      ! Set up the ray physical sampling locations
      myRay % emissionWeight = 1.0_prec/REAL( myRay % nSamples-1, prec )
      DO i = 1, myRay % nSamples
         t = REAL(i-1,prec)*myRay % emissionWeight
         myRay % physicalEmissionPoints(1:3,i) = myRay % entryPoint + &
                                                 (myRay % exitPoint - myRay % entryPoint )*t
      ENDDO
      
      ! Set the starting and terminating element IDs
      i = MINLOC( coords,1 )
      bFace = bfaces(i)
      iFace = extComm % boundaryIDs(bFace) ! Get the global face ID for this boundary face  
      e1    = mesh % faces(iFace) % elementIDs(1) ! Obtain the element ID associated with this face  
      myRay % elementIDs(1) = e1
      
      ! Terminating element
      i = MAXLOC( coords,1 )
      bFace = bfaces(i)
      iFace = extComm % boundaryIDs(bFace) ! Get the global face ID for this boundary face  
      e1    = mesh % faces(iFace) % elementIDs(1) ! Obtain the element ID associated with this face  
      myRay % elementIDs(myRay % nSamples) = e1
      

 END SUBROUTINE FindMeshEndPoints_Ray
! 
 SUBROUTINE FindComputationalCoordinates_Ray( myRay, mesh, interp )
   ! This routine traverses the "physicalEmissionPoints" of the ray and determines the element ID
   ! that each point is located in and the computational coordinate. This knowledge is necessary
   ! for interpolating a scalar field onto the ray for sampling and compositing.
   
   IMPLICIT NONE
   CLASS( Ray ), INTENT(inout)  :: myRay
   TYPE( HexMesh ), INTENT(in)  :: mesh
   TYPE( Lagrange ), INTENT(in) :: interp
   ! Local
   INTEGER :: i, j, k, e1, e2, e3
   LOGICAL :: successful
   REAL(prec) :: s(1:3)
   
   
      
      DO i = 1, myRay % nSamples-1
      
         successful = .FALSE.
         ! For the first emission point on the ray, FindMeshEndPoints gave us the element ID
         ! where the intersection occurred. We can use this element ID to find the computational
         ! coordinate
         e1 = myRay % elementIDs(i)
         IF( e1 > 0 )THEN 
        
            CALL mesh % geom(e1) % CalculateComputationalCoordinates( interp, &
                                  myRay % physicalEmissionPoints(1:3,i), s, successful )

            IF( successful )THEN
           
               myRay % elementIDs(i) = e1
               myRay % computationalEmissionPoints(1:3,i) = s
               ! Set the initial element guess for the next ray point to be this element ID
               myRay % elementIDs(i+1) = e1
              
            ELSE
           
            ! Search through this elements neighbors (one levels)
               DO j = 1, 6 ! Loop over the faces of the element
                  e2 = mesh % elements(e1) % neighbors(j)
                 
                  IF( e2 > 0 )THEN
                     CALL mesh % geom(e2) % CalculateComputationalCoordinates( interp, &
                                     myRay % physicalEmissionPoints(1:3,i), s, successful )
                                     
                     IF( successful )THEN
                        myRay % elementIDs(i) = e2
                        myRay % computationalEmissionPoints(1:3,i) = s
                        ! Set the initial element guess for the next ray point to be this element ID
                        myRay % elementIDs(i+1) = e2
                        EXIT
                     ENDIF
                  ENDIF
                 
               ENDDO
               
               IF( .NOT. successful ) THEN
                  ! Search through two levels of neighbors
                  DO j = 1, 6 ! Loop over the faces of the element
                     e2 = mesh % elements(e1) % neighbors(j)
                 
                     IF( e2 > 0 )THEN
                        DO k = 1, 6
                           e3 = mesh % elements(e2) % neighbors(k)
                           IF( e3 > 0 )THEN
                              CALL mesh % geom(e2) % CalculateComputationalCoordinates( interp, &
                                           myRay % physicalEmissionPoints(1:3,i), s, successful )
                                        
                              IF( successful )THEN
                                 myRay % elementIDs(i) = e2
                                 myRay % computationalEmissionPoints(1:3,i) = s
                                 ! Set the initial element guess for the next ray point to be this element ID
                                 myRay % elementIDs(i+1) = e2
                                 EXIT
                              ENDIF
                              
                           ENDIF ! e3 > 0
                        ENDDO ! k
                        
                        IF( successful )THEN
                           EXIT
                        ENDIF
                        
                     ENDIF ! e2 > 0
                  ENDDO ! j-loop over neighbors
                 
               ENDIF
              
            ENDIF
            
         ENDIF ! e1 > 0
         
      ENDDO ! Loop over physicalEmissionPoints
      
      e1 = myRay % elementIDs( myRay % nSamples )
      ! Obtain the computational coordinates of the termination point
      CALL mesh % geom(e1) % CalculateComputationalCoordinates( interp, &
                     myRay % physicalEmissionPoints(1:3,myRay % nSamples), s )
   
 END SUBROUTINE FindComputationalCoordinates_Ray
!
 SUBROUTINE SampleScalar_Ray( myRay, scalar, interp, nElems, minscalar, maxscalar )
   ! This routine traverses the emission points on the ray and interpolates the
   ! solution from the "solStorage" structure onto the array.
   IMPLICIT NONE
   CLASS( Ray ), INTENT(inout)         :: myRay
   INTEGER, INTENT(in)                 :: nElems
   TYPE( Lagrange ), INTENT(in)        :: interp
   REAL(prec), INTENT(in)              :: scalar(0:interp % N, &
                                                 0:interp % N, &
                                                 0:interp % N, &
                                                 1:nElems)
   REAL(prec), INTENT(in) :: minScalar, maxScalar
   INTEGER :: e1, i, j
   
      DO i = 1, myRay % nSamples
      
         e1 = myRay % elementIDs(i) ! Obtain the element ID for this emission point
         
         myRay % samples(i) = interp % Interpolate_3D( scalar(:,:,:,e1), &
                                       myRay % computationalEmissionPoints(1:3,i) )
                                       
         ! Convert scalar to RGB
         myRay % rgba(1:4,i) = TransferFunction( myRay % samples(i), minScalar, maxScalar )
         
      ENDDO
      
 END SUBROUTINE SampleScalar_Ray
!
 FUNCTION TransferFunction( scalarvalue, maxScalar, minScalar ) RESULT( rgba )
    IMPLICIT NONE
    REAL(prec) :: scalarvalue, maxScalar, minScalar
    REAL(prec) :: rgba(1:4)
    REAL(prec) :: nrml
    
       nrml = (scalarvalue - minScalar)/(maxScalar - minScalar)
       
       ! For now we set the red and blue channels to zero, but set the green 
       ! and opacity channels....
       rgba(1) = 0.0_SP ! red !
       rgba(2) = REAL( nrml**2, SP )  ! green !
       rgba(3) = 0.0_SP ! blue !
       rgba(4) = REAL( nrml**2 , SP ) ! opacity
 
       ! Saturate if above 1.0
       IF( nrml > 1.0_prec )THEN
          rgba(2) = 1.0_SP
          rgba(4) = 1.0_SP
       ENDIF
       
 END FUNCTION TransferFunction
 
END MODULE Ray_Class
