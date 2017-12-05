! CompressibleNavierStokes_InitialConditions.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
PROGRAM CompressibleNavierStokes_InitialConditions

! src/common/
USE ModelPrecision
! src/geom/
!USE HexMesh_Class
! src/highend/euler/
USE CompressibleNavierStokesParams_Class
! src/highend/euler/
USE CompressibleNavierStokes_Class

 IMPLICIT NONE

 TYPE( CompressibleNavierStokes ) :: myeu
 INTEGER       :: myRank, mpiErr, nProcs
 CHARACTER(4)  :: rankChar
 
#ifdef HAVE_MPI
      ! MPI Initialization
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
      ! Sanity check
      PRINT*, 'CompressibleNavierStokes_InitialConditions_MPI : Greetings from Process ', myRank, ' of ',nProcs
#else
      myRank = 0
      nProcs = 1
#endif

      CALL myeu % Build( myRank, nProcs )
      
      CALL InitialCondition( myeu )
      
      CALL ResetBoundaryConditions( myeu, myRank )
      
      CALL myeu % WritePickup( 0, myRank ) 
      
      CALL myeu % WriteTecplot( 0, myeu % params % nPlot, myRank )
      
      ! Before we write the mesh to file again, we need to "unscale" the mesh so that, upon running the 
      ! integrator, the mesh scaling is not applied a second time 
      CALL myeu % mesh % ScaleTheMesh( myeu % dgStorage % interp, &
                                          1.0_prec/myeu % params % xScale, &
                                          1.0_prec/myeu % params % yScale, &
                                          1.0_prec/myeu % params % zScale )
      WRITE( rankChar, '(I4.4)' )myRank
      CALL myeu % mesh % WritePeaceMeshFile( TRIM(myeu % params % PeaceMeshFile)//'.'//rankChar )
      
      CALL myeu % Trash( )
      
#ifdef HAVE_MPI
      CALL MPI_FINALIZE( mpiErr )
#endif

 CONTAINS
 
 SUBROUTINE ResetBoundaryConditions( myDGSEM, myRank )
  
   TYPE( CompressibleNavierStokes ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)          :: myRank
   ! Local
   INTEGER :: iFace, iFace2, e1, e2, s1, p2, i, j
   REAL(prec) :: Lx, Ly, x, y, T, Tbar, rho
   
      Lx = myDGSEM % params % xScale
      Ly = myDGSEM % params % yScale
      
      DO iFace = 1, myDGSEM % nBoundaryFaces

         iFace2 = myDGSEM % extComm % boundaryIDs( iFace )
         e1    = myDGSEM % mesh % Faces(iFace2) % elementIDs(1)
         s1    = myDGSEM % mesh % Faces(iFace2) % elementSides(1)
         e2    = myDGSEM % mesh % Faces(iFace2) % elementIDs(2)
         p2    = myDGSEM % extComm % extProcIDs( iFace )
         
         IF( e2 < 0 .AND. p2 == myRank )THEN
         
            IF( ABS(myDGSEM % mesh % geom(e1) % nHat(3, 0, 0, s1) ) >&
                SQRT( myDGSEM % mesh % geom(e1) % nHat(1, 0, 0, s1)**2 + &
                      myDGSEM % mesh % geom(e1) % nHat(2, 0, 0, s1)**2 ) ) THEN ! Top or bottom surface
                      
               IF( myDGSEM % mesh % geom(e1) % nHat(3, 0, 0, s1) < 0.0_prec ) THEN ! Bottom
                  myDGSEM % mesh % faces(iFace2) % elementIDs(2) = PRESCRIBED
                  
                  DO j = 0, myDGSEM % N
                     DO i = 0, myDGSEM % N
                     
                        x    = myDGSEM % mesh % geom(e1) % xBound(i,j,s1)
                        y    = myDGSEM % mesh % geom(e1) % yBound(i,j,s1)
                        
                        T    = 10.0_prec*exp( -0.5_prec*( (x-0.5_prec*Lx)**2 + (y-0.5_prec*Ly)**2 )/(50.0_prec)**2 ) ! Prescribe a warm temperature plume in the center of the domain.
                        Tbar = myDGSEM % static % boundarysolution(i,j,5,s1,e1)/myDGSEM % static % boundarysolution(i,j,4,s1,e1)
                        
                        ! Set the density anomaly
                        myDGSEM % prescribedState(i,j,4,iFace) = -myDGSEM % static % boundarysolution(i,j,4,s1,e1)*T/(Tbar + T)
                        
                        ! Prescribe a vertical velocity associated with the plume
                        rho =  myDGSEM % prescribedState(i,j,4,iFace) + myDGSEM % static % boundarysolution(i,j,4,s1,e1)
                        myDGSEM % prescribedState(i,j,3,iFace) = rho*10.0_prec*&
                               exp( -0.5_prec*( (x-0.5_prec*Lx)**2 + (y-0.5_prec*Ly)**2 )/(50.0_prec)**2 )
                     ENDDO
                  ENDDO
                  
               ELSE ! Top boundary
                  myDGSEM % mesh % faces(iFace2) % elementIDs(2) = NO_NORMAL_FLOW
               ENDIF
               
            ELSE
               myDGSEM % mesh % faces(iFace2) % elementIDs(2) = RADIATION
            ENDIF
               
         ENDIF
 
      ENDDO
      
      
 END SUBROUTINE ResetBoundaryConditions
 
 SUBROUTINE InitialCondition( myDGSEM )

   TYPE( CompressibleNavierStokes ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: i, j, k, iEl, iFace
   REAL(prec) :: x, y, z, hScale, Lx, Ly, r

   
      Lx = myDGSEM % params % xScale   
      Lx = myDGSEM % params % yScale   
      CALL myDGSEM % CalculateStaticState( ) !! CPU Kernel

      ! ////////////////////////////////////////////////////////////////////////////////// !
      DO iEl = 1, myDGSEM % mesh % nElems
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
            
                  x = myDGSEM % mesh % geom(iEl) % x(i,j,k)
                  y = myDGSEM % mesh % geom(iEl) % y(i,j,k)
                  z = myDGSEM % mesh % geom(iEl) % z(i,j,k)
                  
                  r = ( x-0.5_prec*Lx )**2 + ( y-0.5_prec*Ly )**2 + ( z )**2
                
                  ! Drag profile
                  myDGSEM % dragProfile(i,j,k,iEl) = myDGSEM % params % Cd*exp( -z/myDGSEM % params % dragScale )*&
                                                    (1.0_prec-exp( -r/(0.1_prec*Lx)**2 ))

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
#ifdef HAVE_CUDA     
     ! Copy the static solution from host to device
     myDGSEM % static % solution_dev = myDGSEM % static % solution
#endif
     CALL myDGSEM % CalculateStaticBoundarySolution( )
#ifdef HAVE_CUDA
     ! Copy the boundary solution from the device to the host
     myDGSEM % static % boundarySolution = myDGSEM % static % boundarySolution_dev
#endif

     
 END SUBROUTINE InitialCondition
!
END PROGRAM CompressibleNavierStokes_InitialConditions

