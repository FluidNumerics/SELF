! Fluid_InitialConditions.f90
! 
! Copyright 2017 Joseph Schoonover <schoonover.numerics@gmail.com>
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
PROGRAM Fluid_InitialConditions

! src/common/
USE ModelPrecision
! src/geom/
USE HexMesh_Class
! src/highend/euler/
USE FluidParams_Class
! src/highend/euler/
USE Fluid_Class

 IMPLICIT NONE

 TYPE( Fluid ) :: myeu
 INTEGER       :: myRank, mpiErr, nProcs
 CHARACTER(4)  :: rankChar

#ifdef HAVE_MPI
      ! MPI Initialization
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
      ! Sanity check
      PRINT*, 'Fluid_InitialConditions_MPI : Greetings from Process ', myRank, ' of ',nProcs
#else
      myRank = 0
      nProcs = 1
#endif

      CALL myeu % Build( myRank, nProcs )
      
   !   CALL SetTopography( myeu, myRank )
      
      CALL InitialCondition( myeu )
      
      CALL ResetBoundaryConditions( myeu, myRank )
      
      CALL myeu % WritePickup( 0, myRank ) 
      
      CALL myeu % WriteTecplot( 0, myeu % params % nPlot, myRank )
      
      CALL myeu % mesh % WriteTecplot( 'mesh.'//rankChar )
      
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
 FUNCTION Topography( x, y, xc, yc, Lscale, H) RESULT( z )
   IMPLICIT NONE
   REAL(prec) :: x, y, xc, yc, Lscale, H
   REAL(prec) :: z
 
      z = H*exp( -( (x-xc)**2 + (y-yc)**2 )/(2.0_prec*Lscale**2) )
      
 END FUNCTION Topography
!
 SUBROUTINE SetTopography( myDGSEM, myRank )
  
   TYPE( Fluid ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)          :: myRank
   ! Local
   INTEGER :: iFace, iFace2, e1, e2, s1, p2
   INTEGER :: iEl, i, j, k
   REAL(prec) :: zbot, zold, xc, yc, Lscale, H, hillsize
   
      xc = myDGSEM % params % xScale/2.0_prec
      yc = myDGSEM % params % yScale/2.0_prec
      LScale = myDGSEM % params % xScale/10.0_prec
      hillsize = 0.35_prec*myDGSEM % params % zScale
      H = myDGSEM % params % zScale

      ! Reshape the mesh to incorporate topography
      DO iEl = 1, myDGSEM % mesh % nElems
      
 
         DO k = 1, 6
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
               
                  zold = myDGSEM % mesh % geom(iEl) % zBound(i,j,k)
                  zbot = Topography( myDGSEM % mesh % geom(iEl) % xBound(i,j,k), &
                                     myDGSEM % mesh % geom(iEl) % yBound(i,j,k), &
                                     xc, yc, Lscale, hillsize )
                  myDGSEM % mesh % geom(iEl) % zBound(i,j,k) = ( (H-zbot)/H )*zold + zbot
                  
               ENDDO
            ENDDO
         ENDDO
      
         CALL myDGSEM % mesh % geom(iEl) % ResetInternalMesh( myDGSEM % dgStorage % interp )
         
      ENDDO
                  
                                     
   
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
                      
               IF( myDGSEM % mesh % geom(e1) % nHat(3, 0, 0, s1) < 0.0_prec ) THEN ! Bottom surface
                  ! Reset the boundary surface so that this element conforms to the bottom topography
                  DO j = 0, myDGSEM % N
                     DO i = 0, myDGSEM % N
                        myDGSEM % mesh % geom(e1) % zBound(i,j,s1) = Topography( &
                                                                        myDGSEM % mesh % geom(e1) % xBound(i,j,s1), &
                                                                        myDGSEM % mesh % geom(e1) % yBound(i,j,s1), &
                                                                        xc, yc, Lscale, hillsize )
                     ENDDO
                  ENDDO
                  
                  CALL myDGSEM % mesh % geom(e1) % ResetInternalMesh( myDGSEM % dgStorage % interp )
                  CALL myDGSEM % mesh % geom(e1) % GenerateMetrics( myDGSEM % dgStorage % interp )
               
               ENDIF
            ENDIF
            
               
         ENDIF
 
      ENDDO
      
      
 END SUBROUTINE SetTopography 
!
 SUBROUTINE ResetBoundaryConditions( myDGSEM, myRank )
  
   TYPE( Fluid ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)          :: myRank
   ! Local
   INTEGER :: iFace, iFace2, e1, e2, s1, p2
   INTEGER :: iEl, i, j, k
   
      DO iFace = 1, myDGSEM % nBoundaryFaces

         iFace2 = myDGSEM % extComm % boundaryIDs( iFace )
         e1    = myDGSEM % mesh % Faces(iFace2) % elementIDs(1)
         s1    = myDGSEM % mesh % Faces(iFace2) % elementSides(1)
         e2    = myDGSEM % mesh % Faces(iFace2) % elementIDs(2)
         p2    = myDGSEM % extComm % extProcIDs( iFace )
         
         IF( e2 < 0 .AND. p2 == myRank )THEN
         
            myDGSEM % mesh % faces(iFace2) % elementIDs(2) = NO_NORMAL_FLOW
               
         ENDIF
 
      ENDDO
      
      
 END SUBROUTINE ResetBoundaryConditions
 
 SUBROUTINE InitialCondition( myDGSEM )

   TYPE( Fluid ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: i, j, k, iEl, iFace
   REAL(prec) :: x, y, z, hScale, Lx, Ly, H, Tbar, T

   
      Lx = myDGSEM % params % xScale   
      !$OMP PARALLEL
      CALL myDGSEM % CalculateStaticState( ) !! CPU Kernel
      !$OMP END PARALLEL
      ! ////////////////////////////////////////////////////////////////////////////////// !
      DO iEl = 1, myDGSEM % mesh % nElems
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
            
                  x = myDGSEM % mesh % geom(iEl) % x(i,j,k)
                  y = myDGSEM % mesh % geom(iEl) % y(i,j,k)
                  z = myDGSEM % mesh % geom(iEl) % z(i,j,k)
                  
                  IF( x <= 0.5_prec*Lx )THEN
                     T = -0.25_prec
                  ELSE
                     T = 0.25_prec
                  ENDIF
                  Tbar = myDGSEM % static % solution(i,j,k,5,iEl)/myDGSEM % static % solution(i,j,k,4,iEl)
                  myDGSEM % state % solution(i,j,k,4,iEl) = -myDGSEM % static % solution(i,j,k,4,iEl)*T/(Tbar + T)
                  
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
#ifdef HAVE_CUDA     
     ! Copy the static solution from host to device
     myDGSEM % static % solution_dev = myDGSEM % static % solution
     myDGSEM % state % solution_dev  = myDGSEM % state % solution
#endif
     !$OMP PARALLEL
     CALL myDGSEM % CalculateStaticBoundarySolution( )
     CALL myDGSEM % CalculateBoundarySolution( )
     !$OMP END PARALLEL
#ifdef HAVE_CUDA
     ! Copy the boundary solution from the device to the host
     myDGSEM % static % boundarySolution = myDGSEM % static % boundarySolution_dev
     myDGSEM % state % boundarySolution  = myDGSEM % state % boundarySolution_dev
#endif
     
 END SUBROUTINE InitialCondition
!
END PROGRAM Fluid_InitialConditions
