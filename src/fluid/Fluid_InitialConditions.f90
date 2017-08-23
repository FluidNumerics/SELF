! Fluid_InitialConditions.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 
PROGRAM Fluid_InitialConditions

! src/common/
USE ModelPrecision
! src/geom/
!USE HexMesh_Class
! src/highend/euler/
USE FluidParams_Class
! src/highend/euler/mpi-cuda
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
 
      CALL InitialCondition( myeu )
      
      PRINT*, "Reset Boundary conditions"
      CALL ResetBoundaryConditions( myeu, myRank )
      PRINT*, "DONE!"
      
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
      
#ifdef HAVE_MPI
      CALL MPI_BARRIER( )
#endif
      
      CALL myeu % Trash( )

#ifdef HAVE_MPI
      CALL MPI_FINALIZE( mpiErr )
#endif

 CONTAINS
 SUBROUTINE ResetBoundaryConditions( myDGSEM, myRank )
  
   TYPE( Fluid ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)          :: myRank
   ! Local
   INTEGER :: iFace, iFace2, e1, e2, s1, p2
   
   
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
!
 SUBROUTINE InitialCondition( myDGSEM )

   TYPE( Fluid ), INTENT(inout) :: myDGSEM
   ! Local
   INTEGER    :: i, j, k, iEl, iFace
   REAL(prec) :: x, y, z, r, hScale, Lx, Ly, H, T, Tbar

   
      Lx = myDGSEM % params % xScale
      Ly = myDGSEM % params % yScale
      H  = myDGSEM % params % zScale
      
      myDGSEM % state % solution = 0.0_prec
      
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
                  
                  r = sqrt( ( x-0.5_prec*Lx )**2 + ( y-0.5_prec*Ly )**2 + ( z-0.25_prec*H )**2 )
                
                  IF( r <= 250.0_prec )THEN
                     T = 0.25_prec*(1.0_prec + cos( pi*r/250.0_prec ) ) ! Potential temperature anomaly
                     Tbar = myDGSEM % static % solution(i,j,k,5,iEl)/myDGSEM % static % solution(i,j,k,4,iEl)
                     myDGSEM % state % solution(i,j,k,4,iEl) = -myDGSEM % static % solution(i,j,k,4,iEl)*T/(Tbar + T)
                  ENDIF
                  
                  ! Drag profile
                  !myDGSEM % dragProfile(i,j,k,iEl) = myDGSEM % params % Cd*exp( -z/myDGSEM % params % dragScale )

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
#ifdef HAVE_CUDA
     myDGSEM % state % solution_dev = myDGSEM % state % solution 
#endif
     !$OMP PARALLEL
     CALL myDGSEM % EquationOfState( )
     !$OMP END PARALLEL
     
#ifdef HAVE_CUDA
     myDGSEM % state % solution = myDGSEM % state % solution_dev
#endif
      
 END SUBROUTINE InitialCondition
!
END PROGRAM Fluid_InitialConditions
