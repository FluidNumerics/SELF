! Fluid_InitialConditions.f90
! 
! Copyright 2017 Joseph Schoonover <schoonover.numerics@gmail.com>
!
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
      
      CALL InitialCondition( myeu, myRank )

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

            IF( ABS(myDGSEM % mesh % geom(e1) % nHat(3, 0, 0, s1) ) >&
                SQRT( myDGSEM % mesh % geom(e1) % nHat(1, 0, 0, s1)**2 + &
                      myDGSEM % mesh % geom(e1) % nHat(2, 0, 0, s1)**2 ) ) THEN ! Top or bottom surface

               IF( myDGSEM % mesh % geom(e1) % nHat(3, 0, 0, s1) < 0.0_prec ) THEN ! Bottom
                  myDGSEM % mesh % faces(iFace2) % elementIDs(2) = NO_NORMAL_FLOW

               ELSE ! Top
                  myDGSEM % mesh % faces(iFace2) % elementIDs(2) = RADIATION
               ENDIF

            ENDIF
               myDGSEM % mesh % faces(iFace2) % elementIDs(2) = RADIATION

         ENDIF

      ENDDO


 END SUBROUTINE ResetBoundaryConditions
!
 SUBROUTINE InitialCondition( myDGSEM, myRank )
   IMPLICIT NONE
   TYPE( Fluid ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)                             :: myRank
   ! Local
   INTEGER    :: i, j, k, iEl, iFace
   REAL(prec) :: x, y, z, Lx, Ly, T

   
      Lx = myDGSEM % params % xScale   
      Ly = myDGSEM % params % yScale   

      !$OMP PARALLEL
      CALL myDGSEM % CalculateStaticState( ) !! CPU Kernel
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      CALL myDGSEM % CalculateStaticBoundarySolution( )
      !$OMP END PARALLEL
     
      ! ////////////////////////////////////////////////////////////////////////////////// !
      DO iEl = 1, myDGSEM % mesh % nElems
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
            
                  x = myDGSEM % mesh % geom(iEl) % x(i,j,k)
                  y = myDGSEM % mesh % geom(iEl) % y(i,j,k)
                  z = myDGSEM % mesh % geom(iEl) % z(i,j,k)
                  
                 
                  T    = 100.0_prec*exp( -( (x-0.5_prec*Lx)**2 + &
                                            (y-0.5_prec*Ly)**2 + &
                                            (z-1000.0_prec)**2 )/&
                                            (2.0_prec*100.0_prec**2) )

                  myDGSEM % state % solution(i,j,k,4,iEl) = 10.0_prec*exp( -( (x-0.5_prec*Lx)**2 + &
                                                                              (y-0.5_prec*Ly)**2 + &
                                                                              (z-1000.0_prec)**2 )/&
                                                                              (2.0_prec*100.0_prec**2) )

                  ! Set the anomalous potential temperature
                  myDGSEM % state % solution(i,j,k,5,iEl) = ( myDGSEM % static % solution(i,j,k,4,iEl) + &
                                                              myDGSEM % state % solution(i,j,k,4,iEl) )*&
                                                            ( T +  myDGSEM % static % solution(i,j,k,5,iEl)/&
                                                                   myDGSEM % static % solution(i,j,k,4,iEl) )-&
                                                            myDGSEM % static % solution(i,j,k,5,iEl)
                                                            

                  ! Drag profile
                  myDGSEM % dragProfile(i,j,k,iEl) = myDGSEM % params % Cd*exp( -z/myDGSEM % params % dragScale )
                  
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
#ifdef HAVE_CUDA     
     myDGSEM % state % solution_dev  = myDGSEM % state % solution
#endif
     CALL myDGSEM % EquationOfState( ) ! Calculate the pressure
#ifdef HAVE_CUDA     
     myDGSEM % state % solution = myDGSEM % state % solution_dev
#endif

     
 END SUBROUTINE InitialCondition
!
 END PROGRAM Fluid_InitialConditions
