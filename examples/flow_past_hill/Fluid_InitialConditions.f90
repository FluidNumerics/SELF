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
 INTEGER       :: mpiErr
 CHARACTER(4)  :: rankChar
 LOGICAL       :: setupSuccess

#ifdef HAVE_MPI
      ! MPI Initialization
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myeu % myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, myeu % nProc, mpiErr )
      ! Sanity check
      PRINT*, 'Fluid_InitialConditions_MPI : Greetings from Process ', myeu % myRank, ' of ',myeu % nProc
#else
      myeu % myRank = 0
      myeu % nProc  = 1
#endif
      CALL myeu % Build( setupSuccess )
 
      IF( SetupSuccess )THEN
         CALL InitialCondition( myeu, myeu % myRank )
         
         PRINT*, "Reset Boundary conditions"
         CALL ResetBoundaryConditions( myeu, myeu % myRank )
         PRINT*, "DONE!"
         
         CALL myeu % WritePickup( ) 
         
         CALL myeu % WriteTecplot( )
         
         ! Before we write the mesh to file again, we need to "unscale" the mesh so that, upon running the 
         ! integrator, the mesh scaling is not applied a second time 
         CALL myeu % mesh % ScaleTheMesh( myeu % dgStorage % interp, &
                                             1.0_prec/myeu % params % xScale, &
                                             1.0_prec/myeu % params % yScale, &
                                             1.0_prec/myeu % params % zScale )
         WRITE( rankChar, '(I4.4)' )myeu % myRank
         CALL myeu % mesh % WritePeaceMeshFile( TRIM(myeu % params % PeaceMeshFile)//'.'//rankChar )
      
#ifdef HAVE_MPI
         CALL MPI_BARRIER( )
#endif
      

      ENDIF

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

               IF( myDGSEM % mesh % geom(e1) % nHat(3, 0, 0, s1) > 0.0_prec ) THEN ! Top
                  myDGSEM % mesh % faces(iFace2) % elementIDs(2) =PRESCRIBED
                  myDGSEM % prescribedState(:,:,1,iFace) = myDGSEM % static % boundarySolution(:,:,4,s1,e1)*myDGSEM % params % v0
               ELSE ! Bottom
                  myDGSEM % mesh % faces(iFace2) % elementIDs(2) = NO_NORMAL_FLOW
               ENDIF

            ENDIF

         ENDIF

      ENDDO


 END SUBROUTINE ResetBoundaryConditions
!
 SUBROUTINE InitialCondition( myDGSEM, myRank )
   IMPLICIT NONE
   TYPE( Fluid ), INTENT(inout) :: myDGSEM
   INTEGER, INTENT(in)          :: myRank
   ! Local
   INTEGER    :: i, j, k, iEl, iFace
   REAL(prec) :: x, y, z, Lx, Ly, r

   
      Lx = myDGSEM % params % xScale   
      Ly = myDGSEM % params % yScale   

      !$OMP PARALLEL
      CALL myDGSEM % CalculateStaticState( ) !! CPU Kernel
      !$OMP END PARALLEL

#ifdef HAVE_CUDA
      myDGSEM % static % solution_dev = myDGSEM % static % solution
#endif
      !$OMP PARALLEL
      CALL myDGSEM % CalculateStaticBoundarySolution( )
      !$OMP END PARALLEL
#ifdef HAVE_CUDA
      ! Copy the boundary solution from the device to the host
      myDGSEM % static % boundarySolution  = myDGSEM % static % boundarySolution_dev
#endif
     
      ! ////////////////////////////////////////////////////////////////////////////////// !
      DO iEl = 1, myDGSEM % mesh % nElems
         DO k = 0, myDGSEM % N
            DO j = 0, myDGSEM % N
               DO i = 0, myDGSEM % N
            
                  x = myDGSEM % mesh % geom(iEl) % x(i,j,k)
                  y = myDGSEM % mesh % geom(iEl) % y(i,j,k)
                  z = myDGSEM % mesh % geom(iEl) % z(i,j,k)
                  
                  myDGSEM % state % solution(i,j,k,1,iEl) = myDGSEM % static % solution(i,j,k,4,iEl)*myDGSEM % params % v0
                  
               ENDDO
            ENDDO
         ENDDO
      ENDDO

#ifdef HAVE_CUDA     
           ! Copy the static solution from host to device
           myDGSEM % state % solution_dev  = myDGSEM % state % solution
#endif
           !$OMP PARALLEL
           CALL myDGSEM % CalculateBoundarySolution( )
           !$OMP END PARALLEL
#ifdef HAVE_CUDA
           ! Copy the boundary solution from the device to the host
           myDGSEM % state % boundarySolution  = myDGSEM % state % boundarySolution_dev
#endif
           
      

     
 END SUBROUTINE InitialCondition
!
END PROGRAM Fluid_InitialConditions
