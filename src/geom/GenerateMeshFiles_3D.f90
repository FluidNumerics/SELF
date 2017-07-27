! GenerateMeshFiles_3D.f90
! 
! Copyright 2017 Joseph Schoonover <schoonover.numerics@gmail.com>
!
! GenerateMeshFiles_3D.f90 is part of the Spectral Element Libraries in Fortran (SELF).
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM GenerateMeshFiles_3D

! src/common/
USE ModelPrecision
! src/spectralops/
USE NodalStorage_Class
! src/geom/
USE HexMesh_Class
USE TopographicShapes
! src/boundary/
USE BoundaryCommunicator_Class
! src/templates
USE Params_Class

 IMPLICIT NONE
 
   TYPE( RunParams )            :: params
   TYPE( HexMesh )              :: mesh
   TYPE( BoundaryCommunicator ) :: extComm
   TYPE( NodalStorage )         :: dgStorage
   INTEGER                      :: N, nPlot, nBe, iFace, e2
   CHARACTER(4)                 :: pIDChar
   
      CALL params % Build( )
      N     = params % polyDeg
      nPlot = params % nPlot
 
      CALL dGStorage % Build( N, nPlot, GAUSS, DG )
      
      IF( params % topographicShape == Gaussian )THEN
         PRINT*, 'Using Gaussian Hill for topography'
         TopographicShape => GaussianHill
      ELSE
         TopographicShape => DefaultTopography
      ENDIF
      
      IF( TRIM( params % UCDMeshFile ) == nada )THEN
         IF( params % MeshType == DoublyPeriodic )THEN
            PRINT*,' Loading doubly periodic mesh.'
            CALL mesh % LoadDoublyPeriodicMesh( dgStorage % interp, &
                                                params % nXelem, &
                                                params % nYelem, &
                                                params % nZelem )
         ELSE
            PRINT*,' Loading default mesh.'
            CALL mesh % LoadDefaultMesh( dgStorage % interp, &
                                         params % nXelem, &
                                         params % nYelem, &
                                         params % nZelem )
         ENDIF
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(params % UCDMeshFile)//'.'
         CALL  mesh % ReadUCDMeshFile( dgStorage % interp, &
                                        params % UCDMeshFile )
      ENDIF
      
      nBe = 0
      DO iFace = 1, mesh % nFaces
         e2 = mesh % faces(iFace) % ElementIDs(2)
         IF( e2 <= 0 )THEN
            nBe = nBe + 1
         ENDIF
      ENDDO
      CALL extComm % Initialize( nBe )
      
      nBe = 0
      DO iFace = 1, mesh % nFaces
         e2 = mesh % faces(iFace) % ElementIDs(2)
         IF( e2 <= 0 )THEN
            nBe = nBe + 1
            extComm % boundaryIDs(nBe) = iFace
            extComm % extElemIDs(nBe)  = e2
            extComm % extProcIDs(nBe)  = 0
            mesh % faces(iFace) % boundaryID = nBe
         ENDIF

      ENDDO
      PRINT*, 'Found ',nBe, 'Boundary Faces'

      CALL mesh % WriteTecplot( )
      WRITE( pIDChar, '(I4.4)' ) 0
      CALL mesh % WriteTecplot( 'mesh.'//pIDChar )
      CALL mesh % WritePeaceMeshFile( TRIM(params % PeaceMeshFile)//'.'//pIDChar )
      CALL extComm % WritePickup( 'ExtComm.'//pIDChar )
    
      CALL mesh % Trash( )
      CALL extComm % Trash( )

 
END PROGRAM GenerateMeshFiles_3D
