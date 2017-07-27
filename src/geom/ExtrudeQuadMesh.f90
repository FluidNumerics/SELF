! ExtrudeQuadMesh.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
! Copyright 2016 Joseph Schoonover <jschoonover@lanl.gov>, Los Alamos National Laboratory
!
! The SELF and accompanying documentation were produced in part under the 
! support of Florida State University and the National Science Foundation 
! through Grant OCE-1049131 during 2015 and in part under the support of the 
! Center for Nonlinear Studies and the Department of Energy through the 
! LANL/LDRD program in 2016.
!
! ExtrudeQuadMesh.f90 is part of the Spectral Element Libraries in Fortran (SELF).
! 
! Licensed under the Apache License, Version 2.0 (the "License"); 
! You may obtain a copy of the License at 
!
! http://www.apache.org/licenses/LICENSE-2.0 
!
! Unless required by applicable law or agreed to in writing, software 
! distributed under the License is distributed on an "AS IS" BASIS, 
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and  
! limitations under the License.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM ExtrudeQuadMesh

! src/common/
USE ModelPrecision
! src/nodal/
USE NodalStorage_Class
! src/geom/
USE QuadMesh_Class
USE HexMesh_Class
! src/boundary/
USE BoundaryCommunicator_Class
! src/highend/shallowwater
USE Params_Class

 IMPLICIT NONE
 
   TYPE( RunParams )            :: params
   TYPE( QuadMesh )             :: qMesh
   TYPE( HexMesh )              :: hMesh
   TYPE( BoundaryCommunicator ) :: extComm
   TYPE( NodalStorage )         :: dgStorage
   INTEGER                      :: N, nPlot, nBe, iFace, e1,s1, e2
   INTEGER                      :: nHNodes, nHexElems, k, i, nID, elID, nIDs(1:4)
   INTEGER                      :: ii, jj, kk
   REAL(prec)                   :: x, y, z
   REAL(prec), ALLOCATABLE      :: s(:), sq(:)
   REAL(prec)                   :: zt(1:4), zb(1:4)
   
      CALL params % Build( )
      N     = params % polyDeg
      nPlot = params % nPlot
 
      ALLOCATE( s(0:params % nZelem), sq(0:N) )
      s = UniformPoints( ZERO, ONE, params % nZelem )

      CALL dGStorage % Build( N, nPlot, GAUSS, DG )
      sq = dgStorage % interp % s

      ! Load in the mesh of quadrilaterals (x,y)
      IF( TRIM( params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default mesh.'
         CALL qmesh % LoadDefaultMesh( dgStorage % interp, &
                                      params % nXelem, &
                                      params % nYelem )
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(params % SpecMeshFile)//'.'
         CALL  qmesh % ReadSpecMeshFile( dgStorage % interp, &
                                        params % SpecMeshFile )
      ENDIF
      
      nHexElems = qmesh % nElems * params % nZelem
      nHNodes = qmesh % nNodes * (params % nZelem + 1)

      CALL hMesh % Initialize( nHNodes, nHexElems, 1, N )

      ! Extrude the nodes
      DO k = 0, params % nZelem 
         DO i = 1, qmesh % nNodes

            nID = i + k*qmesh % nNodes

            hMesh % nodes(nID) % nodeId = nID

            x = qMesh % nodes(i) % x
            y = qMesh % nodes(i) % y
            hMesh % nodes(nID) % x = x
            hMesh % nodes(nID) % y = y

            ! Linear interpolation between the top and bottom surfaces
            hMesh % nodes(nID) % z = TopSurface( x, y )*s(k) + &
                                     BottomSurface( x, y )*(ONE-s(k))

            IF( k == 0 )THEN ! this is the bottom surface
               hMesh % nodes(nID) % nodeType = NO_NORMAL_FLOW
            ELSEIF( k == params % nZelem )THEN ! This is the top surface  
               hMesh % nodes(nID) % nodeType = RADIATION
            ELSE
               hMesh % nodes(nID) % nodeType = qMesh % nodes(i) % nodeType
            ENDIF

         ENDDO
      ENDDO

      ! Construct the elements
      DO k = 1, params % nZelem 
         DO i = 1, qmesh % nElems

            elID = i + (k-1)*qmesh % nElems
            nIDs    = qmesh % elements(i) % nodeIDs
            
            hmesh % elements(elID) % elementID    = elID
            hmesh % elements(elID) % nodeIDs(1:4) = nIDs + (k-1)*qmesh % nNodes ! bottom
            hmesh % elements(elID) % nodeIDs(5:8) = nIDs + (k)*qmesh % nNodes   ! top

            DO kk = 1, 4
               zb(kk) = hmesh % nodes(hmesh % elements(elID) % nodeIDs(kk)) % z
               zt(kk) = hmesh % nodes(hmesh % elements(elID) % nodeIDs(kk+4)) % z
            ENDDO


            ! Fill in some geometry
            DO kk = 0, N
               DO jj = 0, N
                  DO ii = 0, N

                     hmesh % geom(elID) % x(ii,jj,kk) = qmesh % geom(i) % x(ii,jj)
                     hmesh % geom(elID) % y(ii,jj,kk) = qmesh % geom(i) % y(ii,jj)

                     hmesh % geom(elID) % z(ii,jj,kk) = HALF*HALF*HALF*( &
                         (ONE-sq(ii))*(ONE-sq(jj))*( (ONE-sq(kk))*zb(1) + (ONE+sq(kk))*zt(1) ) + & 
                         (ONE+sq(ii))*(ONE-sq(jj))*( (ONE-sq(kk))*zb(2) + (ONE+sq(kk))*zt(2) ) + & 
                         (ONE+sq(ii))*(ONE+sq(jj))*( (ONE-sq(kk))*zb(3) + (ONE+sq(kk))*zt(3) ) + &
                         (ONE-sq(ii))*(ONE+sq(jj))*( (ONE-sq(kk))*zb(4) + (ONE+sq(kk))*zt(4) ) )

                  ENDDO
               ENDDO
            ENDDO

         ENDDO
      ENDDO

      CALL hmesh % ScaleTheMesh( dgStorage % interp, ONE, ONE, ONE  )
      ! Construct the element faces
      CALL hmesh % ConstructFaces( )
 
      ! Set up the boundary communicator
      nBe = 0
      DO iFace = 1, hmesh % nFaces
         e2 = hmesh % faces(iFace) % ElementIDs(2)
         IF( e2 <= 0 )THEN
            nBe = nBe + 1
         ENDIF
      ENDDO
      CALL extComm % Initialize( nBe )
      
      nBe = 0
      DO iFace = 1, hmesh % nFaces
         e2 = hmesh % faces(iFace) % ElementIDs(2)
         IF( e2 <= 0 )THEN
            nBe = nBe + 1
            extComm % boundaryIDs(nBe) = iFace
            hmesh % faces(iFace) % boundaryID = nBe

            ! Assign the appropriate boundary conditions
            e1 = hmesh % faces(iFace) % elementIDs(1) ! obtain the primary element ID
            s1 = hmesh % faces(iFace) % elementSides(1) ! and the side ID

           ! DO i = 1,4 ! Obtain the global node types
           !    k   = hmesh % faceMap(i,s1)  ! Local node ID
           !    nID = hmesh % elements(e1) % nodeIDs( k ) ! Global node ID
           !    e2 = hmesh % nodes(nID) % nodeType
           !    IF( hmesh % nodes(nID) % nodeType == NO_NORMAL_FLOW ) THEN
           !       EXIT
           !    ENDIF
           ! ENDDO

            extComm % extElemIDs(nBe)  =  NO_NORMAL_FLOW !e2
            hmesh % faces(iFace) % ElementIDs(2) = NO_NORMAL_FLOW !e2

         ENDIF

      ENDDO

      CALL hmesh % WriteTecplot( )
      CALL hmesh % WritePeaceMeshFile( params % PeaceMeshFile )
      CALL extComm % WritePickup( 'ExtComm' )
    
      CALL qmesh % Trash( )
      CALL hmesh % Trash( )
      CALL extComm % Trash( )

      DEALLOCATE( s, sq )

 CONTAINS

 FUNCTION TopSurface( x, y ) RESULT( z )
   IMPLICIT NONE
   REAL(prec) :: x, y, z

      z = ONE

 END FUNCTION TopSurface

 FUNCTION BottomSurface( x, y ) RESULT( z )
   IMPLICIT NONE
   REAL(prec) :: x, y, z

      z = ZERO

 END FUNCTION BottomSurface
END PROGRAM ExtrudeQuadMesh
