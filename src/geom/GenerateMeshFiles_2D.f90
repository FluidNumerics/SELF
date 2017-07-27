! GenerateMeshFiles_2D.f90
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
! GenerateMeshFiles_2D.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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


PROGRAM GenerateMeshFiles_2D

! src/common/
USE ModelPrecision
! src/nodal/
USE NodalStorage_Class
! src/geom/
USE QuadMesh_Class
! src/boundary/
USE BoundaryCommunicator_Class
! src/highend/shallowwater
USE Params_Class

 IMPLICIT NONE
 
   TYPE( RunParams )            :: params
   TYPE( QuadMesh )             :: mesh
   TYPE( BoundaryCommunicator ) :: extComm
   TYPE( NodalStorage )         :: dgStorage
   INTEGER                      :: N, nPlot, nBe, iEdge, e2, quadType
   
      CALL params % Build( )
      N     = params % polyDeg
      nPlot = params % nPlot
      quadType = params % QuadType
      
      CALL dGStorage % Build( N, nPlot, quadType, DG )
      
      IF( TRIM( params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default mesh.'
         CALL mesh % LoadDefaultMesh( dgStorage % interp, &
                                      params % nXelem, &
                                      params % nYelem )
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(params % SpecMeshFile)//'.'
         CALL  mesh % ReadSpecMeshFile( dgStorage % interp, &
                                        params % SpecMeshFile )
      ENDIF
      
      nBe = 0
      DO iEdge = 1, mesh % nEdges
         e2 = mesh % edges(iEdge) % ElementIDs(2)
         IF( e2 <= 0 )THEN
            nBe = nBe + 1
         ENDIF
      ENDDO
      CALL extComm % Initialize( nBe )
      
      nBe = 0
      DO iEdge = 1, mesh % nEdges
         e2 = mesh % edges(iEdge) % ElementIDs(2)
         IF( e2 <= 0 )THEN
            nBe = nBe + 1
            extComm % boundaryIDs(nBe) = iEdge
            extComm % extElemIDs(nBe)  = e2
            mesh % edges(iEdge) % boundaryID = nBe
         ENDIF

      ENDDO

      CALL mesh % WriteTecplot( )
      CALL mesh % WritePeaceMeshFile( params % PeaceMeshFile )
      CALL extComm % WritePickup( 'ExtComm' )
    
      CALL mesh % Trash( )
      CALL extComm % Trash( )

 
END PROGRAM GenerateMeshFiles_2D
