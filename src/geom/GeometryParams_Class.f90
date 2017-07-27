! GeometryParams_Class.f90
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
! GeometryParams_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

 
MODULE GeometryParams_Class


 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 

 IMPLICIT NONE


    TYPE GeometryParams
       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       CHARACTER(50) :: PeaceMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: geomPolyDeg
       INTEGER       :: nPlot
       INTEGER       :: nProc
       INTEGER       :: nXElems
       INTEGER       :: nYElems
       INTEGER       :: nZElems
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale
       REAL(prec)    :: zScale
      
       CONTAINS

       PROCEDURE :: Build => Build_GeometryParams

    END TYPE GeometryParams 
 

 CONTAINS


 SUBROUTINE Build_GeometryParams( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( GeometryParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       CHARACTER(50) :: PeaceMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nPlot
       INTEGER       :: nProc
       INTEGER       :: nXElems
       INTEGER       :: nYElems
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale


      NAMELIST / SpaceManagement / SpecMeshFile, PeaceMeshFile, polyDeg, nPlot, nProc, nXElems, nYElems, xScale, yScale
 
      SpecMeshFile = nada
      PeaceMeshFile = nada
      polyDeg = 5
      nPlot = 10
      nProc = 2
      nXElems = 4
      nYElems = 4
      xScale = ONE
      yScale = ONE 

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'geometry.params')
         READ( UNIT = nUnit, NML = SpaceManagement )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SpaceManagement )
      
      ! SpaceManagement (Default to isoparametric elements - if overintegration is USEd, default is to double the number of points )
      ! Plotting is defaulted to 10 evenly spaced points
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % PeaceMeshFile = PeaceMeshFile
      thisParam % polyDeg = polyDeg
      thisParam % nPlot = nPlot
      thisParam % nProc = nProc
      thisParam % nXElems = nXElems
      thisParam % nYElems = nYElems
      thisParam % dxPlot = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale = xScale
      thisParam % yScale = yScale
      


 END SUBROUTINE Build_GeometryParams

END MODULE GeometryParams_Class
