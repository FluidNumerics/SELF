! BoundaryCommunicator_Class.f90
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
! BoundaryCommunicator_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file BoundaryCommunicator_Class.f90
!! Contains the \ref BoundaryCommunicator_Class module, and <BR>
!! defines the BoundaryCommunicator data-structure.

!> \defgroup BoundaryCommunicator_Class BoundaryCommunicator_Class 
!! This module defines the BoundaryCommunicator data-structure and its associated routines.
MODULE BoundaryCommunicator_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines

IMPLICIT NONE


!> \addtogroup BoundaryCommunicator_Class 
!! @{

!> \struct BoundaryCommunicator
!! The BoundaryCommunicator class provides a convenient package of attributes for implementing
!! boundary conditions.
!!
!! This structure was motivated by the need for a robust means of implementing boundary conditions
!! on an unstructured mesh. This class makes it trivial to implement message-passing for
!! MPI parallelism.
!!
!! <H2> BoundaryCommunicator </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> nBoundaries <td> INTEGER  <td> The number of boundary edges obtained from a 
!!                                                   search through mesh edges
!!       <tr> <th> extElemIDs(1:nBoundaries) <td> INTEGER  <td> The local element ID's of external
!!                          elements. These are either given as boundary condition flags or as the 
!!                          local element ID of a neighboring processes mesh
!!       <tr> <th> extProcIDs(1:nBoundaries) <td> INTEGER  <td> Process ID for the neighboring element.
!!       <tr> <th> boundaryIDs(1:nBoundaries) <td> INTEGER  <td> Edge ID corresponding to a given
!!                                                                      boundary edge ID.
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref BoundaryCommunicator_Class for more information. The first column lists the "call-name" and 
!!    the second column lists the name of routine that is aliased onto the call-name.
!!    <table>
!!       <tr> <th> Initialize <td> Initialize_BoundaryCommunicator
!!       <tr> <th> Trash <td> Trash_BoundaryCommunicator
!!       <tr> <th> ReadPickup <td> ReadPickup_BoundaryCommunicator
!!       <tr> <th> WritePickup <td> WritePickup_BoundaryCommunicator
!!    </table>
!!

!>@}
   TYPE BoundaryCommunicator
      INTEGER                               :: nBoundaries
      INTEGER, ALLOCATABLE                  :: extProcIDs(:)
      INTEGER, ALLOCATABLE                  :: boundaryIDs(:)
      INTEGER, ALLOCATABLE                  :: unPackMap(:)

      CONTAINS

      PROCEDURE :: Initialize => Initialize_BoundaryCommunicator
      PROCEDURE :: Trash      => Trash_BoundaryCommunicator
      
      PROCEDURE :: ReadPickup  => ReadPickup_BoundaryCommunicator
      PROCEDURE :: WritePickup => WritePickup_BoundaryCommunicator

   END TYPE BoundaryCommunicator

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup BoundaryCommunicator_Class
!! @{ 
! ================================================================================================ !
! S/R Initialize
! 
!> \fn Initialize_BoundaryCommunicator 
!! Allocates space for the BoundaryCommunicator structure and initializes all array values to zero.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(BoundaryCommunicator) :: this <BR>
!! <B>INTEGER</B>                    :: nBe <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Initialize( nBe ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myBC <td> BoundaryCommunicator <td> 
!!   <tr> <td> in <th> N <td> INTEGER <td> Polynomial degree for the solution storage
!!   <tr> <td> in <th> nBe <td> INTEGER <td> The number of boundary edges in the mesh
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Initialize_BoundaryCommunicator( myBC, nBe )

   IMPLICIT NONE
   CLASS(BoundaryCommunicator), INTENT(inout) :: myBC
   INTEGER, INTENT(in)                      :: nBe

      myBC % nBoundaries = nBe

      ALLOCATE( myBC % extProcIDs(1:nBe) )
      ALLOCATE( myBC % boundaryIDs(1:nBe) )
      ALLOCATE( myBC % unPackMap(1:nBe) )

      myBC % extProcIDs  = 0
      myBC % boundaryIDs = 0
      myBC % unPackMap   = 0

 END SUBROUTINE Initialize_BoundaryCommunicator
!
!> \addtogroup BoundaryCommunicator_Class 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_BoundaryCommunicator  
!! Frees memory associated with the attributes of the BoundaryCommunicator data structure. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(BoundaryCommunicator) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myBC <td> BoundaryCommunicator <td> 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_BoundaryCommunicator( myBC )

   IMPLICIT NONE
   CLASS(BoundaryCommunicator), INTENT(inout) :: myBC

      DEALLOCATE( myBC % unPackMap, myBC % extProcIDs, myBC % boundaryIDs )


 END SUBROUTINE Trash_BoundaryCommunicator
!
!
!==================================================================================================!
!-------------------------------------- FILE I/O ROUTINES -----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup BoundaryCommunicator_Class
!! @{ 
! ================================================================================================ !
! S/R WritePickup 
! 
!> \fn WritePickup_BoundaryCommunicator  
!! Writes pickup files for the BoundaryCommunicator data structure. 
!! 
!! Given a file-name base (e.g. "foo"), this routine generates "foo.bcm" (2)
!! The .bcm file contains the boundary edge, external element, and external process information.
!! The file is an ASCII file. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(BoundaryCommunicator) :: this <BR>
!! <B>CHARACTER</B>                  :: filename
!!         .... <BR>
!!     <B>CALL</B> this % WritePickup( filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myBC <td> BoundaryCommunicator <td> Previously constructed boundary-
!!                               communicator data structure
!!   <tr> <td> in <th> filename <td> CHARACTER <td> File base-name for the pickup files  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE WritePickup_BoundaryCommunicator( myBC, filename )

   IMPLICIT NONE
   CLASS( BoundaryCommunicator ), INTENT(in) :: myBC
   CHARACTER(*), INTENT(in)                     :: filename
  ! LOCAL
   INTEGER       :: i, fUnit
     

      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = TRIM(filename)//'.bcm', &
            FORM   ='FORMATTED',&
            ACCESS ='SEQUENTIAL',&
            STATUS ='REPLACE',&
            ACTION ='WRITE' )

      WRITE( fUnit, * ) myBC % nBoundaries
      
      DO i = 1, myBC % nBoundaries

         WRITE( fUnit, * ) myBC % boundaryIDs(i), &
                           myBC % extProcIDs(i), &
                           myBC % unPackMap(i)

      ENDDO 

      CLOSE(fUnit)
      

 END SUBROUTINE WritePickup_BoundaryCommunicator
!
!> \addtogroup BoundaryCommunicator_Class
!! @{ 
! ================================================================================================ !
! S/R ReadPickup 
! 
!> \fn ReadPickup_BoundaryCommunicator  
!! Reads pickup files for the BoundaryCommunicator data structure. 
!! 
!! Given a file-name base (e.g. "foo"), this routine reads "foo.bcm"
!! The .bcm file contains the boundary edge, external element, and external process information.
!! The file is an ASCII file. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(BoundaryCommunicator) :: this <BR>
!! <B>CHARACTER</B>                     :: filename
!!         .... <BR>
!!     <B>CALL</B> this % ReadPickup( filename ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myBC <td> BoundaryCommunicator <td> Previously constructed boundary-
!!                               communicator data structure
!!   <tr> <td> in <th> filename <td> CHARACTER <td> File base-name for the pickup files  
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE ReadPickup_BoundaryCommunicator( myBC, filename )

   IMPLICIT NONE
   CLASS( BoundaryCommunicator ), INTENT(inout) :: myBC
   CHARACTER(*), INTENT(in)                     :: filename
  ! LOCAL
   INTEGER       :: i
   INTEGER       :: fUnit
   INTEGER       :: nBe


      !PRINT *, 'S/R ReadPickup : Reading "'//TRIM(filename)//'.bcm"'
      
      OPEN( UNIT   = NEWUNIT(fUnit), &
            FILE   = TRIM(filename)//'.bcm', &
            FORM   ='FORMATTED',&
            ACCESS ='SEQUENTIAL',&
            STATUS ='OLD',&
            ACTION ='READ' )

      READ( fUnit, * ) nBe

      CALL myBC % Initialize( nBe )

      DO i = 1, myBC % nBoundaries

         READ( fUnit, * ) myBC % boundaryIDs(i), &
                          myBC % extProcIDs(i), &
                          myBC % unPackMap(i)

      ENDDO 

      CLOSE(fUnit)

 END SUBROUTINE ReadPickup_BoundaryCommunicator
!
 END MODULE BoundaryCommunicator_Class



