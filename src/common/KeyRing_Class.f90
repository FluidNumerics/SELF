! KeyRing_Class.f90
! 
! Copyright 2015 Joseph Schoonover <schoonover.numerics@gmail.com>, The Florida State University
! Copyright 2016 Joseph Schoonover <jschoonover@lanl.gov>, Los Alamos National Laboratory
!
! The SELF and accompanying documentation were produced in part under the 
! support of Florida State University and the National Science Foundation 
! through Grant OCE-1049131 during 2015 and in part  the support of the 
! Center for Nonlinear Studies and the Department of Energy through the 
! LANL/LDRD program in 2016.
!
! KeyRing_Class.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

 
MODULE KeyRing_Class
! ========================================= Logs ================================================= !
! 2016-02-02  Joseph Schoonover : schoonover.numerics@gmail.com
!
! The KeyRing class was designed out of the necessity to have a structure that requires multiple
! "key information" in order to gain access to the record's data. In this class, the Record is now
!  the "NotchedKey", and the NotchedKey has "notches" (just as a physical key has notches that give
!  it access to a lock). By having the correct notches (not necessarily in the correct order), one can 
!  gain access to the data held by the key. A "KeyRing" is just a linked-list of keys.
!
!  The application which motivated the development of this module is the construction of the unique 
!  faces in a 3-D mesh given the element-to-corner-node connectivity. Each Face in the mesh is 
!  identifiable uniquely by its four corner nodes. The four corner nodes (in order according to the
!  primary element) make up the notches on the key. The data that is stored is the face ID number.
!  Quick searching through the set of keys is needed in order to create an efficient face-detection
!  algorithm. For this, we use a "KeyCabinet", which is an array of 
!  KeyRings. The index used to assign a KeyRing is the minimum node ID in the set of four nodes that 
!  make up a face. This way, when we query our list to see if a face has already been created, we
!  focus on a single KeyRing which limits our search list to a size less than or equal to the 
!  maximum nodal valence of the mesh. 
!   
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !

 USE CommonRoutines
 
 IMPLICIT NONE

    TYPE NotchedKey
      INTEGER, PRIVATE                     :: listData
      INTEGER, PRIVATE                     :: nNotches
      INTEGER, ALLOCATABLE, PRIVATE        :: notches(:)
      TYPE( NotchedKey ), POINTER          :: next

      CONTAINS
      
      PROCEDURE :: Build => Build_NotchedKey
      PROCEDURE :: Trash => Trash_NotchedKey

      PROCEDURE :: SetNotchesAndData => SetNotchesAndData_NotchedKey
      PROCEDURE :: GetNotchesAndData => GetNotchesAndData_NotchedKey
      PROCEDURE :: SetNotches => SetNotches_NotchedKey
      PROCEDURE :: GetNotches => GetNotches_NotchedKey
      PROCEDURE :: SetData => SetData_NotchedKey
      PROCEDURE :: GetData => GetData_NotchedKey

      PROCEDURE :: IsThisTheNotchedKey

    END TYPE NotchedKey

    TYPE KeyRing      
      TYPE( NotchedKey ), POINTER :: head, tail, current

      CONTAINS

      PROCEDURE :: Build => Build_KeyRing
      PROCEDURE :: Trash => Trash_KeyRing
      
      
      PROCEDURE :: GetData => GetCurrentData_KeyRing
      PROCEDURE :: SetData => SetCurrentData_KeyRing
      PROCEDURE :: GetNotches => GetCurrentNotches_KeyRing
      PROCEDURE :: SetNotches => SetCurrentNotches_KeyRing
      
      PROCEDURE :: ListIsEmpty => ListIsEmpty_KeyRing
      PROCEDURE :: AddToList => AddToList_KeyRing
      PROCEDURE :: FindDataForNotches => FindDataForNotches_KeyRing

      PROCEDURE :: MoveToNext => MoveToNext_KeyRing
      PROCEDURE :: MoveToHead => MoveToHead_KeyRing
      PROCEDURE :: MoveToTail => MoveToTail_KeyRing
      PROCEDURE :: PrintToScreen => PrintToScreen_KeyRing

    END TYPE KeyRing


! Setting up some parameters pertaining to this module
 INTEGER, PARAMETER, PRIVATE :: keyInc   = 1 ! The default increment in the NotchedKey Key
 INTEGER, PARAMETER, PRIVATE :: keyStart = 1 ! The default starting NotchedKey key
 
 CONTAINS
 
! 
! ================================================================================================ !
! ============================== Constructors and Destructors ==================================== !
! ================================================================================================ !
!
 SUBROUTINE Build_NotchedKey( myKey, nNotches, notches )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( NotchedKey ), INTENT(out) :: myKey 
   INTEGER, INTENT(in)              :: nNotches
   INTEGER, INTENT(in)              :: notches(1:nNotches)

    myKey % nNotches =  nNotches
  
    ALLOCATE( myKey % notches(1:nNotches) )
    
    myKey % notches = 0

   ! IF( PRESENT(notches) ) THEN
       myKey % notches = notches
   ! ENDIF

    myKey % next => NULL( )
  
 END SUBROUTINE Build_NotchedKey
!
!
!
 SUBROUTINE Trash_NotchedKey( myKey )
 ! S/R Trash
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( NotchedKey ), INTENT(inout) :: myKey 

    DEALLOCATE( myKey % notches )
  
 END SUBROUTINE Trash_NotchedKey
!
!
!
 SUBROUTINE Build_KeyRing( myList )
 ! S/R Build
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList

    
     myList % head => NULL( )

     myList % tail => NULL()

     myList % current => NULL( )
  
 END SUBROUTINE Build_KeyRing
!
!
!  
 SUBROUTINE Trash_KeyRing( myList )
 ! S/R Trash
 ! 
 !  This subroutine traverses through the linked list and frees the memory
 !  that was allocated for the entries in the linked list.  
 !
 ! ========================================================================== !
  IMPLICIT NONE
  CLASS( KeyRing ) :: myList
  ! LOCAL
  TYPE( NotchedKey ), POINTER :: pNext

     ! Set the current position of the list to the head
     myList % current => myList % head
     
     ! Scroll through the list until the current position is nullified
     DO WHILE ( ASSOCIATED( myList % current ) )

        ! temporarily point to the next in the list
        pNext => myList % current % next 

        CALL myList % current % Trash( )
        ! Deallocate memory pointed to by current position
        DEALLOCATE( myList % current ) 

        ! Update current position
        myList % current => pNext 

     ENDDO
  
 END SUBROUTINE Trash_KeyRing
! 
! ================================================================================================ !
! ======================================== Accessors ============================================= !
! ================================================================================================ !
!
! ------------------------------------------------------------------------------------------------ !
! -------------------------------------- Notched Key --------------------------------------------- !
! ------------------------------------------------------------------------------------------------ !
!
 SUBROUTINE SetNotchesAndData_NotchedKey( myKey, inNotches, inData )
 ! S/R SetNotchesAndData
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( NotchedKey ), INTENT(inout) :: myKey
   INTEGER, INTENT(in)                :: inNotches(1:myKey % nNotches)
   INTEGER, INTENT(in)                :: inData
   ! LOCAL
   LOGICAL :: keyMatches

    CALL myKey % SetNotches( inNotches )

    CALL myKey % SetData( inNotches, inData, keyMatches  )
    
 END SUBROUTINE SetNotchesAndData_NotchedKey
!
!
!
 SUBROUTINE GetNotchesAndData_NotchedKey( myKey, outNotches, outData )
 ! S/R GetNotchesAndData
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( NotchedKey ), INTENT(in) :: myKey
   INTEGER, INTENT(out)            :: outNotches(1:myKey % nNotches)
   INTEGER, INTENT(out)            :: outData
   ! LOCAL
   LOGICAL :: keyMatches
   INTEGER :: theseNotches(1:myKey % nNotches)

    CALL myKey % GetNotches( theseNotches )
    CALL myKey % GetData( theseNotches, outData, keyMatches )
    outNotches = theseNotches
    
 END SUBROUTINE GetNotchesAndData_NotchedKey
!
!
!
 SUBROUTINE SetData_NotchedKey( myKey, testNotches, inData, dataSet )
 ! S/R SetData
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( NotchedKey ), INTENT(inout) :: myKey
   INTEGER, INTENT(in)                :: testNotches(1:myKey % nNotches)
   INTEGER, INTENT(in)                :: inData
   LOGICAL, INTENT(out)               :: dataSet
   ! LOCAL
   LOGICAL :: keyMatches

    CALL myKey % IsThisTheNotchedKey( testNotches, keyMatches )

    dataSet = .FALSE.

    IF( keyMatches )THEN
       myKey % listData = inData
       dataSet          = .TRUE.
    ENDIF
    
 END SUBROUTINE SetData_NotchedKey
!
!
!
 SUBROUTINE GetData_NotchedKey( myKey, testNotches, outData, dataRetrieved )
 ! S/R GetData
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( NotchedKey ), INTENT(in) :: myKey
   INTEGER, INTENT(in)             :: testNotches(1:myKey % nNotches)
   INTEGER, INTENT(out)            :: outData
   LOGICAL, INTENT(out)            :: dataRetrieved
   ! LOCAL
   LOGICAL :: keyMatches

    CALL myKey % IsThisTheNotchedKey( testNotches, keyMatches )

    outData = 0
    dataRetrieved = .FALSE.

    IF( keyMatches )THEN
       outData       = myKey % listData
       dataRetrieved = .TRUE.
    ENDIF
    
 END SUBROUTINE GetData_NotchedKey
!
!
!
 SUBROUTINE SetNotches_NotchedKey( myKey, notches )
 ! S/R SetNotches
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( NotchedKey ), INTENT(inout) :: myKey
   INTEGER, INTENT(in)                :: notches(1:myKey % nNotches)

    myKey % notches = notches
    
 END SUBROUTINE SetNotches_NotchedKey
!
!
!
 SUBROUTINE GetNotches_NotchedKey( myKey, notches )
 ! S/R GetNotches
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( NotchedKey ), INTENT(in) :: myKey
   INTEGER, INTENT(out)            :: notches(1:myKey % nNotches)

    notches = myKey % notches
    
 END SUBROUTINE GetNotches_NotchedKey
!
! ------------------------------------------------------------------------------------------------ !
! ---------------------------------------- KeyRing ----------------------------------------------- !
! ------------------------------------------------------------------------------------------------ !
!
 SUBROUTINE SetCurrentData_KeyRing( myList, inData )
 ! S/R SetCurrentData
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList
   INTEGER             :: inData

    myList % current % listData = inData

 END SUBROUTINE SetCurrentData_KeyRing
!
!
!
 SUBROUTINE GetCurrentData_KeyRing( myList, outData )
 ! S/R GetCurrentData
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList
   INTEGER          :: outData

    outData = myList % current % listData

 END SUBROUTINE GetCurrentData_KeyRing
!
!
!
 SUBROUTINE SetCurrentNotches_KeyRing( myList, inNotches, nNotches )
 ! S/R SetCurrentKey
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList
   INTEGER          :: nNotches
   INTEGER          :: inNotches(1:nNotches)

    CALL myList % current % SetNotches( inNotches )

 END SUBROUTINE SetCurrentNotches_KeyRing
!
!
!
 SUBROUTINE GetCurrentNotches_KeyRing( myList, outNotches, nNotches )
 ! S/R GetCurrentNotches
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList
   INTEGER          :: nNotches
   INTEGER          :: outNotches(1:nNotches)

    CALL myList % current % GetNotches( outNotches )

 END SUBROUTINE GetCurrentNotches_KeyRing
! 
! ================================================================================================ !
! ================================ Data Structure Operations ===================================== !
! ================================================================================================ !
!
!
! ------------------------------------------------------------------------------------------------ !
! -------------------------------------- Notched Key --------------------------------------------- !
! ------------------------------------------------------------------------------------------------ !
!
 SUBROUTINE IsThisTheNotchedKey( myKey, testNotches, keyMatches )
 ! S/R IsThisTheNotchedKey
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   CLASS( NotchedKey ), INTENT(in) :: myKey
   INTEGER, INTENT(in)             :: testNotches(1:myKey % nNotches)
   LOGICAL, INTENT(out)            :: keyMatches
   ! LOCAL
   INTEGER :: theseNotches(1:myKey % nNotches)
   INTEGER :: thoseNotches(1:myKey % nNotches)
   INTEGER :: i

    CALL InsertionSort( myKey % notches, theseNotches, myKey % nNotches )
    CALL InsertionSort( testNotches, thoseNotches, myKey % nNotches )

    keyMatches = .TRUE.

    DO i = 1, myKey % nNotches
       IF( .NOT.( theseNotches(i) == thoseNotches(i) ) )THEN
          keyMatches = .FALSE.
          EXIT
       ENDIF
    ENDDO
 
 END SUBROUTINE IsThisTheNotchedKey
!
! ------------------------------------------------------------------------------------------------ !
! ---------------------------------------- KeyRing ----------------------------------------------- !
! ------------------------------------------------------------------------------------------------ !
!
 FUNCTION ListIsEmpty_KeyRing( myList ) RESULT( TorF)
 ! Function ListIsEmpty
 !
 !  This function checks if "myList" is an empty list. A logical is returned
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList
   LOGICAL          :: TorF

     TorF = .NOT.( ASSOCIATED( myList % head  ) )
     
 END FUNCTION ListIsEmpty_KeyRing
!
!
!
 SUBROUTINE AddToList_KeyRing( myList, inData, inNotches, nNotches )
 ! S/R AddToList
 !
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList
   INTEGER          :: inData
   INTEGER          :: nNotches
   INTEGER          :: inNotches(1:nNotches)
   ! LOCAL
   TYPE( NotchedKey ), POINTER :: previous
   INTEGER :: allocationStatus

    ! Check to see if this list is empty
    IF( myList % ListIsEmpty() )THEN
    
       ALLOCATE( myList % head, STAT = allocationStatus )
       IF( allocationStatus /=0 )THEN
          PRINT*, 'MODULE KeyRing_Class.f90 : S/R AddToList : Memory not allocated for next entry in list.'
          ! An exception handler should be built to handle these problems
       ENDIF      
      
       ! Point the current position to the head
       myList % current => myList % head
       ! Set the data
       CALL myList % current % Build( nNotches, inNotches )
       !CALL myList % SetData( inData )
       myList % current % listData = inData
       
       ! Point the next to null and the tail to current
       myList % current % next => NULL( )
       myList % tail => myList % current
       
    ELSE ! the list is not empty
    
       ! Then we allocate space for the next item in the list    
       ALLOCATE( myList % tail % next, STAT = allocationStatus )
       IF( allocationStatus /=0 )THEN
          PRINT*, 'MODULE KeyRing_Class.f90 : S/R AddToList : Memory not allocated for next entry in list.'
          ! An exception handler should be built to handle these problems
       ENDIF      
       
       !  Temporarily hold onto the tail in case we need the key (if inKey is not present)  
       previous => myList % tail
       ! Reassign the tail
       myList % tail => myList % tail % next
       
       ! Set the current to the tail
       myList % current => myList % tail
          ! Fill in the data
       CALL myList % current % Build( nNotches, inNotches )
      ! CALL myList % SetData( inData )
       myList % current % listData = inData
       
       ! Point the next to null and the tail to current
       myList % current % next => NULL( )
       
    ENDIF
 
 END SUBROUTINE AddToList_KeyRing
!
!
!
 SUBROUTINE FindDataForNotches_KeyRing( myList, inNotches, nNotches, outData, dataRetrieved )
 ! S/R FindDataForNotches
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ), INTENT(inout) :: myList
   INTEGER, INTENT(in)             :: nNotches 
   INTEGER, INTENT(in)             :: inNotches(1:nNotches)
   INTEGER, INTENT(out)            :: outData
   LOGICAL, INTENT(out)            :: dataRetrieved
   ! LOCAL 
   
    CALL myList % MoveToHead( )
    dataRetrieved = .FALSE.

    DO WHILE( ASSOCIATED(myList % current) )

       CALL myList % current % GetData( inNotches, outData, dataRetrieved )
       IF( .NOT.dataRetrieved )THEN

         CALL myList % MoveToNext( )
       ELSE
         EXIT
       ENDIF
    ENDDO
  
 END SUBROUTINE FindDataForNotches_KeyRing
!
!
!
 SUBROUTINE MoveToNext_KeyRing( myList )
 ! S/R MoveToNext
 !
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList

    myList % current => myList % current % next

 END SUBROUTINE MoveToNext_KeyRing
!
!
!
 SUBROUTINE MoveToHead_KeyRing( myList )
 ! S/R MoveToHead
 !
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList

    myList % current => myList % head

 END SUBROUTINE MoveToHead_KeyRing
!
!
!
 SUBROUTINE MoveToTail_KeyRing( myList )
 ! S/R MoveToTail
 !
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList

    myList % current => myList % tail

 END SUBROUTINE MoveToTail_KeyRing
!
!
!
 
 SUBROUTINE GetCount_KeyRing( myList, numberOfListItems )
 ! S/R GetCount
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing )  :: myList
   INTEGER, INTENT(out) :: numberOfListItems

    numberOfListItems = 0 ! Initialize the number of list items
    ! Check to see if this list is empty
    IF( myList % ListIsEmpty() )THEN
    
       PRINT*, 'Module KeyRing_Class.f90 : S/R ListCount : List is empty.'
       RETURN
       
    ELSE ! the list is not empty
   
       CALL myList % MoveToHead( ) ! Rewind the list
       
       DO WHILE( ASSOCIATED( myList % current ) )    
          numberOfListItems = numberOfListItems + 1
          CALL myList % moveToNext( )
        ENDDO

     ENDIF

 END SUBROUTINE GetCount_KeyRing
!
!
!
 SUBROUTINE PrintToScreen_KeyRing( myList )
 ! S/R PrintToScreen
 !
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( KeyRing ) :: myList

    myList % current => myList % head

    PRINT*, '          Data    |    Notches'
    DO WHILE( ASSOCIATED( myList % current ) )
 
       PRINT*, myList % current % listData,' | ', myList % current % notches
       CALL myList % MoveToNext()

    ENDDO

 END SUBROUTINE PrintToScreen_KeyRing
!
!
!
END MODULE KeyRing_Class
