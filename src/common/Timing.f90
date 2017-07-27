! Timing.f90
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
! Timing.f90 is part of the Spectral Element Libraries in Fortran (SELF).
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

!> \file Timing.f90
!! Contains the \ref Timing module, and <BR>
!! defines the \ref RoutineTimer and \ref MultiTimers data-structures.


!> \defgroup Timing Timing
!! This module defines the RoutineTimer and \ref MultiTimers data-structure and associated routines.
MODULE Timing

!src/common/
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
#ifdef HAVE_OPENMP
USE OMP_LIB
#endif
 IMPLICIT NONE

    INTEGER, PARAMETER :: timerNameLength = 30

!> \addtogroup Timing
!! @{

!> \struct RoutineTimer
!! A PRIVATE data-structure for handling timing of Fortran code.
!!
!! A RoutineTimer keeps track of a single start and stop time following calls to type-bound routines
!! that call the intrinsic "CPU_TIME( )". Additionally, upon issuing a call to the type-bound 
!! \ref StopTimer, the elapsed time (difference between start and stop time) is recorded 
!! is accumulated into an attribute of the RoutineTimer. The timer can be assigned a name to allow
!! for human-readable I/O. <BR>
!!
!! The RoutineTimer is meant to be used in a LinkedList of RoutineTimer's so that the end user
!! does not need to know a'priori how many timers they would like to instrument there code with.
!! (See \ref MultiTimers for more details)
!!
!! <H2> RoutineTimer </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> started <td> LOGICAL  <td> Flag to determine if the timer has been started
!!       <tr> <th> stopped <td> LOGICAL  <td> Flag to determine if the timer has been stopped
!!       <tr> <th> startTime <td> REAL(prec) <td> Time that CPU_TIME( ) was called via a call to 
!!                                                \ref StartTimer
!!       <tr> <th> stopTime <td> REAL(prec) <td> Time that CPU_TIME( ) was called via a call to 
!!                                                \ref StopTimer
!!       <tr> <th> accumulatedTime <td> REAL(prec) <td> Total time spent between calls to 
!!                                                      \ref StartTimer and \ref StopTimer
!!       <tr> <th> nObs <td> INTEGER <td> Number of times that \ref StopTimer has been called.
!!                                           This is the number of timings that have contributed
!!                                           to the accumulated time. 
!!       <tr> <th> timerID <td> INTEGER <td> A unique identifying integer for the timer. The timer-ID
!!                                           serves as a book-keeping device when implemented
!!                                           in a LinkedList (as in \ref MultiTimers)
!!                                           
!!       <tr> <th> whatYourTiming <td> CHARACTER() <td> A "human-readable" name for the timer.
!!       <tr> <th> next <td> TYPE(RoutineTimer) <td> A POINTER to another RoutineTimer so that a
!!                                                   LinkedList of RoutineTimers can be traversed.
!!    </table>
!!
!> @}
   TYPE RoutineTimer
      LOGICAL                      :: started = .FALSE., stopped = .FALSE.
      REAL(prec)                   :: startTime = ZERO
      REAL(prec)                   :: stopTime = ZERO
      REAL(prec)                   :: accumulatedTime = ZERO
      REAL(prec)                   :: nObs = ZERO
      INTEGER                      :: timerID
      CHARACTER(timerNameLength)   :: whatYourTiming
      TYPE( RoutineTimer ),POINTER :: next

      CONTAINS

      PROCEDURE :: SetName            => SetName_RoutineTimer
      PROCEDURE :: GetName            => GetName_RoutineTimer
      PROCEDURE :: SetAccumulatedTime => SetAccumulatedTime_RoutineTimer
      PROCEDURE :: GetAccumulatedTime => SetAccumulatedTime_RoutineTimer
  
      PROCEDURE :: StartTimer  => StartTimer_RoutineTimer
      PROCEDURE :: StopTimer   => StopTimer_RoutineTimer
      PROCEDURE :: ElapsedTime => ElapsedTime_RoutineTimer

   END TYPE RoutineTimer
!> \addtogroup Timing
!! @{

!> \struct MultiTimers
!! A Linked-List of \ref RoutineTimer 's.
!!
!! When instrumenting code with timers, it may not be known a'priori how many timers may be needed.
!! A Linked-List structure allows for a programmer to allocate space needed for each timer "on-the-
!! fly". The Linked-List structure provides pointers to the list "head" (start of the list), "tail"
!! (end of the list) and "current" position. Each entry of the list contains a pointer to the next
!! item in the list so that it can be traversed.  Additionally, through the use of the "timerID"
!! attribute of the RoutineTimer data-structure, one can also search for a specified timer.
!!
!! <H2> MultiTimers </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> head <td> RoutineTimer (POINTER)  <td> 
!!                           A pointer that points to the first RoutineTimer in the Linked-List
!!                           of timers
!!       <tr> <th> current <td> RoutineTimer (POINTER)  <td> 
!!                              A pointer that points to the last RoutineTimer in the Linked-List 
!!                              of timers
!!       <tr> <th> tail <td> RoutineTimer (POINTER) <td> 
!!                           A pointer that points to the last RoutineTimer in the Linked-List 
!!                           of timers
!!    </table>
!!
!! <H3> Procedures </H3>
!!    See \ref Timing for more information. The first column lists the "call-name" and the second
!!    column lists the name of routine that is aliased onto the call-name. This is the list of 
!!    PUBLIC type-bound procedures for the MultiTimers data-structure.
!!
!!    <table>
!!       <tr> <th> Build <td> Build_MultiTimers
!!       <tr> <th> Trash <td> Trash_MultiTimers
!!       <tr> <th> AddTimer <td> AddTimer_MultiTimers
!!       <tr> <th> StartTimer <td> StartTimer_MultiTimers
!!       <tr> <th> StopTimer <td>  StopTimer_MultiTimers
!!       <tr> <th> Write_MultiTimers <td>
!!    </table>
!!
!! <H3> Examples </H3>
!!    To use the MultiTimers in your code, first issue a call to "Build" and add a timer using
!!    "AddTimer". <BR>
!!     \verbatim
!!         TYPE(MultiTimers) :: timers
!!         ...
!!         ... -Some Code -
!!         ...
!!            CALL timers % Build( )
!!         ...
!!            ! Add a timer
!!            CALL timers % AddTimer( timername = 'Routine1', timerID = 1)
!!         ...
!!            ! Time a routine
!!            CALL timers % StartTimer( timerID = 1 )
!!               CALL Routine1( inputs, outputs )
!!            CALL timers % StopTimer( timerID = 1 )
!!         ...
!!            ! Write the results to the screen and to the file "Timing.stats"
!!            CALL timers % Write_MultiTimers( )
!!            ! Clear memory held by the timers
!!            CALL timers % Trash( )
!!
!!     \endverbatim
!> @}
   TYPE MultiTimers
      TYPE( RoutineTimer ), POINTER :: head, current, tail

      CONTAINS

      PROCEDURE :: Build      => Build_MultiTimers
      PROCEDURE :: Trash      => Trash_MultiTimers 
      PROCEDURE :: AddTimer   => AddTimer_MultiTimers 
      PROCEDURE :: StartTimer => StartTimer_MultiTimers
      PROCEDURE :: StopTimer  => StopTimer_MultiTimers
      PROCEDURE :: Write_MultiTimers

      PROCEDURE, PRIVATE :: SetName          => SetName_MultiTimers
      PROCEDURE, PRIVATE :: GetName          => GetName_MultiTimers
      PROCEDURE, PRIVATE :: ThereAreNoTimers
      PROCEDURE :: PointToTimer     => PointToTimer_MultiTimers
      PROCEDURE, PRIVATE :: MoveToNext       => MoveToNext_MultiTimers
      PROCEDURE, PRIVATE :: AccumulateTimings
        
    END TYPE MultiTimers

 INTEGER, PRIVATE, PARAMETER :: defaultNameLength = 40
 CONTAINS
!
!
!==================================================================================================!
!---------------------------- CONSTRUCTOR/DESTRUCTOR ROUTINES -------------------------------------!
!==================================================================================================!
!
!> \addtogroup Timing 
!! @{ 
! ================================================================================================ !
! S/R subroutine-name Build
!> \fn Build_MultiTimers  
!! Initializes the Linked-List of RoutineTimer's
!! 
!! The "head", "tail", and "current" RoutineTimer pointers initially point to NULL()
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MultiTimers) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Build( Inputs/Outputs ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> theTimers <td> TYPE( MultiTimers) <td> 
!!                         On output, the "head", "tail", and "current" attributes are nullified.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Build_MultiTimers( theTimers  )

   IMPLICIT NONE
   CLASS( MultiTimers ) :: theTimers
   
     theTimers % head => NULL( )
     ! Point the tail to null
     theTimers % tail => NULL()
     ! Set the current position to Null
     theTimers % current => NULL( )


 END SUBROUTINE Build_MultiTimers
!
!> \addtogroup Timing 
!! @{ 
! ================================================================================================ !
! S/R Trash 
! 
!> \fn Trash_MultiTimers  
!! Deallocates space for each RoutineTimer in the MultiTimers Linked-List.
!!  
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MultiTimers) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Trash( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> theTimers <td> TYPE(MultiTimers) <td> 
!!                         On <B>input</B>, the linked-list of RoutineTimers,<BR>
!!                         On <B>output</B>, the memory held by each RoutineTimer is freed and the
!!                         MultiTimer attributes are nullified.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Trash_MultiTimers( theTimers )

   IMPLICIT NONE
   CLASS( MultiTimers ) :: theTimers
   ! LOCAL 
   TYPE( RoutineTimer ), POINTER :: pNext
   
     ! Set the current position of the list to the head
     theTimers % current => theTimers % head
     
     ! Scroll through the list until the current position is nullified
     DO WHILE ( ASSOCIATED( theTimers % current ) )

        ! temporarily point to the next in the list
        pNext => theTimers % current % next 

        ! Deallocate memory pointed to by current position
        DEALLOCATE( theTimers % current ) 

        ! Update current position
        theTimers % current => pNext 

     ENDDO
      

 END SUBROUTINE Trash_MultiTimers
!
!
!==================================================================================================!
!---------------------------------- ACCESSOR ROUTINES ---------------------------------------------!
!==================================================================================================!
!
! 
!> \addtogroup Timing 
!! @{ 
! ================================================================================================ !
! S/R SetName
! 
!> \fn SetName_MultiTimers  
!! Sets the name of a routine timer given a "timerID".
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DataType) :: this <BR>
!! <B>CHARACTER</B>(len) :: timerName <BR>
!! <B>INTEGER</B>        :: timerID
!!         .... <BR>
!!     <B>CALL</B> this % SetName( timerName, timerID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> theTimers <td> TYPE(MultTimers) <td> Linked-list of RoutineTimer
!!   <tr> <td> in <th> timerName <td> CHARACTER(*) <td> 
!!                     Name of the timer that you want associated with the timerID.
!!   <tr> <td> in <th> timerID <td> INTEGER <td> 
!!                     Unique identifier for the RoutineTimer that you want to name
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE SetName_MultiTimers( theTimers, timerName, timerID )

   IMPLICIT NONE
   CLASS( MultiTimers ), INTENT(inout) :: theTimers
   CHARACTER(*), INTENT(IN)            :: timerName
   INTEGER, OPTIONAL, INTENT(IN)       :: timerID
   ! LOCAL
   
      IF( PRESENT(timerID) )THEN
      
         CALL theTimers % PointToTimer( timerID ) ! Point the "current" pointer to the RoutineTimer with 
                                                  ! timerID = timerID.
      ENDIF
      
      CALL theTimers % current % SetName( timerName )

 END SUBROUTINE SetName_MultiTimers
!
!> \addtogroup Timing 
!! @{ 
! ================================================================================================ !
! S/R GetName 
! 
!> \fn GetName_MultiTimers  
!! Returns the name of a RoutineTimer with the given timerID. 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DataType) :: this <BR>
!! <B>CHARACTER</B>(len) :: timerName <BR>
!! <B>INTEGER</B>        :: timerID
!!         .... <BR>
!!     <B>CALL</B> this % GetName( timerName, timerID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> theTimers <td> TYPE(MultTimers) <td> Linked-list of RoutineTimer
!!   <tr> <td> out <th> timerName <td> CHARACTER(*) <td> 
!!                      Name of the timer that is associated with the timerID.
!!   <tr> <td> in <th> timerID <td> INTEGER <td> 
!!                     Unique identifier for the RoutineTimer that you want to return the name of.
!!  </table> 
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE GetName_MultiTimers( theTimers, timerName, timerID )

   IMPLICIT NONE
   CLASS( MultiTimers )          :: theTimers
   CHARACTER(*), INTENT(OUT)     :: timerName
   INTEGER, OPTIONAL, INTENT(IN) :: timerID
   ! LOCAL
   
      IF( PRESENT(timerID) )THEN
         CALL theTimers % PointToTimer( timerID ) 
      ENDIF
      
      CALL theTimers % current % GetName( timerName )

 END SUBROUTINE GetName_MultiTimers
!
 SUBROUTINE SetName_RoutineTimer( theTimer, timerName )

   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(INOUT) :: theTimer
   CHARACTER(*), INTENT(IN)             :: timerName


      theTimer % whatYourTiming = timerName
      theTimer % nObs = ZERO

 END SUBROUTINE SetName_RoutineTimer
!
 SUBROUTINE GetName_RoutineTimer( theTimer, timerName )

   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(IN) :: theTimer
   CHARACTER(*), INTENT(OUT)         :: timerName

      timerName = theTimer % whatYourTiming

 END SUBROUTINE GetName_RoutineTimer
!
 SUBROUTINE SetAccumulatedTime_RoutineTimer( theTimer, accTime )

   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(INOUT)  :: theTimer
   REAL(PREC), INTENT(IN)                :: accTime


      theTimer % accumulatedTime = accTime

 END SUBROUTINE SetAccumulatedTime_RoutineTimer
!
 FUNCTION GetAccumulatedTime_RoutineTimer( theTimer ) RESULT( accTime )

   IMPLICIT NONE
   CLASS( RoutineTimer ) :: theTimer
   REAL(prec)            :: accTime

      accTime = theTimer % accumulatedTime

 END FUNCTION GetAccumulatedTime_RoutineTimer
!
!
!==================================================================================================!
!-------------------------------- Linked-List Type Operations -------------------------------------!
!==================================================================================================!
!
! 
 FUNCTION ThereAreNoTimers( theTimers ) RESULT( TorF )
  IMPLICIT NONE
  CLASS( MultiTimers ) :: theTimers
  LOGICAL              :: TorF

     TorF = .NOT.( ASSOCIATED( theTimers % head  ) )
     
 END FUNCTION ThereAreNoTimers
!
!> \addtogroup Timing 
!! @{ 
! ================================================================================================ !
! S/R AddTimer 
! 
!> \fn AddTimer_MultiTimers
!! Adds a RoutineTimer to the MultiTimers linked-list with the specified name and timerID. 
!!
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DataType) :: this <BR>
!! <B>CHARACTER</B>(len) :: timerName <BR>
!! <B>INTEGER</B>        :: timerID
!!         .... <BR>
!!     <B>CALL</B> this % AddTimer( timerName, timerID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> theTimers <td> TYPE(MultTimers) <td> Linked-list of RoutineTimer
!!   <tr> <td> out <th> timerName <td> CHARACTER(*) <td> 
!!                      Name of the timer that you want associated with the timerID.
!!   <tr> <td> in <th> timerID <td> INTEGER <td> 
!!                     Unique identifier for the new RoutineTimer.
!!  </table>   
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE AddTimer_MultiTimers( theTimers, timername, timerID )
 
   IMPLICIT NONE
   CLASS( MultiTimers ), INTENT(inout) :: theTimers
   CHARACTER(*)                        :: timername
   INTEGER                             :: timerID
   ! LOCAL
   INTEGER :: allocationStatus

     ! Check to see if this list is empty
     IF( theTimers % ThereAreNoTimers() )THEN
     
        ALLOCATE( theTimers % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE Timing.f90 : S/R AddTimer : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
      
        ! Point the current position to the head
        theTimers % current => theTimers % head
        ! Set the data
        CALL theTimers % SetName( timername )
        
        theTimers % current % timerID = timerID
        theTimers % current % nObs = ZERO
        
        ! Point the next to null and the tail to current
        theTimers % current % next => NULL( )
        theTimers % tail => theTimers % current
        
     ELSE ! the list is not empty
    
        ! Then we allocate space for the next item in the list    
        ALLOCATE( theTimers % tail % next, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE Timing.f90 : S/R AddTimer : Memory not allocated for next entry in list.'
           ! An exception handler should be built to handle these problems
        ENDIF      
        
        ! Reassign the tail
        theTimers % tail => theTimers % tail % next
        
        ! Set the current to the tail
        theTimers % current => theTimers % tail
  
        ! Fill in the data
        CALL theTimers % SetName( timername )
        
        ! Fill in the key information
        theTimers % current % timerID = timerID

        ! Point the next to null and the tail to current
        theTimers % current % next => NULL( )
        
     ENDIF

 END SUBROUTINE AddTimer_MultiTimers
!
 SUBROUTINE PointToTimer_MultiTimers( theTimers, timerID )
 
   IMPLICIT NONE
   CLASS( MultiTimers ) :: theTimers
   INTEGER, INTENT(IN)  :: timerID
   
   
      theTimers % current => theTimers % head ! Point to the head of the list
      
      DO WHILE(ASSOCIATED(theTimers % current))
      
         IF( theTimers % current % timerID == timerID )EXIT
         
         CALL theTimers % MoveToNext( )
      
      ENDDO
      
      
 END SUBROUTINE PointToTimer_MultiTimers
!
 SUBROUTINE MoveToNext_MultiTimers( theTimers )
 
  IMPLICIT NONE
  CLASS( MultiTimers ) :: theTimers

     theTimers % current => theTimers % current % next

 END SUBROUTINE MoveToNext_MultiTimers
!
!> \addtogroup Timing 
!! @{ 
! ================================================================================================ !
! S/R StartTimer 
! 
!> \fn StartTimer_MultiTimers  
!! Starts a RoutineTimer with the specified timerID. 
!!
!! This routine issues a call to CPU_TIME( ) "under the hood", assigning the "startTime" attribute
!! of the RoutineTimer. Additionally, the "started" attribute is set to .TRUE. and the "stopped"
!! attribute is set to .FALSE.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DataType) :: this <BR>
!! <B>INTEGER</B>        :: timerID
!!         .... <BR>
!!     <B>CALL</B> this % StartTimer( timerID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> theTimers <td> TYPE(MultTimers) <td> 
!!                         Linked-list of RoutineTimer. On output, the "current" attribute points
!!                         to the timer with the given timerID.
!!   <tr> <td> in <th> timerID <td> INTEGER <td> 
!!                     Unique identifier for the RoutineTimer that you want to start.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE StartTimer_MultiTimers( theTimers, timerID )

   IMPLICIT NONE
   CLASS( MultiTimers ), INTENT(inout) :: theTimers
   INTEGER, INTENT(in)                 :: timerID

      CALL theTimers % PointToTimer( timerID )
      CALL theTimers % current % StartTimer( )

 END SUBROUTINE StartTimer_MultiTimers
!
!> \addtogroup Timing 
!! @{ 
! ================================================================================================ !
! S/R StopTimer 
! 
!> \fn StopTimer_MultiTimers  
!! Stops a RoutineTimer with the specified timerID. 
!!
!! This routine issues a call to CPU_TIME( ) "under the hood", assigning the "stopime" attribute
!! of the RoutineTimer. Additionally, the  "stopped" attribute is set to .FALSE.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(DataType) :: this <BR>
!! <B>INTEGER</B>        :: timerID
!!         .... <BR>
!!     <B>CALL</B> this % StartTimer( timerID ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> theTimers <td> TYPE(MultTimers) <td> 
!!                         Linked-list of RoutineTimer. On output, the "current" attribute points
!!                         to the timer with the given timerID.
!!   <tr> <td> in <th> timerID <td> INTEGER <td> 
!!                     Unique identifier for the RoutineTimer that you want to start.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE StopTimer_MultiTimers( theTimers, timerID )

   IMPLICIT NONE
   CLASS( MultiTimers ), INTENT(inout) :: theTimers
   INTEGER, INTENT(in)                 :: timerID

      CALL theTimers % PointToTimer( timerID )
      CALL theTimers % current % StopTimer( )
      CALL theTimers % AccumulateTimings( )

 END SUBROUTINE StopTimer_MultiTimers
!
 SUBROUTINE StartTimer_RoutineTimer( thetimer )

   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(inout) :: theTimer

      theTimer % started = .TRUE.
      theTimer % stopped = .FALSE.
      
#ifdef HAVE_OPENMP
      theTimer % startTime = OMP_GET_WTIME( )
#else
      CALL CPU_TIME( theTimer % startTime )
#endif     
 
 END SUBROUTINE StartTimer_RoutineTimer
!
 SUBROUTINE StopTimer_RoutineTimer( thetimer )

   IMPLICIT NONE
   CLASS( RoutineTimer ), INTENT(inout) :: theTimer
   
#ifdef HAVE_OPENMP
      theTimer % stopTime = OMP_GET_WTIME( )
#else
      CALL CPU_TIME( theTimer % stopTime )
#endif
      theTimer % stopped = .TRUE.

 END SUBROUTINE StopTimer_RoutineTimer
!
 FUNCTION ElapsedTime_RoutineTimer( theTimer ) RESULT( eTime )

   IMPLICIT NONE
   CLASS( RoutineTimer ) :: theTimer
   REAL(prec)            :: eTime

      IF( theTimer % stopped )THEN
         eTime = theTimer % stopTime - theTimer % startTime 
      ELSE
         PRINT*, 'Module Timing.f90 : S/R ElapsedTime : Warning! Timer "', TRIM(theTimer % whatYourTiming),'" is not stopped'
         eTime = ZERO 
      ENDIF
 END FUNCTION ElapsedTime_RoutineTimer
!
!
!
 SUBROUTINE AccumulateTimings( theTimers, timerID )

   IMPLICIT NONE
   CLASS( MultiTimers )          :: theTimers
   INTEGER, INTENT(IN), OPTIONAL :: timerID
   ! LOCAL
   
    
      IF( PRESENT(timerID) )THEN
      
         CALL theTimers % PointToTimer( timerID ) ! Point the "current" pointer to the RoutineTimer with 
                                              ! timerID = timerID.
      ENDIF
  
  
  
      theTimers % current % accumulatedTime = theTimers % current % accumulatedTime + &
                                              theTimers % current % ElapsedTime( )

      theTimers % current % nObs = theTimers % current % nObs + ONE 


 END SUBROUTINE AccumulateTimings
!
!> \addtogroup Timing 
!! @{ 
! ================================================================================================ !
! S/R Write_MultiTimers 
! 
!> \fn Write_MultiTimers  
!! Writes all of the RoutineTimer attributes to the screen and to the file "Timing.stats" 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(MuliTimers) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Write_MultiTimers( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> theTimers <td> MultiTimers <td> Linked-list of RoutineTimer. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Write_MultiTimers( theTimers )

   IMPLICIT NONE
   CLASS( MultiTimers ), INTENT(INOUT) :: theTimers 
   ! LOCAL
   INTEGER :: k, fUnit
   CHARACTER(defaultNameLength) :: tName

      OPEN( UNIT = NewUnit(fUnit), FILE = 'Timing.stats' ) 

      WRITE(fUnit,*) '====================== Timing Results ======================'
      WRITE(fUnit,*) ' '
      
      theTimers % current => theTimers % head
      k = 0
      DO WHILE( ASSOCIATED(theTimers % current) )
         k = k+1

         CALL theTimers % GetName( tName )

         WRITE(fUnit,*) tName
         WRITE(fUnit,*) 'Number of Measurements : ', theTimers % current % nObs
         WRITE(fUnit,*) 'Accumulated Time       : ', theTimers % current % accumulatedTime
         WRITE(fUnit,*) 'Average Time           : ', theTimers % current % accumulatedTime/theTimers % current % nObs
         WRITE(fUnit,*) '------------------------------------------------------------'

         CALL theTimers % MoveToNext( )

      ENDDO
      CLOSE(fUnit)


      theTimers % current => theTimers % head
      
      PRINT*, '====================== Timing Results ======================'
      PRINT*, ' '
      
      DO WHILE( ASSOCIATED(theTimers % current) )

         CALL theTimers % GetName( tName ) 

         PRINT*, tName
         PRINT*, 'Number of Measurements : ', theTimers % current % nObs
         PRINT*, 'Accumulated Time       : ', theTimers % current % accumulatedTime
         PRINT*, 'Average Time           : ', theTimers % current % accumulatedTime/theTimers % current % nObs
         PRINT*, '------------------------------------------------------------'

         CALL theTimers % MoveToNext( )

      ENDDO

      PRINT*, '------------------------------------------------------------'
      PRINT*,' Timing Results saved to  "Timing.Stats" '

 END SUBROUTINE Write_MultiTimers

END MODULE Timing
