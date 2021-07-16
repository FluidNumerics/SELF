! SELF_HashTable.f90
! 
! Copyright 2021 Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file SELF_HashTable.f90
!! Contains the \ref SELF_HashTable module, and <BR>
!! defines the \ref HashTable data-structures.


!> \defgroup SELF_HashTable SELF_HashTable
!! This module defines the \ref HashTable data-structure and associated routines.

MODULE SELF_HashTable

 USE SELF_Constants
 USE SELF_LinkedList

 IMPLICIT NONE
!> \addtogroup SELF_HashTable
!! @{

!> \struct HashTable
!! A template hash-table data structure that is comprised of an array of linked lists ( see 
!! SELF_LinkedList ).
!!
!! <H2> HashTable </H2>
!! <H3> Attributes </H3>
!!    <table>
!!       <tr> <th> list(:) <td> LinkedList  <td> Single Linked-List
!!    </table>
!!
!> @}

   TYPE HashTable      
      TYPE( LinkedList ), ALLOCATABLE :: list(:)

      CONTAINS

      PROCEDURE :: Init => Init_HashTable
      PROCEDURE :: Free => Free_HashTable
      PROCEDURE :: AddDataForKeys => AddDataForKeys_HashTable
      PROCEDURE :: ContainsKeys => ContainsKeys_HashTable
      PROCEDURE :: GetDataForKeys => GetDataForKeys_HashTable

   END TYPE HashTable

   

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
!> \addtogroup SELF_HashTable
!! @{ 
! ================================================================================================ !
! S/R Init
! 
!> \fn Init_HashTable  
!! Allocates space for the hash table and initializes each linked list within the hash-table 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HashTable) :: this <BR>
!! <B>INTEGER</B>         :: N <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Init( N ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> out <th> myTable <td> HashTable <td> Initialized hash-table
!!   <tr> <td> in  <th> N <td> INTEGER <td> The number of linked-lists in the hash table
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Init_HashTable( myTable, N )
 
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  INTEGER            :: N
  ! LOCAL
  INTEGER :: i


     ALLOCATE( myTable % list(1:N) )

     DO i = 1,N
        CALL myTable % list(i) % Init( )
     ENDDO

 END SUBROUTINE Init_HashTable
!
!> \addtogroup SELF_HashTable
!! @{ 
! ================================================================================================ !
! S/R Free
! 
!> \fn Free_HashTable  
!!  Deallocates memory associated with the Hash-Table
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HashTable) :: this <BR>
!!         .... <BR>
!!     <B>CALL</B> this % Free( ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myTable <td> HashTable <td> On input, initialized hash-table <br>
!!                                                     On output, memory associated with hash-table
!!                                                     has been deallocated. 
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE Free_HashTable( myTable )
 
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  ! LOCAL
  INTEGER :: i, N


     N = SIZE( myTable % list )

     DO i = 1, N
        CALL myTable % list(i) % Free( )
     ENDDO

     DEALLOCATE( myTable % list )
  
 END SUBROUTINE Free_HashTable
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines  ----------------------------------------!
!==================================================================================================!
!
!
!> \addtogroup SELF_HashTable
!! @{ 
! ================================================================================================ !
! S/R AddDataForKeys
! 
!> \fn AddDataForKeys_HashTable
!! Adds data to the i-th linked list with a key set to "j" 
!! 
!! To access data in a hash-table, two keys need to be given. The first key points to the linked-list
!! in the hash table. The second key gives access to a particular record in the linked-list.
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HashTable) :: this <BR>
!! <B>INTEGER</B>         :: data, i, j <BR>
!!         .... <BR>
!!     <B>CALL</B> this % AddDataForKeys( inData, i, j ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in/out <th> myTable <td> HashTable <td> Initialized hash-table
!!   <tr> <td> in <th> inData <td> INTEGER <td> Data for the j-th record of the i-th linked-list in 
!!                                              the hash-table.
!!   <tr> <td> in <th> i <td> INTEGER <td> Key that points to a particular linked-list in the hash-table
!!   <tr> <td> in <th> j <td> INTEGER <td> Key for a particular record in the i-th linked-list.
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 SUBROUTINE AddDataForKeys_HashTable( myTable, inData, i, j )
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  INTEGER            :: inData, i, j
  
     CALL myTable % list(i) % AddToList( inData, j ) 

 END SUBROUTINE AddDataForKeys_HashTable
!
!> \addtogroup SELF_HashTable
!! @{ 
! ================================================================================================ !
! S/R ContainsKeys 
! 
!> \fn ContainsKeys_HashTable 
!! Checks to see if an entry in the hash-table exists for the keys "i" and "j". 
!! 
!! This function checks the i-th linked-list for a record with a key set to "j". If the keys exist,
!! the function returns a logical set to "TRUE". Otherwise, the function returns "FALSE". 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HashTable) :: this <BR>
!! <B>INTEGER</B>         :: i,j <BR>
!! <B>LOGICAL</B>         :: keysExist <BR>
!!         .... <BR>
!!     keysExist = this % ContainsKeys( i, j ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myTable <td> HashTable <td> Initialized hash-table
!!   <tr> <td> in <th> i <td> INTEGER <td> Key that points to a particular linked-list in the hash-table
!!   <tr> <td> in <th> j <td> INTEGER <td> Key for a particular record in the i-th linked-list.
!!   <tr> <td> out <th> doesContain <td> LOGICAL <td> Logical indicating if the keys exist in the hash-table
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
 FUNCTION ContainsKeys_HashTable( myTable, i, j ) RESULT( doesContain )
 
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  INTEGER            :: i, j
  LOGICAL            :: doesContain
  ! LOCAL
  INTEGER :: thisKey

     doesContain = .FALSE.
     IF(  myTable % list(i) % ListIsEmpty() )THEN ! this list hasn't been started
        RETURN
     ENDIF
 
     ! Rewind the list
     myTable % list(i) % current => myTable % list(i) % head

     DO WHILE( ASSOCIATED( myTable % list(i) % current ) )
    
        CALL myTable % list(i) % GetKey( thisKey )

        IF( thisKey == j )THEN ! This list already has this key
           doesContain =.TRUE.
           RETURN
        ENDIF

        ! other wise we move to the next element in the list
        CALL myTable % list(i) % MoveToNext( )

     ENDDO
 

 END FUNCTION ContainsKeys_HashTable
!
!> \addtogroup SELF_HashTable
!! @{ 
! ================================================================================================ !
! S/R GetDataForKeys 
! 
!> \fn GetDataForKeys_HashTable 
!! Returns data from the hash table associated with the keys "i" and "j". 
!! 
!! <H2> Usage : </H2> 
!! <B>TYPE</B>(HashTable) :: this <BR>
!! <B>INTEGER</B>         :: i, j, data <BR>
!!         .... <BR>
!!     <B>CALL<B> this % GetDataForKeys( data, i, j ) <BR>
!! 
!!  <H2> Parameters : </H2>
!!  <table> 
!!   <tr> <td> in <th> myTable <td> HashTable <td> Initialized hash-table
!!   <tr> <td> in <th> i <td> INTEGER <td> Key that points to a particular linked-list in the hash-table
!!   <tr> <td> in <th> j <td> INTEGER <td> Key for a particular record in the i-th linked-list.
!!   <tr> <td> out <th> outData <td> INTEGER <td> Data associated with keys "i" and "j"
!!  </table>  
!!   
! ================================================================================================ ! 
!>@}
  SUBROUTINE GetDataForKeys_HashTable( myTable, outData, i, j )
  
  IMPLICIT NONE
  CLASS( HashTable ) :: myTable
  INTEGER            :: outData, i, j
  ! LOCAL
  INTEGER :: thisData, thisKey
  
     IF(  myTable % list(i) % ListIsEmpty() )THEN ! this table entry is not pointing to a linked list
        PRINT*, 'MODULE HASHTABLE_CLASS : S/R GetDataForKeys :'
        PRINT*, 'List ', i,' is empty.'
        outData = fillValueInt
        RETURN
     ENDIF

     myTable % list(i) % current => myTable % list(i) % head

     DO WHILE( ASSOCIATED( myTable % list(i) % current ) ) ! this table entry does not contain the keys i,j

        CALL myTable % list(i) % GetData( thisData )
        CALL myTable % list(i) % GetKey( thisKey )

        IF( thisKey == j )THEN
           outData = thisData
           RETURN
        ENDIF
        
        CALL myTable % list(i) % MoveToNext( )
      
     ENDDO
   
 END SUBROUTINE GetDataForKeys_HashTable

END MODULE SELF_HashTable
