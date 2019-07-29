! EquationParser.f90
! 
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! EquationParser defines a public class that can be used to parse and evaluate strings
! representative of equations. An equation, written in infix form, is converted to
! postfix form and evaluated using a postfix calculator.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE EquationParser_Class


USE ModelPrecision
USE ConstantsDictionary

IMPLICIT NONE

  INTEGER, PARAMETER :: Error_Message_Length = 256
  INTEGER, PARAMETER :: Max_Equation_Length  = 1024 
  INTEGER, PARAMETER :: Max_Function_Length  = 5
  INTEGER, PARAMETER :: Max_Variable_Length  = 12 
  INTEGER, PARAMETER :: Token_Length         = 48
  INTEGER, PARAMETER :: Stack_Length         = 128

  ! Token types 
  INTEGER, PARAMETER, PRIVATE :: None_Token               = 0
  INTEGER, PARAMETER, PRIVATE :: Number_Token             = 1
  INTEGER, PARAMETER, PRIVATE :: Variable_Token           = 2
  INTEGER, PARAMETER, PRIVATE :: Operator_Token           = 3
  INTEGER, PARAMETER, PRIVATE :: Function_Token           = 4
  INTEGER, PARAMETER, PRIVATE :: OpeningParentheses_Token = 5
  INTEGER, PARAMETER, PRIVATE :: ClosingParentheses_Token = 6
  INTEGER, PARAMETER, PRIVATE :: Monadic_Token            = 7

  INTEGER, PARAMETER, PRIVATE :: nFunctions = 14
  INTEGER, PARAMETER, PRIVATE :: nSeparators = 7

  TYPE String
    CHARACTER(10) :: str
  END TYPE String

  CHARACTER(1), DIMENSION(7), PRIVATE  :: separators = (/ "+", "-", "*", "/", "(", ")", "^" /) 
  CHARACTER(1), DIMENSION(5), PRIVATE  :: operators  = (/ "+", "-", "*", "/", "^" /) 
  CHARACTER(1), DIMENSION(10), PRIVATE :: numbers    = (/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" /) 
  TYPE(String), DIMENSION(14), PRIVATE :: functions 
  CHARACTER(1), DIMENSION(4), PRIVATE  :: variables  = (/ "x", "y", "z", "t" /)

  ! Private Types !

  TYPE Token  
    CHARACTER(Token_Length) :: tokenString
    INTEGER                 :: tokenType

    CONTAINS  
      PROCEDURE :: Equals_Token

  END TYPE Token


  TYPE TokenStack
    TYPE(Token), ALLOCATABLE :: tokens(:)
    INTEGER                  :: top_index = 0

    CONTAINS
  
      PROCEDURE :: Construct => Construct_TokenStack
      PROCEDURE :: Destruct  => Destruct_TokenStack
      
      PROCEDURE :: Push      => Push_TokenStack
      PROCEDURE :: Pop       => Pop_TokenStack
      PROCEDURE :: Peek      => Peek_TokenStack

      PROCEDURE :: IsEmpty   => IsEmpty_TokenStack
      PROCEDURE :: TopToken
  
  END TYPE TokenStack



  TYPE NumberStack
    REAL(prec), ALLOCATABLE :: tokens(:)
    INTEGER                 :: top_index

    CONTAINS
  
      PROCEDURE :: Construct => Construct_NumberStack
      PROCEDURE :: Destruct  => Destruct_NumberStack
      
      PROCEDURE :: Push      => Push_NumberStack
      PROCEDURE :: Pop       => Pop_NumberStack
      PROCEDURE :: Peek      => Peek_NumberStack

      PROCEDURE :: IsEmpty   => IsEmpty_NumberStack
  
  END TYPE NumberStack

  PRIVATE :: Token, TokenStack, NumberStack
  PRIVATE :: Construct_TokenStack, Destruct_TokenStack, Push_TokenStack, Pop_TokenStack, Peek_TokenStack, IsEmpty_TokenStack
  PRIVATE :: Construct_NumberStack, Destruct_NumberStack, Push_NumberStack, Pop_NumberStack, Peek_NumberStack, IsEmpty_NumberStack
  PRIVATE :: IsNumber, IsVariable, IsFunction, IsOperator, IsSeparator, FindLastFunctionIndex, F_of_X, Priority


  TYPE EquationParser
    CHARACTER(Max_Equation_Length)     :: equation
    CHARACTER(Max_Variable_Length)     :: variableName
    CHARACTER(Max_Equation_Length)     :: inFixFormula
    TYPE( TokenStack )                 :: inFix
    TYPE( TokenStack )                 :: postFix

    CONTAINS

      PROCEDURE :: CleanEquation
      PROCEDURE :: Tokenize
      PROCEDURE :: ConvertToPostfix

      PROCEDURE :: Evaluate

      PROCEDURE :: Print_InFixTokens
      PROCEDURE :: Print_PostFixTokens

  END TYPE EquationParser

  INTERFACE EquationParser
    PROCEDURE Construct_EquationParser
  END INTERFACE EquationParser
 

!  INTERFACE ASSIGNMENT (=) 
!    FUNCTION Equals_Token( tok1 ) RESULT( tok2 )
!      TYPE(Token), INTENT(in)  :: tok1
!      TYPE(Token), INTENT(out) :: tok2
!    END FUNCTION Equals_Token
!  END FUNCTION INTERFACE

CONTAINS

   FUNCTION Construct_EquationParser( equation ) RESULT( parser )
    TYPE( EquationParser ) :: parser
    CHARACTER(*)           :: equation
    ! Local
    CHARACTER(Error_Message_Length) :: errorMsg
    LOGICAL                         :: equationIsClean, tokenized, success

      functions(1) % str = "cos"
      functions(2) % str = "sin"
      functions(3) % str = "tan"
      functions(4) % str = "tanh"
      functions(5) % str = "sqrt"
      functions(6) % str = "abs"
      functions(7) % str = "exp"
      functions(8) % str = "ln"
      functions(9) % str = "log"
      functions(10) % str = "acos"
      functions(11) % str = "asin"
      functions(12) % str = "atan"
      functions(13) % str = "sech"
      functions(14) % str = "rand"

      parser % inFixFormula = " "
      parser % equation = equation
      errorMsg = " "

      CALL parser % CleanEquation( equationIsClean, errorMsg )

      IF( equationIsClean )THEN

        CALL parser % Tokenize( tokenized, errorMsg )

        IF( tokenized )THEN

          CALL parser % ConvertToPostFix( )
        
        ELSE

           PRINT*, TRIM( errorMsg )
           success = .false.

        ENDIF

      END IF
         
  END FUNCTION Construct_EquationParser

  SUBROUTINE CleanEquation( parser, equationCleaned, errorMsg )
    CLASS( EquationParser ), INTENT(inout)       :: parser
    LOGICAL, INTENT(out)                         :: equationCleaned
    CHARACTER(Error_Message_Length), INTENT(out) :: errorMsg
    ! Local
    INTEGER :: nChar, equalSignLoc, j, i
   

      equationCleaned = .FALSE.
      parser % variableName    = '#noname'

      nChar = LEN_TRIM( parser % equation )
      equalSignLoc = INDEX( parser % equation, "=" )

      IF( equalSignLoc == 0 )THEN
        errorMsg = "No equal sign found"
        RETURN
      ENDIF

      parser % variableName = TRIM( parser % equation(1:equalSignLoc-1) )

      ! Grab the formula to the right of the equal sign and left adjust the formula
      parser % inFixFormula = parser % equation(equalSignLoc+1:)
      parser % inFixFormula = ADJUSTL(parser % inFixFormula)
  

      ! Remove any spaces
      j = 1
      DO i = 1, LEN_TRIM(parser % inFixFormula)
        IF( parser % inFixFormula(i:i) /= " " )THEN
          parser % inFixFormula(j:j) = parser % inFixFormula(i:i)
          j = j + 1
        ENDIF
      ENDDO

      parser % inFixFormula(j:Max_Equation_Length) = " "
  
      equationCleaned = .TRUE.

  END SUBROUTINE CleanEquation

  SUBROUTINE Tokenize( parser, tokenized, errorMsg )
    CLASS( EquationParser ), INTENT(inout) :: parser
    LOGICAL, INTENT(out)                   :: tokenized
    CHARACTER(Error_Message_Length)        :: errorMsg
    ! Local
    INTEGER :: i, j 
 

      tokenized = .FALSE.
      errorMsg  = " "

      CALL parser % infix % Construct( Stack_Length )

      i = 1
      DO WHILE( parser % inFixFormula(i:i) /= " " )

        IF( IsVariable( parser % inFixFormula(i:i) ) )THEN
          
          parser % inFix % top_index = parser % inFix % top_index + 1
          parser % inFix % tokens( parser % inFix % top_index ) % tokenString = parser % inFixFormula(i:i)
          parser % inFix % tokens( parser % inFix % top_index ) % tokenType   = Variable_Token 
          i = i+1

          ! Next item must be an operator, closing parentheses, or end of equation

          IF( .NOT. IsOperator( parser % infixFormula(i:i) ) .AND. &
              parser % inFixFormula(i:i) /= ")" .AND. parser % inFixFormula(i:i) /= " "  )THEN

            errorMsg = "Missing operator or closing parentheses after token : "//&
                       TRIM( parser % inFix % tokens( parser % inFix % top_index ) % tokenString )
            RETURN

          ENDIF

        ELSEIF( IsNumber( parser % inFixFormula(i:i) ) )THEN

          parser % inFix % top_index = parser % inFix % top_index + 1
          parser % inFix % tokens( parser % inFix % top_index ) % tokenString = ''


          IF( parser % inFixFormula(i:i) == 'p' .OR. parser % inFixFormula(i:i) == 'P' )THEN

            ! Conditional for using built in "pi" definition
            parser % inFix % tokens( parser % inFix % top_index ) % tokenString(1:2) = parser % inFixFormula(i:i+1)
            j = 2

          ELSE

            j = 0
            DO WHILE( IsNumber( parser % inFixFormula(i+j:i+j) ) )

              parser % inFix % tokens( parser % inFix % top_index ) % tokenString(j+1:j+1) = parser % inFixFormula(i+j:i+j) 
              j = j+1

            ENDDO

          ENDIF

          parser % inFix % tokens( parser % inFix % top_index ) % tokenType = Number_Token
        
          i = i + j

          ! Next item must be an operator or a closing parentheses
          IF( .NOT. IsOperator( parser % infixFormula(i:i) ) .AND. &
              parser % inFixFormula(i:i) /= ")" .AND. parser % inFixFormula(i:i) /= " " )THEN

            errorMsg = "Missing operator or closing parentheses after token : "//&
                       TRIM( parser % inFix % tokens( parser % inFix % top_index ) % tokenString )
            RETURN

          ENDIF

        ELSEIF( IsSeparator( parser % inFixFormula(i:i) ) )THEN


          parser % inFix % top_index = parser % inFix % top_index + 1 
          parser % inFix % tokens( parser % inFix % top_index ) % tokenString = parser % inFixFormula(i:i)

          IF( parser % inFixFormula(i:i) == "(" )THEN
            parser % inFix % tokens( parser % inFix % top_index ) % tokenType   = OpeningParentheses_Token 
          ELSEIF( parser % inFixFormula(i:i) == ")" )THEN
            parser % inFix % tokens( parser % inFix % top_index ) % tokenType   = ClosingParentheses_Token 
          ELSE
            parser % inFix % tokens( parser % inFix % top_index ) % tokenType   = Operator_Token 
          ENDIF

          i = i + 1


        ELSEIF( IsFunction( parser % inFixFormula(i:i) ) )THEN

          parser % inFix % top_index = parser % inFix % top_index + 1
          parser % inFix % tokens( parser % inFix % top_index ) % tokenString = ''

          j = FindLastFunctionIndex( parser % inFixFormula(i:i+Max_Function_Length-1) )

          parser % inFix % tokens( parser % inFix % top_index ) % tokenString = parser % inFixFormula(i:i+j)
          parser % inFix % tokens( parser % inFix % top_index ) % tokenType   = Function_Token 
          i = i+j+1

          ! Check to see if the next string
          IF( parser % inFixFormula(i:i) /= "(" )THEN
            errorMsg = "Missing opening parentheses after token : "//&
                       TRIM( parser % inFix % tokens( parser % inFix % top_index ) % tokenString )

            RETURN
          ENDIF


        ELSE

          errorMsg = "Invalid Token : "//&
                     TRIM( parser % inFixFormula(i:i) )

          RETURN

        ENDIF

      ENDDO

      
      IF( parser % inFix % tokens(1) % tokenType == Operator_Token )THEN
         IF( TRIM( parser % inFix % tokens(1) % tokenString ) == "+" .OR. TRIM( parser % inFix % tokens(1) % tokenString ) == "-" )     THEN
            parser % inFix % tokens(1) % tokenType = Monadic_Token
         END IF
      END IF
      
      DO i = 2, parser % inFix % top_index
         IF( parser % inFix % tokens(i) % tokenType == Operator_Token .AND. parser % inFix % tokens(i-1) % tokenType == OpeningParentheses_Token )     THEN
            parser % inFix % tokens(i) % tokenType = Monadic_Token
         END IF
      END DO

      
      tokenized = .TRUE.

  END SUBROUTINE Tokenize


  SUBROUTINE ConvertToPostFix( parser )
    CLASS( EquationParser ), INTENT(inout) :: parser
    ! Local
    CHARACTER(Error_Message_Length) :: errorMsg
    TYPE( TokenStack )              :: operator_stack
    TYPE( Token )                   :: tok
    INTEGER                         :: i
    
      !success = .FALSE. 

      CALL parser % postfix % Construct( Stack_Length )
      CALL operator_stack % Construct( Stack_Length )
  
      DO i = 1, parser % infix % top_index
     
        IF( parser % inFix % tokens(i) % tokenType == Variable_Token .OR. &
            parser % inFix % tokens(i) % tokenType == Number_Token )THEN

          
          CALL parser % postFix % push( parser % inFix % tokens(i) )

  
        ELSEIF( parser % inFix % tokens(i) % tokenType == Function_Token )THEN

          CALL operator_stack % push( parser % inFix % tokens(i) )

        ELSEIF( parser % inFix % tokens(i) % tokenType == Operator_Token .OR. parser % inFix % tokens(i) % tokenType == Monadic_Token )THEN


          IF( .NOT. operator_stack % IsEmpty( ) )THEN

            tok = operator_stack % TopToken( )
              
            DO WHILE( TRIM(tok % tokenString) /= "(" .AND. &
                      Priority( TRIM(tok % tokenString) ) >  Priority( TRIM(parser % inFix % tokens(i) % tokenString) ) )
       
              CALL parser % postFix % push( tok )
              CALL operator_stack % pop( tok )
              tok = operator_stack % TopToken( )

            ENDDO

          ENDIF

          CALL operator_stack % push( parser % inFix % tokens(i) )

        ELSEIF( parser % inFix % tokens(i) % tokenType == OpeningParentheses_Token )THEN

          CALL operator_stack % push( parser % inFix % tokens(i) )


        ELSEIF( parser % inFix % tokens(i) % tokenType == ClosingParentheses_Token )THEN

          tok = operator_stack % TopToken( )

          DO WHILE( .NOT.( operator_stack % IsEmpty( ) ) .AND. TRIM(tok % tokenString) /= "(" )
            
            CALL parser % postFix % push( tok )
            CALL operator_stack % pop( tok )
            tok = operator_stack % TopToken( )

          ENDDO

          ! Pop the opening parenthesis
          CALL operator_stack % pop( tok )

        ENDIF

      ENDDO

      ! Pop the remaining operators
      DO WHILE( .NOT.( operator_stack % IsEmpty( ) ) )
        
        tok = operator_stack % TopToken( )
        CALL parser % postFix % push( tok )
        CALL operator_stack % pop( tok )
   
      ENDDO
      
  END SUBROUTINE ConvertToPostFix

  FUNCTION Evaluate( parser, x ) RESULT( f )
    CLASS(EquationParser) :: parser
    REAL(prec)            :: x(1:3)
    REAL(prec)            :: f
    ! Local
    INTEGER           :: k
    TYPE(Token)       :: t
    TYPE(NumberStack) :: stack
    REAL(prec)        :: v, a, b, c
         
      CALL stack % Construct( Stack_Length )

      IF( .NOT.( ALLOCATED( parser % postfix % tokens ) ) )THEN

        f = 0.0_prec

      ELSE

        DO k = 1, parser % postfix % top_index 
  
          t = parser % postfix % tokens(k) % Equals_Token( )
  
          SELECT CASE ( t % tokenType )
           
            CASE( Number_Token )
  
              IF( t % tokenString == "pi" .OR. t % tokenString == "PI" )     THEN
                 v = pi
              ELSE
                READ( t % tokenString, * ) v
              END IF
  
              CALL stack % Push( v )
                 
            CASE ( Variable_Token )
  
              IF( TRIM( t % tokenString ) == "x" )THEN
  
                 CALL stack % Push( x(1) )
  
              ELSEIF( TRIM( t % tokenString ) == "y" )THEN
  
                 CALL stack % Push( x(2) )
  
              ELSEIF( TRIM( t % tokenString ) == "z" )THEN
  
                 CALL stack % Push( x(3) )
  
              ENDIF
  
            CASE ( Operator_Token )
  
              CALL stack % Pop( a )
              CALL stack % Pop( b )
  
              SELECT CASE ( TRIM(t % tokenString) )
  
                 CASE ( "+" )
  
                    c = a + b
  
                 CASE ( "-" )
  
                    c = b - a
  
                 CASE ( "*" )
  
                    c = a*b
  
                 CASE ( "/" )
  
                    c = b/a
  
                 CASE ( "^" )
  
                    c = b**a
                 CASE DEFAULT
  
              END SELECT
  
              CALL stack % Push( c )
              
           CASE ( Function_Token )
  
              CALL stack % Pop( a )
  
              b = F_of_X( TRIM(t % tokenString), a )
  
              CALL stack % Push( b )
              
           CASE ( Monadic_Token )
  
             IF( TRIM(t % tokenString) == "-" )     THEN
  
                CALL stack % Pop( a )
                a = -a
                CALL stack % Push( a )
  
             END IF
             
           CASE DEFAULT
  
         END SELECT
  
       END DO
  
       CALL stack % Pop( a )
       f = a
  
       CALL stack % Destruct( )

     ENDIF
         
  END FUNCTION Evaluate

  SUBROUTINE Print_InfixTokens( parser )
    CLASS( EquationParser ), INTENT(in) :: parser
    ! Local
    INTEGER :: i

      DO i = 1, parser % inFix % top_index
        PRINT*, TRIM( parser % inFix % tokens(i) % tokenString )
      ENDDO


  END SUBROUTINE Print_InfixTokens

  SUBROUTINE Print_PostfixTokens( parser )
    CLASS( EquationParser ), INTENT(in) :: parser
    ! Local
    INTEGER :: i

      DO i = 1, parser % postFix % top_index
        PRINT*, TRIM( parser % postFix % tokens(i) % tokenString )
      ENDDO


  END SUBROUTINE Print_PostfixTokens

  ! TokenStack and NumberStack 

  SUBROUTINE Construct_TokenStack( stack, N )
   CLASS(TokenStack), INTENT(out) :: stack
   INTEGER, INTENT(in)            :: N

     ALLOCATE( stack % tokens(1:N) )
     stack % top_index = 0

  END SUBROUTINE Construct_TokenStack

  SUBROUTINE Destruct_TokenStack( stack )
    CLASS(TokenStack), INTENT(inout) :: stack

      IF( ALLOCATED( stack % tokens ) ) DEALLOCATE( stack % tokens )
      stack % top_index = 0

  END SUBROUTINE Destruct_TokenStack

  SUBROUTINE Push_TokenStack( stack, tok ) 
    CLASS(TokenStack), INTENT(inout) :: stack
    TYPE(Token), INTENT(in)         :: tok

      stack % top_index                  = stack % top_index + 1
      stack % tokens(stack % top_index)  % tokenString = tok % tokenString
      stack % tokens(stack % top_index)  % tokenType   = tok % tokenType
 
  END SUBROUTINE Push_TokenStack

  SUBROUTINE Pop_TokenStack( stack, tok ) 
    CLASS(TokenStack), INTENT(inout) :: stack
    TYPE(Token), INTENT(out)        :: tok
    
      IF( stack % top_index <= 0 ) THEN
        PRINT *, "Attempt to pop from empty token stack"
      ELSE 
        tok % tokenString         = stack % tokens( stack % top_index ) % tokenString
        tok % tokenType           = stack % tokens( stack % top_index ) % tokenType
        stack % top_index = stack % top_index - 1
      END IF


  END SUBROUTINE Pop_TokenStack

  SUBROUTINE Peek_TokenStack( stack, tok ) 
    CLASS(TokenStack), INTENT(in) :: stack
    TYPE(Token), INTENT(out)     :: tok
    
      IF( stack % top_index <= 0 ) THEN
        PRINT *, "Attempt to peek from empty token stack"
      ELSE 
        tok % tokenString = stack % tokens( stack % top_index ) % tokenString
        tok % tokenType   = stack % tokens( stack % top_index ) % tokenType
      END IF
  END SUBROUTINE Peek_TokenStack

  LOGICAL FUNCTION IsEmpty_TokenStack( stack )
    CLASS( TokenStack ) :: stack

      IsEmpty_TokenStack = .FALSE.

      IF( stack % top_index <= 0 )THEN
        IsEmpty_TokenStack = .TRUE.
      ENDIF

  END FUNCTION IsEmpty_TokenStack

  TYPE( Token ) FUNCTION TopToken( stack )
    CLASS( TokenStack ) :: stack

      IF( stack % top_index > 0 )THEN
        TopToken % tokenString = stack % tokens( stack % top_index ) % tokenString
        TopToken % tokenType   = stack % tokens( stack % top_index ) % tokenType
      ELSE
        TopToken % tokenString = ''
      ENDIF

  END FUNCTION TopToken 

  FUNCTION Equals_Token( tok1 ) RESULT( tok2 )
    CLASS(Token) :: tok1
    TYPE(Token)  :: tok2

      tok2 % tokenString = tok1 % tokenString 
      tok2 % tokenType   = tok1 % tokenType

  END FUNCTION Equals_Token

  SUBROUTINE Construct_NumberStack( stack, N )
   CLASS(NumberStack), INTENT(out) :: stack
   INTEGER, INTENT(in)            :: N

     ALLOCATE( stack % tokens(1:N) )
     stack % top_index = 0

  END SUBROUTINE Construct_NumberStack

  SUBROUTINE Destruct_NumberStack( stack )
    CLASS(NumberStack), INTENT(inout) :: stack

      IF( ALLOCATED( stack % tokens) ) DEALLOCATE( stack % tokens )
      stack % top_index = 0

  END SUBROUTINE Destruct_NumberStack

  SUBROUTINE Push_NumberStack( stack, tok ) 
    CLASS(NumberStack), INTENT(inout) :: stack
    REAL(prec), INTENT(in)         :: tok

      stack % top_index                  = stack % top_index + 1
      stack % tokens(stack % top_index) = tok 
 
  END SUBROUTINE Push_NumberStack

  SUBROUTINE Pop_NumberStack( stack, tok ) 
    CLASS(NumberStack), INTENT(inout) :: stack
    REAL(prec), INTENT(out)        :: tok
    
      IF( stack % top_index <= 0 ) THEN
        PRINT *, "Attempt to pop from empty token stack"
      ELSE 
        tok               = stack % tokens( stack % top_index )
        stack % top_index = stack % top_index - 1
      END IF


  END SUBROUTINE Pop_NumberStack

  SUBROUTINE Peek_NumberStack( stack, tok ) 
    CLASS(NumberStack), INTENT(in) :: stack
    REAL(prec), INTENT(out)        :: tok
    
      IF( stack % top_index <= 0 ) THEN
        PRINT *, "Attempt to peek from empty token stack"
      ELSE 
        tok = stack % tokens( stack % top_index )
      END IF
  END SUBROUTINE Peek_NumberStack

  LOGICAL FUNCTION IsEmpty_NumberStack( stack )
    CLASS( NumberStack ) :: stack

      IsEmpty_NumberStack = .FALSE.

      IF( stack % top_index <= 0 )THEN
        IsEmpty_NumberStack = .TRUE.
      ENDIF

  END FUNCTION IsEmpty_NumberStack

  ! Support Functions !

  LOGICAL FUNCTION IsSeparator( eqChar )
    CHARACTER(1) :: eqChar
    ! Local
    INTEGER :: i

      IsSeparator = .FALSE.
      DO i = 1, nSeparators 

        IF( eqChar == separators(i) )THEN
          IsSeparator = .TRUE.
        ENDIF

      ENDDO

  END FUNCTION IsSeparator

  LOGICAL FUNCTION IsNumber( eqChar )
    CHARACTER(1) :: eqChar
    ! Local
    INTEGER :: i

      IsNumber = .FALSE.

      IF( eqChar == '.' .OR. eqChar == 'p' .OR. eqChar == 'P' )THEN
        IsNumber = .TRUE.
        RETURN
      ENDIF
         
      DO i = 1, 10

        IF( eqChar == numbers(i) )THEN
          IsNumber = .TRUE.
        ENDIF

      ENDDO

  END FUNCTION IsNumber

  LOGICAL FUNCTION IsVariable( eqChar )
    CHARACTER(1) :: eqChar
    ! Local
    INTEGER :: i

      IsVariable = .FALSE.
      DO i = 1, 3

        IF( eqChar == variables(i) )THEN
          IsVariable = .TRUE.
        ENDIF

      ENDDO

  END FUNCTION IsVariable

  LOGICAL FUNCTION IsOperator( eqChar )
    CHARACTER(1) :: eqChar
    ! Local
    INTEGER :: i

      IsOperator = .FALSE.
      DO i = 1, 5

        IF( eqChar == operators(i) )THEN
          IsOperator = .TRUE.
        ENDIF

      ENDDO

  END FUNCTION IsOperator

  LOGICAL FUNCTION IsFunction( eqChar )
    CHARACTER(1) :: eqChar
    ! Local
    INTEGER :: i

      IsFunction = .FALSE.
      DO i = 1, nFunctions

        IF( eqChar == functions(i) % str(1:1) )THEN
          IsFunction = .TRUE.
        ENDIF

      ENDDO

  END FUNCTION IsFunction

  FUNCTION FindLastFunctionIndex( eqChar ) RESULT( j )
    CHARACTER(Max_Function_Length) :: eqChar
    INTEGER                        :: i, j

      DO i = 1, Max_Function_Length

        IF( eqChar(i:i) == "(" )THEN
          j = i-2
          EXIT
        ENDIF

      ENDDO
         
  END FUNCTION FindLastFunctionIndex

  REAL(prec) FUNCTION F_of_X( func, x ) 
    CHARACTER(*) :: func
    REAL(prec)   :: x
    ! Local
    REAL(prec)   :: r

      IF( TRIM( func ) == "cos" .OR. TRIM( func ) == "COS" )THEN

        F_of_X = cos( x )

      ELSEIF( TRIM( func ) == "sin" .OR. TRIM( func ) == "SIN" )THEN

        F_of_X = sin( x )

      ELSEIF( TRIM( func ) == "tan" .OR. TRIM( func ) == "TAN" )THEN

        F_of_X = tan( x )

      ELSEIF( TRIM( func ) == "tanh" .OR. TRIM( func ) == "TANH" )THEN

        F_of_X = tanh( x )

      ELSEIF( TRIM( func ) == "sech" .OR. TRIM( func ) == "SECH" )THEN

        F_of_X = 2.0_prec/( exp(x) + exp(-x) )

      ELSEIF( TRIM( func ) == "sqrt" .OR. TRIM( func ) == "SQRT" )THEN

        F_of_X = sqrt( x )

      ELSEIF( TRIM( func ) == "abs" .OR. TRIM( func ) == "ABS" )THEN

        F_of_X = abs( x )

      ELSEIF( TRIM( func ) == "exp" .OR. TRIM( func ) == "EXP" )THEN

        F_of_X = exp( x )

      ELSEIF( TRIM( func ) == "ln" .OR. TRIM( func ) == "LN" )THEN

        F_of_X = log( x )

      ELSEIF( TRIM( func ) == "log" .OR. TRIM( func ) == "LOG" )THEN

        F_of_X = log10( x )

      ELSEIF( TRIM( func ) == "acos" .OR. TRIM( func ) == "ACOS" )THEN

        F_of_X = acos( x )

      ELSEIF( TRIM( func ) == "asin" .OR. TRIM( func ) == "ASIN" )THEN

        F_of_X = asin( x )

      ELSEIF( TRIM( func ) == "atan" .OR. TRIM( func ) == "ATAN" )THEN

        F_of_X = atan( x )

      ELSEIF( TRIM( func ) == "rand" .OR. TRIM( func ) == "RAND" )THEN

        CALL RANDOM_NUMBER( r )
        F_of_X = r*x

      ELSE
 
        F_of_X = 0.0_prec

      ENDIF


  END FUNCTION F_of_X

  INTEGER FUNCTION Priority( operatorString )
    CHARACTER(1) :: operatorString


      IF( IsFunction( operatorString ) )THEN

        Priority = 4

      ELSEIF( operatorString == '^' )THEN

        Priority = 3

      ELSEIF( operatorString == '*' .OR. operatorString == '/' )THEN

        Priority = 2

      ELSEIF( operatorString == '+' .OR. operatorString == '-' )THEN
   
        Priority = 1
 
      ELSE

        Priority = 0
      
      ENDIF

  END FUNCTION Priority

END MODULE EquationParser_Class
