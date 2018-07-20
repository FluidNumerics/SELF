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
  INTEGER, PARAMETER :: Max_Equation_Length  = 512 
  INTEGER, PARAMETER :: Max_Function_Length  = 5
  INTEGER, PARAMETER :: Max_Variable_Length  = 12 
  INTEGER, PARAMETER :: Token_Length         = 48
  INTEGER, PARAMETER :: Stack_Length         = 64

  ! Token types 
  INTEGER, PARAMETER, PRIVATE :: None_Token               = 0
  INTEGER, PARAMETER, PRIVATE :: Number_Token             = 1
  INTEGER, PARAMETER, PRIVATE :: Variable_Token           = 2
  INTEGER, PARAMETER, PRIVATE :: Operator_Token           = 3
  INTEGER, PARAMETER, PRIVATE :: Function_Token           = 4
  INTEGER, PARAMETER, PRIVATE :: OpeningParentheses_Token = 5
  INTEGER, PARAMETER, PRIVATE :: ClosingParentheses_Token = 6
  INTEGER, PARAMETER, PRIVATE :: Monadic_Token            = 7

  INTEGER, PARAMETER, PRIVATE :: nSeparators = 7

  TYPE String
    CHARACTER(10) :: str
  END TYPE String

  CHARACTER(1), DIMENSION(7), PRIVATE  :: separators = (/ "+", "-", "*", "/", "(", ")", "^" /) 
  CHARACTER(1), DIMENSION(5), PRIVATE  :: operators  = (/ "+", "-", "*", "/", "^" /) 
  CHARACTER(1), DIMENSION(10), PRIVATE :: numbers    = (/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" /) 
  TYPE(String), DIMENSION(12), PRIVATE :: functions 
  CHARACTER(1), DIMENSION(4), PRIVATE  :: variables  = (/ "x", "y", "z", "t" /)

  ! Private Types !

  TYPE Token  
    CHARACTER(Token_Length) :: tokenString
    INTEGER                 :: tokenType
  END TYPE Token


  TYPE TokenStack
    TYPE(Token), ALLOCATABLE :: tokens(:)
    INTEGER                  :: top_index

    CONTAINS
  
      PROCEDURE :: Construct => Construct_TokenStack
      PROCEDURE :: Destruct  => Destruct_TokenStack
      
      PROCEDURE :: Push      => Push_TokenStack
      PROCEDURE :: Pop       => Pop_TokenStack
      PROCEDURE :: Peek      => Peek_TokenStack

      PROCEDURE :: IsEmpty   => IsEmpty_TokenStack
      PROCEDURE :: TopToken
  
  END TYPE TokenStack

  PRIVATE :: Token, TokenStack
  PRIVATE :: Construct_TokenStack, Destruct_TokenStack, Push_TokenStack, Pop_TokenStack, Peek_TokenStack, IsEmpty_TokenStack
  PRIVATE :: IsNumber, IsVariable, IsFunction, IsOperator, IsSeparator


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

      parser % equation = equation
      errorMsg = " "

      CALL parser % CleanEquation( equationIsClean, errorMsg )

      IF( equationIsClean )THEN

        CALL parser % Tokenize( tokenized, errorMsg )

!        CALL parser % Print_InfixTokens( )
!STOP 
        IF( tokenized )THEN

          CALL parser % ConvertToPostFix( )
        
        CALL parser % Print_PostfixTokens( )
STOP
!           IF( FinalSyntaxCheckOK( self ) )     THEN
!              success = .true.
!           ELSE
!              success = .false.
!           END IF

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

      PRINT*, TRIM(parser % infixformula )

  END SUBROUTINE CleanEquation

  SUBROUTINE Tokenize( parser, tokenized, errorMsg )
    CLASS( EquationParser ), INTENT(inout) :: parser
    LOGICAL, INTENT(out)                   :: tokenized
    CHARACTER(Error_Message_Length)        :: errorMsg
    ! Local
    INTEGER                         :: i, j 
 

      PRINT*, 'Tokenizing...'
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

          ! Next item must be an operator or a closing parentheses
          IF( .NOT. IsOperator( parser % infixFormula(i:i) ) .AND. &
              parser % inFixFormula(i:i) /= ")"  )THEN

            errorMsg = "Missing operator or closing parentheses after token : "//&
                       TRIM( parser % inFix % tokens( parser % inFix % top_index ) % tokenString )
            RETURN

          ENDIF

        ELSEIF( IsNumber( parser % inFixFormula(i:i) ) )THEN

          parser % inFix % top_index = parser % inFix % top_index + 1
          parser % inFix % tokens( parser % inFix % top_index ) % tokenString = ''

          j = 0
          DO WHILE( IsNumber( parser % inFixFormula(i+j:i+j) ) )

            parser % inFix % tokens( parser % inFix % top_index ) % tokenString(j+1:j+1) = parser % inFixFormula(i+j:i+j) 
            j = j+1

          ENDDO

          parser % inFix % tokens( parser % inFix % top_index ) % tokenType = Number_Token
        
          i = i + j

          ! Next item must be an operator or a closing parentheses
          IF( .NOT. IsOperator( parser % infixFormula(i:i) ) .AND. &
              parser % inFixFormula(i:i) /= ")"  )THEN

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

  SUBROUTINE ConvertToPostFix( parser )
    CLASS( EquationParser ), INTENT(inout) :: parser
    ! Local
    CHARACTER(Error_Message_Length) :: errorMsg
    TYPE( TokenStack )              :: operator_stack
    TYPE( Token )                   :: tok
    CHARACTER(Max_Equation_Length)  :: postFix
    INTEGER                         :: i, j
    
      !success = .FALSE. 

      CALL parser % postfix % Construct( Stack_Length )
      CALL operator_stack % Construct( Stack_Length )
  
      DO i = 1, parser % infix % top_index
     
         
        IF( parser % inFix % tokens(i) % tokenType == Variable_Token .OR. &
            parser % inFix % tokens(i) % tokenType == Number_Token )THEN

          
          CALL parser % postFix % push( parser % inFix % tokens(i) )

  
        ELSEIF( parser % inFix % tokens(i) % tokenType == Function_Token )THEN

          CALL operator_stack % push( parser % inFix % tokens(i) )

        ELSEIF( parser % inFix % tokens(i) % tokenType == Operator_Token )THEN

          tok = operator_stack % TopToken( )
              
          DO WHILE( .NOT.( operator_stack % IsEmpty( ) ) .AND. TRIM(tok % tokenString) /= "(" .AND. &
                    Priority( TRIM(tok % tokenString) ) >  Priority( TRIM(parser % inFix % tokens(i) % tokenString) ) )
     
            CALL parser % postFix % push( tok )
            CALL operator_stack % pop( tok )
            tok = operator_stack % TopToken( )

          ENDDO

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

  SUBROUTINE Construct_TokenStack( stack, N )
   CLASS(TokenStack), INTENT(out) :: stack
   INTEGER, INTENT(in)            :: N

     ALLOCATE( stack % tokens(1:N) )
     stack % top_index = 0

  END SUBROUTINE Construct_TokenStack

  SUBROUTINE Destruct_TokenStack( stack )
    CLASS(TokenStack), INTENT(inout) :: stack

      DEALLOCATE( stack % tokens )
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

      !IF( stack % top_index <= 0 ) THEN
      !  PRINT *, "  Empty Stack"
      !ELSE 
        TopToken = stack % tokens( stack % top_index )
      !ENDIF

  END FUNCTION TopToken 

  FUNCTION Equals_Token( tok1 ) RESULT( tok2 )
    CLASS(Token) :: tok1
    TYPE(Token)  :: tok2

      tok2 % tokenString = tok1 % tokenString 
      tok2 % tokenType   = tok1 % tokenType

  END FUNCTION Equals_Token

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
      DO i = 1, 10

        IF( eqChar == numbers(i) .OR. eqChar == '.' )THEN
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
      DO i = 1, 12

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

END MODULE EquationParser_Class

!
!
         
         
!
!!////////////////////////////////////////////////////////////////////////
!!
!      FUNCTION EvaluateEquation_At_( self, x ) RESULT(y)
!         IMPLICIT NONE
!!
!!        ---------
!!        Arguments
!!        ---------
!!
!         TYPE(EquationParser) :: self
!         REAL(KIND=EP)           :: x, y
!!
!!        ---------------
!!        Local variables
!!        ---------------
!!
!         INTEGER           :: k, N
!         TYPE(Token)       :: t
!         TYPE(NumberStack) :: stack
!         REAL(KIND=EP)     :: v, a, b, c
!         
!         N = SIZE(self%postfix)
!         CALL ConstructNumberStack( stack, N )
!         
!         DO k = 1, N 
!            t = self%postfix(k)
!            SELECT CASE ( t%tokenType )
!            
!               CASE( TYPE_NUMBER )
!                  IF( t%token == "pi" .OR. t%token == "PI" )     THEN
!                     v = PI
!                  ELSE
!                     READ( t%token, * ) v
!                  END IF
!                  CALL NumberStackPush( stack, v )
!                  
!               CASE ( TYPE_VARIABLE )
!                  CALL NumberStackPush( stack, x )
!
!               CASE ( TYPE_OPERATOR )
!                  CALL NumberStackPop( stack, a )
!                  CALL NumberStackPop( stack, b )
!                  SELECT CASE ( t%token )
!                     CASE ( "+" )
!                        c = a + b
!                     CASE ( "-" )
!                        c = b - a
!                     CASE ( "*" )
!                        c = a*b
!                     CASE ( "/" )
!                        c = b/a
!                     CASE ( "^" )
!                        IF( MOD(a,2.0_EP) == 0.0 )     THEN
!                           c = ABS(b)**a
!                        ELSE
!                           c = b**a
!                        END IF
!                     CASE DEFAULT
!                  END SELECT
!                  CALL NumberStackPush( stack, c )
!               
!               CASE ( TYPE_FUNCTION )
!                  CALL NumberStackPop( stack, a )
!                  CALL FunOfx( t%token, a, b )
!                  CALL NumberStackPush( stack, b )
!                  
!               CASE (TYPE_MONO_OPERATOR )
!                 IF( t%token == "-" )     THEN
!                    CALL NumberStackPop( stack, a )
!                    a = -a
!                    CALL NumberStackPush( stack, a )
!                 END IF
!                 
!               CASE DEFAULT
!            END SELECT
!         END DO
!!
!!        ----
!!        Done
!!        ----
!!
!         CALL NumberStackPop( stack, a )
!         y = a
!!
!!        --------
!!        Clean up
!!        --------
!!
!         CALL DestructNumberStack( stack )
!         
!      END FUNCTION EvaluateEquation_At_
!!
!!///////////////////////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
      
!!///////////////////////////////////////////////////////////////////////
!!                                                                       
!      SUBROUTINE FunOfx(fun,a,result) 
!!                                                                       
!      REAL(KIND=EP)    :: a,result 
!      CHARACTER(LEN=*) :: fun 
!      INTRINSIC        ::  abs,acos,asin,atan,cos,exp,log,log10,sin,sqrt,tan,tanh 
!!---
!!     ..                                                                
!      if ( fun == "cos" .OR. fun == "COS") then 
!          result = cos(a) 
!                                                                        
!      else if ( fun == "sin".OR. fun == "SIN") then 
!          result = sin(a) 
!                                                                        
!      else if ( fun == "exp".OR. fun == "EXP") then 
!          result = exp(a) 
!                                                                        
!      else if ( fun == "sqrt".OR. fun == "SQRT") then 
!          result = sqrt(a) 
!                                                                        
!      else if ( fun == "ln".OR. fun == "LN") then 
!          result = log(a) 
!                                                                        
!      else if ( fun == "log".OR. fun == "LOG") then 
!          result = log10(a) 
!                                                                        
!      else if ( fun == "abs".OR. fun == "ABS") then 
!          result = abs(a) 
!                                                                        
!      else if ( fun == "acos".OR. fun == "ACOS") then 
!          result = acos(a) 
!                                                                        
!      else if ( fun == "asin".OR. fun == "ASIN") then 
!          result = asin(a) 
!                                                                        
!      else if ( fun == "tan".OR. fun == "TAN") then 
!          result = tan(a) 
!                                                                        
!      else if ( fun == "atan".OR. fun == "ATAN") then 
!          result = atan(a) 
!                                                                        
!      else if ( fun == "tanh".OR. fun == "TANH") then 
!          result = tanh(a) 
!                                                                        
!      else 
!          write (6,fmt=*) "unknown function" 
!          result = 0.0d0 
!      end if 
!      END SUBROUTINE funofx                                          
!
!////////////////////////////////////////////////////////////////////////
!
