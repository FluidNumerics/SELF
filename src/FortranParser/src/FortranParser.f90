!
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
MODULE FortranParser
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  ! Fortran 2008 function parser
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  ! 
  ! This is an OOP Fortran 2008 version of the original fparser by Roland Schmehl. This simple class
  ! wrapping of the original fparser has been developed by Jacopo Chevallard, and it is available on
  ! the GitHub repository https://github.com/jacopo-chevallard/FortranParser.
  ! 
  ! For comments and bug reports, please open an issue on
  ! https://github.com/jacopo-chevallard/FortranParser/issues
  !
  ! This function parser module is intended for applications where a set of mathematical
  ! fortran-style expressions is specified at runtime and is then evaluated for a large 
  ! number of variable values. This is done by compiling the set of function strings 
  ! into byte code, which is interpreted efficiently for the various variable values. 
  !
  ! The source code of the original fparser is available from http://fparser.sourceforge.net
  !
  ! Please send comments, corrections or questions realtive to the original fparser to its author:
  ! Roland Schmehl <roland.schmehl@alumni.uni-karlsruhe.de>
  !
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  ! The function parser concept is based on a C++ class library written by  Juha 
  ! Nieminen <warp@iki.fi> available from http://warp.povusers.org/FunctionParser/
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  USE FortranParser_parameters, ONLY: rn,is               ! Import KIND parameters

  IMPLICIT NONE

  PUBLIC                     :: EquationParser

  !------- -------- --------- --------- --------- --------- --------- --------- -------
  PRIVATE

  INTEGER(is),                              PARAMETER :: cImmed   = 1,          &
                                                         cNeg     = 2,          &
                                                         cAdd     = 3,          & 
                                                         cSub     = 4,          & 
                                                         cMul     = 5,          & 
                                                         cDiv     = 6,          & 
                                                         cPow     = 7,          & 
                                                         cAbs     = 8,          &
                                                         cExp     = 9,          &
                                                         cLog10   = 10,         &
                                                         cLog     = 11,         &
                                                         cSqrt    = 12,         &
                                                         cSinh    = 13,         &
                                                         cCosh    = 14,         &
                                                         cTanh    = 15,         &
                                                         cSin     = 16,         &
                                                         cCos     = 17,         &
                                                         cTan     = 18,         &
                                                         cAsin    = 19,         &
                                                         cAcos    = 20,         &
                                                         cAtan    = 21,         &
                                                         VarBegin = 22

  CHARACTER (LEN=1), DIMENSION(cAdd:cPow),  PARAMETER :: Ops      = (/ '+',     &
                                                                       '-',     &
                                                                       '*',     &
                                                                       '/',     &
                                                                       '^' /)

  CHARACTER (LEN=5), DIMENSION(cAbs:cAtan), PARAMETER :: Funcs    = (/ 'abs  ', &
                                                                       'exp  ', &
                                                                       'log10', &
                                                                       'log  ', &
                                                                       'sqrt ', &
                                                                       'sinh ', &
                                                                       'cosh ', &
                                                                       'tanh ', &
                                                                       'sin  ', &
                                                                       'cos  ', &
                                                                       'tan  ', &
                                                                       'asin ', &
                                                                       'acos ', &
                                                                       'atan ' /)

  INTEGER, parameter  :: MAX_FUN_LENGTH = 1024

  TYPE EquationParser

    INTEGER(is), POINTER :: ByteCode(:) => null()
    INTEGER              :: ByteCodeSize = 0
    REAL(rn),    POINTER :: Immed(:) => null()
    INTEGER              :: ImmedSize = 0
    REAL(rn),    POINTER :: Stack(:) => null()
    INTEGER              :: StackSize = 0
    INTEGER              :: StackPtr = 0

    character(len=MAX_FUN_LENGTH) :: funcString = ''
    character(len=MAX_FUN_LENGTH) :: funcStringOrig = ''
    character(len=MAX_FUN_LENGTH), allocatable :: variableNames(:) 
    contains

      private

      procedure, public :: evaluate
      procedure:: parse
      procedure :: Compile
      procedure :: AddCompiledByte
      procedure :: CompileSubstr
      procedure :: MathItemIndex
      procedure :: CheckSyntax

      final :: finalize

  END TYPE EquationParser

  ! Class constructor
  interface EquationParser
    procedure constructor
  end interface EquationParser

CONTAINS

!*****************************************************************************************
  type (EquationParser) function constructor(FuncStr, Var)

    CHARACTER (LEN=*),               INTENT(in) :: FuncStr   ! Function string
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names

    constructor%ByteCode => null()
    constructor%Immed => null()
    constructor%Stack => null()

    constructor%ByteCodeSize = 0
    constructor%ImmedSize = 0
    constructor%StackSize = 0
    constructor%StackPtr = 0

    constructor%funcString = FuncStr
    constructor%funcStringOrig = FuncStr

    allocate(constructor%variableNames(size(Var)))

    constructor%variableNames(:) = Var(:)

    call constructor%parse()

  end function constructor

!*****************************************************************************************
  subroutine finalize(this)

    type(EquationParser) :: this

    if (associated(this%ByteCode))  nullify(this%ByteCode)
    if (associated(this%Immed))     nullify(this%Immed)
    if (associated(this%Stack))     nullify(this%Stack)

  end subroutine finalize

!*****************************************************************************************
  SUBROUTINE parse(this)
    ! Parse ith function string FuncStr and compile it into bytecode
    class(EquationParser) :: this

    CALL Replace('**','^ ', this%funcString)                            ! Exponent into 1-Char. format

    CALL RemoveSpaces(this%funcString)                                 ! Condense function string

    CALL this%CheckSyntax()

    CALL this%Compile()                                ! Compile into bytecode

  END SUBROUTINE parse

!*****************************************************************************************
  FUNCTION evaluate(this, Val) RESULT (res)
    ! Evaluate bytecode of ith function for the values passed in array Val(:)
    class(EquationParser) :: this
    REAL(rn), DIMENSION(:), INTENT(in) :: Val                ! Variable values

    REAL(rn)                           :: res                ! Result
    INTEGER                            :: IP,              & ! Instruction pointer
                                          DP,              & ! Data pointer
                                          SP                 ! Stack pointer
    REAL(rn),                PARAMETER :: zero = 0._rn
    integer       :: EvalErrType

    DP = 1
    SP = 0
    EvalErrType=0

    DO IP=1,this%ByteCodeSize

       SELECT CASE (this%ByteCode(IP))

       CASE (cImmed); SP=SP+1; this%Stack(SP)=this%Immed(DP); DP=DP+1

       CASE   (cNeg); this%Stack(SP)=-this%Stack(SP)

       CASE   (cAdd); this%Stack(SP-1)=this%Stack(SP-1)+this%Stack(SP); SP=SP-1

       CASE   (cSub); this%Stack(SP-1)=this%Stack(SP-1)-this%Stack(SP); SP=SP-1

       CASE   (cMul); this%Stack(SP-1)=this%Stack(SP-1)*this%Stack(SP); SP=SP-1

       CASE   (cDiv)

         IF (this%Stack(SP)==0._rn) THEN 
           EvalErrType=1
           res=zero
           exit
         ENDIF
         this%Stack(SP-1)=this%Stack(SP-1)/this%Stack(SP); SP=SP-1

       CASE   (cPow); this%Stack(SP-1)=this%Stack(SP-1)**this%Stack(SP); SP=SP-1

       CASE   (cAbs); this%Stack(SP)=ABS(this%Stack(SP))

       CASE   (cExp); this%Stack(SP)=EXP(this%Stack(SP))

       CASE (cLog10)

         IF (this%Stack(SP)<=0._rn) THEN
           EvalErrType=3
           res=zero
           exit
         ENDIF
         this%Stack(SP)=LOG10(this%Stack(SP))

       CASE   (cLog)

         IF (this%Stack(SP)<=0._rn) THEN
           EvalErrType=3
           res=zero
           exit
         ENDIF
         this%Stack(SP)=LOG(this%Stack(SP))

       CASE  (cSqrt)

         IF (this%Stack(SP)<0._rn) THEN
           EvalErrType=3
           res=zero
           exit
         ENDIF
         this%Stack(SP)=SQRT(this%Stack(SP))

       CASE  (cSinh); this%Stack(SP)=SINH(this%Stack(SP))

       CASE  (cCosh); this%Stack(SP)=COSH(this%Stack(SP))

       CASE  (cTanh); this%Stack(SP)=TANH(this%Stack(SP))

       CASE   (cSin); this%Stack(SP)=SIN(this%Stack(SP))

       CASE   (cCos); this%Stack(SP)=COS(this%Stack(SP))

       CASE   (cTan); this%Stack(SP)=TAN(this%Stack(SP))

       CASE  (cAsin) 

         IF ((this%Stack(SP)<-1._rn) .OR. (this%Stack(SP)>1._rn)) THEN
           EvalErrType=4
           res=zero
           exit
         ENDIF
         this%Stack(SP)=ASIN(this%Stack(SP))

       CASE  (cAcos); 
         IF ((this%Stack(SP)<-1._rn).OR.(this%Stack(SP)>1._rn)) THEN
           EvalErrType=4
           res=zero
           exit
         ENDIF
         this%Stack(SP)=ACOS(this%Stack(SP))

       CASE  (cAtan); this%Stack(SP)=ATAN(this%Stack(SP))

       CASE  DEFAULT; SP=SP+1; this%Stack(SP)=Val(this%ByteCode(IP)-VarBegin+1)

       END SELECT

    END DO

    IF (EvalErrType > 0) then
      WRITE(*,*)'*** Error: ',EvalErrMsg(EvalErrType)
    else
      res = this%Stack(1)
    endif

  END FUNCTION evaluate

!*****************************************************************************************
  SUBROUTINE CheckSyntax(this)
    ! Check syntax of function string,  returns 0 if syntax is ok
    class(EquationParser) :: this
    INTEGER(is)                                 :: n
    CHARACTER (LEN=1)                           :: c
    REAL(rn)                                    :: r
    LOGICAL                                     :: err
    INTEGER                                     :: ParCnt, & ! Parenthesis counter
                                                   j,ib,in,lFunc

    j = 1
    ParCnt = 0
    lFunc = LEN_TRIM(this%funcString)
    step: DO
       IF (j > lFunc) CALL ParseErrMsg (j, this%funcStringOrig)
       c = this%funcString(j:j)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Check for valid operand (must appear)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       IF (c == '-' .OR. c == '+') THEN                      ! Check for leading - or +
          j = j+1
          IF (j > lFunc) CALL ParseErrMsg (j, this%funcStringOrig, 'Missing operand')
          c = this%funcString(j:j)
          IF (ANY(c == Ops)) CALL ParseErrMsg (j, this%funcStringOrig, 'Multiple operators')
       END IF
       n = MathFunctionIndex (this%funcString(j:))
       IF (n > 0) THEN                                       ! Check for math function
          j = j+LEN_TRIM(Funcs(n))
          IF (j > lFunc) CALL ParseErrMsg (j, this%funcStringOrig, 'Missing function argument')
          c = this%funcString(j:j)
          IF (c /= '(') CALL ParseErrMsg (j, this%funcStringOrig, 'Missing opening parenthesis')
       END IF
       IF (c == '(') THEN                                    ! Check for opening parenthesis
          ParCnt = ParCnt+1
          j = j+1
          CYCLE step
       END IF
       IF (SCAN(c,'0123456789.') > 0) THEN                   ! Check for number
          r = RealNum (this%funcString(j:),ib,in,err)
          IF (err) CALL ParseErrMsg (j, this%funcStringOrig, 'Invalid number format:  '//this%funcString(j+ib-1:j+in-2))
          j = j+in-1
          IF (j > lFunc) EXIT
          c = this%funcString(j:j)
       ELSE                                                  ! Check for variable
          n = VariableIndex (this%funcString(j:),this%variableNames,ib,in)
          IF (n == 0) CALL ParseErrMsg (j, this%funcStringOrig, 'Invalid element: '//this%funcString(j+ib-1:j+in-2))
          j = j+in-1
          IF (j > lFunc) EXIT
          c = this%funcString(j:j)
       END IF
       DO WHILE (c == ')')                                   ! Check for closing parenthesis
          ParCnt = ParCnt-1
          IF (ParCnt < 0) CALL ParseErrMsg (j, this%funcStringOrig, 'Mismatched parenthesis')
          IF (this%funcString(j-1:j-1) == '(') CALL ParseErrMsg (j-1, this%funcStringOrig, 'Empty parentheses')
          j = j+1
          IF (j > lFunc) EXIT
          c = this%funcString(j:j)
       END DO
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Now, we have a legal operand: A legal operator or end of string must follow
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       IF (j > lFunc) EXIT
       IF (ANY(c == Ops)) THEN                               ! Check for multiple operators
          IF (j+1 > lFunc) CALL ParseErrMsg (j, this%funcStringOrig)
          IF (ANY(this%funcString(j+1:j+1) == Ops)) CALL ParseErrMsg (j+1, this%funcStringOrig, 'Multiple operators')
       ELSE                                                  ! Check for next operand
          CALL ParseErrMsg (j, this%funcStringOrig, 'Missing operator')
       END IF
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Now, we have an operand and an operator: the next loop will check for another 
       ! operand (must appear)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       j = j+1
    END DO step
    IF (ParCnt > 0) CALL ParseErrMsg (j, this%funcStringOrig, 'Missing )')
  END SUBROUTINE CheckSyntax

!*****************************************************************************************
  FUNCTION EvalErrMsg(EvalErrType) RESULT (msg)
    ! Return error message
    integer, intent(in) :: EvalErrType
    CHARACTER (LEN=*), DIMENSION(4), PARAMETER :: m = (/ 'Division by zero                ', &
                                                         'Argument of SQRT negative       ', &
                                                         'Argument of LOG negative        ', &
                                                         'Argument of ASIN or ACOS illegal' /)
    CHARACTER (LEN=LEN(m))                     :: msg
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (EvalErrType < 1 .OR. EvalErrType > SIZE(m)) THEN
       msg = ''
    ELSE
       msg = m(EvalErrType)
    ENDIF

  END FUNCTION EvalErrMsg

!*****************************************************************************************
  SUBROUTINE ParseErrMsg (j, FuncStr, Msg)
    ! Print error message and terminate program
    INTEGER,                     INTENT(in) :: j
    CHARACTER (LEN=*),           INTENT(in) :: FuncStr       ! Original function string
    CHARACTER (LEN=*), OPTIONAL, INTENT(in) :: Msg

    INTEGER                                 :: k

    IF (PRESENT(Msg)) THEN
       WRITE(*,*) '*** Error in syntax of function string: '//Msg
    ELSE
       WRITE(*,*) '*** Error in syntax of function string:'
    ENDIF

    WRITE(*,*)
    WRITE(*,'(A)') ' '//FuncStr

    WRITE(*,'(A)') '?'
    STOP
  END SUBROUTINE ParseErrMsg

!*****************************************************************************************
  FUNCTION OperatorIndex (c) RESULT (n)
    ! Return operator index
    CHARACTER (LEN=1), INTENT(in) :: c
    INTEGER(is)                   :: n,j

    n = 0

    DO j=cAdd,cPow
       IF (c == Ops(j)) THEN
          n = j
          EXIT
       END IF
    END DO

  END FUNCTION OperatorIndex

!*****************************************************************************************
  FUNCTION MathFunctionIndex (str) RESULT (n)
    ! Return index of math function beginnig at 1st position of string str
    CHARACTER (LEN=*), INTENT(in) :: str

    INTEGER(is)                   :: n,j
    INTEGER                       :: k
    CHARACTER (LEN=LEN(Funcs))    :: fun

    n = 0

    DO j=cAbs,cAtan                                          ! Check all math functions
       k = MIN(LEN_TRIM(Funcs(j)), LEN(str))   
       CALL LowCase (str(1:k), fun)
       IF (fun == Funcs(j)) THEN                             ! Compare lower case letters
          n = j                                              ! Found a matching function
          EXIT
       END IF
    END DO

  END FUNCTION MathFunctionIndex

!*****************************************************************************************
  FUNCTION VariableIndex (str, Var, ibegin, inext) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return index of variable at begin of string str (returns 0 if no variable found)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),               INTENT(in) :: str       ! String
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n         ! Index of variable
    INTEGER, OPTIONAL,              INTENT(out) :: ibegin, & ! Start position of variable name
                                                   inext     ! Position of character after name
    INTEGER                                     :: j,ib,in,lstr
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    lstr = LEN_TRIM(str)
    IF (lstr > 0) THEN
       DO ib=1,lstr                                          ! Search for first character in str
          IF (str(ib:ib) /= ' ') EXIT                        ! When lstr>0 at least 1 char in str
       END DO                        
       DO in=ib,lstr                                         ! Search for name terminators
          IF (SCAN(str(in:in),'+-*/^) ') > 0) EXIT
       END DO
       DO j=1,SIZE(Var)
          IF (str(ib:in-1) == Var(j)) THEN                     
             n = j                                           ! Variable name found
             EXIT
          END IF
       END DO
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
  END FUNCTION VariableIndex

!*****************************************************************************************
  SUBROUTINE RemoveSpaces (str)
    ! Remove Spaces from string, remember positions of characters in old string
    CHARACTER (LEN=*), INTENT(inout) :: str

    INTEGER                          :: k,lstr

    lstr = LEN_TRIM(str)

    k = 1

    DO WHILE (str(k:lstr) /= ' ')                             
       IF (str(k:k) == ' ') THEN
          str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
          k = k-1
       END IF
       k = k+1
    END DO

  END SUBROUTINE RemoveSpaces

!*****************************************************************************************
  SUBROUTINE Replace (ca,cb,str)
    ! Replace ALL appearances of character set ca in string str by character set cb
    CHARACTER (LEN=*),       INTENT(in) :: ca
    CHARACTER (LEN=LEN(ca)), INTENT(in) :: cb                ! LEN(ca) must be LEN(cb)
    CHARACTER (LEN=*),    INTENT(inout) :: str

    INTEGER                             :: j,lca

    lca = LEN(ca)

    DO j=1,LEN_TRIM(str)-lca+1
       IF (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
    END DO

  END SUBROUTINE Replace

!*****************************************************************************************
  SUBROUTINE Compile(this)
    ! Compile i-th function string F into bytecode
    class(EquationParser) :: this
    INTEGER                                     :: istat

    IF (ASSOCIATED(this%ByteCode)) DEALLOCATE ( this%ByteCode, &
                                                   this%Immed,    &
                                                   this%Stack     )
    this%ByteCodeSize = 0
    this%ImmedSize    = 0
    this%StackSize    = 0
    this%StackPtr     = 0

    CALL this%CompileSubstr(1,LEN_TRIM(this%funcString))               ! Compile string to determine size

    ALLOCATE ( this%ByteCode(this%ByteCodeSize), & 
               this%Immed(this%ImmedSize),       &
               this%Stack(this%StackSize),       &
               STAT = istat                            )
    IF (istat /= 0) THEN
       WRITE(*,*) '*** Parser error: Memmory allocation for byte code failed'
       STOP
    ELSE
       this%ByteCodeSize = 0
       this%ImmedSize    = 0
       this%StackSize    = 0
       this%StackPtr     = 0
       CALL this%CompileSubstr(1,LEN_TRIM(this%funcString))            ! Compile string into bytecode
    END IF

  END SUBROUTINE Compile

!*****************************************************************************************
  SUBROUTINE AddCompiledByte(this, b)
    ! Add compiled byte to bytecode
    class(EquationParser) :: this
    INTEGER(is), INTENT(in) :: b                             ! Value of byte to be added

    this%ByteCodeSize = this%ByteCodeSize + 1

    IF (ASSOCIATED(this%ByteCode)) then
      this%ByteCode(this%ByteCodeSize) = b
    endif

  END SUBROUTINE AddCompiledByte

!*****************************************************************************************
  FUNCTION MathItemIndex(this, b, e) RESULT (n)
    ! Return math item index, if item is real number, enter it into Comp-structure
    class(EquationParser) :: this

    INTEGER,           INTENT(in) :: b,e                     ! First and last pos. of substring
    INTEGER(is)                                 :: n         ! Byte value of math item

    n = 0

    IF (SCAN(this%funcString(b:b),'0123456789.') > 0) THEN                 ! Check for begin of a number
       this%ImmedSize = this%ImmedSize + 1
       IF (ASSOCIATED(this%Immed)) this%Immed(this%ImmedSize) = RealNum(this%funcString(b:e))
       n = cImmed
    ELSE                                                     ! Check for a variable
       n = VariableIndex(this%funcString(b:e), this%variableNames)
       IF (n > 0) n = VarBegin+n-1
    END IF

  END FUNCTION MathItemIndex

!*****************************************************************************************
  FUNCTION CompletelyEnclosed (F, b, e) RESULT (res)
    ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
    CHARACTER (LEN=*), INTENT(in) :: F                       ! Function substring
    INTEGER,           INTENT(in) :: b,e                     ! First and last pos. of substring

    LOGICAL                       :: res
    INTEGER                       :: j,k

    res=.false.

    IF (F(b:b) == '(' .AND. F(e:e) == ')') THEN
       k = 0
       DO j=b+1,e-1
          IF     (F(j:j) == '(') THEN
             k = k+1
          ELSEIF (F(j:j) == ')') THEN
             k = k-1
          END IF
          IF (k < 0) EXIT
       END DO
       IF (k == 0) res=.true.                                ! All opened parenthesis closed
    END IF

  END FUNCTION CompletelyEnclosed

!*****************************************************************************************
  RECURSIVE SUBROUTINE CompileSubstr(this, b, e)
    ! Compile i-th function string funcString into bytecode
    class(EquationParser) :: this
    INTEGER,                         INTENT(in) :: b,e       ! Begin and end position substring

    INTEGER(is)                                 :: n
    INTEGER                                     :: b2,j,k,io
    CHARACTER (LEN=*),                PARAMETER :: calpha = 'abcdefghijklmnopqrstuvwxyz'// &
                                                            'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    ! Check for special cases of substring

    IF     (this%funcString(b:b) == '+') THEN                              ! Case 1: funcString(b:e) = '+...'
!      WRITE(*,*)'1. funcString(b:e) = "+..."'
       CALL this%CompileSubstr(b+1, e)
       RETURN
    ELSEIF (CompletelyEnclosed (this%funcString, b, e)) THEN               ! Case 2: funcString(b:e) = '(...)'
!      WRITE(*,*)'2. funcString(b:e) = "(...)"'
       CALL this%CompileSubstr(b+1, e-1)
       RETURN
    ELSEIF (SCAN(this%funcString(b:b), calpha) > 0) THEN        
       n = MathFunctionIndex(this%funcString(b:e))
       IF (n > 0) THEN
          b2 = b+INDEX(this%funcString(b:e),'(')-1
          IF (CompletelyEnclosed(this%funcString, b2, e)) THEN             ! Case 3: funcString(b:e) = 'fcn(...)'
!            WRITE(*,*)'3. funcString(b:e) = "fcn(...)"'
             CALL this%CompileSubstr(b2+1, e-1)
             CALL this%AddCompiledByte(n)
             RETURN
          END IF
       END IF

    ELSEIF (this%funcString(b:b) == '-') THEN
       IF (CompletelyEnclosed(this%funcString, b+1, e)) THEN              ! Case 4: this%funcString(b:e) = '-(...)'
!         WRITE(*,*)'4. this%funcString(b:e) = "-(...)"'
          CALL this%CompileSubstr(b+2, e-1)
          CALL this%AddCompiledByte(cNeg)
          RETURN
       ELSEIF (SCAN(this%funcString(b+1:b+1),calpha) > 0) THEN
          n = MathFunctionIndex(this%funcString(b+1:e))
          IF (n > 0) THEN
             b2 = b+INDEX(this%funcString(b+1:e),'(')
             IF (CompletelyEnclosed(this%funcString, b2, e)) THEN          ! Case 5: this%funcString(b:e) = '-fcn(...)'
!               WRITE(*,*)'5. this%funcString(b:e) = "-fcn(...)"'
                CALL this%CompileSubstr(b2+1, e-1);
                CALL this%AddCompiledByte(n)
                CALL this%AddCompiledByte(cNeg)
                RETURN
             END IF
          END IF
       ENDIF
    END IF

    ! Check for operator in substring: check only base level (k=0), exclude expr. in ()

    DO io=cAdd,cPow                                          ! Increasing priority +-*/^
       k = 0
       DO j=e,b,-1
          IF     (this%funcString(j:j) == ')') THEN
             k = k+1
          ELSEIF (this%funcString(j:j) == '(') THEN
             k = k-1
          END IF
          IF (k == 0 .AND. this%funcString(j:j) == Ops(io) .AND. IsBinaryOp (j, this%funcString)) THEN
             IF (ANY(this%funcString(j:j) == Ops(cMul:cPow)) .AND. this%funcString(b:b) == '-') THEN ! Case 6: this%funcString(b:e) = '-...Op...' with Op > -
!               WRITE(*,*)'6. this%funcString(b:e) = "-...Op..." with Op > -'
                CALL this%CompileSubstr(b+1, e)
                CALL this%AddCompiledByte(cNeg)
                RETURN                 
             ELSE                                                        ! Case 7: this%funcString(b:e) = '...BinOp...'
!               WRITE(*,*)'7. Binary operator',this%funcString(j:j)
                CALL this%CompileSubstr(b, j-1)
                CALL this%CompileSubstr(j+1, e)
                CALL this%AddCompiledByte(OperatorIndex(Ops(io)))
                this%StackPtr = this%StackPtr - 1
                RETURN
             END IF
          END IF
       END DO
    END DO

    ! Check for remaining items, i.e. variables or explicit numbers

    b2 = b

    IF (this%funcString(b:b) == '-') b2 = b2+1

    n = this%MathItemIndex(b2, e)

!   WRITE(*,*)'8. AddCompiledByte ',n
    CALL this%AddCompiledByte(n)

    this%StackPtr = this%StackPtr + 1
    IF (this%StackPtr > this%StackSize) this%StackSize = this%StackSize + 1

    IF (b2 > b) CALL this%AddCompiledByte(cNeg)

  END SUBROUTINE CompileSubstr

!*****************************************************************************************
  FUNCTION IsBinaryOp (j, F) RESULT (res)
    ! Check if operator F(j:j) in string F is binary operator
    ! Special cases already covered elsewhere:              (that is corrected in v1.1)
    ! - operator character F(j:j) is first character of string (j=1)
    INTEGER,           INTENT(in) :: j                       ! Position of Operator
    CHARACTER (LEN=*), INTENT(in) :: F                       ! String

    LOGICAL                       :: res                     ! Result
    INTEGER                       :: k
    LOGICAL                       :: Dflag,Pflag

    res=.true.

    IF (F(j:j) == '+' .OR. F(j:j) == '-') THEN               ! Plus or minus sign:
       IF (j == 1) THEN                                      ! - leading unary operator ?
          res = .false.
       ELSEIF (SCAN(F(j-1:j-1),'+-*/^(') > 0) THEN           ! - other unary operator ?
          res = .false.
       ELSEIF (SCAN(F(j+1:j+1),'0123456789') > 0 .AND. &     ! - in exponent of real number ?
               SCAN(F(j-1:j-1),'eEdD')       > 0) THEN
          Dflag=.false.; Pflag=.false.
          k = j-1
          DO WHILE (k > 1)                                   !   step to the left in mantissa 
             k = k-1
             IF     (SCAN(F(k:k),'0123456789') > 0) THEN
                Dflag=.true.
             ELSEIF (F(k:k) == '.') THEN
                IF (Pflag) THEN
                   EXIT                                      !   * EXIT: 2nd appearance of '.'
                ELSE
                   Pflag=.true.                              !   * mark 1st appearance of '.'
                ENDIF
             ELSE
                EXIT                                         !   * all other characters
             END IF
          END DO
          IF (Dflag .AND. (k == 1 .OR. SCAN(F(k:k),'+-*/^(') > 0)) res = .false.
       END IF
    END IF
  END FUNCTION IsBinaryOp

!*****************************************************************************************
  FUNCTION RealNum(str, ibegin, inext, error) RESULT (res)
    ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
    CHARACTER (LEN=*),  INTENT(in) :: str                    ! String
    REAL(rn)                       :: res                    ! Real number
    INTEGER, OPTIONAL, INTENT(out) :: ibegin,              & ! Start position of real number
                                      inext                  ! 1st character after real number
    LOGICAL, OPTIONAL, INTENT(out) :: error                  ! Error flag

    INTEGER                        :: ib,in,istat
    LOGICAL                        :: Bflag,               & ! .T. at begin of number in str
                                      InMan,               & ! .T. in mantissa of number
                                      Pflag,               & ! .T. after 1st '.' encountered
                                      Eflag,               & ! .T. at exponent identifier 'eEdD'
                                      InExp,               & ! .T. in exponent of number
                                      DInMan,              & ! .T. if at least 1 digit in mant.
                                      DInExp,              & ! .T. if at least 1 digit in exp.
                                      err                    ! Local error flag
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    Bflag=.true.; InMan=.false.; Pflag=.false.; Eflag=.false.; InExp=.false.
    DInMan=.false.; DInExp=.false.
    ib   = 1
    in   = 1
    DO WHILE (in <= LEN_TRIM(str))
       SELECT CASE (str(in:in))
       CASE (' ')                                            ! Only leading blanks permitted
          ib = ib+1
          IF (InMan .OR. Eflag .OR. InExp) EXIT
       CASE ('+','-')                                        ! Permitted only
          IF     (Bflag) THEN           
             InMan=.true.; Bflag=.false.                     ! - at beginning of mantissa
          ELSEIF (Eflag) THEN               
             InExp=.true.; Eflag=.false.                     ! - at beginning of exponent
          ELSE
             EXIT                                            ! - otherwise STOP
          ENDIF
       CASE ('0':'9')                                        ! Mark
          IF     (Bflag) THEN           
             InMan=.true.; Bflag=.false.                     ! - beginning of mantissa
          ELSEIF (Eflag) THEN               
             InExp=.true.; Eflag=.false.                     ! - beginning of exponent
          ENDIF
          IF (InMan) DInMan=.true.                           ! Mantissa contains digit
          IF (InExp) DInExp=.true.                           ! Exponent contains digit
       CASE ('.')
          IF     (Bflag) THEN
             Pflag=.true.                                    ! - mark 1st appearance of '.'
             InMan=.true.; Bflag=.false.                     !   mark beginning of mantissa
          ELSEIF (InMan .AND..NOT.Pflag) THEN
             Pflag=.true.                                    ! - mark 1st appearance of '.'
          ELSE
             EXIT                                            ! - otherwise STOP
          END IF
       CASE ('e','E','d','D')                                ! Permitted only
          IF (InMan) THEN
             Eflag=.true.; InMan=.false.                     ! - following mantissa
          ELSE
             EXIT                                            ! - otherwise STOP
          ENDIF
       CASE DEFAULT
          EXIT                                               ! STOP at all other characters
       END SELECT
       in = in+1
    END DO
    err = (ib > in-1) .OR. (.NOT.DInMan) .OR. ((Eflag.OR.InExp).AND..NOT.DInExp)
    IF (err) THEN
       res = 0.0_rn
    ELSE
       READ(str(ib:in-1),*,IOSTAT=istat) res
       err = istat /= 0
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
    IF (PRESENT(error))  error  = err
  END FUNCTION RealNum

!*****************************************************************************************
  SUBROUTINE LowCase (str1, str2)
    ! Transform upper case letters in str1 into lower case letters, result is str2
    IMPLICIT NONE
    CHARACTER (LEN=*),  INTENT(in) :: str1
    CHARACTER (LEN=*), INTENT(out) :: str2
    INTEGER                        :: j,k
    CHARACTER (LEN=*),   PARAMETER :: lc = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),   PARAMETER :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    str2 = str1

    DO j=1,LEN_TRIM(str1)
       k = INDEX(uc,str1(j:j))
       IF (k > 0) str2(j:j) = lc(k:k)
    END DO

  END SUBROUTINE LowCase

  END MODULE FortranParser
