# FortranParser

Fortran 2008 parser of mathematical expressions, based on Roland Schmehl [fparser](http://fparser.sourceforge.net)

## Table of contents

- [Changes](#changes)
- [Compilation](#compilation)
- [Basic usage](#basic-usage)
- [Usage in combination with JSON-Fortran](#usage-in-combination-with-json-fortran)
- [Error handling](#error-handling)
- [Function string syntax](#function-string-syntax)
- [Notes](#notes)
- [Credits](#credits)

## Changes
- What's new in version 2.0:

  * Renamed fparser to FortranParser
  * Changed approach to OOP, now everything happens by using the ``EquationParser`` class

- What's new in version 1.1:       (thanks to Wilton P. Silva and Juha Mäkipelto for the bug reports)

  * EXP failed: Corrected typo in alphabet string in subroutine LowCase.
  * Expression containing binary operator with precedence over unary 
  minus (e.g. "-x^2") failed and has been corrected in subroutines
  IsBinaryOp and CompileSubstr.
  * Leading plus (e.g. "+x") is now understood by correcting subroutines
  CompileSubstr and CheckSyntax
  * Multiple operators produce error message in subroutine CheckSyntax

## Compilation

To compile FortranParser you need [CMake](http://cmake.org). If CMake is installed on your machine, you can compile and install FortranParser with the following commands

- Clone the GitHub repository 
 ```bash
 git clone https://github.com/jacopo-chevallard/FortranParser.git
 ```
 this command will clone into a new directory ``FortranParser`` the current master branch.

- Move into the newly created ``FortranParser`` directory, create a ``build`` directory and move into it
 ```bash
 mkdir build ; cd build
 ```

- Run CMake
 ```bash
  cmake -DCMAKE_INSTALL_PREFIX=<install_dir> ..
 ```
 where ``<install_dir>`` is your installation directory. The FortranParser
 library will be installed in ``<install_dir>/lib``, the *mod files in
 ``<install_dir>/include``. The above CMake command include the compilation of
 some tests. This can be avoided by passing the option ``-DENABLE_TESTING=OFF``
 to the cmake command.

- Compile and install the FortranParser library
 ```bash
 make install
 ```

## Basic usage

### Step 0 - Module Import
In all program units where you want to use the function parser procedures 
and variables you must import the module by:

```fortran
USE FortranParser, only : EquationParser
```

This command imports only the public class ``EquationParser``, which has only
two public methods, the class ``constructor``, and the method ``evaluate``

### Step 1 - Constructor and function parsing

An instance of the ``EquationParser`` class is created with the following syntax
```fortran
  use FortranParser, only : EquationParser

  implicit none

  type(EquationParser) :: eqParser
  character(len=100)    :: stringEquation
  character(len=10)    :: variables(3)

  stringEquation = '10 + 3*x - 5*x*y + exp(-z**2)'
  variables = ['x', 'y', 'z']

  eqParser = EquationParser(stringEquation, variables)

```

The constructor deals with the parsing (checking and compilation) into the
bytecode. 

### Step 2 - Function evaluation
The function value is evaluated for a specific set of variable values 
by calling the method
```fortran
  value = eqParser%evaluate(varValues)
```
where ``varValues`` is 1-dimensional array containing the variable values.

## Usage in combination with JSON-Fortran

FortranParser in combination with
[JSON-Fortran](https://github.com/jacobwilliams/json-fortran) opens an
easy-to-use way of reading and evaluating mathematical expressions at runtime.
If you have both packages installed, than you can use:
```fortran
  use FortranParser, only : EquationParser
  use json_module

  implicit none

  type(json_file) :: json
  type(EquationParser) :: eqParser

  jsonString = '{"variables" : ["x", "y"], "function" : "2*x+sin(y**2)"}'
  call json%load_from_string(jsonString)

  call json%get('function', func, found)
  call json%get('variables', vars, found)

  eqParser = EquationParser(func, vars)

  value = eqParser%evaluate([2., 5.])
  value = eqParser%evaluate([-2., -5.])
```

## Error handling

An error in the function parsing step leads to a detailed error message 
(Type and position of error) and program termination.

An error during function evaluation returns a function value of 0.0 and
trigger an error message from the bytecode-interpreter.

## Function string syntax

Although they have to be passed as array elements of the same declared 
length (Fortran 90 restriction), the variable names can be of arbitrary 
actual length for the parser. Parsing for variables is case sensitive. 

The syntax of the function string is similar to the Fortran convention. 
Mathematical Operators recognized are +, -, *, /, ** or alternatively ^, 
whereas symbols for brackets must be (). 

The function parser recognizes the (single argument) Fortran 90 intrinsic 
functions
 - abs 
 - exp 
 - log10
 - log
 - sqrt
 - sinh
 - cosh
 - tanh
 - sin
 - cos
 - tan
 - asin
 - acos
 - atan

Parsing for intrinsic functions is case INsensitive.

Operations are evaluated in the correct order:

 - ()          expressions in brackets first
 - -A          unary minus (or plus)
 - A**B A^B    exponentiation (A raised to the power B)
 - A*B  A/B    multiplication and division
 - A+B  A-B    addition and subtraction

The function string can contain integer or real constants. To be recognized
as explicit constants these must conform to the format
```
[+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
```

where nnn means any number of digits. The mantissa must contain at least
one digit before or following an optional decimal point. Valid exponent 
identifiers are 'e', 'E', 'd' or 'D'. If they appear they must be followed 
by a valid exponent!

## Notes

* The precision of real numbers can be adapted to the calling program by 
  adjusting the KIND parameter rn in the external module parameters.

* The package compilation is based on CMake 

* The package contains some test programs to demonstrate implementation and
  performance of the function parser.

## Credits

The original fparser, by Roland Schmehl can be found at http://fparser.sourceforge.net.

The function parser concept is based on a C++ class library written by 
Juha Nieminen <warp@iki.fi> available from:

http://warp.povusers.org/FunctionParser/
