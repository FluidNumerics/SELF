program test

  implicit none
  integer :: exit_code

  exit_code = vectordivergence_2d_cpu_constant()
  stop exit_code

contains
  integer function vectordivergence_2d_cpu_constant() result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Scalar_2D
    use SELF_Vector_2D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
    integer,parameter :: nelem = 100
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(Vector2D) :: f
    type(Scalar2D) :: df
    type(Lagrange),target :: interp

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Initialize vectors
    call f%Init(interp,nvar,nelem)

    call df%Init(interp,nvar,nelem)

    ! Set the source vector (on the control grid) to a non-zero constant
    f%interior = 1.0_prec

    call f%Divergence(df%interior)

    ! Calculate diff from exact
    df%interior = abs(df%interior-0.0_prec)

    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    call f%free()
    call df%free()
    call interp%free()

  endfunction vectordivergence_2d_cpu_constant
endprogram test
