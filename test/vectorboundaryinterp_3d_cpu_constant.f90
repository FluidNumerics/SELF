program test

  implicit none
  integer :: exit_code

  exit_code = vectorboundaryinterp_3d_cpu_constant()
  stop exit_code

contains
  integer function vectorboundaryinterp_3d_cpu_constant() result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Data

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
    type(Vector3D) :: f
    type(Lagrange),target :: interp

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Initialize vectors
    call f%Init(interp,nvar,nelem)

    ! Set the source vector (on the control grid) to a non-zero constant
    f%interior = 1.0_prec

    call f%BoundaryInterp()

    ! Calculate diff from exact
    f%boundary = abs(f%boundary-1.0_prec)

    if(maxval(f%boundary) <= tolerance) then
      r = 0
    else
      r = 1
      print*,maxval(f%boundary)
    endif

    call f%free()
    call interp%free()

  endfunction vectorboundaryinterp_3d_cpu_constant
endprogram test
