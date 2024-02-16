program test

  implicit none
  integer :: exit_code
  
  exit_code = scalargradient_3d_cpu_constant()
  stop exit_code

contains
integer function scalargradient_3d_cpu_constant() result(r)
  use SELF_Constants
  use SELF_Memory
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
  type(Scalar3D) :: f
  type(Vector3D) :: df
  type(Lagrange),target :: interp

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  ! Initialize scalars
  call f % Init(interp,nvar,nelem)

  call df % Init(interp,nvar,nelem)

  ! Set the source scalar (on the control grid) to a non-zero constant
  f % interior  = 1.0_prec

  ! Interpolate with gpuAccel = .FALSE.
  call f % Gradient(df, .false.)

  ! Calculate diff from exact
  df % interior  = abs(df % interior  - 0.0_prec)

  if (maxval(df % interior ) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  call f % free()
  call df % free()
  call interp % free()

end function scalargradient_3d_cpu_constantend program test
