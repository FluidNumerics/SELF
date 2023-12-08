integer function scalargradient_2d_gpu_constant() result(r)
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
  type(Scalar2D) :: f
  type(Vector2D) :: df
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
  f % interior % hostdata = 1.0_prec

  call f % interior % updatedevice()

  ! Interpolate with gpuAccel = .true.
  call f % Gradient(df, .true.)

  call df % interior % updatehost()

  ! Calculate diff from exact
  df % interior % hostdata = abs(df % interior % hostdata - 0.0_prec)

  if (maxval(df % interior % hostdata) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  call f % free()
  call df % free()
  call interp % free()

end function scalargradient_2d_gpu_constant
