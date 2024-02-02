integer function scalargridinterp_1d_gpu_constant() result(r)
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
  type(Scalar1D) :: f
  type(Scalar1D) :: fTarget
  type(Lagrange),target :: interp
  type(Lagrange),target :: interpTarget

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  call interpTarget % Init(N=targetDegree, &
                           controlNodeType=UNIFORM, &
                           M=targetDegree, &
                           targetNodeType=UNIFORM)

  ! Initialize scalars
  call f % Init(interp,nvar,nelem)
  call fTarget % Init(interpTarget,nvar,nelem)

  ! Set the source scalar (on the control grid) to a non-zero constant
  f % interior % hostdata = 1.0_prec

  call f % interior % updatedevice()

  ! Interpolate with gpuAccel = .true.
  call f % GridInterp(fTarget,.true.)

  call fTarget % interior % updatehost()

  ! Calculate diff from exact
  fTarget % interior % hostdata = abs(fTarget % interior % hostdata - 1.0_prec)

  if (maxval(fTarget % interior % hostdata) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  call f % free()
  call fTarget % free()
  call interp % free()
  call interpTarget % free()

end function scalargridinterp_1d_gpu_constant
