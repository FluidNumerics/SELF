integer function vectorboundaryinterp_3d_gpu_constant() result(r)
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
  type(Vector3D) :: f
  type(Lagrange),target :: interp

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  ! Initialize vectors
  call f % Init(interp,nvar,nelem)

  ! Set the source vector (on the control grid) to a non-zero constant
  f % interior % hostdata = 1.0_prec

  ! copy data from host to device
  call f % interior % updatedevice()

  ! Interpolate with gpuAccel = .true.
  call f % BoundaryInterp(.true.)

  ! copy data from device to host
  call f % boundary % updatehost()

  ! Calculate diff from exact
  f % boundary % hostdata = abs(f % boundary % hostdata - 1.0_prec)

  if (maxval(f % boundary % hostdata) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  call f % free()
  call interp % free()

end function vectorboundaryinterp_3d_gpu_constant
