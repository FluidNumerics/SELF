program test

  implicit none
  integer :: exit_code
  
  exit_code = mappedscalardgderivative_1d_gpu_constant()
  stop exit_code

contains
integer function mappedscalardgderivative_1d_gpu_constant() result(r)
  use SELF_Constants
  use SELF_Lagrange
  use SELF_MappedData
  use SELF_Mesh
  use SELF_Geometry
  use iso_c_binding
  use hipfort_hipblas


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
  type(MappedScalar1D) :: f
  type(MappedScalar1D) :: df
  type(Lagrange),target :: interp
  type(Mesh1D),TARGET :: mesh
  type(Geometry1D),TARGET :: geometry
  type(MPILayer),TARGET :: decomp
  type(c_ptr) :: handle

  call hipblasCheck(hipblasCreate(handle))

  call decomp % Init(enableMPI=.false.)
  call mesh % UniformBlockMesh(nGeo=1,&
                               nElem=nelem,&
                               x=(/0.0_prec,10.0_prec/))

  CALL decomp % GenerateDecomposition(nelem,nelem+1)

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry % Init(interp,mesh % nElem)
  call geometry % GenerateFromMesh(mesh)

  call geometry % updatedevice()
  call mesh % updatedevice()

  ! Initialize scalars
  call f % Init(interp,nvar,nelem)
  call df % Init(interp,nvar,nelem)

  call f % SetEquation( 1, 'f = 1.0')
  call f % SetInteriorFromEquation( geometry, 0.0_prec ) 
  print*, "min, max (interior)", minval(f % interior ), maxval(f % interior )

  call f % updatedevice()

  call f % BoundaryInterp(handle)

  ! Set boundary conditions

  f % boundary(1,1,1) = 1.0_prec ! Left most
  f % boundary(2,nelem,1) = 1.0_prec ! Right most
  
  ! Adjust for -\hat{x} direction on left element boundaries
  f % boundary(1,:,1) = -f % boundary(1,:,1)

  print*, "min, max (boundary)", minval(f % boundary ), maxval(f % boundary )

  call f % DGDerivative(geometry, df, handle)
  call hipcheck(hipdevicesynchronize())

  ! Calculate diff from exact
  df % interior  = abs(df % interior  - 0.0_prec)

  if (maxval(df % interior ) <= tolerance) then
    r = 0
  else
    print*, "AbsMax error: ", maxval(df % interior )
    r = 1
  end if

  ! Clean up
  call decomp % Free()
  call mesh % Free()
  call geometry % Free()
  call interp % free()
  call f % free()
  call df % free()
  call hipblasCheck(hipblasDestroy(handle))

end function mappedscalardgderivative_1d_gpu_constant

end program test
