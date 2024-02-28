program test

  implicit none
  integer :: exit_code
  
  exit_code = mappedscalarderivative_1d_gpu_constant()
  stop exit_code

contains
integer function mappedscalarderivative_1d_gpu_constant() result(r)
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh
  use SELF_Geometry
  use SELF_MappedData
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
  type(c_ptr) :: handle

  call hipblasCheck(hipblasCreate(handle))

  call mesh % UniformBlockMesh(nGeo=1,&
                               nElem=nelem,&
                               x=(/0.0_prec,10.0_prec/))

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

  call f % updatedevice()

  call f % Derivative(geometry, df, handle)
  call hipcheck(hipdevicesynchronize())

  ! Calculate diff from exact
  df % interior  = abs(df % interior  - 0.0_prec)

  if (maxval(df % interior ) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  ! Clean up
  call mesh % Free()
  call geometry % Free()
  call interp % free()
  call f % free()
  call df % free()
  call hipblasCheck(hipblasDestroy(handle))


end function mappedscalarderivative_1d_gpu_constant
end program test
