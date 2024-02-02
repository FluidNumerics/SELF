integer function mappedscalarbrderivative_1d_cpu_constant() result(r)
  use SELF_Constants
  use SELF_Memory
  use SELF_Lagrange
  use SELF_MappedData
  use SELF_Mesh
  use SELF_Geometry

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


  ! Initialize scalars
  call f % Init(interp,nvar,nelem)
  call df % Init(interp,nvar,nelem)

  call f % SetEquation( 1, 'f = 1.0')

  call f % SetInteriorFromEquation( geometry, 0.0_prec ) 
  print*, "min, max (interior)", minval(f % interior % hostdata), maxval(f % interior % hostdata)

  call f % BoundaryInterp(.false.)
  print*, "min, max (boundary)", minval(f % boundary % hostdata), maxval(f % boundary % hostdata)

  call f % SideExchange( mesh, decomp, .false.)
  ! Set boundary conditions
  f % extBoundary % hostData(1,1,1) = 1.0_prec ! Left most
  f % extBoundary % hostData(1,2,nelem) = 1.0_prec ! Right most
  print*, "min, max (extboundary)", minval(f % extBoundary % hostdata), maxval(f % extBoundary % hostdata)

  call f % BassiRebaySides(.false.)
  print*, "min, max (avgboundary)", minval(f % avgBoundary % hostdata), maxval(f % avgBoundary % hostdata)

  call f % Derivative(geometry, df, selfWeakBRForm, .false.)

  ! Calculate diff from exact
  df % interior % hostdata = abs(df % interior % hostdata - 0.0_prec)

  if (maxval(df % interior % hostdata) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  ! Clean up
  call decomp % Free()
  call mesh % Free()
  call geometry % Free()
  call interp % free()
  call f % free()
  call df % free()

end function mappedscalarbrderivative_1d_cpu_constant
