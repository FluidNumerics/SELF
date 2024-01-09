integer function mappedscalargradient_2d_gpu_constant() result(r)

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh
  use SELF_Geometry
  use SELF_MappedData

  implicit none

  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
  real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
  type(Lagrange),target :: interp
  type(Mesh2D),TARGET :: mesh
  type(SEMQuad),TARGET :: geometry
  type(MappedScalar2D) :: f
  type(MappedVector2D) :: df
  type(MPILayer),TARGET :: decomp
  CHARACTER(LEN=255) :: WORKSPACE


  ! Initialize a domain decomposition
  ! Here MPI is disabled, since scaling is currently
  ! atrocious with the uniform block mesh
  call decomp % Init(enableMPI=.false.)

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  ! Create a uniform block mesh
  call get_environment_variable("WORKSPACE",WORKSPACE)
  call mesh % Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5",decomp)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry % Init(interp,mesh % nElem)
  call geometry % GenerateFromMesh(mesh)
  
  call f % Init(interp,nvar,mesh % nelem)
  call df % Init(interp,nvar,mesh % nelem)

  call f % SetEquation( 1, 'f = 1.0')

  call f % SetInteriorFromEquation( geometry, 0.0_prec ) 
  print*, "min, max (interior)", minval(f % interior % hostdata), maxval(f % interior % hostdata)

!  call f % BoundaryInterp(.false.)
!  print*, "min, max (boundary)", minval(f % boundary % hostdata), maxval(f % boundary % hostdata)

!  call f % SideExchange( mesh, decomp, .false.)
  ! Set boundary conditions
!  f % extBoundary % hostData(1,1,1) = 1.0_prec ! Left most
!  f % extBoundary % hostData(1,2,nelem) = 1.0_prec ! Right most
!  print*, "min, max (extboundary)", minval(f % extBoundary % hostdata), maxval(f % extBoundary % hostdata)

!  call f % BassiRebaySides(.false.)
!  print*, "min, max (avgboundary)", minval(f % avgBoundary % hostdata), maxval(f % avgBoundary % hostdata)

  call f % interior % updatedevice()

  call f % Gradient( geometry, df, selfStrongForm, .true. ) 

  call df % interior % updatehost()

  ! Calculate diff from exact
  df % interior % hostdata = abs(df % interior % hostdata - 0.0_prec)

  if (maxval(df % interior % hostdata) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  ! Clean up
  call decomp % Free()
  call geometry % Free()
  call mesh % Free()
  call interp % Free()
  call f % free()
  call df % free()
    
  r = 0

end function mappedscalargradient_2d_gpu_constant
