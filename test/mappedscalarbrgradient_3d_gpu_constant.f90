integer function mappedscalarbrgradient_3d_gpu_constant() result(r)

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh
  use SELF_Geometry
  use SELF_MappedData

  implicit none

  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  integer,parameter :: nvar = 1
#ifdef doUBLE_PRECISION
  real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
  real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
  type(Lagrange),target :: interp
  type(Mesh3D),TARGET :: mesh
  type(SEMHex),TARGET :: geometry
  type(MappedScalar3D) :: f
  type(MappedVector3D) :: df
  type(MPILayer),TARGET :: decomp
  integer :: iel
  integer :: iside
  integer :: i
  integer :: j
  integer :: e2
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
  call mesh % Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D_mesh.h5",decomp)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry % Init(interp,mesh % nElem)
  call geometry % GenerateFromMesh(mesh)
  
  call f % Init(interp,nvar,mesh % nelem)
  call df % Init(interp,nvar,mesh % nelem)

  call f % SetEquation( 1, 'f = 1.0')

  call f % SetInteriorFromEquation( geometry, 0.0_prec ) 
  print*, "min, max (interior)", minval(f % interior % hostdata), maxval(f % interior % hostdata)

  call f % interior % updatedevice()

  call f % BoundaryInterp(.true.)
  print*, "min, max (boundary)", minval(f % boundary % hostdata), maxval(f % boundary % hostdata)

  call f % SideExchange( mesh, decomp, .true.)

  call f % boundary % updatehost()
  call f % extboundary % updatehost()

  ! Set boundary conditions by prolonging the "boundary" attribute to the domain boundaries
  do iel = 1,f % nElem
    do iside = 1,4
      e2 = mesh % sideInfo % hostData(3,iside,iel) ! Neighboring Element ID
      if (e2 == 0)then
        do j = 0,f % interp % N
          do i = 0,f % interp % N
            f % extBoundary % hostData(i,j,1,iside,iel) = f % boundary % hostdata(i,j,1,iside,iel) 
          end do
        end do
      end if
    end do
  end do

  call f % extboundary % updatedevice()

  print*, "min, max (extboundary)", minval(f % extBoundary % hostdata), maxval(f % extBoundary % hostdata)

  call f % Gradient( geometry, df, selfWeakBRForm, .true. ) 

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

end function mappedscalarbrgradient_3d_gpu_constant
