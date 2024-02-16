program test

  implicit none
  integer :: exit_code
  
  exit_code = mappedscalarbrgradient_2d_cpu_constant()
  stop exit_code

contains
integer function mappedscalarbrgradient_2d_cpu_constant() result(r)

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
  type(Mesh2D),TARGET :: mesh
  type(SEMQuad),TARGET :: geometry
  type(MappedScalar2D) :: f
  type(MappedVector2D) :: df
  type(MPILayer),TARGET :: decomp
  integer :: iel
  integer :: iside
  integer :: i
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
  call mesh % Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5",decomp)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry % Init(interp,mesh % nElem)
  call geometry % GenerateFromMesh(mesh)
  
  call f % Init(interp,nvar,mesh % nelem)
  call df % Init(interp,nvar,mesh % nelem)

  call f % SetEquation( 1, 'f = 1.0')

  call f % SetInteriorFromEquation( geometry, 0.0_prec ) 
  print*, "min, max (interior)", minval(f % interior ), maxval(f % interior )

  call f % BoundaryInterp(.false.)
  print*, "min, max (boundary)", minval(f % boundary ), maxval(f % boundary )

  call f % SideExchange( mesh, decomp, .false.)

  ! Set boundary conditions by prolonging the "boundary" attribute to the domain boundaries
  do iel = 1,f % nElem
    do iside = 1,4
      e2 = mesh % sideInfo % hostData(3,iside,iel) ! Neighboring Element ID
      if (e2 == 0)then
        do i = 0,f % interp % N
          f % extBoundary % hostData(i,1,iside,iel) = f % boundary (i,1,iside,iel) 
        end do
      end if
    end do
  end do

  print*, "min, max (extboundary)", minval(f % extBoundary ), maxval(f % extBoundary )

  call f % Gradient( geometry, df, selfWeakBRForm, .false. ) 

  ! Calculate diff from exact
  df % interior  = abs(df % interior  - 0.0_prec)

  if (maxval(df % interior ) <= tolerance) then
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

end function mappedscalarbrgradient_2d_cpu_constant
end program test
