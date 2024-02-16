program test

  implicit none
  integer :: exit_code
  
  exit_code = mesh3d_setup()
  stop exit_code

contains
integer function mesh3d_setup() result(r)

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh
  use SELF_Geometry

  implicit none

  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
  real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
  type(Lagrange),target :: interp
  type(Mesh3D),TARGET :: mesh
  type(SEMHex),TARGET :: geometry
  type(MPILayer),TARGET :: decomp
  CHARACTER(LEN=255) :: WORKSPACE


  call get_environment_variable("WORKSPACE",WORKSPACE)

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

  ! Clean up
  call decomp % Free()
  call geometry % Free()
  call mesh % Free()
  call interp % Free()
    
  r = 0

end function mesh3d_setup
end program test
