program test

  implicit none
  integer :: exit_code
  
  exit_code = mappedscalargradient_3d_cpu_linear()
  stop exit_code

contains
integer function mappedscalargradient_3d_cpu_linear() result(r)

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
  type(Mesh3D),TARGET :: mesh
  type(SEMHex),TARGET :: geometry
  type(MappedScalar3D) :: f
  type(MappedVector3D) :: df
  type(MPILayer),TARGET :: decomp
  CHARACTER(LEN=255) :: WORKSPACE
  integer :: i,j,k,iel


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

  call f % SetEquation( 1, 'f = x*y*z')

  call f % SetInteriorFromEquation( geometry, 0.0_prec ) 
  print*, "min, max (interior)", minval(f % interior ), maxval(f % interior )

  call f % Gradient( geometry, df ) 

  ! Calculate diff from exact
  ! Calculate diff from exact
  do iel = 1,mesh % nelem
    do k = 1,controlDegree + 1
      do j = 1,controlDegree + 1
        do i = 1,controlDegree + 1
          df % interior(i,j,k,iel,1,1) = abs(df % interior(i,j,k,iel,1,1) - geometry % x % interior(i,j,k,iel,1,2)*geometry % x % interior(i,j,k,iel,1,3)) ! df/dx = y*z
          df % interior(i,j,k,iel,1,2) = abs(df % interior(i,j,k,iel,1,2) - geometry % x % interior(i,j,k,iel,1,1)*geometry % x % interior(i,j,k,iel,1,3)) ! df/dy = x*z
          df % interior(i,j,k,iel,1,3) = abs(df % interior(i,j,k,iel,1,3) - geometry % x % interior(i,j,k,iel,1,1)*geometry % x % interior(i,j,k,iel,1,2)) ! df/dy = x*y
        end do
      end do
    end do
  end do
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

end function mappedscalargradient_3d_cpu_linear
end program test
