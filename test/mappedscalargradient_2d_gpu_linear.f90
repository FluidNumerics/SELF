program test

  implicit none
  integer :: exit_code
  
  exit_code = mappedscalargradient_2d_gpu_linear()
  stop exit_code

contains
integer function mappedscalargradient_2d_gpu_linear() result(r)

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Mesh
  use SELF_Geometry
  use SELF_MappedData
  use iso_c_binding
  use hipfort_hipblas
  use HDF5
  use SELF_HDF5

  implicit none

  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
  real(prec),parameter :: tolerance = 10.0_prec**(-2)
#endif
  type(Lagrange),target :: interp
  type(Mesh2D),TARGET :: mesh
  type(SEMQuad),TARGET :: geometry
  type(MappedScalar2D) :: f
  type(MappedVector2D) :: df
  type(MPILayer),TARGET :: decomp
  CHARACTER(LEN=255) :: WORKSPACE
  type(c_ptr) :: handle
  integer :: iel, j, i
  INTEGER(HID_T) :: fileId

  call hipblasCheck(hipblasCreate(handle))

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

  call f % SetEquation( 1, 'f = x*y')

  call f % SetInteriorFromEquation( geometry, 0.0_prec ) 
  print*, "min, max (interior)", minval(f % interior ), maxval(f % interior )

  call f % updatedevice()

  call f % Gradient( geometry, df, handle ) 
  call hipcheck(hipdevicesynchronize())


  CALL f % SetName(1, "f")
  CALL f % SetUnits(1,"[null]")

  CALL Open_HDF5('output.h5',H5F_ACC_TRUNC_F,fileId)

  ! Write the interpolant to the file
  print*, "Writing interpolant data to file"
  CALL f % interp % WriteHDF5( fileId )

  ! Write the model state to file
  print*, "Writing control grid solution to file"
  CALL CreateGroup_HDF5(fileId,'/controlgrid')
  CALL f % WriteHDF5( fileId, '/controlgrid/solution' )

  print*, "Writing control grid solution gradient to file"
  CALL CreateGroup_HDF5(fileId,'/controlgrid')
  CALL df % WriteHDF5( fileId, '/controlgrid/solution_gradient' )

  ! Write the geometry to file
  print*, "Writing control grid  geometry to file"
  CALL CreateGroup_HDF5(fileId,'/controlgrid/geometry')
  CALL geometry % x % WriteHDF5(fileId,'/controlgrid/geometry/x')

  CALL Close_HDF5(fileId)

    ! Calculate diff from exact
  do iel = 1, mesh % nelem
    do j = 1, controlDegree+1
      do i = 1, controlDegree+1
        df % interior(i,j,iel,1,1) = abs(df % interior(i,j,iel,1,1) - geometry % x % interior(i,j,iel,1,2)) ! df/dx = y
        df % interior(i,j,iel,1,2) = abs(df % interior(i,j,iel,1,2) - geometry % x % interior(i,j,iel,1,1)) ! df/dy = x
      end do
    end do
  end do

  print*, maxval(df % interior ), tolerance 


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
  call hipblasCheck(hipblasDestroy(handle))

end function mappedscalargradient_2d_gpu_linear
end program test
