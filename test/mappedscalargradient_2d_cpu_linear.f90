program test

  implicit none
  integer :: exit_code

  exit_code = mappedscalargradient_2d_cpu_linear()
  stop exit_code

contains
  integer function mappedscalargradient_2d_cpu_linear() result(r)

    use SELF_constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_MappedVector_2D
    use SELF_MPI

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
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    type(MappedVector2D) :: df
    type(MPILayer),target :: decomp
    character(LEN=255) :: WORKSPACE
    integer :: iel,j,i
    integer(HID_T) :: fileId

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    call decomp%Init(enableMPI=.false.)

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Create a uniform block mesh
    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5",decomp)

    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call df%Init(interp,nvar,mesh%nelem)

    call f%SetEquation(1,'f = x*y')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    print*,"min, max (interior)",minval(f%interior),maxval(f%interior)

    call f%Gradient(geometry,df%interior)

    call f%SetName(1,"f")
    call f%SetUnits(1,"[null]")

    call Open_HDF5('output.h5',H5F_ACC_TRUNC_F,fileId)

    ! Write the interpolant to the file
    print*,"Writing interpolant data to file"
    call f%interp%WriteHDF5(fileId)

    ! Write the model state to file
    print*,"Writing control grid solution to file"
    call CreateGroup_HDF5(fileId,'/controlgrid')
    call f%WriteHDF5(fileId,'/controlgrid/solution')

    print*,"Writing control grid solution gradient to file"
    call CreateGroup_HDF5(fileId,'/controlgrid')
    call df%WriteHDF5(fileId,'/controlgrid/solution_gradient')

    ! Write the geometry to file
    print*,"Writing control grid  geometry to file"
    call CreateGroup_HDF5(fileId,'/controlgrid/geometry')
    call geometry%x%WriteHDF5(fileId,'/controlgrid/geometry/x')

    call Close_HDF5(fileId)

    ! Calculate diff from exact
    do iel = 1,mesh%nelem
      do j = 1,controlDegree+1
        do i = 1,controlDegree+1
          df%interior(i,j,iel,1,1) = abs(df%interior(i,j,iel,1,1)-geometry%x%interior(i,j,iel,1,2)) ! df/dx = y
          df%interior(i,j,iel,1,2) = abs(df%interior(i,j,iel,1,2)-geometry%x%interior(i,j,iel,1,1)) ! df/dy = x
        enddo
      enddo
    enddo

    print*,maxval(df%interior),tolerance
    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    ! Clean up
    call decomp%Free()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()
    call f%free()
    call df%free()

    r = 0

  endfunction mappedscalargradient_2d_cpu_linear
endprogram test
