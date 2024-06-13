program test

  implicit none
  integer :: exit_code

  exit_code = mappedscalargradient_2d_cpu_constant()
  stop exit_code

contains
  integer function mappedscalargradient_2d_cpu_constant() result(r)

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
    real(prec),parameter :: tolerance = 10.0_prec**(-2)
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    type(MappedVector2D) :: df
    type(MPILayer),target :: decomp
    character(LEN=255) :: WORKSPACE

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

    call f%SetEquation(1,'f = 1.0')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    print*,"min, max (interior)",minval(f%interior),maxval(f%interior)

    call f%Gradient(geometry,df)

    ! Calculate diff from exact
    df%interior = abs(df%interior-0.0_prec)
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

  endfunction mappedscalargradient_2d_cpu_constant
endprogram test
