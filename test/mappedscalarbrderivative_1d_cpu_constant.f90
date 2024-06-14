program test

  implicit none
  integer :: exit_code

  exit_code = mappedscalarbrderivative_1d_cpu_constant()
  stop exit_code

contains
  integer function mappedscalarbrderivative_1d_cpu_constant() result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_MappedScalar_1D
    use SELF_Mesh_1D
    use SELF_Geometry_1D
    use SELF_MPI

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
    type(Mesh1D),target :: mesh
    type(Geometry1D),target :: geometry
    type(MPILayer),target :: decomp

    call decomp%Init(enableMPI=.false.)
    call mesh%UniformBlockMesh(nGeo=1, &
                               nElem=nelem, &
                               x=(/0.0_prec,10.0_prec/))

    call decomp%GenerateDecomposition(nelem,nelem+1)

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    ! Initialize scalars
    call f%Init(interp,nvar,nelem)
    call df%Init(interp,nvar,nelem)

    call f%SetEquation(1,'f = 1.0')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    print*,"min, max (interior)",minval(f%interior),maxval(f%interior)

    call f%BoundaryInterp()
    print*,"min, max (boundary)",minval(f%boundary),maxval(f%boundary)

    call f%SideExchange(mesh,decomp)
    ! Set boundary conditions
    f%extBoundary(1,1,1) = 1.0_prec ! Left most
    f%extBoundary(2,nelem,1) = 1.0_prec ! Right most
    print*,"min, max (extboundary)",minval(f%extBoundary),maxval(f%extBoundary)

    call f%BassiRebaySides()
    print*,"min, max (avgboundary)",minval(f%avgBoundary),maxval(f%avgBoundary)

    call f%BRDerivative(geometry,df%interior)
    ! Calculate diff from exact
    df%interior = abs(df%interior-0.0_prec)

    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    ! Clean up
    call decomp%Free()
    call mesh%Free()
    call geometry%Free()
    call interp%free()
    call f%free()
    call df%free()

  endfunction mappedscalarbrderivative_1d_cpu_constant
endprogram test
