program test

  implicit none
  integer :: exit_code

  exit_code = mappedscalarbrgradient_3d_cpu_linear()
  stop exit_code

contains
  integer function mappedscalarbrgradient_3d_cpu_linear() result(r)

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
    type(Mesh3D),target :: mesh
    type(SEMHex),target :: geometry
    type(MappedScalar3D) :: f
    type(MappedVector3D) :: df
    type(MPILayer),target :: decomp
    integer :: iel
    integer :: iside
    integer :: i
    integer :: j
    integer :: k
    integer :: e2,s2,bcid
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
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block3D/Block3D_mesh.h5",decomp)

    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call df%Init(interp,nvar,mesh%nelem)

    call f%SetEquation(1,'f = x*y*z')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    print*,"min, max (interior)",minval(f%interior),maxval(f%interior)

    call f%BoundaryInterp()
    print*,"min, max (boundary)",minval(f%boundary),maxval(f%boundary)

    call f%SideExchange(mesh,decomp)

    ! Set boundary conditions by prolonging the "boundary" attribute to the domain boundaries
    do iel = 1,f%nElem
      do iside = 1,6
        e2 = mesh%sideInfo(3,iside,iel) ! Neighboring Element ID
        s2 = mesh%sideInfo(4,iside,iel)/10
        bcid = mesh%sideInfo(5,iside,iel)
        if(s2 == 0 .or. bcid /= 0) then
          do j = 1,f%interp%N+1
            do i = 1,f%interp%N+1
              f%extBoundary(i,j,iside,iel,1) = f%boundary(i,j,iside,iel,1)
            enddo
          enddo
        endif
      enddo
    enddo

    print*,"min, max (extboundary)",minval(f%extBoundary),maxval(f%extBoundary)

    call f%BRGradient(geometry,df)

    ! Calculate diff from exact
    do iel = 1,mesh%nelem
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            df%interior(i,j,k,iel,1,1) = abs(df%interior(i,j,k,iel,1,1)- &
                                             geometry%x%interior(i,j,k,iel,1,2)* &
                                             geometry%x%interior(i,j,k,iel,1,3)) ! df/dx = y*z
            df%interior(i,j,k,iel,1,2) = abs(df%interior(i,j,k,iel,1,2)- &
                                             geometry%x%interior(i,j,k,iel,1,1)* &
                                             geometry%x%interior(i,j,k,iel,1,3)) ! df/dy = x*z
            df%interior(i,j,k,iel,1,3) = abs(df%interior(i,j,k,iel,1,3)- &
                                             geometry%x%interior(i,j,k,iel,1,1)* &
                                             geometry%x%interior(i,j,k,iel,1,2)) ! df/dy = x*y
          enddo
        enddo
      enddo
    enddo
    print*,"maxval(df_error)",maxval(df%interior),tolerance

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

  endfunction mappedscalarbrgradient_3d_cpu_linear
endprogram test
