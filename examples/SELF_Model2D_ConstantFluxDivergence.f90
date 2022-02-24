PROGRAM SELF_Model2D_ConstantFluxDivergence

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_Model

  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 7 ! Polynomial degree of solution
  INTEGER, PARAMETER :: quadrature = GAUSS ! Quadrature
  INTEGER, PARAMETER :: M = 15 ! Number of points in the uniform plotting mesh
  INTEGER, PARAMETER :: nXe = 1 ! Number of elements in the x-direction
  INTEGER, PARAMETER :: nYe = 1 ! Number of elements in the y-direction
  INTEGER, PARAMETER :: nvar = 1 ! The number of tracer fields
  REAL(prec), PARAMETER :: Lx = 1.0_prec ! Length of the domain in the x-direction 
  REAL(prec), PARAMETER :: Ly = 1.0_prec ! Length of the domain in the y-direction 
  REAL(prec), PARAMETER :: dt = 0.001_prec ! Time step size
  REAL(prec), PARAMETER :: tn = 0.5_prec ! File time 
  REAL(prec), PARAMETER :: tolerance=5000.0_prec*epsilon(1.0_prec) ! Error tolerance


  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh2D),TARGET :: mesh
  TYPE(SEMQuad),TARGET :: geometry
  TYPE(MPILayer),TARGET :: decomp
  TYPE(Model2D) :: semModel
  INTEGER :: i,j,iSide,iEl
  REAL(prec) :: nhat(1:2),nmag,fn
  REAL(prec) :: dFError(1:nvar)
  LOGICAL :: fail 
  CHARACTER(LEN=255) :: SELF_PREFIX

    fail = .FALSE.
    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=.FALSE.)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
    CALL mesh % Read_HOPr(TRIM(SELF_PREFIX)//"/etc/mesh/Block2D/Block2D_mesh.h5")

    ! Generate a decomposition
    CALL decomp % GenerateDecomposition(mesh)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % Init(interp,mesh % nElem)
    CALL geometry % GenerateFromMesh(mesh)

    ! Initialize the model
    CALL semModel % Init(nvar,mesh,geometry,decomp)

    ! Set the flux vector field
    DO iEl = 1, mesh % nElem
      DO j = 0, interp % N
        DO i = 0, interp % N
          semModel % flux % interior % hostData(1,i,j,1,iEl) = 1.0_prec
          semModel % flux % interior % hostData(2,i,j,1,iEl) = 1.0_prec
        ENDDO
      ENDDO
    ENDDO

    CALL semModel % flux % boundaryInterp(semModel % gpuAccel)

    ! Compute \vec{F} \dot \hat{n}
    ! Similar to "RiemannSolver"
    DO iEl = 1, mesh % nElem
      DO iSide = 1, 4
        DO i = 0, interp % N

           ! Get the boundary normals on cell edges from the mesh geometry
           nhat(1:2) = geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)

           fn = semModel % flux % boundary % hostData(1,i,1,iSide,iEl)*nhat(1)+&
                      semModel % flux % boundary % hostData(1,i,1,iSide,iEl)*nhat(2)

           nmag = semModel % geometry % nScale % boundary % hostData(i,1,iSide,iEl)
           ! Calculate the flux
           semModel % flux % boundaryNormal % hostData(i,1,iSide,iEl) = fn*nmag

        ENDDO
      ENDDO
    ENDDO

    CALL semModel % ReprojectFlux()
    CALL semModel % CalculateFluxDivergence() 

    ! TO DO !
    ! Compute max error in the flux divergence !
    1003 FORMAT('Flux divergence error max(|divF-exact|) greater than tolerance : ',E12.5, ' > ', E12.5)
    dfError = semModel % fluxDivergence % AbsMaxInterior()

    IF( dfError(1) > tolerance )THEN
      WRITE(*,1003) dfError(1), tolerance 
      fail=.TRUE.
    ENDIF

    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()
    CALL semModel % Free()

END PROGRAM SELF_Model2D_ConstantFluxDivergence
