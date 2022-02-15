PROGRAM SELF_Model2D_FluxDivergence

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry
USE SELF_Model

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


  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh2D),TARGET :: mesh
  TYPE(SEMQuad),TARGET :: geometry
  TYPE(MPILayer),TARGET :: decomp
  TYPE(Model2D) :: semModel

    ! Initialize a domain decomposition
    ! Here MPI is disabled, since scaling is currently
    ! atrocious with the uniform block mesh
    CALL decomp % Init(enableMPI=.FALSE.)

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL mesh % UniformBlockMesh(N,(/nXe,nYe/),(/0.0_prec,Lx,0.0_prec,Ly/))

    ! Generate a decomposition
     CALL decomp % GenerateDecomposition(mesh)

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % GenerateFromMesh(mesh,interp)

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

    CALL semModel % CalculateFluxDivergence() 

    PRINT*, semModel % fluxDivergence % AbsMaxInterior()

    CALL decomp % Free()
    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()
    CALL semModel % Free()

END PROGRAM SELF_Model2D_FluxDivergence
