PROGRAM SELF_Geometry1D_IO

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry

  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 7 ! Polynomial degree of solution
  INTEGER, PARAMETER :: quadrature = GAUSS ! Quadrature
  INTEGER, PARAMETER :: M = 15 ! Number of points in the uniform plotting mesh
  INTEGER, PARAMETER :: nXe = 1 ! Number of elements in the x-direction
  INTEGER, PARAMETER :: nvar = 1 ! The number of tracer fields
  REAL(prec), PARAMETER :: Lx = 1.0_prec ! Length of the domain in the x-direction 
  REAL(prec), PARAMETER :: tolerance=10.0_prec*epsilon(1.0_prec) ! Error tolerance

  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh1D),TARGET :: mesh
  TYPE(Geometry1D),TARGET :: geometry
  REAL(prec) :: dxdsExpect
  REAL(prec) :: error
  INTEGER :: i, iEl
  LOGICAL :: fail


    fail = .FALSE.

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL mesh % UniformBlockMesh(N,nXe,(/0.0_prec,Lx/))

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % Init(interp,mesh % nElem)
    CALL geometry % GenerateFromMesh(mesh)

    dxdsExpect = 0.5_prec*Lx/REAL(nXe,prec)
    error = 0.0_prec
    DO iEl = 1, mesh % nElem
      DO i = 0, interp % N
        error = MAX(ABS(geometry % dxds % interior % hostData(i,1,iEl)-dxdsExpect),error)
      ENDDO
    ENDDO

    IF( error > tolerance )THEN
      1000 FORMAT('Element interior metric term error max(|dxds-exact|) greater than tolerance : ',E12.5, ' > ', E12.5)
      WRITE(*,1000) error, tolerance 
      fail=.TRUE.
    ENDIF

    error = 0.0_prec
    DO iEl = 1, mesh % nElem
      error = MAX(ABS(geometry % dxds % boundary % hostData(1,1,iEl)-dxdsExpect),error)
      error = MAX(ABS(geometry % dxds % boundary % hostData(1,2,iEl)-dxdsExpect),error)
    ENDDO

    IF( error > tolerance )THEN
      1001 FORMAT('Element boundary metric term error max(|dxds-exact|) greater than tolerance : ',E12.5, ' > ', E12.5)
      WRITE(*,1001) error, tolerance 
      fail=.TRUE.
    ENDIF

    CALL geometry % Write()
    CALL geometry % WriteTecplot()

    CALL geometry % Free()
    CALL mesh % Free()
    CALL interp % Free()

    IF( fail )THEN
      STOP 1
    ENDIF

END PROGRAM SELF_Geometry1D_IO
