PROGRAM SELF_SEMHex_IO

USE SELF_Constants
USE SELF_Lagrange
USE SELF_Mesh
USE SELF_Geometry

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N = 7 ! Polynomial degree of solution
  INTEGER, PARAMETER :: quadrature = GAUSS ! Hexrature
  INTEGER, PARAMETER :: M = 15 ! Number of points in the uniform plotting mesh
  INTEGER, PARAMETER :: nXe = 5 ! Number of elements in the x-direction
  INTEGER, PARAMETER :: nYe = 5 ! Number of elements in the x-direction
  INTEGER, PARAMETER :: nZe = 5 ! Number of elements in the x-direction
  INTEGER, PARAMETER :: nvar = 1 ! The number of tracer fields
  REAL(prec), PARAMETER :: Lx = 1.0_prec ! Length of the domain in the x-direction 
  REAL(prec), PARAMETER :: Ly = 1.0_prec ! Length of the domain in the x-direction 
  REAL(prec), PARAMETER :: Lz = 1.0_prec ! Length of the domain in the x-direction 
  REAL(prec), PARAMETER :: tolerance=500.0_prec*epsilon(1.0_prec) ! Error tolerance

  TYPE(Lagrange),TARGET :: interp
  TYPE(Mesh3D),TARGET :: mesh
  TYPE(SEMHex),TARGET :: geometry
  REAL(prec) :: dxdsExpect(1:3,1:3)
  REAL(prec) :: JExpect
  REAL(prec) :: dxdserror(1:3,1:3)
  REAL(prec) :: Jerror
  INTEGER :: i,j,k,iSide,row,col,iEl
  LOGICAL :: fail
  CHARACTER(LEN=255) :: SELF_PREFIX

    fail = .FALSE.

    ! Create an interpolant
    CALL interp % Init(N,quadrature,M,UNIFORM)

    ! Create a uniform block mesh
    CALL get_environment_variable("SELF_PREFIX", SELF_PREFIX)
    CALL mesh % Read_HOPr(TRIM(SELF_PREFIX)//"Block3D/Block3D_mesh.h5")

    ! Generate geometry (metric terms) from the mesh elements
    CALL geometry % GenerateFromMesh(mesh,interp)

    dxdsExpect(1,1) = 0.5_prec*Lx/REAL(nXe,prec)
    dxdsExpect(2,1) = 0.0_prec
    dxdsExpect(3,1) = 0.0_prec
    dxdsExpect(1,2) = 0.0_prec
    dxdsExpect(2,2) = 0.5_prec*Ly/REAL(nYe,prec)
    dxdsExpect(3,2) = 0.0_prec
    dxdsExpect(1,3) = 0.0_prec
    dxdsExpect(2,3) = 0.0_prec
    dxdsExpect(3,3) = 0.5_prec*Lz/REAL(nZe,prec)
    dxdserror = 0.0_prec
    DO iEl = 1, mesh % nElem
      DO k = 0, interp % N
        DO j = 0, interp % N
          DO i = 0, interp % N
            DO col = 1, 3
              DO row = 1, 3
                dxdserror(row,col) = MAX(ABS(geometry % dxds % interior % hostData(row,col,i,j,k,1,iEl)-&
                        dxdsExpect(row,col)),dxdserror(row,col))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    1000 FORMAT('Element interior metric term error max(|dxds-exact|) greater than tolerance (',I1,',',I1,'): ',E12.5, ' > ', E12.5)
    DO col = 1, 3
      DO row = 1, 3
        IF( dxdserror(row,col) > tolerance )THEN
          WRITE(*,1000) row, col, dxdserror(row,col), tolerance 
          fail=.TRUE.
        ENDIF
      ENDDO
    ENDDO

    dxdserror = 0.0_prec
    DO iEl = 1, mesh % nElem
      DO iSide = 1, 6
        DO j = 0, interp % N
          DO i = 0, interp % N
            DO col = 1, 3
              DO row = 1, 3
                dxdserror(row,col) = MAX(ABS(geometry % dxds % boundary % hostData(row,col,i,j,1,iSide,iEl)-&
                        dxdsExpect(row,col)),dxdserror(row,col))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    1001 FORMAT('Element boundary metric term error max(|dxds-exact|) greater than tolerance (',I1,',',I1,'): ',E12.5, ' > ', E12.5)
    DO col = 1, 3
      DO row = 1, 3
        IF( dxdserror(row,col) > tolerance )THEN
          WRITE(*,1001) row, col, dxdserror(row,col), tolerance 
          fail=.TRUE.
        ENDIF
      ENDDO
    ENDDO

    Jexpect = 0.125_prec*Lx*Ly/REAL(nXe,prec)/REAL(nYe,prec)/REAL(nZe,prec)
    Jerror = 0.0_prec
    DO iEl = 1, mesh % nElem
      DO k = 0, interp % N
        DO j = 0, interp % N
          DO i = 0, interp % N
            Jerror = MAX(ABS(geometry % J % interior % hostData(i,j,k,1,iEl)-&
                        JExpect),Jerror)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    1002 FORMAT('Element interior Jacobian error max(|J-exact|) greater than tolerance : ',E12.5, ' > ', E12.5)
    IF( Jerror > tolerance )THEN
      WRITE(*,1002) Jerror, tolerance 
      fail=.TRUE.
    ENDIF

    Jerror = 0.0_prec
    DO iEl = 1, mesh % nElem
      DO iSide = 1,6
        DO j = 0, interp % N
          DO i = 0, interp % N
            Jerror = MAX(ABS(geometry % J % boundary % hostData(i,j,1,iSide,iEl)-&
                        JExpect),Jerror)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    1003 FORMAT('Element boundary Jacobian error max(|J-exact|) greater than tolerance : ',E12.5, ' > ', E12.5)
    IF( Jerror > tolerance )THEN
      WRITE(*,1003) Jerror, tolerance 
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

END PROGRAM SELF_SEMHex_IO
