!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_LinearShallowWater

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model


  TYPE,EXTENDS(Model2D) :: LinearShallowWater
    !! iVar = 1 ~> u velocity component
    !! iVar = 2 ~> v velocity component
    !! iVar = 3 ~> free surface height
    REAL(prec) :: fCori ! coriolis parameter ( 1/s )
    REAL(prec) :: g     ! gravity ( m/s^2) 
    REAL(prec) :: H     ! fluid thickness ( m ) 

    CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_LinearShallowWater
    PROCEDURE :: CalculateEntropy => CalculateEntropy_LinearShallowWater

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_LinearShallowWater
    PROCEDURE :: FluxMethod => Flux_LinearShallowWater
    PROCEDURE :: RiemannSolver => RiemannSolver_LinearShallowWater
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_LinearShallowWater

  END TYPE LinearShallowWater

  INTERFACE
    SUBROUTINE SetBoundaryCondition_LinearShallowWater_gpu_wrapper(solution, extBoundary, nHat, sideInfo, N, nVar, nEl) &
      bind(c,name="SetBoundaryCondition_LinearShallowWater_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: solution, extBoundary, nHat, sideInfo
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SetBoundaryCondition_LinearShallowWater_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Source_LinearShallowWater_gpu_wrapper(source, solution, f, N, nVar, nEl) &
      bind(c,name="Source_LinearShallowWater_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: source, solution
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: f
    END SUBROUTINE Source_LinearShallowWater_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Flux_LinearShallowWater_gpu_wrapper(flux, solution, g, H, N, nVar, nEl) &
      bind(c,name="Flux_LinearShallowWater_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, solution
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: g, H
    END SUBROUTINE Flux_LinearShallowWater_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE RiemannSolver_LinearShallowWater_gpu_wrapper(flux, solution, extBoundary, nHat, nScale, g, H, N, nVar, nEl) &
      bind(c,name="RiemannSolver_LinearShallowWater_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, solution, extBoundary, nHat, nScale
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: g, H
    END SUBROUTINE RiemannSolver_LinearShallowWater_gpu_wrapper
  END INTERFACE

CONTAINS

  SUBROUTINE Init_LinearShallowWater(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(LinearShallowWater),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh2D),INTENT(in),TARGET :: mesh
    TYPE(SEMQuad),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp
    ! Local
    INTEGER :: ivar
    CHARACTER(LEN=3) :: ivarChar
    CHARACTER(LEN=25) :: varname
    INTEGER :: nvarloc

    ! Ensure that the number of variables is 3
    ! nvar is unused in this class extension
    nvarloc = 3

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.
    this % fluxDivMethod = SELF_CONSERVATIVE_FLUX 
    this % g = 1.0_prec
    this % H = 1.0_prec
    this % fCori = 0.0_prec

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % velocity % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % compVelocity % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    ! First three variables are treated as u, v, eta
    ! Any additional are treated as passive tracers 
    CALL this % solution % SetName(1,"u")
    CALL this % solution % SetUnits(1,"m/s")
    CALL this % solution % SetDescription(1,"x-component of the barotropic velocity field")

    CALL this % solution % SetName(2,"v")
    CALL this % solution % SetUnits(2,"m/s")
    CALL this % solution % SetDescription(2,"y-component of the barotropic velocity field")

    CALL this % solution % SetName(3,"eta")
    CALL this % solution % SetUnits(3,"m")
    CALL this % solution % SetDescription(3,"Free surface height anomaly")

  END SUBROUTINE Init_LinearShallowWater

  SUBROUTINE CalculateEntropy_LinearShallowWater(this)
  !! Base method for calculating entropy of a model
  !! Calculates the entropy as the integration of the 
  !! squared tracer over the domain
    IMPLICIT NONE
    CLASS(LinearShallowWater), INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, iVar, iEl
    REAL(prec) :: Jacobian, u, v, eta
    REAL(prec) :: wi,wj

    ! TO DO : GPU reduction

    this % entropy = 0.0_prec

    DO iEl = 1, this % geometry % x % nElem
        DO j = 0, this % geometry % x % interp % N
          DO i = 0, this % geometry % x % interp % N

            ! Coordinate mapping Jacobian
            Jacobian = this % geometry % J % interior % hostData(i,j,1,iEl)

            ! Quadrature weights
            wi = this % geometry % x % interp % qWeights % hostData(i) 
            wj = this % geometry % x % interp % qWeights % hostData(j) 

            ! Solution
            u = this % solution % interior % hostData(i,j,1,iEl)
            v = this % solution % interior % hostData(i,j,2,iEl)
            eta = this % solution % interior % hostData(i,j,3,iEl)

            this % entropy = this % entropy + 0.5_prec*( this % H*(u*u + v*v) + this % g*eta*eta )*Jacobian*wi*wj

          ENDDO
        ENDDO
    ENDDO

  END SUBROUTINE CalculateEntropy_LinearShallowWater

  SUBROUTINE SetBoundaryCondition_LinearShallowWater(this)
    IMPLICIT NONE
    CLASS(LinearShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, i
    INTEGER :: bcid, e2
    REAL(prec) :: u, v, nhat(1:2)


    IF( this % gpuAccel )THEN

      CALL SetBoundaryCondition_LinearShallowWater_gpu_wrapper( this % solution % boundary % deviceData,&
            this % solution % extBoundary % deviceData, &
            this % geometry % nHat % boundary % deviceData, &
            this % mesh % sideInfo % deviceData, &
            this % solution % interp % N, &
            this % solution % nVar, &
            this % solution % nElem)

    ELSE

      DO iEl = 1, this % solution % nElem
        DO iSide = 1, 4
            DO i = 0, this % solution % interp % N

              bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
              e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
              IF( e2 == 0 )THEN
                IF( bcid == SELF_BC_RADIATION )THEN

                  this % solution % extBoundary % hostData(i,1,iSide,iEl) = 0.0_prec
                  this % solution % extBoundary % hostData(i,2,iSide,iEl) = 0.0_prec
                  this % solution % extBoundary % hostData(i,3,iSide,iEl) = 0.0_prec

                ELSEIF( bcid == SELF_BC_NONORMALFLOW )THEN

                  nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
                  u = this % solution % boundary % hostData(i,1,iSide,iEl) 
                  v = this % solution % boundary % hostData(i,2,iSide,iEl) 
                  this % solution % extBoundary % hostData(i,1,iSide,iEl) = (nhat(2)**2 - nhat(1)**2)*u - 2.0_prec*nhat(1)*nhat(2)*v
                  this % solution % extBoundary % hostData(i,2,iSide,iEl) = (nhat(1)**2 - nhat(2)**2)*v - 2.0_prec*nhat(1)*nhat(2)*u
                  this % solution % extBoundary % hostData(i,3,iSide,iEl) = this % solution % boundary % hostData(i,3,iSide,iEl)

                ELSE ! Default boundary condition is radiation

                  this % solution % extBoundary % hostData(i,1,iSide,iEl) = 0.0_prec
                  this % solution % extBoundary % hostData(i,2,iSide,iEl) = 0.0_prec
                  this % solution % extBoundary % hostData(i,3,iSide,iEl) = 0.0_prec

                ENDIF
              ENDIF

            ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE SetBoundaryCondition_LinearShallowWater 

  SUBROUTINE Source_LinearShallowWater(this)
    IMPLICIT NONE
    CLASS(LinearShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    IF( this % gpuAccel )THEN

      CALL Source_LinearShallowWater_gpu_wrapper(this % source % interior % deviceData,&
             this % solution % interior % deviceData, &
             this % fCori,&
             this % source % interp % N, &
             this % solution % nVar, &
             this % solution % nElem )

    ELSE

      DO iEl = 1, this % source % nElem
        DO j = 0, this % source % interp % N
          DO i = 0, this % source % interp % N

            this % source % interior % hostData(i,j,1,iEl) = -this % fCori*this % solution % interior % hostData(i,j,2,iEl)
            this % source % interior % hostData(i,j,2,iEl) = this % fCori*this % solution % interior % hostData(i,j,1,iEl)
            this % source % interior % hostData(i,j,3,iEl) = 0.0_prec

          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE Source_LinearShallowWater

  SUBROUTINE Flux_LinearShallowWater(this)
    IMPLICIT NONE
    CLASS(LinearShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    IF( this % gpuAccel )THEN

      CALL Flux_LinearShallowWater_gpu_wrapper(this % flux % interior % deviceData,&
                                               this % solution % interior % deviceData, &
                                               this % g, this % H, this % solution % interp % N, &
                                               this % solution % nVar, this % solution % nElem)

    ELSE
      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              IF ( iVar == 1 )THEN ! u-velocity
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % g*this % solution % interior % hostData(i,j,3,iEl)

                this % flux % interior % hostData(2,i,j,iVar,iEl) = 0.0_prec

              ELSEIF ( iVar == 2 )THEN ! v-velocity

                this % flux % interior % hostData(1,i,j,iVar,iEl) = 0.0_prec

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % g*this % solution % interior % hostData(i,j,3,iEl)


              ELSEIF ( iVar == 3 )THEN ! free surface height
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % H*this % solution % interior % hostData(i,j,1,iEl)

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % H*this % solution % interior % hostData(i,j,2,iEl)
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE Flux_LinearShallowWater

  SUBROUTINE RiemannSolver_LinearShallowWater(this)
    IMPLICIT NONE
    CLASS(LinearShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: nhat(1:2), nmag
    REAL(prec) :: c, unL, unR, etaL, etaR, wL, wR


    IF( this % gpuAccel )THEN

      CALL RiemannSolver_LinearShallowWater_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
             this % solution % boundary % deviceData, &
             this % solution % extBoundary % deviceData, &
             this % geometry % nHat % boundary % deviceData, &
             this % geometry % nScale % boundary % deviceData, &
             this % g, this % H, &
             this % solution % interp % N, &
             this % solution % nVar, &
             this % solution % nElem)

    ELSE

      DO iEl = 1, this % solution % nElem
        DO iSide = 1, 4
          DO i = 0, this % solution % interp % N

             ! Get the boundary normals on cell edges from the mesh geometry
             nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
             nmag = this % geometry % nScale % boundary % hostData(i,1,iSide,iEl)
             c = sqrt( this % g * this % H )

             ! Calculate the normal velocity at the cell edges
             unL = this % solution % boundary % hostData(i,1,iSide,iEl)*nHat(1)+&
                   this % solution % boundary % hostData(i,2,iSide,iEl)*nHat(2)

             unR = this % solution % extBoundary % hostData(i,1,iSide,iEl)*nHat(1)+&
                   this % solution % extBoundary % hostData(i,2,iSide,iEl)*nHat(2)

             etaL = this % solution % boundary % hostData(i,3,iSide,iEl)
             etaR = this % solution % extBoundary % hostData(i,3,iSide,iEl)

             ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
             wL = 0.5_prec*(unL/this % g + etaL/c)
             wR = 0.5_prec*(unR/this % g - etaR/c)

             this % flux % boundaryNormal % hostData(i,1,iSide,iEl) = this % g*c*( wL - wR )*nHat(1)*nmag
             this % flux % boundaryNormal % hostData(i,2,iSide,iEl) = this % g*c*( wL - wR )*nHat(2)*nmag
             this % flux % boundaryNormal % hostData(i,3,iSide,iEl) = c*c*( wL + wR )*nmag

          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE RiemannSolver_LinearShallowWater

END MODULE SELF_LinearShallowWater
