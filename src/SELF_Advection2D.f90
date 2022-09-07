!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Advection2D

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model2D


  TYPE,EXTENDS(Model2D) :: Advection2D

    CONTAINS

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_Advection2D
    PROCEDURE :: FluxMethod => Flux_Advection2D
    PROCEDURE :: RiemannSolver => RiemannSolver_Advection2D
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_Advection2D
    
    ! Overridden Methods
    PROCEDURE :: CalculateEntropy => CalculateEntropy_Advection2D

  END TYPE Advection2D

CONTAINS

  SUBROUTINE CalculateEntropy_Advection2D(this)
  !! Base method for calculating entropy of a model
  !! Calculates the entropy as the integration of the 
  !! squared tracer over the domain
    IMPLICIT NONE
    CLASS(Advection2D), INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, iVar, iEl
    REAL(prec) :: Jacobian, s
    REAL(prec) :: wi,wj

    ! TO DO : GPU reduction

    this % entropy = 0.0_prec

    DO iEl = 1, this % geometry % x % nElem
      DO iVar = 1, this % geometry % x % nVar
        DO j = 0, this % geometry % x % interp % N
          DO i = 0, this % geometry % x % interp % N

            ! Coordinate mapping Jacobian
            Jacobian = this % geometry % J % interior % hostData(i,j,1,iEl)

            ! Quadrature weights
            wi = this % geometry % x % interp % qWeights % hostData(i) 
            wj = this % geometry % x % interp % qWeights % hostData(j) 

            ! Solution
            s = this % solution % interior % hostData(i,j,iVar,iEl)

            this % entropy = this % entropy + s*s*Jacobian*wi*wj

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE CalculateEntropy_Advection2D

  SUBROUTINE SetBoundaryCondition_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl,iVar,e2      

    DO iEl = 1, this % solution % nElem
      DO iSide = 1, 4
        DO iVar = 1, this % solution % nVar
          DO i = 0, this % solution % interp % N

            e2 = this % mesh % sideInfo % hostData(3,iSide,iEl)

            IF( e2 == 0 )THEN
              this % solution % extBoundary % hostData(i,iVar,iSide,iEl) = 0.0_prec
            ENDIF

          ENDDO
        ENDDO
      ENDDO
    ENDDO


  END SUBROUTINE SetBoundaryCondition_Advection2D 

  SUBROUTINE Source_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    DO iEl = 1, this % source % nElem
      DO iVar = 1, this % source % nVar
        DO j = 0, this % source % interp % N
          DO i = 0, this % source % interp % N

            this % source % interior % hostData(i,j,iVar,iEl) = 0.0_prec

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE Source_Advection2D

  SUBROUTINE Flux_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    DO iEl = 1, this % source % nElem
      DO iVar = 1, this % source % nVar
        DO j = 0, this % source % interp % N
          DO i = 0, this % source % interp % N

            this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                  this % velocity % interior % hostData(1,i,j,iVar,iEl)*&
                  this % solution % interior % hostData(i,j,iVar,iEl)

            this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                  this % velocity % interior % hostData(2,i,j,iVar,iEl)*&
                  this % solution % interior % hostData(i,j,iVar,iEl)

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE Flux_Advection2D

  SUBROUTINE RiemannSolver_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl,iVar
    REAL(prec) :: extState, intState
    REAL(prec) :: nhat(1:2), nmag, un


    DO iEl = 1, this % solution % nElem
      DO iSide = 1, 4
        DO iVar = 1, this % solution % nVar
          DO i = 0, this % solution % interp % N

             ! Get the boundary normals on cell edges from the mesh geometry
             nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)

             ! Calculate the normal velocity at the cell edges
             un = this % velocity % boundary % hostData(1,i,1,iSide,iEl)*nHat(1)+&
                  this % velocity % boundary % hostData(2,i,1,iSide,iEl)*nHat(2)

             ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
             extState = this % solution % extBoundary % hostData(i,iVar,iSide,iEl)
             intState = this % solution % boundary % hostData(i,iVar,iSide,iEl)
             nmag = this % geometry % nScale % boundary % hostData(i,1,iSide,iEl)

             ! Calculate the flux
             this % flux % boundaryNormal % hostData(i,iVar,iSide,iEl) = 0.5_prec*&
                 ( un*(intState + extState) - abs(un)*(extState - intState) )*nmag

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE RiemannSolver_Advection2D

END MODULE SELF_Advection2D
