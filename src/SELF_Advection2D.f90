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
  USE SELF_Model


  TYPE,EXTENDS(Model2D) :: Advection2D
    TYPE(MappedVector2D) :: velocity

    CONTAINS

    ! Extended methods to handle velocity attribute
    PROCEDURE :: Init => Init_Advection2D
    PROCEDURE :: Free => Free_Advection2D

    ! Concretized Methods
    PROCEDURE :: Source2D => Source_Advection2D
    PROCEDURE :: Flux2D => Flux_Advection2D
    PROCEDURE :: RiemannSolver2D => RiemannSolver_Advection2D

  END TYPE Advection2D

CONTAINS

  SUBROUTINE Init_Advection2D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh2D),INTENT(in),TARGET :: mesh
    TYPE(SEMQuad),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp
    
    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry

    CALL this % velocity % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % solution % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nVar,this % mesh % nElem)

  END SUBROUTINE Init_Advection2D

  SUBROUTINE Free_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this

    CALL this % velocity % Free()
    CALL this % solution % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_Advection2D

  SUBROUTINE Source_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    DO iEl = 1, this % source % nElem
      DO iVar = 1, this % source % nElem
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
      DO iVar = 1, this % source % nElem
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
