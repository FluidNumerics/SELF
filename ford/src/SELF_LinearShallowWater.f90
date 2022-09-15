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

    ! Concretized Methods
    PROCEDURE :: Source2D => Source_LinearShallowWater
    PROCEDURE :: Flux2D => Flux_LinearShallowWater
    PROCEDURE :: RiemannSolver2D => RiemannSolver_LinearShallowWater
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_LinearShallowWater

  END TYPE LinearShallowWater

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

CONTAINS

  SUBROUTINE SetBoundaryCondition_LinearShallowWater(this)
    IMPLICIT NONE
    CLASS(LinearShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, i
    INTEGER :: bcid
    REAL(prec) :: u, v, nhat(1:2)


    DO iEl = 1, this % solution % nElem
      DO iSide = 1, 4
          DO i = 0, this % solution % interp % N

            bcid = this % mesh % self_sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
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

!            ELSEIF( bcid == SELF_BC_PRESCRIBED )THEN
            ELSE

              PRINT*, "WARNING : Unrecognized boundary condition"

            ENDIF
          ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE SetBoundaryCondition_LinearShallowWater 

  SUBROUTINE Source_LinearShallowWater(this)
    IMPLICIT NONE
    CLASS(LinearShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

    DO iEl = 1, this % source % nElem
      DO j = 0, this % source % interp % N
        DO i = 0, this % source % interp % N

          this % source % interior % hostData(i,j,1,iEl) = -this % fCori*this % source % interior % hostData(i,j,2,iEl)
          this % source % interior % hostData(i,j,2,iEl) = this % fCori*this % source % interior % hostData(i,j,1,iEl)
          this % source % interior % hostData(i,j,3,iEl) = 0.0_prec

        ENDDO
      ENDDO
    ENDDO

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
    INTEGER :: i,iSide,iEl,iVar
    REAL(prec) :: nhat(1:2), nmag
    REAL(prec) :: c, unL, unR, etaL, etaR, wL, wR



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
             wL = 0.5_prec*(etaL/c + unL/this % g)
             wR = 0.5_prec*(etaR/c - unR/this % g)

             this % flux % boundaryNormal % hostData(i,1,iSide,iEl) = this % g*c*( wL - wR )*nHat(1)*nmag
             this % flux % boundaryNormal % hostData(i,2,iSide,iEl) = this % g*c*( wL - wR )*nHat(2)*nmag
             this % flux % boundaryNormal % hostData(i,1,iSide,iEl) = c*c*( wL + wR )*nmag

          ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE RiemannSolver_LinearShallowWater

END MODULE SELF_LinearShallowWater
