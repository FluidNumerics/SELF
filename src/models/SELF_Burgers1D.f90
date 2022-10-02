!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Burgers1D

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model1D


  TYPE,EXTENDS(Model1D) :: Burgers1D

    CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_Burgers1D
    PROCEDURE :: CalculateEntropy => CalculateEntropy_Burgers1D

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_Burgers1D
    PROCEDURE :: FluxMethod => Flux_Burgers1D
    PROCEDURE :: RiemannSolver => RiemannSolver_Burgers1D
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_Burgers1D

  END TYPE Burgers1D

!  INTERFACE
!    SUBROUTINE SetBoundaryCondition_Burgers1D_gpu_wrapper(solution, extBoundary, nHat, sideInfo, N, nVar, nEl) &
!      bind(c,name="SetBoundaryCondition_Burgers1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: solution, extBoundary, nHat, sideInfo
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!    END SUBROUTINE SetBoundaryCondition_Burgers1D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE Source_Burgers1D_gpu_wrapper(source, solution, gradH, g, N, nVar, nEl) &
!      bind(c,name="Source_Burgers1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: source, solution, gradH
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!      REAL(c_prec),VALUE :: g
!    END SUBROUTINE Source_Burgers1D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE Flux_Burgers1D_gpu_wrapper(flux, solution, g, H, N, nVar, nEl) &
!      bind(c,name="Flux_Burgers1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: flux, solution
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!      REAL(c_prec),VALUE :: g, H
!    END SUBROUTINE Flux_Burgers1D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE RiemannSolver_Burgers1D_gpu_wrapper(flux, solution, extBoundary, nHat, nScale, g, H, N, nVar, nEl) &
!      bind(c,name="RiemannSolver_Burgers1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: flux, solution, extBoundary, nHat, nScale
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!      REAL(c_prec),VALUE :: g, H
!    END SUBROUTINE RiemannSolver_Burgers1D_gpu_wrapper
!  END INTERFACE

CONTAINS

  SUBROUTINE Init_Burgers1D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh1D),INTENT(in),TARGET :: mesh
    TYPE(Geometry1D),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp
    ! Local
    INTEGER :: ivar
    CHARACTER(LEN=3) :: ivarChar
    CHARACTER(LEN=25) :: varname
    INTEGER :: nvarloc

    ! Ensure that the number of variables is 3
    ! nvar is unused in this class extension
    nvarloc = 1

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.
    this % fluxDivMethod = SELF_CONSERVATIVE_FLUX 

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % velocity % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    CALL this % solution % SetName(1,"u")
    CALL this % solution % SetUnits(1,"[]")
    CALL this % solution % SetDescription(1,"null")


  END SUBROUTINE Init_Burgers1D

  SUBROUTINE CalculateEntropy_Burgers1D(this)
  !! Base method for calculating entropy of a model
  !! Calculates the entropy as the integration of the 
  !! squared tracer over the domain
    IMPLICIT NONE
    CLASS(Burgers1D), INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, iVar, iEl
    REAL(prec) :: Jacobian, u
    REAL(prec) :: wi

    this % entropy = 0.0_prec

    DO iEl = 1, this % geometry % x % nElem
      DO i = 0, this % geometry % x % interp % N

        ! Coordinate mapping Jacobian
        Jacobian = this % geometry % dxds % interior % hostData(i,1,iEl)

        ! Quadrature weights
        wi = this % geometry % x % interp % qWeights % hostData(i) 

        ! Solution
        u = this % solution % interior % hostData(i,1,iEl)

        this % entropy = this % entropy + ( 0.5_prec*u*u )

      ENDDO
    ENDDO

  END SUBROUTINE CalculateEntropy_Burgers1D

  SUBROUTINE SetBoundaryCondition_Burgers1D(this)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, i
    INTEGER :: bcid, e2
    REAL(prec) :: u, v, nhat(1:2)

     ! Left most boundary
     this % solution % extBoundary % hostData(1,1,1) = 1.0_prec

     ! Right most boundary
     this % solution % extBoundary % hostData(1,2,this % solution % nElem) = -1.0_prec


  END SUBROUTINE SetBoundaryCondition_Burgers1D 

  SUBROUTINE Source_Burgers1D(this)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      !IF( this % gpuAccel )THEN

      ! CALL Source_Burgers1D_gpu_wrapper(this % source % interior % deviceData, &
      !                                        this % solution % interior % deviceData, &
      !                                        this % gradH % interior % deviceData, &
      !                                        this % g, &
      !                                        this % source % interp % N, &
      !                                        this % solution % nVar, &
      !                                        this % solution % nElem )

      !ELSE
        DO iEl = 1, this % source % nElem
          DO i = 0, this % source % interp % N

            this % source % interior % hostData(i,1,iEl) = 0.0_prec

          ENDDO
        ENDDO
      !ENDIF

  END SUBROUTINE Source_Burgers1D

  SUBROUTINE Flux_Burgers1D(this)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl

      DO iEl = 1, this % solution % nElem
        DO i = 0, this % solution % interp % N

          ! u^2
          this % flux % interior % hostData(i,1,iEl) = &
                0.5_prec*this % solution % interior % hostData(i,1,iEl)*&
                this % solution % interior % hostData(i,1,iEl)


        ENDDO
      ENDDO

  END SUBROUTINE Flux_Burgers1D

  SUBROUTINE RiemannSolver_Burgers1D(this)
    IMPLICIT NONE
    CLASS(Burgers1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: unL, unR
    REAL(prec) :: alpha, nhat
    REAL(prec) :: fluxL
    REAL(prec) :: fluxR

      DO iEl = 1, this % solution % nElem
        DO iSide = 1, 2

           ! Get the boundary normals on cell edges from the mesh geometry
           IF( iSide ==  1)THEN
             nhat = -1.0_prec
           ELSE
             nhat = 1.0_prec
           ENDIF

           ! Calculate the normal velocity at the cell edges
           unL = this % velocity % boundary % hostData(1,iSide,iEl)
           unR = this % velocity % extBoundary % hostData(1,iSide,iEl)

           fluxL = 0.5_preC*(unL*unL)*nHat

           fluxR = 0.5_prec*(unR*unR)*nHat

           alpha = MAX( ABS(unL), ABS(unR) )

           ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
           this % flux % boundary % hostData(1,iSide,iEl) = 0.5_prec*( &
                   fluxL + fluxR + alpha*( unL-unR ) )

        ENDDO
      ENDDO

  END SUBROUTINE RiemannSolver_Burgers1D

END MODULE SELF_Burgers1D
