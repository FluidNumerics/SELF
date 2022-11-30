!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_ShallowWater1D

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model1D


  TYPE,EXTENDS(Model1D) :: ShallowWater1D
    !! iVar = 1 ~> u velocity component
    !! iVar = 2 ~> free surface height
    REAL(prec) :: g     ! gravity ( m/s^2) 
    TYPE(MappedScalar1D) :: H ! bottom topography ( m )
    TYPE(MappedScalar1D) :: gradH ! bottom topography gradient ( m/m )


    CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_ShallowWater1D
    PROCEDURE :: CalculateEntropy => CalculateEntropy_ShallowWater1D

    PROCEDURE :: PreTendency => PreTendency_ShallowWater1D
!    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_ShallowWater1D
    PROCEDURE :: ProductRuleCorrection => ProductRuleCorrection_ShallowWater1D

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_ShallowWater1D
    PROCEDURE :: FluxMethod => Flux_ShallowWater1D
    PROCEDURE :: RiemannSolver => RiemannSolver_ShallowWater1D
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_ShallowWater1D

    ! New Methods
    GENERIC :: SetTopography => SetTopographyFromChar_ShallowWater1D,&
                              SetTopographyFromEqn_ShallowWater1D
    PROCEDURE,PRIVATE :: SetTopographyFromChar_ShallowWater1D
    PROCEDURE,PRIVATE :: SetTopographyFromEqn_ShallowWater1D

    PROCEDURE :: SetLakeAtRest => SetLakeAtRest_ShallowWater1D

  END TYPE ShallowWater1D

!  INTERFACE
!    SUBROUTINE SetBoundaryCondition_ShallowWater1D_gpu_wrapper(solution, extBoundary, nHat, sideInfo, N, nVar, nEl) &
!      bind(c,name="SetBoundaryCondition_ShallowWater1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: solution, extBoundary, nHat, sideInfo
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!    END SUBROUTINE SetBoundaryCondition_ShallowWater1D_gpu_wrapper
!  END INTERFACE
!
  INTERFACE
    SUBROUTINE Source_ShallowWater1D_gpu_wrapper(source, solution, gradH, g, N, nVar, nEl) &
      bind(c,name="Source_ShallowWater1D_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: source, solution, gradH
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: g
    END SUBROUTINE Source_ShallowWater1D_gpu_wrapper
  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE Flux_ShallowWater1D_gpu_wrapper(flux, solution, g, H, N, nVar, nEl) &
!      bind(c,name="Flux_ShallowWater1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: flux, solution
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!      REAL(c_prec),VALUE :: g, H
!    END SUBROUTINE Flux_ShallowWater1D_gpu_wrapper
!  END INTERFACE
!
!  INTERFACE
!    SUBROUTINE RiemannSolver_ShallowWater1D_gpu_wrapper(flux, solution, extBoundary, nHat, nScale, g, H, N, nVar, nEl) &
!      bind(c,name="RiemannSolver_ShallowWater1D_gpu_wrapper")
!      USE iso_c_binding
!      USE SELF_Constants
!      IMPLICIT NONE
!      TYPE(c_ptr) :: flux, solution, extBoundary, nHat, nScale
!      INTEGER(C_INT),VALUE :: N,nVar,nEl
!      REAL(c_prec),VALUE :: g, H
!    END SUBROUTINE RiemannSolver_ShallowWater1D_gpu_wrapper
!  END INTERFACE

CONTAINS

  SUBROUTINE Init_ShallowWater1D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(out) :: this
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
    nvarloc = 2

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.
    this % g = 1.0_prec

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % H % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % gradH % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % prevSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % velocity % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nvarloc,this % mesh % nElem)

    ! First three variables are treated as u, v, eta
    ! Any additional are treated as passive tracers 
    CALL this % solution % SetName(1,"Hu")
    CALL this % solution % SetUnits(1,"m^2/s")
    CALL this % solution % SetDescription(1,"x-component of the volume flux")

    CALL this % solution % SetName(2,"H")
    CALL this % solution % SetUnits(2,"m")
    CALL this % solution % SetDescription(2,"Total fluid thickness")

  END SUBROUTINE Init_ShallowWater1D

  SUBROUTINE CalculateEntropy_ShallowWater1D(this)
  !! Base method for calculating entropy of a model
  !! Calculates the entropy as the integration of the 
  !! squared tracer over the domain
    IMPLICIT NONE
    CLASS(ShallowWater1D), INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, iVar, iEl
    REAL(prec) :: Jacobian, Hu, H, b
    REAL(prec) :: wi

    this % entropy = 0.0_prec

    DO iEl = 1, this % geometry % x % nElem
      DO i = 0, this % geometry % x % interp % N

        ! Coordinate mapping Jacobian
        Jacobian = this % geometry % dxds % interior % hostData(i,1,iEl)

        ! Quadrature weights
        wi = this % geometry % x % interp % qWeights % hostData(i) 

        ! Solution
        Hu = this % solution % interior % hostData(i,1,iEl)
        H = this % solution % interior % hostData(i,2,iEl)
        b = this % H % interior % hostData(i,1,iEl)

        this % entropy = this % entropy + ( 0.5_prec*Hu*Hu/H + &
                                            0.5_prec*this % g*H*H - this % g*H*b )*Jacobian*wi

      ENDDO
    ENDDO

  END SUBROUTINE CalculateEntropy_ShallowWater1D

  SUBROUTINE PreTendency_ShallowWater1D(this)
    !! Calculate the velocity at element interior and element boundaries
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    !
    INTEGER :: iEl, i
    REAL(prec) :: H


      DO iEl = 1, this % solution % nElem
        DO i = 0, this % solution % interp % N

          H = this % solution % interior % hostData(i,2,iEl)

          this % velocity % interior % hostData(i,1,iEl) = &
                  this % solution % interior % hostData(i,1,iEl)/H 

        ENDDO 
      ENDDO 

      CALL this % velocity % BoundaryInterp( gpuAccel = this % gpuAccel )
      CALL this % velocity % SideExchange(this % mesh, this % decomp, this % gpuAccel)

  END SUBROUTINE PreTendency_ShallowWater1D


  SUBROUTINE SetLakeAtRest_ShallowWater1D(this) 
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, i


      DO iEl = 1, this % solution % nElem
        DO i = 0, this % solution % interp % N
          this % solution % interior % hostData(i,1,iEl) = 0.0_prec
          this % solution % interior % hostData(i,2,iEl) = this % H % interior % hostData(i,1,iEl)
        ENDDO
      ENDDO


  END SUBROUTINE SetLakeAtRest_ShallowWater1D

  SUBROUTINE SetTopographyFromChar_ShallowWater1D(this, eqnChar) 
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: eqnChar

      CALL this % H % SetEquation(1, eqnChar)

      CALL this % H % SetInteriorFromEquation( this % geometry, this % t )

      CALL this % H % BoundaryInterp( gpuAccel = .FALSE. )

      CALL this % H % Derivative( this % geometry, &
                                this % gradH, &
                                selfStrongForm, &
                                .FALSE.)

      IF( this % gpuAccel )THEN
        CALL this % H % UpdateDevice()
      ENDIF

  END SUBROUTINE SetTopographyFromChar_ShallowWater1D

  SUBROUTINE SetTopographyFromEqn_ShallowWater1D(this, eqn) 
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn

      ! Copy the equation parser
      CALL this % H % SetEquation(1, eqn % equation)

      CALL this % H % SetInteriorFromEquation( this % geometry, this % t )

      CALL this % H % BoundaryInterp( gpuAccel = .FALSE. )

      CALL this % H % Derivative( this % geometry, &
                                this % gradH, &
                                selfStrongForm, &
                                .FALSE.)

      IF( this % gpuAccel )THEN
        CALL this % H % UpdateDevice()
      ENDIF

  END SUBROUTINE SetTopographyFromEqn_ShallowWater1D 

  SUBROUTINE SetBoundaryCondition_ShallowWater1D(this)
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, i
    INTEGER :: bcid, e2
    REAL(prec) :: u, v, nhat(1:2)


     ! Only do radiation condition (for now)
     ! Left most boundary
     this % solution % extBoundary % hostData(1,1,1) = 0.0_prec
     this % solution % extBoundary % hostData(2,1,1) = this % H % boundary % hostData(1,1,1)

     ! Right most boundary
     this % solution % extBoundary % hostData(1,2,this % solution % nElem) = 0.0_prec
     this % solution % extBoundary % hostData(2,2,this % solution % nElem) = &
             this % H % boundary % hostData(1,2,this % solution % nElem)


  END SUBROUTINE SetBoundaryCondition_ShallowWater1D 

  SUBROUTINE Source_ShallowWater1D(this)
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      IF( this % gpuAccel )THEN

       CALL Source_ShallowWater1D_gpu_wrapper(this % source % interior % deviceData, &
                                              this % solution % interior % deviceData, &
                                              this % gradH % interior % deviceData, &
                                              this % g, &
                                              this % source % interp % N, &
                                              this % solution % nVar, &
                                              this % solution % nElem )

      ELSE

        DO iEl = 1, this % source % nElem
          DO i = 0, this % source % interp % N

            this % source % interior % hostData(i,1,iEl) = &
                    this % g*this % solution % interior % hostData(i,2,iEl)*&
                    this % gradH % interior % hostData(i,1,iEl)

            this % source % interior % hostData(i,2,iEl) = 0.0_prec

          ENDDO
        ENDDO

        !CALL this % ProductRuleCorrection()

      ENDIF

  END SUBROUTINE Source_ShallowWater1D

  SUBROUTINE ProductRuleCorrection_ShallowWater1D(this)
    !!
    !! For this class, we use the PostfluxDivergence to add
    !! split form correction terms to the momentum equation
    !!
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    INTEGER :: iEl, i, ii
    REAL(prec) :: Dv, Dhv2, Dhv, Dh, Dh2, dmat
    REAL(prec) :: hv, h, v
    REAL(prec) :: hvi, hi, vi


      DO iEl = 1, this % solution % nElem
        DO i = 0, this % solution % interp % N

          Dv = 0.0_prec
          Dhv2 = 0.0_prec
          Dhv = 0.0_prec
          Dh2 = 0.0_prec

          DO ii = 0, this % solution % interp % N 
            hvi = this % solution % interior % hostData(ii,1,iEl)
            hi = this % solution % interior % hostData(ii,2,iEl)
            vi = this % velocity % interior % hostData(ii,1,iEl)
            dmat = this % solution % interp % dMatrix % hostData(ii,i)

            Dv = Dv + dmat*vi
            Dhv2 = Dhv2 + dmat*hvi*vi
            Dhv = Dhv + dmat*hvi
            Dh2 = Dh2 + dmat*hi*hi
            Dh = Dh + dmat*hi
          ENDDO

          hv = this % solution % interior % hostData(i,1,iEl)
          h = this % solution % interior % hostData(i,2,iEl)
          v = this % velocity % interior % hostData(i,1,iEl)

          this % source % interior % hostData(i,1,iEl) = &
                  this % source % interior % hostData(i,1,iEl) - &
                             ( 0.5_prec*this % g*(2.0_prec*h*Dh - Dh2) ) / &
                  this % geometry % dxds % interior % hostData(i,1,iEl)
         
        ENDDO
      ENDDO

  END SUBROUTINE ProductRuleCorrection_ShallowWater1D

  SUBROUTINE Flux_ShallowWater1D(this)
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO i = 0, this % solution % interp % N

            IF ( iVar == 1 )THEN ! Hu

              ! Hu^2 + (gH^2)/2
              this % flux % interior % hostData(i,iVar,iEl) = &
                    this % velocity % interior % hostData(i,1,iEl)*& ! u
                    this % solution % interior % hostData(i,1,iEl)+& ! Hu
                    0.5_prec*this % g*this % solution % interior % hostData(i,2,iEl)*&
                    this % solution % interior % hostData(i,2,iEl)


            ELSEIF ( iVar == 2 )THEN ! H - total fluid thickness
              this % flux % interior % hostData(i,iVar,iEl) = &
                    this % solution % interior % hostData(i,1,iEl)

            ENDIF

          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE Flux_ShallowWater1D

  SUBROUTINE RiemannSolver_ShallowWater1D(this)
    IMPLICIT NONE
    CLASS(ShallowWater1D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: cL, cR, unL, unR, HL, HR, HuL, HuR, HvL, HvR
    REAL(prec) :: alpha, nhat
    REAL(prec) :: uAvg, hAvg, h2Avg
    REAL(prec) :: fluxL(1:2)
    REAL(prec) :: fluxR(1:2)

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

           HL = this % solution % boundary % hostData(2,iSide,iEl)
           HR = this % solution % extBoundary % hostData(2,iSide,iEl)

           uAvg = (unL+unR)*0.5_prec
           hAvg = (hL+hR)*0.5_prec
           h2Avg = (hL*hL + hR*hR)*0.5_prec

           !cL = sqrt( this % g * this % solution % boundary % hostData(2,iSide,iEl))
           !cR = sqrt( this % g * this % solution % extBoundary % hostData(2,iSide,iEl))

           ! Entropy conserving flux from Gassner et al. (2016)
           this % flux % boundary % hostData(1,iSide,iEl) = ( &
                   uAvg*uAvg*hAvg + 0.5_prec*this % g*h2Avg )*nHat

           this % flux % boundary % hostData(2,iSide,iEl) = ( uAvg*hAvg )*nHat

           !cL = sqrt( this % g * this % solution % boundary % hostData(2,iSide,iEl))
           !cR = sqrt( this % g * this % solution % extBoundary % hostData(2,iSide,iEl))
!           fluxL(1) = (HuL*unL + this % g*HL*HL*0.5_prec)*nHat
!           fluxL(2) = HuL*nHat
!
!           fluxR(1) = (HuR*unR + this % g*HR*HR*0.5_prec)*nHat
!           fluxR(2) = HuR*nHat
!
!           alpha = MAX( ABS(unL+cL), ABS(unR+cR), &
!                        ABS(unL-cL), ABS(unR-cR) )

           ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
!           this % flux % boundary % hostData(1,iSide,iEl) = 0.5_prec*( &
!                   fluxL(1) + fluxR(1) + alpha*( HuL - HuR ) )
!           this % flux % boundary % hostData(2,iSide,iEl) = 0.5_prec*( &
!                   fluxL(2) + fluxR(2) + alpha*( HL - HR ) )

        ENDDO
      ENDDO

!    ENDIF

  END SUBROUTINE RiemannSolver_ShallowWater1D

END MODULE SELF_ShallowWater1D
