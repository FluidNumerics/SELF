!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_ShallowWater

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model


  TYPE,EXTENDS(Model2D) :: ShallowWater
    !! iVar = 1 ~> u velocity component
    !! iVar = 2 ~> v velocity component
    !! iVar = 3 ~> free surface height
    REAL(prec) :: fCori ! coriolis parameter ( 1/s )
    REAL(prec) :: g     ! gravity ( m/s^2) 
    TYPE(MappedScalar2D) :: H ! bottom topography ( m )
    TYPE(MappedVector2D) :: gradH ! bottom topography gradient ( m/m )


    CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_ShallowWater
    PROCEDURE :: CalculateEntropy => CalculateEntropy_ShallowWater

    PROCEDURE :: PreTendency => PreTendency_ShallowWater

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_ShallowWater
    PROCEDURE :: FluxMethod => Flux_ShallowWater
    PROCEDURE :: RiemannSolver => RiemannSolver_ShallowWater
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_ShallowWater

    ! New Methods
    GENERIC :: SetTopography => SetTopographyFromChar_ShallowWater,&
                              SetTopographyFromEqn_ShallowWater
    PROCEDURE,PRIVATE :: SetTopographyFromChar_ShallowWater
    PROCEDURE,PRIVATE :: SetTopographyFromEqn_ShallowWater

  END TYPE ShallowWater

  INTERFACE
    SUBROUTINE SetBoundaryCondition_ShallowWater_gpu_wrapper(solution, extBoundary, nHat, sideInfo, N, nVar, nEl) &
      bind(c,name="SetBoundaryCondition_ShallowWater_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: solution, extBoundary, nHat, sideInfo
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SetBoundaryCondition_ShallowWater_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Source_ShallowWater_gpu_wrapper(source, solution, f, N, nVar, nEl) &
      bind(c,name="Source_ShallowWater_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: source, solution
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: f
    END SUBROUTINE Source_ShallowWater_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE Flux_ShallowWater_gpu_wrapper(flux, solution, g, H, N, nVar, nEl) &
      bind(c,name="Flux_ShallowWater_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, solution
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: g, H
    END SUBROUTINE Flux_ShallowWater_gpu_wrapper
  END INTERFACE

  INTERFACE
    SUBROUTINE RiemannSolver_ShallowWater_gpu_wrapper(flux, solution, extBoundary, nHat, nScale, g, H, N, nVar, nEl) &
      bind(c,name="RiemannSolver_ShallowWater_gpu_wrapper")
      USE iso_c_binding
      USE SELF_Constants
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, solution, extBoundary, nHat, nScale
      INTEGER(C_INT),VALUE :: N,nVar,nEl
      REAL(c_prec),VALUE :: g, H
    END SUBROUTINE RiemannSolver_ShallowWater_gpu_wrapper
  END INTERFACE

CONTAINS

  SUBROUTINE Init_ShallowWater(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(ShallowWater),INTENT(out) :: this
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
    this % fCori = 0.0_prec

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % H % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % gradH % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % workSol % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % velocity % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % compVelocity % Init(geometry % x % interp,1,this % mesh % nElem)
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

    CALL this % solution % SetName(2,"Hv")
    CALL this % solution % SetUnits(2,"m^2/s")
    CALL this % solution % SetDescription(2,"y-component of the volume flux")

    CALL this % solution % SetName(3,"H")
    CALL this % solution % SetUnits(3,"m")
    CALL this % solution % SetDescription(3,"Total fluid thickness")

  END SUBROUTINE Init_ShallowWater

  SUBROUTINE CalculateEntropy_ShallowWater(this)
  !! Base method for calculating entropy of a model
  !! Calculates the entropy as the integration of the 
  !! squared tracer over the domain
    IMPLICIT NONE
    CLASS(ShallowWater), INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, iVar, iEl
    REAL(prec) :: Jacobian, Hu, Hv, H, b
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
            Hu = this % solution % interior % hostData(i,j,1,iEl)
            Hv = this % solution % interior % hostData(i,j,2,iEl)
            H = this % solution % interior % hostData(i,j,3,iEl)
            b = this % H % interior % hostData(i,j,1,iEl)

            this % entropy = this % entropy + ( 0.5_prec*( Hu*Hu/H + Hv*Hv/H ) + &
                                                0.5_prec*this % g*H*H - this % g*H*b )*Jacobian*wi*wj

          ENDDO
        ENDDO
    ENDDO

  END SUBROUTINE CalculateEntropy_ShallowWater

  SUBROUTINE PreTendency_ShallowWater(this)
    !! Calculate the velocity at element interior and element boundaries
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    IMPLICIT NONE
    CLASS(ShallowWater),INTENT(inout) :: this
    !
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: H


      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            H = this % solution % interior % hostData(i,j,3,iEl)

            this % velocity % interior % hostData(1,i,j,1,iEl) = &
                    this % solution % interior % hostData(i,j,1,iEl)/H 

            this % velocity % interior % hostData(2,i,j,1,iEl) = &
                    this % solution % interior % hostData(i,j,2,iEl)/H 

          ENDDO 
        ENDDO 
      ENDDO 

      DO iEl = 1, this % solution % nElem
        DO iSide = 1, 4
          DO i = 0, this % solution % interp % N

            H = this % solution % boundary % hostData(i,3,iSide,iEl)

            this % velocity % boundary % hostData(1,i,1,iSide,iEl) = &
                    this % solution % boundary % hostData(i,1,iSide,iEl)/H 

            this % velocity % boundary % hostData(2,i,1,iSide,iEl) = &
                    this % solution % boundary % hostData(i,2,iSide,iEl)/H 
          ENDDO
        ENDDO
      ENDDO

      CALL this % velocity % SideExchange(this % mesh, this % decomp, this % gpuAccel)

  END SUBROUTINE PreTendency_ShallowWater

  SUBROUTINE SetTopographyFromChar_ShallowWater(this, eqnChar) 
    IMPLICIT NONE
    CLASS(ShallowWater),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: eqnChar

      CALL this % H % SetEquation(1, eqnChar)

      CALL this % H % SetInteriorFromEquation( this % geometry, this % t )

      CALL this % H % BoundaryInterp( gpuAccel = .FALSE. )

      CALL this % H % Gradient( this % geometry, &
                                this % gradH, &
                                selfStrongForm, &
                                .FALSE.)

      IF( this % gpuAccel )THEN
        CALL this % H % UpdateDevice()
      ENDIF

  END SUBROUTINE SetTopographyFromChar_ShallowWater

  SUBROUTINE SetTopographyFromEqn_ShallowWater(this, eqn) 
    IMPLICIT NONE
    CLASS(ShallowWater),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn

      ! Copy the equation parser
      CALL this % H % SetEquation(1, eqn % equation)

      CALL this % H % SetInteriorFromEquation( this % geometry, this % t )

      CALL this % H % BoundaryInterp( gpuAccel = .FALSE. )

      CALL this % H % Gradient( this % geometry, &
                                this % gradH, &
                                selfStrongForm, &
                                .FALSE.)

      IF( this % gpuAccel )THEN
        CALL this % H % UpdateDevice()
      ENDIF

  END SUBROUTINE SetTopographyFromEqn_ShallowWater 

  SUBROUTINE SetBoundaryCondition_ShallowWater(this)
    IMPLICIT NONE
    CLASS(ShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, i
    INTEGER :: bcid, e2
    REAL(prec) :: u, v, nhat(1:2)


!    IF( this % gpuAccel )THEN
!
!      CALL SetBoundaryCondition_ShallowWater_gpu_wrapper( this % solution % boundary % deviceData,&
!            this % solution % extBoundary % deviceData, &
!            this % geometry % nHat % boundary % deviceData, &
!            this % mesh % sideInfo % deviceData, &
!            this % solution % interp % N, &
!            this % solution % nVar, &
!            this % solution % nElem)
!
!    ELSE

      DO iEl = 1, this % solution % nElem
        DO iSide = 1, 4
            DO i = 0, this % solution % interp % N

              bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
              e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
              IF( e2 == 0 )THEN
                IF( bcid == SELF_BC_RADIATION )THEN

                  this % solution % extBoundary % hostData(i,1,iSide,iEl) = 0.0_prec
                  this % solution % extBoundary % hostData(i,2,iSide,iEl) = 0.0_prec
                  this % solution % extBoundary % hostData(i,3,iSide,iEl) = this % H % boundary % hostData(i,1,iSide,iEl)
                  !PRINT*, "Radiation Boundary, H :",this % solution % extBoundary % hostData(i,3,iSide,iEl) 

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
                  this % solution % extBoundary % hostData(i,3,iSide,iEl) = this % H % boundary % hostData(i,1,iSide,iEl)

                ENDIF
              ENDIF

            ENDDO
        ENDDO
      ENDDO

!    ENDIF

  END SUBROUTINE SetBoundaryCondition_ShallowWater 

  SUBROUTINE Source_ShallowWater(this)
    IMPLICIT NONE
    CLASS(ShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

!    IF( this % gpuAccel )THEN
!
!      CALL Source_ShallowWater_gpu_wrapper(this % source % interior % deviceData,&
!             this % solution % interior % deviceData, &
!             this % fCori,&
!             this % source % interp % N, &
!             this % solution % nVar, &
!             this % solution % nElem )
!
!    ELSE

      DO iEl = 1, this % source % nElem
        DO j = 0, this % source % interp % N
          DO i = 0, this % source % interp % N

            this % source % interior % hostData(i,j,1,iEl) = &
                    -this % fCori*this % solution % interior % hostData(i,j,2,iEl) +&
                    this % g*this % solution % interior % hostData(i,j,3,iEl)*&
                    this % gradH % interior % hostData(1,i,j,1,iEl)

            this % source % interior % hostData(i,j,2,iEl) = &
                    this % fCori*this % solution % interior % hostData(i,j,1,iEl) +&
                    this % g*this % solution % interior % hostData(i,j,3,iEl)*&
                    this % gradH % interior % hostData(2,i,j,1,iEl)

            this % source % interior % hostData(i,j,3,iEl) = 0.0_prec

          ENDDO
        ENDDO
      ENDDO

!    ENDIF

  END SUBROUTINE Source_ShallowWater

  SUBROUTINE Flux_ShallowWater(this)
    IMPLICIT NONE
    CLASS(ShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

!    IF( this % gpuAccel )THEN
!
!      CALL Flux_ShallowWater_gpu_wrapper(this % flux % interior % deviceData,&
!                                               this % solution % interior % deviceData, &
!                                               this % g, this % H, this % solution % interp % N, &
!                                               this % solution % nVar, this % solution % nElem)
!
!    ELSE
      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              IF ( iVar == 1 )THEN ! Hu

                ! Hu^2 + (gH^2)/2
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(1,i,j,1,iEl)*& ! u
                      this % solution % interior % hostData(i,j,1,iEl)+& ! Hu
                      0.5_prec*this % g*this % solution % interior % hostData(i,j,3,iEl)*&
                      this % solution % interior % hostData(i,j,3,iEl)

                ! Huv
                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(2,i,j,1,iEl)*& ! v
                      this % solution % interior % hostData(i,j,1,iEl) ! Hu

              ELSEIF ( iVar == 2 )THEN !Hv

                ! Huv
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(1,i,j,1,iEl)*& ! u
                      this % solution % interior % hostData(i,j,2,iEl) ! Hv

                ! (Hv^2 + (gH^2)/2)
                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(2,i,j,1,iEl)*& ! v
                      this % solution % interior % hostData(i,j,2,iEl)+& ! Hv
                      0.5_prec*this % g*this % solution % interior % hostData(i,j,3,iEl)*&
                      this % solution % interior % hostData(i,j,3,iEl)



              ELSEIF ( iVar == 3 )THEN ! H - total fluid thickness
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % solution % interior % hostData(i,j,1,iEl)

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % solution % interior % hostData(i,j,2,iEl)
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO
!    ENDIF

  END SUBROUTINE Flux_ShallowWater

  SUBROUTINE RiemannSolver_ShallowWater(this)
    IMPLICIT NONE
    CLASS(ShallowWater),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: nhat(1:2), nmag
    REAL(prec) :: cL, cR, unL, unR, HL, HR, HuL, HuR, HvL, HvR
    REAL(prec) :: alpha
    REAL(prec) :: fluxL(1:3)
    REAL(prec) :: fluxR(1:3)


!    IF( this % gpuAccel )THEN
!
!      CALL RiemannSolver_ShallowWater_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
!             this % solution % boundary % deviceData, &
!             this % solution % extBoundary % deviceData, &
!             this % geometry % nHat % boundary % deviceData, &
!             this % geometry % nScale % boundary % deviceData, &
!             this % g, this % H, &
!             this % solution % interp % N, &
!             this % solution % nVar, &
!             this % solution % nElem)
!
!    ELSE

      DO iEl = 1, this % solution % nElem
        DO iSide = 1, 4
          DO i = 0, this % solution % interp % N

             ! Get the boundary normals on cell edges from the mesh geometry
             nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
             nmag = this % geometry % nScale % boundary % hostData(i,1,iSide,iEl)

             ! Calculate the normal velocity at the cell edges
             unL = this % velocity % boundary % hostData(1,i,1,iSide,iEl)*nHat(1)+&
                   this % velocity % boundary % hostData(2,i,1,iSide,iEl)*nHat(2)

             unR = this % velocity % extBoundary % hostData(1,i,1,iSide,iEl)*nHat(1)+&
                   this % velocity % extBoundary % hostData(2,i,1,iSide,iEl)*nHat(2)

             cL = sqrt( this % g * this % solution % boundary % hostData(i,3,iSide,iEl))
             cR = sqrt( this % g * this % solution % extBoundary % hostData(i,3,iSide,iEl))

             HuL = this % solution % boundary % hostData(i,1,iSide,iEl)
             HuR = this % solution % extBoundary % hostData(i,1,iSide,iEl)

             HvL = this % solution % boundary % hostData(i,2,iSide,iEl)
             HvR = this % solution % extBoundary % hostData(i,2,iSide,iEl)

             HL = this % solution % boundary % hostData(i,3,iSide,iEl)
             HR = this % solution % extBoundary % hostData(i,3,iSide,iEl)

             fluxL(1) = HuL*unL + this % g*HL*HL*0.5_prec*nHat(1)
             fluxL(2) = HvL*unL + this % g*HL*HL*0.5_prec*nHat(2)
             fluxL(3) = HuL*nHat(1) + HvL*nHat(2)

             fluxR(1) = HuR*unR + this % g*HR*HR*0.5_prec*nHat(1)
             fluxR(2) = HvR*unR + this % g*HR*HR*0.5_prec*nHat(2)
             fluxR(3) = HuR*nHat(1) + HvR*nHat(2)

             alpha = MAX( ABS(unL+cL), ABS(unR+cR), &
                          ABS(unL-cL), ABS(unR-cR) )

             ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
             this % flux % boundaryNormal % hostData(i,1,iSide,iEl) = 0.5_prec*( &
                     fluxL(1) + fluxR(1) + alpha*( HuL - HuR ) )*nmag
             this % flux % boundaryNormal % hostData(i,2,iSide,iEl) = 0.5_prec*( &
                     fluxL(2) + fluxR(2) + alpha*( HvL - HvR ) )*nmag
             this % flux % boundaryNormal % hostData(i,3,iSide,iEl) = 0.5_prec*( &
                     fluxL(3) + fluxR(3) + alpha*( HL - HR ) )*nmag

          ENDDO
        ENDDO
      ENDDO

!    ENDIF

  END SUBROUTINE RiemannSolver_ShallowWater

END MODULE SELF_ShallowWater
