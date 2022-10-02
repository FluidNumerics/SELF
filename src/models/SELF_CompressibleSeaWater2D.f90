!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_CompressibleSeaWater2D

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model2D


  TYPE,EXTENDS(Model2D) :: CompressibleSeaWater2D
    !! iVar = 1 ~> rho*u (x-momentum)
    !! iVar = 2 ~> rho*v (y-momentum)
    !! iVar = 3 ~> rho (density)
    !! iVar = 4 ~> rho*E (Total Energy)
    !! iVar = 5 ~> rho*S (Salinity)
    TYPE(MappedScalar2D) :: kineticEnergy
    TYPE(MappedScalar2D) :: pressure
    TYPE(MappedScalar2D) :: enthalpy
    TYPE(MappedScalar2D) :: soundSpeed


    CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_CompressibleSeaWater2D
    PROCEDURE :: PreTendency => PreTendency_CompressibleSeaWater2D

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_CompressibleSeaWater2D
    PROCEDURE :: FluxMethod => Flux_CompressibleSeaWater2D
    PROCEDURE :: RiemannSolver => RiemannSolver_CompressibleSeaWater2D
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_CompressibleSeaWater2D

  END TYPE CompressibleSeaWater2D

 ! INTERFACE
 !   SUBROUTINE SetBoundaryCondition_CompressibleSeaWater2D_gpu_wrapper(solution, extBoundary, nHat, sideInfo, N, nVar, nEl) &
 !     bind(c,name="SetBoundaryCondition_CompressibleSeaWater2D_gpu_wrapper")
 !     USE iso_c_binding
 !     USE SELF_Constants
 !     IMPLICIT NONE
 !     TYPE(c_ptr) :: solution, extBoundary, nHat, sideInfo
 !     INTEGER(C_INT),VALUE :: N,nVar,nEl
 !   END SUBROUTINE SetBoundaryCondition_CompressibleSeaWater2D_gpu_wrapper
 ! END INTERFACE

 ! INTERFACE
 !   SUBROUTINE Source_CompressibleSeaWater2D_gpu_wrapper(source, solution, f, N, nVar, nEl) &
 !     bind(c,name="Source_CompressibleSeaWater2D_gpu_wrapper")
 !     USE iso_c_binding
 !     USE SELF_Constants
 !     IMPLICIT NONE
 !     TYPE(c_ptr) :: source, solution
 !     INTEGER(C_INT),VALUE :: N,nVar,nEl
 !     REAL(c_prec),VALUE :: f
 !   END SUBROUTINE Source_CompressibleSeaWater2D_gpu_wrapper
 ! END INTERFACE

 ! INTERFACE
 !   SUBROUTINE Flux_CompressibleSeaWater2D_gpu_wrapper(flux, solution, g, H, N, nVar, nEl) &
 !     bind(c,name="Flux_CompressibleSeaWater2D_gpu_wrapper")
 !     USE iso_c_binding
 !     USE SELF_Constants
 !     IMPLICIT NONE
 !     TYPE(c_ptr) :: flux, solution
 !     INTEGER(C_INT),VALUE :: N,nVar,nEl
 !     REAL(c_prec),VALUE :: g, H
 !   END SUBROUTINE Flux_CompressibleSeaWater2D_gpu_wrapper
 ! END INTERFACE

 ! INTERFACE
 !   SUBROUTINE RiemannSolver_CompressibleSeaWater2D_gpu_wrapper(flux, solution, extBoundary, nHat, nScale, g, H, N, nVar, nEl) &
 !     bind(c,name="RiemannSolver_CompressibleSeaWater2D_gpu_wrapper")
 !     USE iso_c_binding
 !     USE SELF_Constants
 !     IMPLICIT NONE
 !     TYPE(c_ptr) :: flux, solution, extBoundary, nHat, nScale
 !     INTEGER(C_INT),VALUE :: N,nVar,nEl
 !     REAL(c_prec),VALUE :: g, H
 !   END SUBROUTINE RiemannSolver_CompressibleSeaWater2D_gpu_wrapper
 ! END INTERFACE

CONTAINS

  SUBROUTINE Init_CompressibleSeaWater2D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(CompressibleSeaWater2D),INTENT(out) :: this
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
    nvarloc = 5

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.
    this % fluxDivMethod = SELF_CONSERVATIVE_FLUX 

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % kineticEnergy % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % pressure % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % enthalpy % Init(geometry % x % interp,1,this % mesh % nElem)
    CALL this % soundSpeed % Init(geometry % x % interp,1,this % mesh % nElem)
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
    CALL this % solution % SetName(1,"rho*u")
    CALL this % solution % SetUnits(1,"kg/m^2/s")
    CALL this % solution % SetDescription(1,"x-component of the momentum per unit volume")

    CALL this % solution % SetName(2,"rho*v")
    CALL this % solution % SetUnits(2,"kg/m^2/s")
    CALL this % solution % SetDescription(2,"y-component of the momentum per unit volume")

    CALL this % solution % SetName(3,"rho")
    CALL this % solution % SetUnits(3,"kg/m^3")
    CALL this % solution % SetDescription(3,"fluid density; mass per unit volume")

    CALL this % solution % SetName(4,"rho*E")
    CALL this % solution % SetUnits(4,"kg/m/s^2")
    CALL this % solution % SetDescription(4,"Density weighted total energy per unit volume")

    CALL this % solution % SetName(5,"rho*S")
    CALL this % solution % SetUnits(5,"kg/m^3 [kg/kg]")
    CALL this % solution % SetDescription(5,"Density weighted salinity per unit volume")

  END SUBROUTINE Init_CompressibleSeaWater2D

  SUBROUTINE EquationOfState_CompressibleSeaWater2D(this)
    !! Calculates the fluid pressure, given density, energy and salinity

    IMPLICIT NONE
    CLASS(CompressibleSeaWater2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            this % pressure % interior % hostData(i,j,1,iEl) = 

          ENDDO
        ENDDO
      ENDDO

      CALL this % pressure % BoundaryInterp(this % gpuAccel)
      CALL this % pressure % SideExchange(this % mesh, this % decomp, this % gpuAccel)

  END SUBROUTINE EquationOfState_CompressibleSeaWater2D

  SUBROUTINE CalculateSoundSpeed_CompressibleSeaWater2D(this)
    !! Calculates the speed of sound throughout the fluid

    IMPLICIT NONE
    CLASS(CompressibleSeaWater2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N



          ENDDO
        ENDDO
      ENDDO

      CALL this % soundSpeed % BoundaryInterp(this % gpuAccel)
      CALL this % soundSpeed % SideExchange(this % mesh, this % decomp, this % gpuAccel)

  END SUBROUTINE CalculateSoundSpeed_CompressibleSeaWater2D

  SUBROUTINE PreTendency_CompressibleSeaWater2D(this)
    !! Calculate the velocity and density weighted enthalpy at element interior and element boundaries
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    IMPLICIT NONE
    CLASS(CompressibleSeaWater2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rho


! TO DO : 
!      CALL this % CalculateKineticEnergy()
!      CALL this % EquationOfState()
!      CALL this % CalculateSoundSpeed()

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            rho = this % solution % interior % hostData(i,j,3,iEl)

            this % velocity % interior % hostData(1,i,j,1,iEl) = &
                    this % solution % interior % hostData(i,j,1,iEl)/rho 

            this % velocity % interior % hostData(2,i,j,1,iEl) = &
                    this % solution % interior % hostData(i,j,2,iEl)/rho

            this % enthalpy % interior % hostData(i,j,1,iEl) = &
                    this % solution % interior % hostData(i,j,4,iEl)+&
                    this % pressure % interior % hostData(i,j,1,iEl)

          ENDDO 
        ENDDO 
      ENDDO 



      CALL this % velocity % BoundaryInterp(this % gpuAccel)
      CALL this % velocity % SideExchange(this % mesh, this % decomp, this % gpuAccel)

      CALL this % enthalpy % BoundaryInterp(this % gpuAccel)
      CALL this % enthalpy % SideExchange(this % mesh, this % decomp, this % gpuAccel)

  END SUBROUTINE PreTendency_CompressibleSeaWater2D

  SUBROUTINE SetBoundaryCondition_CompressibleSeaWater2D(this)
    IMPLICIT NONE
    CLASS(CompressibleSeaWater2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, i
    INTEGER :: bcid, e2
    REAL(prec) :: u, v, nhat(1:2)


      DO iEl = 1, this % solution % nElem
        DO iSide = 1, 4
            DO i = 0, this % solution % interp % N

              bcid = this % mesh % sideInfo % hostData(5,iSide,iEl) ! Boundary Condition ID
              e2 = this % mesh % sideInfo % hostData(3,iSide,iEl) ! Neighboring Element ID
              IF( e2 == 0 )THEN
                IF( bcid == SELF_BC_NONORMALFLOW )THEN

                  nhat(1:2) = this % geometry % nHat % boundary % hostData(1:2,i,1,iSide,iEl)
                  u = this % solution % boundary % hostData(i,1,iSide,iEl) 
                  v = this % solution % boundary % hostData(i,2,iSide,iEl) 
                  this % solution % extBoundary % hostData(i,1,iSide,iEl) = (nhat(2)**2 - nhat(1)**2)*u - 2.0_prec*nhat(1)*nhat(2)*v
                  this % solution % extBoundary % hostData(i,2,iSide,iEl) = (nhat(1)**2 - nhat(2)**2)*v - 2.0_prec*nhat(1)*nhat(2)*u
                  this % solution % extBoundary % hostData(i,3,iSide,iEl) = this % solution % boundary % hostData(i,3,iSide,iEl)
                  this % solution % extBoundary % hostData(i,4,iSide,iEl) = this % solution % boundary % hostData(i,4,iSide,iEl)
                  this % solution % extBoundary % hostData(i,5,iSide,iEl) = this % solution % boundary % hostData(i,5,iSide,iEl)

                ENDIF
              ENDIF

            ENDDO
        ENDDO
      ENDDO

!    ENDIF

  END SUBROUTINE SetBoundaryCondition_CompressibleSeaWater2D 

  SUBROUTINE Source_CompressibleSeaWater2D(this)
    IMPLICIT NONE
    CLASS(CompressibleSeaWater2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      DO iEl = 1, this % source % nElem
        DO iVar = 1, this % source % nvar
          DO j = 0, this % source % interp % N
            DO i = 0, this % source % interp % N

              this % source % interior % hostData(i,j,iVar,iEl) = 0.0_prec

            ENDDO
          ENDDO
        ENDDO
      ENDDO

!    ENDIF

  END SUBROUTINE Source_CompressibleSeaWater2D

  SUBROUTINE Flux_CompressibleSeaWater2D(this)
    IMPLICIT NONE
    CLASS(CompressibleSeaWater2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              IF ( iVar == 1 )THEN !

                ! rho*u*u + p 
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(1,i,j,1,iEl)*& ! u
                      this % solution % interior % hostData(i,j,1,iEl)+& ! rho*u
                      this % pressure % interior % hostData(i,j,1,iEl)

                ! rho*u*v
                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(2,i,j,1,iEl)*& ! v
                      this % solution % interior % hostData(i,j,1,iEl) ! rho*u

              ELSEIF ( iVar == 2 )THEN !

                ! rho*v*u
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(1,i,j,1,iEl)*& ! u
                      this % solution % interior % hostData(i,j,2,iEl) ! rho*v

                ! rho*v*v + p
                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(2,i,j,1,iEl)*& ! v
                      this % solution % interior % hostData(i,j,2,iEl)+& ! rho*v
                      this % pressure % interior % hostData(i,j,1,iEl)



              ELSEIF ( iVar == 3 )THEN ! density
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % solution % interior % hostData(i,j,1,iEl) !rho*u

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % solution % interior % hostData(i,j,2,iEl) !rho*v

              ELSEIF ( iVar == 4 )THEN ! total energy (rho*u*H)

                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(1,i,j,1,iEl)*&
                      this % enthalpy % interior % hostData(i,j,1,iEl) !rho*u*H

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(2,i,j,2,iEl)*&
                      this % enthalpy % interior % hostData(i,j,1,iEl) !rho*v*H


              ELSEIF ( iVar == 5 )THEN ! density weighted salinity

                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(1,i,j,1,iEl)*&
                      this % solution % interior % hostData(i,j,iVar,iEl) !rho*u*S

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(2,i,j,1,iEl)*&
                      this % solution % interior % hostData(i,j,iVar,iEl) !rho*v*S

              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO
!    ENDIF

  END SUBROUTINE Flux_CompressibleSeaWater2D

  SUBROUTINE RiemannSolver_CompressibleSeaWater2D(this)
    IMPLICIT NONE
    CLASS(CompressibleSeaWater2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: nhat(1:2), nmag
    REAL(prec) :: cL, cR, unL, unR, HL, HR, HuL, HuR, HvL, HvR
    REAL(prec) :: alpha
    REAL(prec) :: fluxL(1:5)
    REAL(prec) :: fluxR(1:5)
    REAL(prec) :: jump(1:5)


!    IF( this % gpuAccel )THEN
!
!      CALL RiemannSolver_CompressibleSeaWater2D_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
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

             cL = this % soundSpeed % boundary % hostData(i,1,iSide,iEl)
             cR = this % soundSpeed % extBoundary % hostData(i,1,iSide,iEl)

             fluxL(1) = unL*this % solution % boundary % hostData(i,1,iSide,iEl) +&
                        this % pressure % boundary % hostData(i,1,iSide,iEl)*nHat(1)

             fluxL(2) = unL*this % solution % boundary % hostData(i,2,iSide,iEl) +&
                        this % pressure % boundary % hostData(i,1,iSide,iEl)*nHat(2)

             fluxL(3) = unL*this % solution % boundary % hostData(i,3,iSide,iEl)
             fluxL(4) = unL*this % enthalpy % boundary % hostData(i,1,iSide,iEl)
             fluxL(5) = unL*this % solution % boundary % hostData(i,5,iSide,iEl)

             fluxR(1) = unL*this % solution % extBoundary % hostData(i,1,iSide,iEl) +&
                        this % pressure % extBoundary % hostData(i,1,iSide,iEl)*nHat(1)

             fluxR(2) = unL*this % solution % extBoundary % hostData(i,2,iSide,iEl) +&
                        this % pressure % extBoundary % hostData(i,1,iSide,iEl)*nHat(2)

             fluxR(3) = unL*this % solution % extBoundary % hostData(i,3,iSide,iEl)
             fluxR(4) = unL*this % enthalpy % extBoundary % hostData(i,1,iSide,iEl)
             fluxR(5) = unL*this % solution % extBoundary % hostData(i,5,iSide,iEl)

             jump(1:5) = this % solution % boundary % hostData(i,1:5,iSide,iEl)-&
                         this % solution % extBoundary % hostData(i,1:5,iSide,iEl)

             alpha = MAX( ABS(unL+cL), ABS(unR+cR), &
                          ABS(unL-cL), ABS(unR-cR) )

             ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
             this % flux % boundaryNormal % hostData(i,1:5,iSide,iEl) = 0.5_prec*( &
                     fluxL(1:5) + fluxR(1:5) + alpha*( jump(1:5) ) )*nmag

          ENDDO
        ENDDO
      ENDDO

!    ENDIF

  END SUBROUTINE RiemannSolver_CompressibleSeaWater2D

END MODULE SELF_CompressibleSeaWater2D
