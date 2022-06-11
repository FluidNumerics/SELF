!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_CompressibleIdealGas2D

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model


  TYPE,EXTENDS(Model2D) :: CompressibleIdealGas2D
    !!
    !! For the solution attribute, we use the following convention
    !! iVar = 1 ~> rho*u (x-momentum)
    !! iVar = 2 ~> rho*v (y-momentum)
    !! iVar = 3 ~> rho (density)
    !! iVar = 4 ~> rho*E (Total Energy)

    !! For the requiredDiagnostics attribute, we use the following convention
    !! iVar = 1 ~> Kinetic Energy
    !! iVar = 2 ~> Pressure
    !! iVar = 3 ~> Enthalpy
    !! iVar = 4 ~> Sound Speed
    !! (Planned) iVar = 5 ~> In-Situ Temperature
    !!
    TYPE(MappedScalar2D) :: prescribedSolution
    TYPE(MappedScalar2D) :: requiredDiagnostics
    REAL(prec) :: expansionFactor
    REAL(prec) :: Cp ! Heat capacity at constant pressure ( J/g/K )
    REAL(prec) :: Cv ! Heat capacity at constant volume ( J/g/K )


    CONTAINS

    ! Overridden Methods
    PROCEDURE :: Init => Init_CompressibleIdealGas2D
    PROCEDURE :: Free => Free_CompressibleIdealGas2D
    PROCEDURE :: PreTendency => PreTendency_CompressibleIdealGas2D

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_CompressibleIdealGas2D
    PROCEDURE :: FluxMethod => Flux_CompressibleIdealGas2D
    PROCEDURE :: RiemannSolver => RiemannSolver_CompressibleIdealGas2D
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_CompressibleIdealGas2D

    ! New Methods
    GENERIC :: SetPrescribedSolution => SetPrescribedSolutionFromChar_CompressibleIdealGas2D,&
                              SetPrescribedSolutionFromEqn_CompressibleIdealGas2D,&
                              SetPrescribedSolutionFromSolution_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetPrescribedSolutionFromChar_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetPrescribedSolutionFromEqn_CompressibleIdealGas2D
    PROCEDURE,PRIVATE :: SetPrescribedSolutionFromSolution_CompressibleIdealGas2D

    PROCEDURE :: CalculateVelocity => CalculateVelocity_CompressibleIdealGas2D
    PROCEDURE :: CalculateKineticEnergy => CalculateKineticEnergy_CompressibleIdealGas2D
    PROCEDURE :: EquationOfState => EquationOfState_CompressibleIdealGas2D
    PROCEDURE :: CalculateSoundSpeed => CalculateSoundSpeed_CompressibleIdealGas2D
    PROCEDURE :: CalculateEnthalpy => CalculateEnthalpy_CompressibleIdealGas2D

  END TYPE CompressibleIdealGas2D

  INTEGER, PARAMETER, PRIVATE :: keIndex = 1 ! Index for kinetic energy
  INTEGER, PARAMETER, PRIVATE :: prIndex = 2 ! Index for pressure
  INTEGER, PARAMETER, PRIVATE :: enIndex = 3 ! Index for enthalpy
  INTEGER, PARAMETER, PRIVATE :: soIndex = 4 ! Index for sound speed
  INTEGER, PARAMETER, PRIVATE :: teIndex = 5 ! Index for in-situ temperature
  INTEGER, PARAMETER, PRIVATE :: nRequiredDiagnostics = 5

 ! INTERFACE
 !   SUBROUTINE SetBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper(solution, extBoundary, nHat, sideInfo, N, nVar, nEl) &
 !     bind(c,name="SetBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper")
 !     USE iso_c_binding
 !     USE SELF_Constants
 !     IMPLICIT NONE
 !     TYPE(c_ptr) :: solution, extBoundary, nHat, sideInfo
 !     INTEGER(C_INT),VALUE :: N,nVar,nEl
 !   END SUBROUTINE SetBoundaryCondition_CompressibleIdealGas2D_gpu_wrapper
 ! END INTERFACE

 ! INTERFACE
 !   SUBROUTINE Source_CompressibleIdealGas2D_gpu_wrapper(source, solution, f, N, nVar, nEl) &
 !     bind(c,name="Source_CompressibleIdealGas2D_gpu_wrapper")
 !     USE iso_c_binding
 !     USE SELF_Constants
 !     IMPLICIT NONE
 !     TYPE(c_ptr) :: source, solution
 !     INTEGER(C_INT),VALUE :: N,nVar,nEl
 !     REAL(c_prec),VALUE :: f
 !   END SUBROUTINE Source_CompressibleIdealGas2D_gpu_wrapper
 ! END INTERFACE

 ! INTERFACE
 !   SUBROUTINE Flux_CompressibleIdealGas2D_gpu_wrapper(flux, solution, g, H, N, nVar, nEl) &
 !     bind(c,name="Flux_CompressibleIdealGas2D_gpu_wrapper")
 !     USE iso_c_binding
 !     USE SELF_Constants
 !     IMPLICIT NONE
 !     TYPE(c_ptr) :: flux, solution
 !     INTEGER(C_INT),VALUE :: N,nVar,nEl
 !     REAL(c_prec),VALUE :: g, H
 !   END SUBROUTINE Flux_CompressibleIdealGas2D_gpu_wrapper
 ! END INTERFACE

 ! INTERFACE
 !   SUBROUTINE RiemannSolver_CompressibleIdealGas2D_gpu_wrapper(flux, solution, extBoundary, nHat, nScale, g, H, N, nVar, nEl) &
 !     bind(c,name="RiemannSolver_CompressibleIdealGas2D_gpu_wrapper")
 !     USE iso_c_binding
 !     USE SELF_Constants
 !     IMPLICIT NONE
 !     TYPE(c_ptr) :: flux, solution, extBoundary, nHat, nScale
 !     INTEGER(C_INT),VALUE :: N,nVar,nEl
 !     REAL(c_prec),VALUE :: g, H
 !   END SUBROUTINE RiemannSolver_CompressibleIdealGas2D_gpu_wrapper
 ! END INTERFACE

CONTAINS

  SUBROUTINE Init_CompressibleIdealGas2D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh2D),INTENT(in),TARGET :: mesh
    TYPE(SEMQuad),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp
    ! Local
    INTEGER :: ivar
    CHARACTER(LEN=3) :: ivarChar
    CHARACTER(LEN=25) :: varname
    INTEGER :: nvarloc

    ! Ensure that the number of variables is 4
    ! nvar is unused in this class extension
    nvarloc = 4

    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry
    this % gpuAccel = .FALSE.
    this % fluxDivMethod = SELF_CONSERVATIVE_FLUX 

    CALL this % solution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % prescribedSolution % Init(geometry % x % interp,nvarloc,this % mesh % nElem)
    CALL this % requiredDiagnostics % Init(geometry % x % interp,nRequiredDiagnostics,this % mesh % nElem)
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
    CALL this % solution % SetDescription(4,"Density weighted total (kinetic + internal) energy per unit volume")

    CALL this % requiredDiagnostics % SetName(1,"KE")
    CALL this % requiredDiagnostics % SetUnits(1,"kg/m/s^2")
    CALL this % requiredDiagnostics % SetDescription(1,"Kinetic energy per unit volume")

    CALL this % requiredDiagnostics % SetName(2,"P")
    CALL this % requiredDiagnostics % SetUnits(2,"kg/m/s^2")
    CALL this % requiredDiagnostics % SetDescription(2,"Fluid pressure")

    CALL this % requiredDiagnostics % SetName(3,"H")
    CALL this % requiredDiagnostics % SetUnits(3,"kg/m/s^2")
    CALL this % requiredDiagnostics % SetDescription(3,"Density weighted enthalpy")

    CALL this % requiredDiagnostics % SetName(4,"c")
    CALL this % requiredDiagnostics % SetUnits(4,"m/s")
    CALL this % requiredDiagnostics % SetDescription(4,"Sound speed")

    CALL this % requiredDiagnostics % SetName(5,"T")
    CALL this % requiredDiagnostics % SetUnits(5,"K")
    CALL this % requiredDiagnostics % SetDescription(5,"In-Situ Temperature")


  END SUBROUTINE Init_CompressibleIdealGas2D

  SUBROUTINE Free_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this

    CALL this % solution % Free()
    CALL this % prescribedSolution % Free()
    CALL this % requiredDiagnostics % Free()
    CALL this % workSol % Free()
    CALL this % velocity % Free()
    CALL this % compVelocity % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_CompressibleIdealGas2D

  SUBROUTINE SetPrescribedSolutionFromEqn_CompressibleIdealGas2D(this, eqn) 
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn(1:this % prescribedSolution % nVar)
    ! Local
    INTEGER :: iVar

      ! Copy the equation parser
      DO iVar = 1, this % prescribedSolution % nVar
        CALL this % prescribedSolution % SetEquation(ivar, eqn(iVar) % equation)
      ENDDO

      CALL this % prescribedSolution % SetInteriorFromEquation( this % geometry, this % t )
      CALL this % prescribedSolution % BoundaryInterp( gpuAccel = .FALSE. )

      ! Store the entropy for this state
      CALL this % CalculateEntropy()
      CALL this % ReportEntropy()

      IF( this % gpuAccel )THEN
        CALL this % prescribedSolution % UpdateDevice()
      ENDIF

  END SUBROUTINE SetPrescribedSolutionFromEqn_CompressibleIdealGas2D 

  SUBROUTINE SetPrescribedSolutionFromChar_CompressibleIdealGas2D(this, eqnChar) 
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % prescribedSolution % nVar)
    ! Local
    INTEGER :: iVar

      DO iVar = 1, this % prescribedSolution % nVar
        CALL this % prescribedSolution % SetEquation(ivar, eqnChar(iVar))
      ENDDO

      CALL this % prescribedSolution % SetInteriorFromEquation( this % geometry, this % t )
      CALL this % prescribedSolution % BoundaryInterp( gpuAccel = .FALSE. )

      ! Store the entropy for this state
      CALL this % CalculateEntropy()
      CALL this % ReportEntropy()

      IF( this % gpuAccel )THEN
        CALL this % prescribedSolution % UpdateDevice()
      ENDIF

  END SUBROUTINE SetPrescribedSolutionFromChar_CompressibleIdealGas2D

  SUBROUTINE SetPrescribedSolutionFromSolution_CompressibleIdealGas2D(this) 
  !! Sets the prescribed solution using the current solution attribute
  !! This can be useful for situations where you want to set the
  !! boundary conditions and initial conditions to be identical
  !!
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iEl,iVar

      DO iEl = 1, this % source % nElem
        DO iVar = 1, this % source % nvar
          DO j = 0, this % source % interp % N
            DO i = 0, this % source % interp % N

              this % prescribedSolution % interior % hostData(i,j,iVar,iEl) = &
                this % solution % interior % hostData(i,j,iVar,iEl)
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL this % prescribedSolution % BoundaryInterp( gpuAccel = .FALSE. )

      ! Store the entropy for this state
      CALL this % CalculateEntropy()
      CALL this % ReportEntropy()

      IF( this % gpuAccel )THEN
        CALL this % prescribedSolution % UpdateDevice()
      ENDIF

  END SUBROUTINE SetPrescribedSolutionFromSolution_CompressibleIdealGas2D

  SUBROUTINE SetCp_CompressibleIdealGas2D(this, Cp)
  !! Accessor routine to set the heat capacity at constant pressure
  !! Also updates the expansionFactor attribute when called.
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    REAL(prec), INTENT(in) :: Cp

      this % Cp = Cp
      this % expansionFactor = this % Cp/this % Cv

  END SUBROUTINE SetCp_CompressibleIdealGas2D

  SUBROUTINE SetCv_CompressibleIdealGas2D(this, Cv)
  !! Accessor routine to set the heat capacity at constant volume
  !! Also updates the expansionFactor attribute when called.
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D), INTENT(inout) :: this
    REAL(prec), INTENT(in) :: Cv

      this % Cv = Cv
      this % expansionFactor = this % Cp/this % Cv

  END SUBROUTINE SetCv_CompressibleIdealGas2D

  SUBROUTINE CalculateKineticEnergy_CompressibleIdealGas2D(this)
    !! Calculates the kinetic energy from momentum and density
    !! and stores the output in the requiredDiagnostics output.
    !!
    !! We recognize there are two ways we can calculate the kinetic
    !! energy with the given prognostic and diagnostic variables 
    !! we track.
    !!
    !! Option 1.
    !!
    !!   Use the velocity field diagnostic with the prognostic density.
    !!   This would result in
    !!
    !!    KE = 0.5_prec*rho*( u*u + v*v )
    !!
    !!   where 
    !!
    !!     u = (rho*u)/rho
    !!     v = (rho*v)/rho
    !!
    !! Option 2. 
    !!
    !!   Use the prognostic momentum and density fields. This would
    !!   result in
    !!
    !!    KE = 0.5_prec*( (rho*u)*(rho*u) + (rho*v)*(rho*v) )/rho
    !!
    !! Analytically, the two options are identical. In floating point
    !! arithmetic, these are different.
    !!
    !! It's currently unclear which option is more advantageous (and when),
    !! and I am arbitrarily implementing Option 2.
    !!
    !! If you find a good reason to use Option 1, or some other approach to 
    !! calculate kinetic energy that is more advantageous, develop an
    !! example that highlights the benefits of your approach and open a
    !! pull request.
    !!

    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            this % requiredDiagnostics % interior % hostData(i,j,keIndex,iEl) = 0.5_prec*&
                    (this % solution % interior % hostData(i,j,1,iEl)*&
                    this % solution % interior % hostData(i,j,1,iEl)-&
                    this % solution % interior % hostData(i,j,2,iEl)*&
                    this % solution % interior % hostData(i,j,2,iEl))/&
                    this % solution % interior % hostData(i,j,4,iEl)

          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE CalculateKineticEnergy_CompressibleIdealGas2D

  SUBROUTINE CalculateVelocity_CompressibleIdealGas2D(this)
    !! Calculates the velocity field from momentum and stores
    !! the output in the velocity attribute.
    !!
    !!     u = (rho*u)/rho
    !!     v = (rho*v)/rho
    !!

    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rho

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            rho = this % solution % interior % hostData(i,j,3,iEl)

            this % velocity % interior % hostData(1,i,j,1,iEl) = &
                    this % solution % interior % hostData(i,j,1,iEl)/rho 

            this % velocity % interior % hostData(2,i,j,1,iEl) = &
                    this % solution % interior % hostData(i,j,2,iEl)/rho

          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE CalculateVelocity_CompressibleIdealGas2D

  SUBROUTINE EquationOfState_CompressibleIdealGas2D(this)
    !! Calculates the fluid pressure given other diagnostic fields.
    !! We use the Ideal Gas Law
    !! 
    !!    p = (\gamma-1)*\rho*e
    !!
    !! where $\gamma = \frac{C_p}{C_v}$ is the expansion coefficient,
    !! $\rho$ is the fluid density, and $e$ is the internal energy.
    !!
    !! We calculate $rho*e$ as
    !!
    !!   rho*e = (rho*E - 0.5_prec*rho*KE)
    !!
    !! where rho*E is the total energy, a prognostic variable, modelled
    !! by this class, and 

    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rho, rhoE, rhoKE

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            rho = this % solution % interior % hostData(i,j,3,iEl)
            rhoE = this % solution % interior % hostData(i,j,4,iEl)
            rhoKE = this % requiredDiagnostics % interior % hostData(i,j,keIndex,iEl)

            this % requiredDiagnostics % interior % hostData(i,j,prIndex,iEl) = &
                    (this % expansionFactor - 1.0_prec)*(rhoE - rhoKE)/rho

          ENDDO
        ENDDO
      ENDDO


  END SUBROUTINE EquationOfState_CompressibleIdealGas2D

  SUBROUTINE CalculateSoundSpeed_CompressibleIdealGas2D(this)
    !! Calculates the speed of sound given other diagnostic fields.
    !! We use the Ideal Gas Law
    !!
    !!    p = (\gamma-1)*\rho*e
    !!
    !! The speed of sound is defined through the relation
    !!
    !!   \frac{\partial P}{\partial \rho} = c^{2}
    !!
    !! Then, we have that
    !!
    !!   c = ((\gamma-1)*e)^{1/2}
    !!
    !! where gamma = Cp/Cv is the expansion coefficient,
    !! rho is the fluid density, and e is the internal energy.
    !!
    !! We calculate $e$ as
    !!
    !!   e = (\rho*E - 0.5_prec*\rho*KE)/\rho
    !!
    !! where rho*E is the total energy, a prognostic variable, modelled
    !! by this class, $\rho*KE$ is the kinetic energy (a required diagnostic)
    !! and, $\rho$ is the density ( a prognostic variable ).
    !!

    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rho, rhoE, rhoKE

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            rho = this % solution % interior % hostData(i,j,3,iEl)
            rhoE = this % solution % interior % hostData(i,j,4,iEl)
            rhoKE = this % requiredDiagnostics % interior % hostData(i,j,keIndex,iEl)

            this % requiredDiagnostics % interior % hostData(i,j,soIndex,iEl) = &
                    sqrt((this % expansionFactor - 1.0_prec)*(rhoE - rhoKE)/rho)


          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE CalculateSoundSpeed_CompressibleIdealGas2D

  SUBROUTINE CalculateEnthalpy_CompressibleIdealGas2D(this)
    !! Calculates the dynamic enthalpy from the total energy
    !! and pressure field and stores the output in the 
    !! requiredDiagnostics attribute
    !!
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rhoE, pressure

      DO iEl = 1, this % solution % nElem
        DO j = 0, this % solution % interp % N
          DO i = 0, this % solution % interp % N

            rhoE = this % solution % interior % hostData(i,j,4,iEl)
            pressure = this % requiredDiagnostics % interior % hostData(i,j,prIndex,iEl)

            this % requiredDiagnostics % interior % hostData(i,j,enIndex,iEl) = &
                    rhoE + pressure

          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE CalculateEnthalpy_CompressibleIdealGas2D

  SUBROUTINE PreTendency_CompressibleIdealGas2D(this)
    !! Calculate the velocity and density weighted enthalpy at element interior and element boundaries
    !! PreTendency is a template routine that is used to house any additional calculations
    !! that you want to execute at the beginning of the tendency calculation routine.
    !! This default PreTendency simply returns back to the caller without executing any instructions
    !!
    !! The intention is to provide a method that can be overridden through type-extension, to handle
    !! any steps that need to be executed before proceeding with the usual tendency calculation methods.
    !!
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iSide, j, i
    REAL(prec) :: rho

      ! Calculate the velocity field from
      ! the momentum and density
      CALL this % CalculateVelocity()

      ! Calculate the kinetic energy from
      ! the momentum and density
      CALL this % CalculateKineticEnergy()

      ! Calculate the pressure using an
      ! equation of state.
      ! Requires knowledge of the fluid kinetic energy
      ! and total energy (to diagnose internal energy)
      ! and therefore depends on the CalculateKineticEnergy
      ! call above
      CALL this % EquationOfState()

      ! Calculate the speed of sound
      ! Requires knowledge of the fluid kinetic energy
      ! and total energy (to diagnose internal energy)
      ! and therefore depends on the CalculateKineticEnergy
      ! call above
      CALL this % CalculateSoundSpeed()

      ! Calculates the fluid enthalpy
      ! Requires knowledge of the fluid pressure
      ! and total energy and therefore depends on
      ! the EquationOfState call above.
      CALL this % CalculateEnthalpy()

      ! Interpolate velocity and required diagnostics to the element boundaries
      CALL this % velocity % BoundaryInterp(this % gpuAccel)
      CALL this % requiredDiagnostics % BoundaryInterp(this % gpuAccel)

      ! Perform any MPI exchanges for the velocity and the required diagnostics
      ! across shared element faces between neighboring processes.
      CALL this % velocity % SideExchange(this % mesh, this % decomp, this % gpuAccel)
      CALL this % requiredDiagnostics % SideExchange(this % mesh, this % decomp, this % gpuAccel)

  END SUBROUTINE PreTendency_CompressibleIdealGas2D

  SUBROUTINE SetBoundaryCondition_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
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

              ELSEIF( bcid == SELF_BC_PRESCRIBED )THEN

                this % solution % extBoundary % hostData(i,1,iSide,iEl) = &
                        this % prescribedSolution % boundary % hostData(i,1,iSide,iEl)
                this % solution % extBoundary % hostData(i,2,iSide,iEl) = &
                        this % prescribedSolution % boundary % hostData(i,2,iSide,iEl)
                this % solution % extBoundary % hostData(i,3,iSide,iEl) = &
                        this % prescribedSolution % boundary % hostData(i,3,iSide,iEl)
                this % solution % extBoundary % hostData(i,4,iSide,iEl) = &
                        this % prescribedSolution % boundary % hostData(i,4,iSide,iEl)
                this % solution % extBoundary % hostData(i,5,iSide,iEl) = &
                        this % prescribedSolution % boundary % hostData(i,5,iSide,iEl)

              ENDIF

            ENDIF

          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE SetBoundaryCondition_CompressibleIdealGas2D 

  SUBROUTINE Source_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
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

  END SUBROUTINE Source_CompressibleIdealGas2D

  SUBROUTINE Flux_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
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
                      this % requiredDiagnostics % interior % hostData(i,j,prIndex,iEl)

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
                      this % requiredDiagnostics % interior % hostData(i,j,prIndex,iEl)



              ELSEIF ( iVar == 3 )THEN ! density
                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % solution % interior % hostData(i,j,1,iEl) !rho*u

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % solution % interior % hostData(i,j,2,iEl) !rho*v

              ELSEIF ( iVar == 4 )THEN ! total energy (rho*u*H)

                this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(1,i,j,1,iEl)*&
                      this % requiredDiagnostics % interior % hostData(i,j,enIndex,iEl) !rho*u*H

                this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(2,i,j,2,iEl)*&
                      this % requiredDiagnostics % interior % hostData(i,j,enIndex,iEl) !rho*v*H


              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO
!    ENDIF

  END SUBROUTINE Flux_CompressibleIdealGas2D

  SUBROUTINE RiemannSolver_CompressibleIdealGas2D(this)
    IMPLICIT NONE
    CLASS(CompressibleIdealGas2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,iSide,iEl
    REAL(prec) :: nhat(1:2), nmag
    REAL(prec) :: cL, cR, unL, unR, HL, HR, HuL, HuR, HvL, HvR
    REAL(prec) :: alpha
    REAL(prec) :: fluxL(1:4)
    REAL(prec) :: fluxR(1:4)
    REAL(prec) :: jump(1:4)


!    IF( this % gpuAccel )THEN
!
!      CALL RiemannSolver_CompressibleIdealGas2D_gpu_wrapper(this % flux % boundaryNormal % deviceData, &
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

             cL = this % requiredDiagnostics % boundary % hostData(i,soIndex,iSide,iEl)
             cR = this % requiredDiagnostics % extBoundary % hostData(i,soIndex,iSide,iEl)

             fluxL(1) = unL*this % solution % boundary % hostData(i,1,iSide,iEl) +&
                        this % requiredDiagnostics % boundary % hostData(i,prIndex,iSide,iEl)*nHat(1)

             fluxL(2) = unL*this % solution % boundary % hostData(i,2,iSide,iEl) +&
                        this % requiredDiagnostics % boundary % hostData(i,prIndex,iSide,iEl)*nHat(2)

             fluxL(3) = unL*this % solution % boundary % hostData(i,3,iSide,iEl)
             fluxL(4) = unL*this % requiredDiagnostics % boundary % hostData(i,enIndex,iSide,iEl)

             fluxR(1) = unL*this % solution % extBoundary % hostData(i,1,iSide,iEl) +&
                        this % requiredDiagnostics % extBoundary % hostData(i,prIndex,iSide,iEl)*nHat(1)

             fluxR(2) = unL*this % solution % extBoundary % hostData(i,2,iSide,iEl) +&
                        this % requiredDiagnostics % extBoundary % hostData(i,prIndex,iSide,iEl)*nHat(2)

             fluxR(3) = unL*this % solution % extBoundary % hostData(i,3,iSide,iEl)
             fluxR(4) = unL*this % requiredDiagnostics % extBoundary % hostData(i,enIndex,iSide,iEl)

             jump(1:4) = this % solution % boundary % hostData(i,1:4,iSide,iEl)-&
                         this % solution % extBoundary % hostData(i,1:4,iSide,iEl)

             alpha = MAX( ABS(unL+cL), ABS(unR+cR), &
                          ABS(unL-cL), ABS(unR-cR) )

             ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
             this % flux % boundaryNormal % hostData(i,1:4,iSide,iEl) = 0.5_prec*( &
                     fluxL(1:4) + fluxR(1:4) + alpha*( jump(1:4) ) )*nmag

          ENDDO
        ENDDO
      ENDDO

!    ENDIF

  END SUBROUTINE RiemannSolver_CompressibleIdealGas2D

END MODULE SELF_CompressibleIdealGas2D
