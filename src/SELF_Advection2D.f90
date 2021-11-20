MODULE SELF_Advection2D

USE SELF_Constants
USE SELF_Mesh
USE SELF_DG
USE FEQParse
USE FLAP

! Needed for Fortran-C interoperability
! Helps expose HIP kernels to Fortran
USE ISO_C_BINDING


  TYPE,EXTENDS(DG2D), PUBLIC :: Advection2D

    TYPE(MappedVector2D),PUBLIC :: velocity

    ! Model Settings !
    REAL(prec) :: Lx, Ly ! Domain lengths
    REAL(prec) :: dt ! Default time step size
    INTEGER :: controlDegree
    INTEGER :: targetDegree
    INTEGER :: controlQuadrature ! ENUMS in SELF_Constants.f90
    INTEGER :: targetQuadrature ! ENUMS in SELF_Constants.f90
    CHARACTER(LEN=self_FileNameLength) :: meshFile
    INTEGER :: nxElements
    INTEGER :: nyElements
    INTEGER :: timeIntegrator ! ENUMS needed in SELF_Constants.f90 !! TO DO !!
    CHARACTER(LEN=self_EquationLength) :: velEqnX ! Velocity Equation (x-direction)
    CHARACTER(LEN=self_EquationLength) :: velEqnY ! Velocity Equation (y-direction)
    CHARACTER(LEN=self_EquationLength) :: icEqn ! Initial condition Equation
    CHARACTER(LEN=self_EquationLength) :: bcEqn ! Boundary condition Equation

    

    CONTAINS

      PROCEDURE,PUBLIC :: Init => Init_Advection2D

      PROCEDURE,PUBLIC :: InitFromCLI => InitFromCLI_Advection2D

      PROCEDURE,PUBLIC :: Free => Free_Advection2D


      GENERIC, PUBLIC :: SetSolution => SetSolutionFromEquation!, &
!                                        SetSolutionFromFile
      PROCEDURE, PRIVATE :: SetSolutionFromEquation!, &
!                            SetSolutionFromFile

      GENERIC, PUBLIC :: SetVelocity => SetVelocityFromEquation
      PROCEDURE, PRIVATE :: SetVelocityFromEquation

      GENERIC, PUBLIC :: SetBoundaryCondition => SetBoundaryConditionFromEquation!, &
!                                        SetBoundaryConditionFromFile
      PROCEDURE, PRIVATE :: SetBoundaryConditionFromEquation!, &

 !     PROCEDURE, PUBLIC :: WriteTecplot => WriteTecplot_Advection2D

 !     PROCEDURE, PUBLIC :: TimeStepRK3 => TimeStepRK3_Advection2D
      PROCEDURE, PUBLIC :: Tendency => Tendency_Advection2D
      PROCEDURE, PUBLIC :: InternalFlux => InternalFlux_Advection2D
      PROCEDURE, PUBLIC :: SideFlux => SideFlux_Advection2D

  END TYPE Advection2D


  ! Interfaces to GPU kernels !
  INTERFACE
    SUBROUTINE InternalFlux_Advection2D_gpu_wrapper(flux, solution, velocity, N, nVar, nEl) &
      BIND(c,name="InternalFlux_Advection2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, solution, velocity
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE InternalFlux_Advection2D_gpu_wrapper 
  END INTERFACE

  INTERFACE
    SUBROUTINE SideFlux_Advection2D_gpu_wrapper(flux, boundarySol, extSol, velocity, nHat, nScale, N, nVar, nEl) &
      BIND(c,name="SideFlux_Advection2D_gpu_wrapper")
      USE ISO_C_BINDING
      IMPLICIT NONE
      TYPE(c_ptr) :: flux, boundarySol, extSol, velocity, nHat, nScale
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE SideFlux_Advection2D_gpu_wrapper 
  END INTERFACE

CONTAINS

  SUBROUTINE Init_Advection2D(this,cqType,tqType,cqDegree,tqDegree,nvar,enableMPI,spec)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(out) :: this
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nvar
    LOGICAL,INTENT(in) :: enableMPI
    TYPE(MeshSpec),INTENT(in) :: spec

    CALL this % decomp % Init(enableMPI)

    ! Load Mesh
    CALL this % mesh % Load(spec,this % decomp)

    CALL this % decomp % SetMaxMsg(this % mesh % nUniqueSides)

    ! Create geometry from mesh
    CALL this % geometry % GenerateFromMesh(&
            this % mesh,cqType,tqType,cqDegree,tqDegree)

    CALL this % solution % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % solutionGradient % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % flux % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % velocity % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % source % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
                this % mesh % nElem)

    CALL this % fluxDivergence % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % workScalar % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % workVector % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % workTensor % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    CALL this % compFlux % Init(&
            cqDegree,cqType,tqDegree,tqType,nVar,&
            this % mesh % nElem)

    ALLOCATE (this % solutionMetaData(1:nvar))


  END SUBROUTINE Init_Advection2D

  SUBROUTINE InitFromCLI_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this

    ! Local
    TYPE(COMMAND_LINE_INTERFACE) :: cli

    CALL cli % init(progname="sadv2d", &
                    version="v0.0.0", &
                    description="SELF Advection in 2-D", &
                    license="ANTI-CAPITALIST SOFTWARE LICENSE (v 1.4)", &
                    authors="Joseph Schoonover (Fluid Numerics LLC)")

    ! Get the control degree
    CALL cli % add(switch="--control-degree", &
                   switch_ab="-c", &
                   help="The polynomial degree of the control points."//NEW_LINE("A"), &
                   def="7", &
                   required=.FALSE.)

    ! Get the target degree (assumed for plotting)
    CALL cli % add(switch="--target-degree", &
                   switch_ab="-t", &
                   help="The polynomial degree for the"//&
                  &" target points for interpolation."//&
                  &" Typically used for plotting"//NEW_LINE("A"), &
                   def="14", &
                   required=.FALSE.)

    ! Get the control quadrature
    ! Everyone know Legendre-Gauss Quadrature is the best...
    CALL cli % add(switch="--control-quadrature", &
                   switch_ab="-cq", &
                   def="gauss", &
                   help="The quadrature type for control points."//NEW_LINE("A"), &
                   choices="gauss,gauss-lobatto,chebyshev-gauss,chebyshev-gauss-lobatto", &
                   required=.FALSE.)


    ! Set the target grid quadrature
    ! Default to uniform (assumed for plotting)
    CALL cli % add(switch="--target-quadrature", &
                   switch_ab="-tq", &
                   def="uniform", &
                   help="The quadrature type for target points."//NEW_LINE("A"), &
                   choices="gauss,gauss-lobatto,uniform", &
                   required=.FALSE.)

    ! (Optional) Provide a file for a mesh
    ! Assumed in HOPR or ISM-v2 format
    CALL cli % add(switch="--mesh", &
                   switch_ab="-m", &
                   help="Path to a mesh file for control mesh."//NEW_LINE("A"), &
                   def="", &
                   required=.FALSE.)

    ! (Optional) If a mesh is not provided, you
    ! can request a structured grid to be generated
    ! just set the nxelement, nyelements..
    CALL cli % add(switch="--nxelements", &
                   switch_ab="-nx", &
                   help="The number of elements in the x-direction for structured mesh generation.", &
                   def="5", &
                   required=.FALSE.)

    CALL cli % add(switch="--nyelements", &
                   switch_ab="-ny", &
                   help="The number of elements in the y-direction for structured mesh generation.", &
                   def="5", &
                   required=.FALSE.)

    ! Alright... now tell me some physical mesh dimensions
    CALL cli % add(switch="--xlength", &
                   switch_ab="-lx", &
                   help="The physical x-scale for structured mesh generation."//&
                   " Ignored if a mesh file is provided", &
                   def="1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--ylength", &
                   switch_ab="-ly", &
                   help="The physical y-scale for structured mesh generation."//&
                   " Ignored if a mesh file is provided", &
                   def="1.0", &
                   required=.FALSE.)

    ! Set the velocity field
    CALL cli % add(switch="--velocity-x", &
                   switch_ab="-vx", &
                   help="Equation for the x-component of the velocity field (x,y dependent only!)",&
                   def="vx=1.0", &
                   required=.FALSE.)

    CALL cli % add(switch="--velocity-y", &
                   switch_ab="-vy", &
                   help="Equation for the y-component of the velocity field (x,y dependent only!)",&
                   def="vy=1.0", &
                   required=.FALSE.)

    ! Set the initial conditions
    ! .. TO DO .. 
    !  > How to handle multiple tracer fields ??
    CALL cli % add(switch="--initial-condition", &
                   switch_ab="-ic", &
                   help="Equation for the initial tracer distributions (x,y dependent only!)",&
                   def="f = exp( -( x^2 + y^2)/0.01 )", &
                   required=.FALSE.)

    CALL cli % add(switch="--boundary-condition", &
                   switch_ab="-ic", &
                   help="Equation for the boundary tracer distributions (can be time dependent!)", &
                   def="f = exp( -( (x-t)^2 + (y-t)^2)/0.01 )", &
                   required=.FALSE.)

    ! Give me a time integrator
    CALL cli % add(switch="--integrator", &
                   switch_ab="-int", &
                   help="Sets the time integration method. Only 'euler' or 'williamson_rk3'", &
                   def="williamson_rk3", &
                   required=.FALSE.)

    ! Set the CLI parameters !
    ! All of these need to go to character variables and then
    ! be converted... TO DO !!
    ! I won't bore  you all with this part...
    ! 
    CALL cli % get(val=this % meshfile,switch='--mesh')


  END SUBROUTINE InitFromCLI_Advection2D

  SUBROUTINE Free_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this

    CALL this % mesh % Free()
    CALL this % geometry % Free()
    CALL this % solution % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()
    CALL this % workScalar % Free()
    CALL this % workVector % Free()
    CALL this % workTensor % Free()
    CALL this % compFlux % Free()
    CALL this % velocity % Free()
    DEALLOCATE (this % solutionMetaData)

  END SUBROUTINE Free_Advection2D

  SUBROUTINE SetSolutionFromEquation( this, eqn )
    IMPLICIT NONE
    CLASS(Advection2D), INTENT(inout) :: this
    TYPE(EquationParser), INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: i, j, iEl, iVar
    REAL(prec) :: x
    REAL(prec) :: y


    DO iEl = 1,this % solution % nElem
      DO iVar = 1, this % solution % nVar
        DO j = 0, this % solution % N
          DO i = 0, this % solution % N

             ! Get the mesh positions
             x = this % geometry % x % interior % hostData(1,i,j,1,iEl)
             y = this % geometry % x % interior % hostData(2,i,j,1,iEl)

             this % solution % interior % hostData(i,j,iVar,iEl) = &
               eqn(iVar) % Evaluate((/x, y/))


          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CALL this % solution % interior % UpdateDevice()

  END SUBROUTINE SetSolutionFromEquation

!  SUBROUTINE SetSolutionFromFile( this, filename )
!    IMPLICIT NONE
!    CLASS(Advection2D), INTENT(inout) :: this
!    CHARACTER(*), INTENT(in) :: filename
!    ! Local
!
!
!
!  END SUBROUTINE SetSolutionFromEquation

  SUBROUTINE SetVelocityFromEquation( this, eqn )
    IMPLICIT NONE
    CLASS(Advection2D), INTENT(inout) :: this
    TYPE(EquationParser), INTENT(in) :: eqn(1:2)
    ! Local
    INTEGER :: i, j, iEl, iVar, iSide
    REAL(prec) :: x
    REAL(prec) :: y


    DO iEl = 1,this % solution % nElem

      ! Set the velocity at the element interiors
      DO j = 0, this % solution % N
        DO i = 0, this % solution % N

           ! Get the mesh positions
           x = this % geometry % x % interior % hostData(1,i,j,1,iEl)
           y = this % geometry % x % interior % hostData(2,i,j,1,iEl)

           ! Set the velocity in the x-direction
           this % velocity % interior % hostData(1,i,j,1,iEl) = &
             eqn(1) % Evaluate((/x, y/))

           ! Set the velocity in the y-direction
           this % velocity % interior % hostData(2,i,j,1,iEl) = &
             eqn(2) % Evaluate((/x, y/))


        ENDDO
      ENDDO

      ! Set the velocity at element edges
      DO iSide = 1, 4
        DO i = 0, this % solution % N

           ! Get the mesh positions
           x = this % geometry % x % boundary % hostData(1,i,1,iSide,iEl)
           y = this % geometry % x % boundary % hostData(2,i,1,iSide,iEl)

           ! Set the velocity in the x-direction
           this % velocity % boundary % hostData(1,i,1,iSide,iEl) = &
             eqn(1) % Evaluate((/x, y/))

           ! Set the velocity in the y-direction
           this % velocity % interior % hostData(2,i,1,iSide,iEl) = &
             eqn(2) % Evaluate((/x, y/))


        ENDDO
      ENDDO

    ENDDO

    CALL this % solution % interior % UpdateDevice()
    CALL this % solution % boundary % UpdateDevice()


  END SUBROUTINE SetVelocityFromEquation

  SUBROUTINE SetBoundaryConditionFromEquation( this, eqn )
    IMPLICIT NONE
    CLASS(Advection2D), INTENT(inout) :: this
    TYPE(EquationParser), INTENT(in) :: eqn(1:this % solution % nVar)
    ! Local
    INTEGER :: i, iEl, iVar, iSide
    REAL(prec) :: x
    REAL(prec) :: y


    DO iEl = 1,this % solution % nElem
      DO iVar = 1, this % solution % nvar
        DO iSide = 1, 4
          DO i = 0, this % solution % N

             ! If this element's side has no neighbor assigned
             ! it is assumed to be a physical boundary.
             ! In this case, we want to assign the external boundary
             ! condition.
             IF( this % mesh % self_sideInfo % hostData(3,iSide,iEl) == 0 )THEN
               ! Get the mesh positions
               x = this % geometry % x % boundary % hostData(1,i,1,iSide,iEl)
               y = this % geometry % x % boundary % hostData(2,i,1,iSide,iEl)

               ! Set the external boundary condition
               ! .... TO DO ... 
               !  > Set time tracking mechanism to set time 
               !    varying boundary conditions... 
               !    wouldn't that be slick ?
               this % solution % extBoundary % hostData(i,1,iSide,iEl) = &
                 eqn(iVar) % Evaluate((/x, y/))
             ENDIF


          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Copy data to the GPU
    CALL this % solution % extBoundary % UpdateDevice()


  END SUBROUTINE SetBoundaryConditionFromEquation

  SUBROUTINE Tendency_Advection2D( this, gpuAccel ) 
    IMPLICIT NONE
    CLASS(Advection2D), INTENT(inout) :: this
    LOGICAL, INTENT(in) :: gpuAccel

      CALL this % InternalFlux( gpuAccel )
      CALL this % SideFlux( gpuAccel )
      CALL this % CalculateFluxDivergence( gpuAccel )
      CALL this % CalculatedSdt( gpuAccel )

  END SUBROUTINE Tendency_Advection2D

  SUBROUTINE SideFlux_Advection2D( this, gpuAccel )
    !! Calculates the Advective Flux on element sides using a Lax-Friedrich's upwind Riemann Solver
    IMPLICIT NONE
    CLASS(Advection2D), INTENT(inout) :: this
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: i,iSide,iVar,iEl
    REAL(prec) :: nhat(1:2)
    REAL(prec) :: nmag
    REAL(prec) :: un
    REAL(prec) :: extState
    REAL(prec) :: intState

      ! Exchange side information between neighboring cells
      ! If GPU acceleration is enabled, SideExchange runs on the GPU
      ! and device pointers are updated.
      CALL this % solution % SideExchange( this % mesh, &
                                           this % decomp, &
                                           gpuAccel )

      IF(gpuAccel)THEN

        CALL SideFlux_Advection2D_gpu_wrapper( this % flux % boundary % deviceData, &
                                               this % solution % boundary % deviceData, &
                                               this % solution % extBoundary % deviceData, &
                                               this % velocity % boundary % deviceData, &
                                               this % geometry % nHat % boundary % deviceData, &
                                               this % geometry % nScale % boundary % deviceData, &
                                               this % solution % N, &
                                               this % solution % nVar, &
                                               this % solution % nElem )

      ELSE

        DO iEl = 1, this % solution % nElem
          DO iSide = 1, 4
            DO iVar = 1, this % solution % nVar
              DO i = 0, this % solution % N

                 ! Get the boundary normals on cell edges from the mesh geometry
                 nhat(1:2) = this % geometry % nHat % boundary % hostData(2,i,1,iSide,iEl)

                 ! Calculate the normal velocity at the cell edges
                 un = this % velocity % boundary % hostData(1,i,1,iSide,iEl)*nHat(1)+&
                      this % velocity % boundary % hostData(2,i,1,iSide,iEl)*nHat(2)

                 ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
                 extState = this % solution % extBoundary % hostData(i,iVar,iSide,iEl)
                 intState = this % solution % boundary % hostData(i,iVar,iSide,iEl)
                 nmag = this % geometry % nScale % boundary % hostData(i,iVar,iSide,iEl)

                 ! Calculate the flux
                 ! Since the flux is a vector, we need to store the normal flux somewhere...
                 ! Although this is not ideal (unused memory for other vector components)
                 ! we'll store the boundary flux in the first dimension 
                 ! This will need to be addressed in the SELF_MappedData class..
                 ! maybe put a "normalFlux" attribute.
                 ! TO DO :: Log this in a ticket..
                 this % flux % boundary % hostData(1,i,iVar,iSide,iEl) = 0.5_prec*&
                     ( un*(intState + extState) - abs(un)*(extState - intState) )*nmag

              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

  END SUBROUTINE SideFlux_Advection2D

  SUBROUTINE InternalFlux_Advection2D( this, gpuAccel )
    !! Calculates the advective flux using the provided velocity
    IMPLICIT NONE
    CLASS(Advection2D), INTENT(inout) :: this
    LOGICAL,INTENT(in) :: gpuAccel
    ! Local
    INTEGER :: i,j,iVar,iEl

    IF( gpuAccel )THEN

      ! When GPU acceleration is enabled (requested by the user)
      ! we call the gpu wrapper interface, which will call the
      ! HIP kernel "under the hood"
      CALL InternalFlux_Advection2D_gpu_wrapper(this % flux % interior % deviceData,&
                                                this % solution % interior % deviceData, &
                                                this % velocity % interior % deviceData, &
                                                this % solution % N, & 
                                                this % solution % nVar, &
                                                this % solution % nElem )

    ELSE

      DO iEl = 1,this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % N
            DO i = 0, this % solution % N

              this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(1,i,j,1,iEl)*&
                      this % solution % interior % hostData(i,j,iVar,iEl)

              this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                      this % velocity % interior % hostData(2,i,j,1,iEl)*&
                      this % solution % interior % hostData(i,j,iVar,iEl)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE InternalFlux_Advection2D

END MODULE SELF_Advection2D
