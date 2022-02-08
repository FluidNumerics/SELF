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

    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_Advection2D

    GENERIC :: SetVelocityField => SetVelocityFieldFromChar_Advection2D,&
                              SetVelocityFieldFromEqn_Advection2D
    PROCEDURE,PRIVATE :: SetVelocityFieldFromChar_Advection2D
    PROCEDURE,PRIVATE :: SetVelocityFieldFromEqn_Advection2D

  END TYPE Advection2D

  INTERFACE
   SUBROUTINE FluxDivergence_Advection2D_SplitForm_gpu_wrapper(dgMatrixT_dev,&
                   dMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,s_dev,v_dev,bf_dev,df_dev,N,nVar,nEl) &
      bind(c,name="FluxDivergence_Advection2D_SplitForm_gpu_wrapper")
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(c_ptr) :: dgMatrixT_dev,dMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,s_dev,v_dev,bf_dev,df_dev
      INTEGER(C_INT),VALUE :: N,nVar,nEl
    END SUBROUTINE FluxDivergence_Advection2D_SplitForm_gpu_wrapper
  END INTERFACE
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

    CALL this % velocity % Init(geometry % x % interp,1,this % mesh % nElem)
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

  SUBROUTINE SetVelocityFieldFromEqn_Advection2D(this, eqn) 
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this
    TYPE(EquationParser),INTENT(in) :: eqn(1:2)

      ! Copy the equation parser
      ! Set the x-component of the velocity
      CALL this % velocity % SetEquation(1,1,eqn(1) % equation)

      ! Set the y-component of the velocity
      CALL this % velocity % SetEquation(2,1,eqn(2) % equation)

      ! Set the velocity values using the equation parser
      CALL this % velocity % SetInteriorFromEquation( this % geometry, this % t )

  END SUBROUTINE SetVelocityFieldFromEqn_Advection2D 

  SUBROUTINE SetVelocityFieldFromChar_Advection2D(this, eqnChar) 
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this
    CHARACTER(LEN=SELF_EQUATION_LENGTH),INTENT(in) :: eqnChar(1:this % velocity % nVar)

      ! Set the x-component of the velocity
      CALL this % velocity % SetEquation(1,1,eqnChar(1))

      ! Set the y-component of the velocity
      CALL this % velocity % SetEquation(2,1,eqnChar(2))

      ! Set the velocity values using the equation parser
      CALL this % velocity % SetInteriorFromEquation( this % geometry, this % t )

  END SUBROUTINE SetVelocityFieldFromChar_Advection2D

  SUBROUTINE CalculateFluxDivergence_Advection2D(this)
    IMPLICIT NONE
    CLASS(Advection2D),INTENT(inout) :: this
    ! Local
    INTEGER    :: i,j,ii,iVar,iEl
    REAL(prec) :: dF

    IF ( this % gpuAccel )THEN
          
      CALL FluxDivergence_Advection2D_SplitForm_gpu_wrapper(this % solution % interp % dgMatrix % deviceData, &
                                             this % solution % interp % dMatrix % deviceData, &
                                             this % solution % interp % bMatrix % deviceData, &
                                             this % solution % interp % qWeights % deviceData, &
                                             this % flux % interior % deviceData, &
                                             this % solution % interior % deviceData, &
                                             this % velocity % interior % deviceData, &
                                             this % flux % boundaryNormal % deviceData, &
                                             this % fluxDivergence % interior % deviceData, &
                                             this % solution % interp % N,&
                                             this % solution % nVar,&
                                             this % solution % nElem)
  
    ELSE 
      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              dF = 0.0_prec
              DO ii = 0,this % solution % interp % N
                dF = dF + 0.5_prec*( this % solution % interp % dgMatrix % hostData(ii,i)*&
                                     this % flux % interior % hostData(1,ii,j,iVar,iEl) + &
                                     this % solution % interp % dgMatrix % hostData(ii,j)*&
                                     this % flux % interior % hostData(2,i,ii,iVar,iEl) + &
                                     this % solution % interior % hostData(i,j,iVar,iEl)*&
                                     ( this % velocity % interior % hostData(1,ii,j,1,iEl) -&
                                       this % velocity % interior % hostData(1,i,j,1,iEl) )*&
                                     this % solution % interp % dMatrix % hostData(ii,i) + &
                                     this % solution % interior % hostData(i,j,iVar,iEl)*&
                                     ( this % velocity % interior % hostData(2,i,ii,1,iEl) -&
                                       this % velocity % interior % hostData(2,i,j,1,iEl) )*&
                                     this % solution % interp % dMatrix % hostData(ii,j) )
              END DO

              ! Ok to stay ! Boundary terms do not change
              dF = dF + (this % solution % interp % bMatrix % hostData(i,1)*&
                         this % flux % boundaryNormal % hostData(j,iVar,2,iEl) + &
                         this % solution % interp % bMatrix % hostData(i,0)*&
                         this % flux % boundaryNormal % hostData(j,iVar,4,iEl))/ &
                         this % solution % interp % qWeights % hostData(i) + &
                         (this % solution % interp % bMatrix % hostData(j,1)*&
                          this % flux % boundaryNormal % hostData(i,iVar,3,iEl) + &
                          this % solution % interp % bMatrix % hostData(j,0)*&
                          this % flux % boundaryNormal % hostData(i,iVar,1,iEl))/ &
                         this % solution % interp % qWeights % hostData(j)

             this % fluxDivergence % interior % hostData(i,j,iVar,iEl) = dF

            END DO
          END DO
        END DO
      END DO
    ENDIF


  END SUBROUTINE CalculateFluxDivergence_Advection2D

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
