!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Model

  USE SELF_Metadata
  USE SELF_MPI
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_HDF5
  USE HDF5

  TYPE,ABSTRACT :: Model
    INTEGER :: timeIntegrator
    LOGICAL :: gpuAccel

    CONTAINS

!    PROCEDURE,PRIVATE :: ForwardStepEuler 
    PROCEDURE(CalculateTendency),DEFERRED :: CalculateTendency

  END TYPE Model


  TYPE,EXTENDS(Model),ABSTRACT :: Model2D
    TYPE(MappedScalar2D) :: solution
    TYPE(MappedVector2D) :: solutionGradient
    TYPE(MappedVector2D) :: flux
    TYPE(MappedScalar2D) :: source
    TYPE(MappedScalar2D) :: fluxDivergence
    TYPE(MappedScalar2D) :: dSdt
    TYPE(MPILayer),POINTER :: decomp
    TYPE(Mesh2D),POINTER :: mesh
    TYPE(SEMQuad),POINTER :: geometry

    CONTAINS

    PROCEDURE :: Init => Init_Model2D
    PROCEDURE :: Free => Free_Model2D

    PROCEDURE :: UpdateHost => UpdateHost_Model2D
    PROCEDURE :: UpdateDevice => UpdateDevice_Model2D

    PROCEDURE :: CalculateTendency => CalculateTendency_Model2D
    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_Model2D

    PROCEDURE(Source2D),DEFERRED :: Source2D
    PROCEDURE(Flux2D),DEFERRED :: Flux2D
    PROCEDURE(RiemannSolver2D),DEFERRED :: RiemannSolver2D

    PROCEDURE :: ReprojectFlux => ReprojectFlux_Model2D

    PROCEDURE :: Read => Read_Model2D
    PROCEDURE :: Write => Write_Model2D

  END TYPE Model2D

  INTERFACE 
    SUBROUTINE CalculateTendency( this )
      IMPORT Model
      IMPLICIT NONE
      CLASS(Model),INTENT(inout) :: this
    END SUBROUTINE CalculateTendency
  END INTERFACE

  INTERFACE 
    SUBROUTINE Flux2D( this )
      IMPORT Model2D
      IMPLICIT NONE
      CLASS(Model2D),INTENT(inout) :: this
    END SUBROUTINE Flux2D
  END INTERFACE

  INTERFACE 
    SUBROUTINE Source2D( this )
      IMPORT Model2D
      IMPLICIT NONE
      CLASS(Model2D),INTENT(inout) :: this
    END SUBROUTINE Source2D
  END INTERFACE

  INTERFACE 
    SUBROUTINE RiemannSolver2D( this )
      IMPORT Model2D
      IMPLICIT NONE
      CLASS(Model2D),INTENT(inout) :: this
    END SUBROUTINE RiemannSolver2D
  END INTERFACE

CONTAINS

  SUBROUTINE Init_Model2D(this,nvar,mesh,geometry,decomp)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(out) :: this
    INTEGER,INTENT(in) :: nvar
    TYPE(Mesh2D),INTENT(in),TARGET :: mesh
    TYPE(SEMQuad),INTENT(in),TARGET :: geometry
    TYPE(MPILayer),INTENT(in),TARGET :: decomp

    
    this % decomp => decomp
    this % mesh => mesh
    this % geometry => geometry

    CALL this % decomp % SetMaxMsg(this % mesh % nUniqueSides)
    CALL this % decomp % setElemToRank(this % mesh % nElem)

    CALL this % solution % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % dSdt % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % solutionGradient % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % flux % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % source % Init(geometry % x % interp,nVar,this % mesh % nElem)
    CALL this % fluxDivergence % Init(geometry % x % interp,nVar,this % mesh % nElem)

  END SUBROUTINE Init_Model2D

  SUBROUTINE Free_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % mesh % Free()
    CALL this % geometry % Free()
    CALL this % solution % Free()
    CALL this % dSdt % Free()
    CALL this % solutionGradient % Free()
    CALL this % flux % Free()
    CALL this % source % Free()
    CALL this % fluxDivergence % Free()

  END SUBROUTINE Free_Model2D

  SUBROUTINE UpdateHost_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % mesh % UpdateHost()
    CALL this % geometry % UpdateHost()
    CALL this % solution % UpdateHost()
    CALL this % dSdt % UpdateHost()
    CALL this % solutionGradient % UpdateHost()
    CALL this % flux % UpdateHost()
    CALL this % source % UpdateHost()
    CALL this % fluxDivergence % UpdateHost()

  END SUBROUTINE UpdateHost_Model2D

  SUBROUTINE UpdateDevice_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % mesh % UpdateDevice()
    CALL this % geometry % UpdateDevice()
    CALL this % dSdt % UpdateDevice()
    CALL this % solutionGradient % UpdateDevice()
    CALL this % flux % UpdateDevice()
    CALL this % source % UpdateDevice()
    CALL this % fluxDivergence % UpdateDevice()

  END SUBROUTINE UpdateDevice_Model2D

  SUBROUTINE ReprojectFlux_Model2D(this) 
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    ! Local
    INTEGER :: iEl, iVar, j, i
    REAL(prec) :: Fx, Fy

      DO iEl = 1,this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              Fx = this % flux % interior % hostData(1,i,j,iVar,iEl)
              Fy = this % flux % interior % hostData(2,i,j,iVar,iEl)

              this % flux % interior % hostData(1,i,j,iVar,iEl) = &
                this % geometry % dsdx % interior % hostData(1,1,i,j,1,iel)*Fx + &
                this % geometry % dsdx % interior % hostData(2,1,i,j,1,iel)*Fy 

              this % flux % interior % hostData(2,i,j,iVar,iEl) = &
                this % geometry % dsdx % interior % hostData(1,2,i,j,1,iel)*Fx + &
                this % geometry % dsdx % interior % hostData(2,2,i,j,1,iel)*Fy 


            ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE ReprojectFlux_Model2D

  SUBROUTINE CalculateFluxDivergence_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this

    CALL this % flux % Divergence(this % geometry, &
                                  this % fluxDivergence, &
                                  selfWeakDGForm,&
                                  this % gpuAccel)

  END SUBROUTINE CalculateFluxDivergence_Model2D

  SUBROUTINE CalculateTendency_Model2D(this)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, iVar, iEl

!      CALL this % solution % AverageSides()
!      CALL this % solution % DiffSides()
      CALL this % Source2D()
      CALL this % RiemannSolver2D()
      CALL this % Flux2D()
      CALL this % ReprojectFlux()
      CALL this % CalculateFluxDivergence()

   ! IF( gpuAccel )THEN

   !   CALL CalculateDSDt_Model2D_gpu_wrapper( this % fluxDivergence % interior % deviceData, &
   !                                   this % source % interior % deviceData, &
   !                                   this % dSdt % interior % deviceData, &
   !                                   this % solution % interp % N, &
   !                                   this % solution % nVar, &
   !                                   this % solution % nElem ) 
   !                                   
   ! ELSE

      DO iEl = 1, this % solution % nElem
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              this % dSdt % interior % hostData(i,j,iVar,iEl) = &
                      this % source % interior % hostData(i,j,iVar,iEl) -&
                      this % fluxDivergence % interior % hostData(i,j,iVar,iEl)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

   ! ENDIF

  END SUBROUTINE CalculateTendency_Model2D

  SUBROUTINE Write_Model2D(this,fileName)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(in) :: this
    CHARACTER(*),INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:4)
    INTEGER(HID_T) :: xOffset(1:5)
    INTEGER(HID_T) :: bOffset(1:4)
    INTEGER(HID_T) :: bxOffset(1:5)
    INTEGER(HID_T) :: solGlobalDims(1:4)
    INTEGER(HID_T) :: xGlobalDims(1:5)
    INTEGER(HID_T) :: bGlobalDims(1:4)
    INTEGER(HID_T) :: bxGlobalDims(1:5)
    INTEGER :: firstElem

    IF (this % decomp % mpiEnabled) THEN

      CALL Open_HDF5(fileName,H5F_ACC_TRUNC_F,fileId,this % decomp % mpiComm)

      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)
      solOffset(1:4) = (/0,0,1,firstElem/)
      solGlobalDims(1:4) = (/this % solution % interp % N, &
                             this % solution % interp % N, &
                             this % solution % nVar, &
                             this % decomp % nElem/)


      xOffset(1:5) = (/1,0,0,1,firstElem/)
      xGlobalDims(1:5) = (/2, &
                           this % solution % interp % N, &
                           this % solution % interp % N, &
                           this % solution % nVar, &
                           this % decomp % nElem/)

      ! Offsets and dimensions for element boundary data
      bOffset(1:4) = (/0,1,1,firstElem/)
      bGlobalDims(1:4) = (/this % solution % interp % N, &
                           this % solution % nVar, &
                           4,&
                           this % decomp % nElem/)

      bxOffset(1:5) = (/1,0,1,1,firstElem/)
      bxGlobalDims(1:5) = (/2,&
                           this % solution % interp % N, &
                           this % solution % nVar, &
                           4,&
                           this % decomp % nElem/)

      
      CALL CreateGroup_HDF5(fileId,'/quadrature')

      IF( this % decomp % rankId == 0 )THEN
        CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
                             this % solution % interp % controlPoints)

        CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
                             this % solution % interp % qWeights)

        CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
                             this % solution % interp % dgMatrix)

        CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
                             this % solution % interp % dMatrix)
      ENDIF

      CALL CreateGroup_HDF5(fileId,'/state')

      CALL CreateGroup_HDF5(fileId,'/state/interior')

      CALL CreateGroup_HDF5(fileId,'/state/boundary')

      CALL CreateGroup_HDF5(fileId,'/mesh')

      CALL CreateGroup_HDF5(fileId,'/mesh/interior')

      CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

      CALL WriteArray_HDF5(fileId,'/state/interior/solution', &
                           this % solution % interior,solOffset,solGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solution', &
                           this % solution % boundary,bOffset,bGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/fluxDivergence', &
                           this % fluxDivergence % interior,solOffset,solGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/flux', &
                           this % flux % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/flux', &
                           this % flux % boundary,bxOffset,bxGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/interior/solutionGradient', &
                           this % solutionGradient % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solutionGradient', &
                           this % solutionGradient % boundary,bxOffset,bxGlobalDims)

      CALL WriteArray_HDF5(fileId,'/mesh/interior/x', &
                           this % geometry % x % interior,xOffset,xGlobalDims)

      CALL WriteArray_HDF5(fileId,'/mesh/boundary/x', &
                           this % geometry % x % boundary,bxOffset,bxGlobalDims)

      CALL Close_HDF5(fileId)

    ELSE

      CALL Open_HDF5(fileName,H5F_ACC_TRUNC_F,fileId)

      CALL CreateGroup_HDF5(fileId,'/quadrature')

      CALL WriteArray_HDF5(fileId,'/quadrature/xi', &
                           this % solution % interp % controlPoints)

      CALL WriteArray_HDF5(fileId,'/quadrature/weights', &
                           this % solution % interp % qWeights)

      CALL WriteArray_HDF5(fileId,'/quadrature/dgmatrix', &
                           this % solution % interp % dgMatrix)

      CALL WriteArray_HDF5(fileId,'/quadrature/dmatrix', &
                           this % solution % interp % dMatrix)

      CALL CreateGroup_HDF5(fileId,'/state')

      CALL CreateGroup_HDF5(fileId,'/state/interior')

      CALL CreateGroup_HDF5(fileId,'/state/boundary')

      CALL CreateGroup_HDF5(fileId,'/mesh')

      CALL CreateGroup_HDF5(fileId,'/mesh/interior')

      CALL CreateGroup_HDF5(fileId,'/mesh/boundary')

      CALL WriteArray_HDF5(fileId,'/state/interior/solution',this % solution % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solution',this % solution % boundary)

      CALL WriteArray_HDF5(fileId,'/state/interior/fluxDivergence',this % fluxDivergence % interior)

      CALL WriteArray_HDF5(fileId,'/state/interior/flux',this % flux % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/flux',this % flux % boundary)

      CALL WriteArray_HDF5(fileId,'/state/interior/solutionGradient',this % solutionGradient % interior)

      CALL WriteArray_HDF5(fileId,'/state/boundary/solutionGradient',this % solutionGradient % boundary)

      CALL WriteArray_HDF5(fileId,'/mesh/interior/x',this % geometry % x % interior)

      CALL WriteArray_HDF5(fileId,'/mesh/boundary/x',this % geometry % x % boundary)

      CALL Close_HDF5(fileId)

    END IF

  END SUBROUTINE Write_Model2D

  SUBROUTINE Read_Model2D(this,fileName)
    IMPLICIT NONE
    CLASS(Model2D),INTENT(inout) :: this
    CHARACTER(*),INTENT(in) :: fileName
    ! Local
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:4)
    INTEGER :: firstElem
    INTEGER :: N

    IF (this % decomp % mpiEnabled) THEN
      CALL Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId, &
                     this % decomp % mpiComm)
    ELSE
      CALL Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId)
    END IF

    CALL ReadAttribute_HDF5(fileId,'N',N)

    IF (this % solution % interp % N /= N) THEN
      STOP 'Error : Solution polynomial degree does not match input file'
    END IF

    IF (this % decomp % mpiEnabled) THEN
      firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId) + 1
      solOffset(1:4) = (/0,0,1,firstElem/)
      CALL ReadArray_HDF5(fileId,'/state/interior/solution', &
                          this % solution % interior,solOffset)
    ELSE
      CALL ReadArray_HDF5(fileId,'/state/interior/solution',this % solution % interior)
    END IF

    CALL Close_HDF5(fileId)

  END SUBROUTINE Read_Model2D

END MODULE SELF_Model
