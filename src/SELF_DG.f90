!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_DG

USE SELF_Metadata
USE SELF_MPI
USE SELF_Mesh
USE SELF_MappedData
USE SELF_HDF5
USE HDF5

!
! User supplies : 
!  * Interior physical flux
!  * Boundary physical flux
!  * Source Terms
!
! How someone should use...
!  > Create the object
!  > Loop over time
!     > Fill in source
!     > Fill in interior flux (Usually a Riemmann Solver); if 2nd order operators are included, bassi-rebay flux routine is available to
!     support gradient calculation.
!     > ForwardStep
!  > 

TYPE, PUBLIC :: DG2D
    TYPE(MappedScalar2D), PUBLIC :: solution
    TYPE(MappedVector2D), PUBLIC :: solutionGradient
    TYPE(MappedVector2D), PUBLIC :: flux
    TYPE(MappedScalar2D), PUBLIC :: source
    TYPE(MappedScalar2D), PUBLIC :: fluxDivergence
    TYPE(MPILayer), PUBLIC :: decomp
    TYPE(Mesh2D), PUBLIC :: mesh
    TYPE(SEMQuad), PUBLIC :: geometry
    TYPE(Metadata), ALLOCATABLE, PUBLIC :: solutionMetaData(:)

    ! Work arrays
    TYPE(MappedScalar2D), PRIVATE :: workScalar
    TYPE(MappedVector2D), PRIVATE :: workVector
    TYPE(MappedVector2D), PRIVATE :: compFlux
    TYPE(MappedTensor2D), PRIVATE :: workTensor

    CONTAINS

      PROCEDURE, PUBLIC :: Init => Init_DG2D
      PROCEDURE, PUBLIC :: Free => Free_DG2D

      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_DG2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_DG2D

      ! Accessors
    !  PROCEDURE, PUBLIC :: GetNumberOfVars => GetNumberOfVars_DG2D
    !  PROCEDURE, PUBLIC :: GetControlDegree => GetControlDegree_DG2D
    !  PROCEDURE, PUBLIC :: GetTargetDegree => GetTargetDegree_DG2D
    !  PROCEDURE, PUBLIC :: GetControlQuadrature => GetControlQuadrature_DG2D
    !  PROCEDURE, PUBLIC :: GetTargetQuadrature => GetTargetQuadrature_DG2D
    !  PROCEDURE, PUBLIC :: GetNumberOfElement => GetNumberOfElements_DG2D
    !  PROCEDURE, PUBLIC :: GetNumberOfGlobalSides => GetNumberOfGlobalSides_DG2D
    !  PROCEDURE, PUBLIC :: GetNumberOfUniqueSides => GetNumberOfUniqueSides_DG2D
    !  PROCEDURE, PUBLIC :: GetNumberOfGlobalNodes => GetNumberOfGlobalNodes_DG2D
    !  PROCEDURE, PUBLIC :: GetNumberOfUniqueNodes => GetNumberOfUniqueNodes_DG2D

      PROCEDURE, PUBLIC :: CalculateSolutionGradient => CalculateSolutionGradient_DG2D 
      PROCEDURE, PUBLIC :: CalculateFluxDivergence => CalculateFluxDivergence_DG2D

      PROCEDURE, PUBLIC :: Read => Read_DG2D 
      PROCEDURE, PUBLIC :: Write => Write_DG2D

  END TYPE DG2D



  TYPE, PUBLIC :: DG3D
    TYPE(MappedScalar3D), PUBLIC :: solution
    TYPE(MappedVector3D), PUBLIC :: solutionGradient
    TYPE(MappedVector3D), PUBLIC :: flux
    TYPE(MappedScalar3D), PUBLIC :: source
    TYPE(MappedScalar3D), PUBLIC :: fluxDivergence
    TYPE(MPILayer), PUBLIC :: decomp
    TYPE(Mesh3D), PUBLIC :: mesh
    TYPE(SEMHex), PUBLIC :: geometry
    TYPE(Metadata), ALLOCATABLE, PUBLIC :: solutionMetaData(:)

    ! Work arrays
    TYPE(MappedScalar3D), PRIVATE :: workScalar
    TYPE(MappedVector3D), PRIVATE :: workVector
    TYPE(MappedVector3D), PRIVATE :: compFlux
    TYPE(MappedTensor3D), PRIVATE :: workTensor

    CONTAINS

      PROCEDURE, PUBLIC :: Init => Init_DG3D
      PROCEDURE, PUBLIC :: Free => Free_DG3D

      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_DG3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_DG3D

      ! Accessors
    !  PROCEDURE, PUBLIC :: GetNumberOfVars => GetNumberOfVars_DG3D
    !  PROCEDURE, PUBLIC :: GetControlDegree => GetControlDegree_DG3D
    !  PROCEDURE, PUBLIC :: GetTargetDegree => GetTargetDegree_DG3D
    !  PROCEDURE, PUBLIC :: GetControlQuadrature => GetControlQuadrature_DG3D
    !  PROCEDURE, PUBLIC :: GetTargetQuadrature => GetTargetQuadrature_DG3D
    !  PROCEDURE, PUBLIC :: GetNumberOfElement => GetNumberOfElements_DG3D
    !  PROCEDURE, PUBLIC :: GetNumberOfGlobalSides => GetNumberOfGlobalSides_DG3D
    !  PROCEDURE, PUBLIC :: GetNumberOfUniqueSides => GetNumberOfUniqueSides_DG3D
    !  PROCEDURE, PUBLIC :: GetNumberOfGlobalNodes => GetNumberOfGlobalNodes_DG3D
    !  PROCEDURE, PUBLIC :: GetNumberOfUniqueNodes => GetNumberOfUniqueNodes_DG3D

      PROCEDURE, PUBLIC :: CalculateSolutionGradient => CalculateSolutionGradient_DG3D 
      PROCEDURE, PUBLIC :: CalculateFluxDivergence => CalculateFluxDivergence_DG3D

      PROCEDURE, PUBLIC :: Read => Read_DG3D 
      PROCEDURE, PUBLIC :: Write => Write_DG3D

  END TYPE DG3D

  INTEGER, PARAMETER :: SELF_DG_BASSIREBAY = 100

CONTAINS
SUBROUTINE Init_DG2D(this,cqType,tqType,cqDegree,tqDegree,nvar,spec)
    IMPLICIT NONE
    CLASS(DG2D), INTENT(out) :: this
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nvar
    TYPE(MeshSpec), INTENT(in) :: spec

      CALL this % decomp % Init()

      ! Load Mesh
      CALL this % mesh % Load(spec,this % decomp)

      ! Create geometry from mesh
      CALL this % geometry % GenerateFromMesh(this % mesh,cqType,tqType,cqDegree,tqDegree)

      CALL this % solution % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % solutionGradient % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % flux % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % source % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % fluxDivergence % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)

      CALL this % workScalar % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % workVector % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % workTensor % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % compFlux % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)

      ALLOCATE(this % solutionMetaData(1:nvar))

  END SUBROUTINE Init_DG2D

  SUBROUTINE Free_DG2D(this)
    IMPLICIT NONE
    CLASS(DG2D), INTENT(inout) :: this

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
      DEALLOCATE(this % solutionMetaData)
      
  END SUBROUTINE Free_DG2D

  SUBROUTINE UpdateHost_DG2D(this)
    IMPLICIT NONE
    CLASS(DG2D), INTENT(inout) :: this

      CALL this % mesh % UpdateHost()
      CALL this % geometry % UpdateHost()
      CALL this % solution % UpdateHost()
      CALL this % solutionGradient % UpdateHost()
      CALL this % flux % UpdateHost()
      CALL this % source % UpdateHost()
      CALL this % fluxDivergence % UpdateHost()
      CALL this % workScalar % UpdateHost()
      CALL this % workVector % UpdateHost()
      CALL this % workTensor % UpdateHost()
      CALL this % compFlux % UpdateHost()

  END SUBROUTINE UpdateHost_DG2D

  SUBROUTINE UpdateDevice_DG2D(this)
    IMPLICIT NONE
    CLASS(DG2D), INTENT(inout) :: this

      CALL this % mesh % UpdateDevice()
      CALL this % geometry % UpdateDevice()
      CALL this % solution % UpdateDevice()
      CALL this % solutionGradient % UpdateDevice()
      CALL this % flux % UpdateDevice()
      CALL this % source % UpdateDevice()
      CALL this % fluxDivergence % UpdateDevice()
      CALL this % workScalar % UpdateDevice()
      CALL this % workVector % UpdateDevice()
      CALL this % workTensor % UpdateDevice()
      CALL this % compFlux % UpdateDevice()
      
  END SUBROUTINE UpdateDevice_DG2D 

  SUBROUTINE CalculateSolutionGradient_DG2D(this,gpuAccel)
    IMPLICIT NONE
    CLASS(DG2D), INTENT(inout) :: this
    LOGICAL, INTENT(in), OPTIONAL :: gpuAccel

    CALL this % solution % SideExchange(this % mesh,&
                                        this % decomp,&
                                        gpuAccel)

    CALL this % solution % BassiRebaySides(gpuAccel)

    CALL this % solution % Gradient(this % workTensor, &
                                    this % geometry, &
                                    this % solutionGradient, &
                                    selfWeakDGForm,gpuAccel)

  END SUBROUTINE CalculateSolutionGradient_DG2D

  SUBROUTINE CalculateFluxDivergence_DG2D(this,gpuAccel)
    IMPLICIT NONE
    CLASS(DG2D), INTENT(inout) :: this
    LOGICAL, INTENT(in), OPTIONAL :: gpuAccel

    CALL this % flux % Divergence(this % compFlux, &
                                  this % geometry, &
                                  this % fluxDivergence, &
                                  selfWeakDGForm,gpuAccel)

  END SUBROUTINE CalculateFluxDivergence_DG2D

  SUBROUTINE Write_DG2D(this,fileName)
    IMPLICIT NONE
    CLASS(DG2D), INTENT(in) :: this
    CHARACTER(*), INTENT(in) :: fileName
    ! Local 
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:4)
    INTEGER(HID_T) :: xOffset(1:5)
    INTEGER :: firstElem, nLocalElems
    INTEGER :: nGeo, nBCs

      IF(this % decomp % mpiEnabled)THEN
        CALL Open_HDF5(fileName, H5F_ACC_TRUNC_F, fileId, this % decomp % mpiComm)
      ELSE
        CALL Open_HDF5(fileName, H5F_ACC_TRUNC_F, fileId)
      ENDIF

      IF(this % decomp % mpiEnabled)THEN
        firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)+1
      ELSE
        firstElem = 1
      ENDIF

      solOffset(1:4) = (/0,0,1,firstElem/)
      CALL WriteArray_HDF5(fileId, '/solution', solOffset, this % solution % interior)

      CALL WriteArray_HDF5(fileId, '/fluxDivergence', solOffset, this % fluxDivergence % interior)

      xOffset(1:5) = (/1,0,0,1,firstElem/)
      CALL WriteArray_HDF5(fileId, '/flux', xOffset, this % flux % interior)

      CALL WriteArray_HDF5(fileId, '/solutionGradient', xOffset, this % solutionGradient % interior)

      CALL WriteArray_HDF5(fileId, '/x', xOffset, this % geometry % x % interior)

      CALL Close_HDF5(fileId)

  END SUBROUTINE Write_DG2D

  SUBROUTINE Read_DG2D(this,fileName)
    IMPLICIT NONE
    CLASS(DG2D), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: fileName
    ! Local 
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:4)
    INTEGER :: firstElem
    INTEGER :: N

      IF(this % decomp % mpiEnabled)THEN
        CALL Open_HDF5(fileName, H5F_ACC_RDWR_F, fileId, this % decomp % mpiComm)
      ELSE
        CALL Open_HDF5(fileName, H5F_ACC_RDWR_F, fileId)
      ENDIF

      CALL ReadAttribute_HDF5(fileId, 'N', N)

      IF(this % solution % N /= N)THEN
        STOP 'Error : Solution polynomial degree does not match input file'
      ENDIF

      IF(this % decomp % mpiEnabled)THEN
        firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)+1
      ELSE
        firstElem = 1
      ENDIF

      solOffset(1:4) = (/0,0,1,firstElem/)
      CALL ReadArray_HDF5(fileId, 'solution', solOffset, this % solution % interior)

      CALL Close_HDF5(fileId)

  END SUBROUTINE Read_DG2D

  SUBROUTINE Init_DG3D(this,cqType,tqType,cqDegree,tqDegree,nvar,spec)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(out) :: this
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nvar
    TYPE(MeshSpec), INTENT(in) :: spec

      CALL this % decomp % Init()

      ! Load Mesh
      CALL this % mesh % Load(spec,this % decomp)

      ! Create geometry from mesh
      CALL this % geometry % GenerateFromMesh(this % mesh,cqType,tqType,cqDegree,tqDegree)

      CALL this % solution % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % solutionGradient % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % flux % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % source % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % fluxDivergence % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)

      CALL this % workScalar % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % workVector % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % workTensor % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % compFlux % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)

      ALLOCATE(this % solutionMetaData(1:nvar))

  END SUBROUTINE Init_DG3D

  SUBROUTINE Free_DG3D(this)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(inout) :: this

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
      DEALLOCATE(this % solutionMetaData)
      
  END SUBROUTINE Free_DG3D

  SUBROUTINE UpdateHost_DG3D(this)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(inout) :: this

      CALL this % mesh % UpdateHost()
      CALL this % geometry % UpdateHost()
      CALL this % solution % UpdateHost()
      CALL this % solutionGradient % UpdateHost()
      CALL this % flux % UpdateHost()
      CALL this % source % UpdateHost()
      CALL this % fluxDivergence % UpdateHost()
      CALL this % workScalar % UpdateHost()
      CALL this % workVector % UpdateHost()
      CALL this % workTensor % UpdateHost()
      CALL this % compFlux % UpdateHost()

  END SUBROUTINE UpdateHost_DG3D

  SUBROUTINE UpdateDevice_DG3D(this)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(inout) :: this

      CALL this % mesh % UpdateDevice()
      CALL this % geometry % UpdateDevice()
      CALL this % solution % UpdateDevice()
      CALL this % solutionGradient % UpdateDevice()
      CALL this % flux % UpdateDevice()
      CALL this % source % UpdateDevice()
      CALL this % fluxDivergence % UpdateDevice()
      CALL this % workScalar % UpdateDevice()
      CALL this % workVector % UpdateDevice()
      CALL this % workTensor % UpdateDevice()
      CALL this % compFlux % UpdateDevice()
      
  END SUBROUTINE UpdateDevice_DG3D 

  SUBROUTINE CalculateSolutionGradient_DG3D(this,gpuAccel)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(inout) :: this
    LOGICAL, INTENT(in), OPTIONAL :: gpuAccel

    CALL this % solution % SideExchange(this % mesh,&
                                        this % decomp,&
                                        gpuAccel)

    CALL this % solution % BassiRebaySides(gpuAccel)

    CALL this % solution % Gradient(this % workTensor, &
                                    this % geometry, &
                                    this % solutionGradient, &
                                    selfWeakDGForm,gpuAccel)

  END SUBROUTINE CalculateSolutionGradient_DG3D

  SUBROUTINE CalculateFluxDivergence_DG3D(this,gpuAccel)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(inout) :: this
    LOGICAL, INTENT(in), OPTIONAL :: gpuAccel

    CALL this % flux % Divergence(this % compFlux, &
                                  this % geometry, &
                                  this % fluxDivergence, &
                                  selfWeakDGForm,gpuAccel)

  END SUBROUTINE CalculateFluxDivergence_DG3D

  SUBROUTINE Write_DG3D(this,fileName)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(in) :: this
    CHARACTER(*), INTENT(in) :: fileName
    ! Local 
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:5)
    INTEGER(HID_T) :: xOffset(1:6)
    INTEGER :: firstElem, nLocalElems
    INTEGER :: nGeo, nBCs

      IF(this % decomp % mpiEnabled)THEN
        CALL Open_HDF5(fileName, H5F_ACC_TRUNC_F, fileId, this % decomp % mpiComm)
      ELSE
        CALL Open_HDF5(fileName, H5F_ACC_TRUNC_F, fileId)
      ENDIF

      IF(this % decomp % mpiEnabled)THEN
        firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)+1
      ELSE
        firstElem = 1
      ENDIF

      solOffset(1:5) = (/0,0,0,1,firstElem/)
      CALL WriteArray_HDF5(fileId, '/solution', solOffset, this % solution % interior)

      CALL WriteArray_HDF5(fileId, '/fluxDivergence', solOffset, this % fluxDivergence % interior)

      xOffset(1:6) = (/1,0,0,0,1,firstElem/)
      CALL WriteArray_HDF5(fileId, '/flux', xOffset, this % flux % interior)

      CALL WriteArray_HDF5(fileId, '/solutionGradient', xOffset, this % solutionGradient % interior)


      CALL WriteArray_HDF5(fileId, '/x', xOffset, this % geometry % x % interior)

      CALL Close_HDF5(fileId)

  END SUBROUTINE Write_DG3D

  SUBROUTINE Read_DG3D(this,fileName)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: fileName
    ! Local 
    INTEGER(HID_T) :: fileId
    INTEGER(HID_T) :: solOffset(1:5)
    INTEGER :: firstElem
    INTEGER :: N

      IF(this % decomp % mpiEnabled)THEN
        CALL Open_HDF5(fileName, H5F_ACC_RDWR_F, fileId, this % decomp % mpiComm)
      ELSE
        CALL Open_HDF5(fileName, H5F_ACC_RDWR_F, fileId)
      ENDIF

      CALL ReadAttribute_HDF5(fileId, 'N', N)

      IF(this % solution % N /= N)THEN
        STOP 'Error : Solution polynomial degree does not match input file'
      ENDIF

      IF(this % decomp % mpiEnabled)THEN
        firstElem = this % decomp % offsetElem % hostData(this % decomp % rankId)+1
      ELSE
        firstElem = 1
      ENDIF

      solOffset(1:5) = (/0,0,0,1,firstElem/)
      CALL ReadArray_HDF5(fileId, 'solution', solOffset, this % solution % interior)

      CALL Close_HDF5(fileId)

  END SUBROUTINE Read_DG3D

END MODULE SELF_DG 
