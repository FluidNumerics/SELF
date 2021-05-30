!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_DG

!USE SELF_Metadata
USE SELF_Mesh
USE SELF_MappedData

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


  TYPE, PUBLIC :: DG3D
    TYPE(MappedScalar3D), PUBLIC :: solution
    TYPE(MappedVector3D), PUBLIC :: solutionGradient
    TYPE(MappedVector3D), PUBLIC :: flux
    TYPE(MappedScalar3D), PUBLIC :: source
    TYPE(MappedScalar3D), PUBLIC :: tendency
    TYPE(Mesh3D), PUBLIC :: mesh
    TYPE(SEMHex), PUBLIC :: geometry

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
    !  PROCEDURE, PUBLIC :: CalculateTendency => CalculateTendency_DG3D
      
    !  PROCEDURE, PUBLIC :: ForwardStepEuler => ForwardStepEuler_DG3D
    !  PROCEDURE, PUBLIC :: ForwardStepRK2 => ForwardStepRK2_DG3D
    !  PROCEDURE, PUBLIC :: ForwardStepRK3 => ForwardStepRK3_DG3D

    !  PROCEDURE, PUBLIC :: Read => Read_DG3D ! Load from file
    !  PROCEDURE, PUBLIC :: Write => Write_DG3D ! Drop to file, Options for restart file,
                                                ! tecplot file, ...

  END TYPE DG3D

  INTEGER, PARAMETER :: SELF_DG_BASSIREBAY = 100

CONTAINS

  SUBROUTINE Init_DG3D(this,cqType,tqType,cqDegree,tqDegree,nvar,spec,nRanks,myRank,mpiComm)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(out) :: this
    INTEGER,INTENT(in) :: cqType
    INTEGER,INTENT(in) :: tqType
    INTEGER,INTENT(in) :: cqDegree
    INTEGER,INTENT(in) :: tqDegree
    INTEGER,INTENT(in) :: nvar
    TYPE(MeshSpec), INTENT(in) :: spec
    INTEGER, INTENT(in) :: nRanks
    INTEGER, INTENT(in) :: myRank
    INTEGER, OPTIONAL, INTENT(in) :: mpiComm

      ! Load Mesh
      IF(PRESENT(mpiComm))THEN
        CALL this % mesh % Load(spec,nRanks,myRank,mpiComm)
      ELSE
        CALL this % mesh % Load(spec,nRanks,myRank)
      ENDIF
      ! Create geometry from mesh
      CALL this % geometry % GenerateFromMesh(this % mesh,cqType,tqType,cqDegree,tqDegree)

      CALL this % solution % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % solutionGradient % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % flux % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % source % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)
      CALL this % tendency % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)

      CALL this % workScalar % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % workVector % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % workTensor % Init(cqDegree,cqType,tqDegree,tqType,3*nVar,this % mesh % nElem)
      CALL this % compFlux % Init(cqDegree,cqType,tqDegree,tqType,nVar,this % mesh % nElem)

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
      CALL this % tendency % Free()
      CALL this % workScalar % Free()
      CALL this % workVector % Free()
      CALL this % workTensor % Free()
      CALL this % compFlux % Free()

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
      CALL this % tendency % UpdateHost()
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
      CALL this % tendency % UpdateDevice()
      CALL this % workScalar % UpdateDevice()
      CALL this % workVector % UpdateDevice()
      CALL this % workTensor % UpdateDevice()
      CALL this % compFlux % UpdateDevice()
      
  END SUBROUTINE UpdateDevice_DG3D 

  SUBROUTINE CalculateSolutionGradient_DG3D(this,gpuAccel)
    IMPLICIT NONE
    CLASS(DG3D), INTENT(inout) :: this
    LOGICAL, INTENT(in), OPTIONAL :: gpuAccel

    CALL this % solution % SideExchange(this % mesh,gpuAccel)

    CALL this % solution % BassiRebaySides(gpuAccel)

    CALL this % solution % Gradient(this % workTensor, &
                                    this % geometry, &
                                    this % solutionGradient, &
                                    selfWeakDGForm,gpuAccel)

  END SUBROUTINE CalculateSolutionGradient_DG3D

END MODULE SELF_DG 
