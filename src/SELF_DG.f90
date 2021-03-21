MODULE SELF_DG

! How someone should use...
!  > Create the object
!  > Loop over time
!     > Fill in source
!     > Fill in flux (Usually a Riemmann Solver); if 2nd order operators are included, bassi-rebay flux routine is available to
!     support gradient calculation.
!     > ForwardStep
!  > 
  TYPE, PUBLIC :: DG1D

    TYPE(Mapped_Scalar1D), PUBLIC :: solution
    TYPE(Mapped_Scalar1D), PUBLIC :: solutionDerivative
    TYPE(Mapped_Scalar1D), PUBLIC :: flux
    TYPE(Mapped_Scalar1D), PUBLIC :: source
    TYPE(Mapped_Scalar1D), PUBLIC :: tendency
    TYPE(Mesh1D), PUBLIC :: mesh
    TYPE(Geometry1D), PUBLIC :: geometry

    CONTAINS

      PROCEDURE, PUBLIC :: Init => Init_DG1D
      PROCEDURE, PUBLIC :: Free => Free_DG1D

      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_DG1D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_DG1D

      ! Accessors
      PROCEDURE, PUBLIC :: GetNumberOfVars => GetNumberOfVars_DG1D
      PROCEDURE, PUBLIC :: GetControlDegree => GetControlDegree_DG1D
      PROCEDURE, PUBLIC :: GetTargetDegree => GetTargetDegree_DG1D
      PROCEDURE, PUBLIC :: GetControlQuadrature => GetControlQuadrature_DG1D
      PROCEDURE, PUBLIC :: GetTargetQuadrature => GetTargetQuadrature_DG1D
      PROCEDURE, PUBLIC :: GetNumberOfElement => GetNumberOfElements_DG1D
      PROCEDURE, PUBLIC :: GetNumberOfGlobalSides => GetNumberOfGlobalSides_DG1D
      PROCEDURE, PUBLIC :: GetNumberOfUniqueSides => GetNumberOfUniqueSides_DG1D
      PROCEDURE, PUBLIC :: GetNumberOfGlobalNodes => GetNumberOfGlobalNodes_DG1D
      PROCEDURE, PUBLIC :: GetNumberOfUniqueNodes => GetNumberOfUniqueNodes_DG1D

      PROCEDURE, PUBLIC :: SolutionDerivativeBR => SolutionDerivativeBR_DG1D ! Use Bassi-Rebay flux to calculate weak derivative

      PROCEDURE, PUBLIC :: CalculateTendency => CalculateTendency_DG1D
      
      PROCEDURE, PUBLIC :: ForwardStepEuler => ForwardStepEuler_DG1D
      PROCEDURE, PUBLIC :: ForwardStepRK3 => ForwardStepRK3_DG1D
      PROCEDURE, PUBLIC :: ForwardStepRK3 => ForwardStepRK3_DG1D

      PROCEDURE, PUBLIC :: Load => Load_DG1D ! Load from file
      PROCEDURE, PUBLIC :: Drop => Drop_DG1D ! Drop to file


  END TYPE DG1D

  TYPE, PUBLIC :: DG2D
    TYPE(Mapped_Scalar2D), PUBLIC :: solution
    TYPE(Mapped_Vector2D), PUBLIC :: solutionGradient
    TYPE(Mapped_Vector2D), PUBLIC :: flux
    TYPE(Mapped_Scalar2D), PUBLIC :: source
    TYPE(Mapped_Scalar2D), PUBLIC :: tendency
    TYPE(Mesh2D), PUBLIC :: mesh
    TYPE(SEMQuad), PUBLIC :: geometry

    CONTAINS

      PROCEDURE, PUBLIC :: Init => Init_DG2D
      PROCEDURE, PUBLIC :: Free => Free_DG2D

      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_DG2D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_DG2D

      ! Accessors
      PROCEDURE, PUBLIC :: GetNumberOfVars => GetNumberOfVars_DG2D
      PROCEDURE, PUBLIC :: GetControlDegree => GetControlDegree_DG2D
      PROCEDURE, PUBLIC :: GetTargetDegree => GetTargetDegree_DG2D
      PROCEDURE, PUBLIC :: GetControlQuadrature => GetControlQuadrature_DG2D
      PROCEDURE, PUBLIC :: GetTargetQuadrature => GetTargetQuadrature_DG2D
      PROCEDURE, PUBLIC :: GetNumberOfElement => GetNumberOfElements_DG2D
      PROCEDURE, PUBLIC :: GetNumberOfGlobalSides => GetNumberOfGlobalSides_DG2D
      PROCEDURE, PUBLIC :: GetNumberOfUniqueSides => GetNumberOfUniqueSides_DG2D
      PROCEDURE, PUBLIC :: GetNumberOfGlobalNodes => GetNumberOfGlobalNodes_DG2D
      PROCEDURE, PUBLIC :: GetNumberOfUniqueNodes => GetNumberOfUniqueNodes_DG2D

      PROCEDURE, PUBLIC :: SolutionDerivativeBR => SolutionDerivativeBR_DG2D ! Use Bassi-Rebay flux to calculate weak derivative

      PROCEDURE, PUBLIC :: CalculateTendency => CalculateTendency_DG2D
      
      PROCEDURE, PUBLIC :: ForwardStepEuler => ForwardStepEuler_DG2D
      PROCEDURE, PUBLIC :: ForwardStepRK3 => ForwardStepRK3_DG2D
      PROCEDURE, PUBLIC :: ForwardStepRK3 => ForwardStepRK3_DG2D

      PROCEDURE, PUBLIC :: Load => Load_DG2D ! Load from file
      PROCEDURE, PUBLIC :: Drop => Drop_DG2D ! Drop to file

  END TYPE DG2D

  TYPE, PUBLIC :: DG3D
    TYPE(Mapped_Scalar3D), PUBLIC :: solution
    TYPE(Mapped_Vector3D), PUBLIC :: solutionGradient
    TYPE(Mapped_Vector3D), PUBLIC :: flux
    TYPE(Mapped_Scalar3D), PUBLIC :: source
    TYPE(Mapped_Scalar3D), PUBLIC :: tendency
    TYPE(Mesh3D), PUBLIC :: mesh
    TYPE(SEMQuad), PUBLIC :: geometry

    CONTAINS

      PROCEDURE, PUBLIC :: Init => Init_DG3D
      PROCEDURE, PUBLIC :: Free => Free_DG3D

      PROCEDURE, PUBLIC :: UpdateHost => UpdateHost_DG3D
      PROCEDURE, PUBLIC :: UpdateDevice => UpdateDevice_DG3D

      ! Accessors
      PROCEDURE, PUBLIC :: GetNumberOfVars => GetNumberOfVars_DG3D
      PROCEDURE, PUBLIC :: GetControlDegree => GetControlDegree_DG3D
      PROCEDURE, PUBLIC :: GetTargetDegree => GetTargetDegree_DG3D
      PROCEDURE, PUBLIC :: GetControlQuadrature => GetControlQuadrature_DG3D
      PROCEDURE, PUBLIC :: GetTargetQuadrature => GetTargetQuadrature_DG3D
      PROCEDURE, PUBLIC :: GetNumberOfElement => GetNumberOfElements_DG3D
      PROCEDURE, PUBLIC :: GetNumberOfGlobalSides => GetNumberOfGlobalSides_DG3D
      PROCEDURE, PUBLIC :: GetNumberOfUniqueSides => GetNumberOfUniqueSides_DG3D
      PROCEDURE, PUBLIC :: GetNumberOfGlobalNodes => GetNumberOfGlobalNodes_DG3D
      PROCEDURE, PUBLIC :: GetNumberOfUniqueNodes => GetNumberOfUniqueNodes_DG3D

      PROCEDURE, PUBLIC :: SolutionDerivativeBR => SolutionDerivativeBR_DG3D ! Use Bassi-Rebay flux to calculate weak derivative

      PROCEDURE, PUBLIC :: CalculateTendency => CalculateTendency_DG3D
      
      PROCEDURE, PUBLIC :: ForwardStepEuler => ForwardStepEuler_DG3D
      PROCEDURE, PUBLIC :: ForwardStepRK3 => ForwardStepRK3_DG3D
      PROCEDURE, PUBLIC :: ForwardStepRK3 => ForwardStepRK3_DG3D

      PROCEDURE, PUBLIC :: Load => Load_DG3D ! Load from file
      PROCEDURE, PUBLIC :: Drop => Drop_DG3D ! Drop to file

  END TYPE DG3D

CONTAINS
END MODULE SELF_DG
