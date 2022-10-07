
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_CGModel2D

  USE SELF_Model
  USE SELF_Model2D

  TYPE, EXTENDS(Model2D) :: CGModel2D
  !! A class that contains additional methods for solving
  !! 
  !!  s_t + div( flux ) = source
  !!
  !! or 
  !!
  !!  div( flux ) = source
  !!
  !! using the Continuous Galerkin Spectral Element Method
  !!
  
    CONTAINS

    ! Overridden methods
    PROCEDURE :: Init => Init_CGModel2D

    PROCEDURE :: CalculateTendency => CalculateTendency_CGModel2D
    PROCEDURE :: CalculateFluxDivergence => CalculateFluxDivergence_CGModel2D

  END TYPE CGModel2D
END MODULE SELF_CGModel2D
