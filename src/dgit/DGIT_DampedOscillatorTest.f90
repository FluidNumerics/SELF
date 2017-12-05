PROGRAM DGIT_DampedOscillatorTest


! src/common/
 USE ModelPrecision
! src/spectralops/
 USE NodalStorage_Class
 

 INTEGER, PARAMETER    :: N = 3
 REAL(prec), PARAMETER :: dampFrequency = 0.1_prec
 REAL(prec), PARAMETER :: oscFrequency  = 0.0_prec
 REAL(prec), PARAMETER :: dt = 1.0_prec


 TYPE( DGITIntegrator ) :: integrator
 REAL(prec)             :: solution(1:2)
 REAL(prec)             :: rhs(1:2,0:N), subStates(1:2,0:N), dsdt(1:2,0:N)

  CALL integrator % Build( N, dt, GAUSS )


  CALL integrator % Trash( )


CONTAINS

 SUBROUTINE CalculateTendency( s, tendency )
   IMPLICIT NONE
   REAL(prec), INTENT(in)  :: s(1:2)
   REAL(prec), INTENT(out) :: tendency(1:2)


      tendency(1) = -dampFrequency*s(1) + oscFrequency*s(2)
      tendency(2) = -dampFrequency*s(2) - oscFrequency*s(1)


 END SUBROUTINE CalculateTendency 


END PROGRAM DGIT_DampedOscillatorTest
