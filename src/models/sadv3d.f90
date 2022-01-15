PROGRAM sadv3d
        
USE SELF_Advection3D
USE SELF_Mesh
USE FEQParse
USE SELF_Constants

  TYPE( Advection3D ) :: model

  CALL model % InitCLI( )

  ! Now that we have CLI variables,
  ! we need to determine what to do
  !
  ! When convergence check is requested, we want to advance the model multiple
  ! times, gradually increasing N. Additionally, we want to calculate the
  ! difference between the exact solution and the numerical solution to obtain
  ! the max(abs(error)).

  CALL model % ModelExecute() 


END PROGRAM sadv3d
