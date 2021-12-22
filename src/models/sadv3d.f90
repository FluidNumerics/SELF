PROGRAM sadv3d
        
USE SELF_Advection3D
USE SELF_Mesh
USE FEQParse
USE SELF_Constants

  TYPE( Advection3D ) :: model
  INTEGER :: nDumps
  INTEGER :: i
  REAL(prec) :: endTime

  CALL model % InitFromCLI( )

  CALL model %  WriteTecplot()

  nDumps = INT(( model % endTime - model % initialTime )/( model % outputInterval ) )
  DO i = 1, nDumps
  
    endTime = model % simulationTime + model % outputInterval
        
    CALL model % ForwardStep( endTime )
    CALL model % WriteTecplot()
    CALL model % WritePickup()

  ENDDO

  CALL model % Free()

END PROGRAM sadv3d
