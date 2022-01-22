PROGRAM sadv2d
        
USE SELF_Advection2D
USE SELF_Mesh
USE FEQParse
USE SELF_Constants

  TYPE( Advection2D ) :: model
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

  ENDDO

  CALL model % Free()

 

END PROGRAM sadv2d
