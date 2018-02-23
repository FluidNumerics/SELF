PROGRAM MappedGradient_Testing

USE ModelPrecision
USE CommonRoutines
USE NodalDG_Class
USE NodalDGSolution_3D_Class
USE HexMesh_Class
USE ModelParameters_Class
USE TopographicShapes


IMPLICIT NONE

  TYPE( ModelParameters )    :: params
  TYPE( NodalDG )            :: dgStorage
  TYPE( NodalDGSolution_3D ) :: state
  TYPE( HexMesh )            :: mesh

  LOGICAL      :: readSuccess
  INTEGER      :: i, j, k, iEq, iEl, iFace, e1, e2, s1, bID, fUnit
  REAL(prec)   :: x, y, z
  CHARACTER(5) :: zoneID

    CALL params % Build( readSuccess )

    CALL dGStorage % Build( UniformPoints(-1.0_prec,1.0_prec,params % nPlot), &
                            params % polyDeg, params % nPlot, GAUSS )

    IF( params % topographicShape == Gaussian )THEN
       TopographicShape => GaussianHill
    ELSE
       TopographicShape => DefaultTopography
    ENDIF

    CALL mesh % ConstructStructuredMesh( dgStorage % interp, &
                                         params % nXelem, &
                                         params % nYelem, &
                                         params % nZelem, &
                                         .TRUE.  )

    CALL mesh % ScaleTheMesh( dgStorage % interp, &
                              params % xScale, &
                              params % yScale, &
                              params % zScale )

    CALL state % Build( params % polyDeg, 1, &
                        mesh % elements % nElements, &
                        mesh % NumberOfBoundaryFaces( ) )
   

    DO iEl = 1, state % nElements
      DO iEq = 1, state % nEquations
        DO k = 0, state % N
          DO j = 0, state % N
            DO i = 0, state % N

              x = mesh % elements % x(i,j,k,1,iEl)
              y = mesh % elements % x(i,j,k,2,iEl)
              z = mesh % elements % x(i,j,k,3,iEl)

              state % solution(i,j,k,iEq,iEl) = exp( -(z-params % zScale)/(0.1_prec*params % zScale) ) 

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CALL state % Calculate_Solution_At_Boundaries( dgStorage ) 

    
!    DO iFace = 1, mesh % faces % nFaces
!
!      e1 = mesh % faces % elementIDs(1,iFace)
!      e2 = mesh % faces % elementIDs(2,iFace)
!      s1 = mesh % faces % elementSides(1,iFace)
!      bID = ABS(mesh % faces % boundaryID(iFace))
!
!      IF( e2 < 0 )THEN
!
!        DO iEq = 1, state % nEquations
!          DO j = 0, state % N
!            DO i = 0, state % N
!          
!              state % externalState(i,j,iEq,bID) = state % boundarySolution(i,j,iEq,s1,e1)
!
!            ENDDO
!          ENDDO
!        ENDDO
!
!      ENDIF
!
!    ENDDO

    CALL state % Mapped_BassiRebay_Gradient( dgStorage, mesh )

    PRINT*, 'Max dp/dx',MAXVAL( state % solutionGradient(1,:,:,:,:,:) )
    PRINT*, 'Max dp/dy',MAXVAL( state % solutionGradient(2,:,:,:,:,:) )
    PRINT*, 'Max dp/dz',MAXVAL( state % solutionGradient(3,:,:,:,:,:) )


    CALL mesh % WriteTecplot( )

    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= 'State.tec', &
      FORM='formatted', &
      STATUS='replace')

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "p", "dpdx", "dpdy", "dpdz"'

    DO iEl = 1, mesh % elements % nElements

      WRITE(zoneID,'(I5.5)') mesh % elements % elementID(iEl)
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',params % polyDeg+1,&
                                                 ', J=',params % polyDeg+1,&
                                                 ', K=',params % polyDeg+1,',F=POINT'

      DO k = 0, params % polyDeg
        DO j = 0, params % polyDeg
          DO i = 0, params % polyDeg

            WRITE(fUnit,'(17(E15.7,1x))') mesh % elements % x(i,j,k,1,iEl),&
                                          mesh % elements % x(i,j,k,2,iEl),&
                                          mesh % elements % x(i,j,k,3,iEl),&
                                          state % solution(i,j,k,1,iEl), &
                                          state % solutionGradient(1,i,j,k,1,iEl), &
                                          state % solutionGradient(2,i,j,k,1,iEl), &
                                          state % solutionGradient(3,i,j,k,1,iEl)
                                          
          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)

    CALL state % Trash( )
    CALL mesh % Trash( )
    CALL dgStorage % Trash( )
    CALL params % Trash( )

END PROGRAM MappedGradient_Testing
