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
  INTEGER      :: i, j, k, iEq, iEl, iFace, e1, e2, s1, bID, fUnit, idir
  REAL(prec)   :: x, y, z
  CHARACTER(5) :: zoneID

    CALL params % Build( readSuccess )

    CALL dGStorage % Build( UniformPoints(-1.0_prec,1.0_prec,params % nPlot), &
                            params % polyDeg, params % nPlot, GAUSS )

    ! This loads in the mesh from the "pc-mesh file" and sets up the device arrays for the mesh
    CALL mesh % ReadSELFMeshFile( TRIM(params % SELFMeshFile)//'.0000' )

    ! Multiply the mesh positions to scale the size of the mesh
    CALL mesh % ScaleTheMesh( dgStorage % interp, &
                              params % xScale, &
                              params % yScale, &
                              params % zScale )

    CALL state % Build( params % polyDeg, 1, &
                        mesh % elements % nElements, &
                        mesh % NumberOfBoundaryFaces( ) )

    DO idir = 1, 3   
    DO iEl = 1, state % nElements
      DO iEq = 1, state % nEquations
        DO k = 0, state % N
          DO j = 0, state % N
            DO i = 0, state % N

              state % flux(1,i,j,k,iEq,iEl) = mesh % elements % Ja(i,j,k,idir,1,iEl)
              state % flux(2,i,j,k,iEq,iEl) = mesh % elements % Ja(i,j,k,idir,2,iEl)
              state % flux(3,i,j,k,iEq,iEl) = mesh % elements % Ja(i,j,k,idir,3,iEl)

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CALL state % Calculate_Strong_Flux_Divergence( dgStorage )
   
    PRINT*, MAXVAL( ABS( state % fluxDivergence ) ) 
    ENDDO


    DO iEl = 1, state % nElements
      DO iEq = 1, state % nEquations
        DO k = 0, state % N
          DO j = 0, state % N
            DO i = 0, state % N

              x = mesh % elements % x(i,j,k,1,iEl)
              y = mesh % elements % x(i,j,k,2,iEl)
              z = mesh % elements % x(i,j,k,3,iEl)

              !state % solution(i,j,k,iEq,iEl) = 0.001*exp( -(z-params % zScale)/(0.1_prec*params % zScale) ) 
              state % solution(i,j,k,iEq,iEl) =  1000.0_prec!z*(10000.0_prec)/(params % zScale ) 

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CALL state % Calculate_Solution_At_Boundaries( dgStorage ) 
    DO iFace = 1, mesh % faces % nFaces

      e1 = mesh % faces % elementIDs(1,iFace)
      e2 = mesh % faces % elementIDs(2,iFace)
      s1 = mesh % faces % elementSides(1,iFace)
      bID = ABS(mesh % faces % boundaryID(iFace))

      IF( e2 < 0 )THEN

        DO iEq = 1, state % nEquations
          DO j = 0, state % N
            DO i = 0, state % N
          
              state % externalState(i,j,iEq,bID) = state % boundarySolution(i,j,iEq,s1,e1)

            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ENDDO

    CALL state % Mapped_BassiRebay_Gradient( dgStorage, mesh )

    PRINT*, 'Max dp/dx',MAXVAL( ABS( state % solutionGradient(1,:,:,:,:,:) ) )
    PRINT*, 'Max dp/dy',MAXVAL( ABS( state % solutionGradient(2,:,:,:,:,:) ) )
    PRINT*, 'Max dp/dz',MAXVAL( ABS( state % solutionGradient(3,:,:,:,:,:) ) )


    CALL mesh % WriteTecplot( )

    CALL WriteToTecplot( mesh, state, dgStorage, params )

    CALL state % Trash( )
    CALL mesh % Trash( )
    CALL dgStorage % Trash( )
    CALL params % Trash( )

CONTAINS

 SUBROUTINE WriteToTecplot( mesh, state, dgStorage, params )
 
    IMPLICIT NONE
    TYPE( HexMesh ), INTENT(in)            :: mesh
    TYPE( NodalDGSolution_3D ), INTENT(in) :: state
    TYPE( NodalDG ), INTENT(in)            :: dgStorage
    TYPE( ModelParameters ), INTENT(in)    :: params
    !Local
    INTEGER     :: i, j, k, iEl
    REAL(prec)  :: x(0:params % nPlot,0:params % nPlot,0:params % nPlot,1:3,1:mesh % elements % nElements)
    REAL(prec)  :: sol(0:params % nPlot,0:params % nPlot,0:params % nPlot,1:state % nEquations,1:mesh % elements % nElements)
    REAL(prec)  :: solg(1:3,0:params % nPlot,0:params % nPlot,0:params % nPlot,1:state % nEquations, 1:mesh % elements % nElements)


    sol = ApplyInterpolationMatrix_3D_Lagrange( dgStorage % interp, &
                                                state % solution, &
                                                state % nEquations, &
                                                mesh % elements % nElements )

    x = ApplyInterpolationMatrix_3D_Lagrange( dgStorage % interp, mesh % elements % x, &
                                              3, mesh % elements % nElements )
    
    DO i = 1, 3 
      solg(i,:,:,:,:,:) = ApplyInterpolationMatrix_3D_Lagrange( dgStorage % interp, &
                                                                state % solutionGradient(i,:,:,:,:,:), &
                                                                state % nEquations, &
                                                                mesh % elements % nElements )
    ENDDO

    OPEN( UNIT=NEWUNIT(fUnit), &
      FILE= 'State.tec', &
      FORM='formatted', &
      STATUS='replace')

    WRITE(fUnit,*) 'VARIABLES = "X", "Y", "Z", "p", "dpdx", "dpdy", "dpdz"'

    DO iEl = 1, mesh % elements % nElements

      WRITE(zoneID,'(I5.5)') mesh % elements % elementID(iEl)
      WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',params % nPlot+1,&
                                                 ', J=',params % nPlot+1,&
                                                 ', K=',params % nPlot+1,',F=POINT'

      DO k = 0, params % nPlot
        DO j = 0, params % nPlot
          DO i = 0, params % nPlot

            WRITE(fUnit,'(17(E15.7,1x))') x(i,j,k,1,iEl),&
                                          x(i,j,k,2,iEl),&
                                          x(i,j,k,3,iEl),&
                                          sol(i,j,k,1,iEl), &
                                          solG(1,i,j,k,1,iEl), &
                                          solG(2,i,j,k,1,iEl), &
                                          solG(3,i,j,k,1,iEl)
                                          
          ENDDO
        ENDDO
      ENDDO

    ENDDO

    CLOSE(UNIT=fUnit)
 END SUBROUTINE WriteToTecplot

END PROGRAM MappedGradient_Testing
