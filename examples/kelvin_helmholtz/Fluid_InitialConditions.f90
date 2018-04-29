! Fluid_InitialConditions.f90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM Fluid_InitialConditions

  USE ModelPrecision
  USE Fluid_Class

  IMPLICIT NONE

  TYPE( Fluid ) :: myeu
  LOGICAL       :: setupSuccess
  CHARACTER(4)  :: rankChar

  CALL myeu % Build( setupSuccess )

  IF( SetupSuccess )THEN

    CALL InitialCondition( myeu )

    CALL ResetBoundaryConditions( myeu )

    CALL myeu % WritePickup( )

    CALL myeu % WriteTecplot( )

    ! Before we write the mesh to file again, we need to "unscale" the mesh so that, upon running the
    ! integrator, the mesh scaling is not applied a second time
    CALL myeu % mesh % ScaleTheMesh( myeu % dgStorage % interp, &
                                     1.0_prec/myeu % params % xScale, &
                                     1.0_prec/myeu % params % yScale, &
                                     1.0_prec/myeu % params % zScale )
    WRITE( rankChar, '(I4.4)' ) myeu % extComm % myRank
    CALL myeu % mesh % WriteSELFMeshFile( TRIM(myeu % params % SELFMeshFile)//'.'//rankChar )

    CALL myeu % Trash( )

  ENDIF


CONTAINS
  SUBROUTINE ResetBoundaryConditions( myDGSEM )

    TYPE( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    INTEGER :: IFace, IFace2, e1, e2, s1, p2, i, j, iEq


    DO IFace = 1, myDGSEM % extComm % nBoundaries

      IFace2 = myDGSEM % extComm % boundaryIDs( IFace )
      e1    = myDGSEM % mesh % Faces % elementIDs(1,iFace2)
      s1    = myDGSEM % mesh % Faces % elementSides(1,iFace2)
      e2    = myDGSEM % mesh % Faces % elementIDs(2,iFace2)
      p2    = myDGSEM % extComm % extProcIDs( iFace )

      IF( e2 < 0 .AND. p2 == myeu % extComm % myRank )THEN

        myDGSEM % mesh % faces % elementIDs(2,iFace2) = PRESCRIBED

        DO iEq = 1, myDGSEM % state % nEquations
          DO i = 0, myDGSEM % params % polyDeg
            DO j = 0, myDGSEM % params % polyDeg

              myDGSEM % state % prescribedState(i,j,iEq,iFace) = myDGSEM % state % boundarySolution(i,j,iEq,s1,e1)

            ENDDO
          ENDDO
        ENDDO

      ENDIF

    ENDDO


  END SUBROUTINE ResetBoundaryConditions
!
  SUBROUTINE InitialCondition( myDGSEM )

    TYPE( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    INTEGER    :: i, j, k, iEl
    REAL(prec) :: x, y, z, zc, Hz, r, Lx, Ly, H, T, Tbar


    Lx = myDGSEM % params % xScale
    Ly = myDGSEM % params % yScale
    H  = myDGSEM % params % zScale

    zc = 0.5_prec*H
    Hz = 0.01_prec*H

    myDGSEM % state % solution = 0.0_prec

    !$OMP PARALLEL
    CALL myDGSEM % CalculateStaticState( ) !! CPU Kernel
    !$OMP END PARALLEL

    ! ////////////////////////////////////////////////////////////////////////////////// !
    DO iEl = 1, myDGSEM % mesh % elements % nElements
      DO k = 0, myDGSEM % params % polyDeg
        DO j = 0, myDGSEM % params % polyDeg
          DO i = 0, myDGSEM % params % polyDeg

            x = myDGSEM % mesh % elements % x(i,j,k,1,iEl)
            y = myDGSEM % mesh % elements % x(i,j,k,2,iEl)
            z = myDGSEM % mesh % elements % x(i,j,k,3,iEl)

            T = myDGSEM % params % dTdz*z

           ! IF( z > 0.1_prec*H .AND. z < 0.9_prec*H ) THEN
           !   CALL RANDOM_NUMBER( r )
           !   T = T + 2.0_prec*(r-0.5_prec)*0.25_prec
           ! ENDIF

            Tbar = myDGSEM % static % solution(i,j,k,5,iEl)/myDGSEM % static % solution(i,j,k,4,iEl)
            myDGSEM % state % solution(i,j,k,4,iEl) = -myDGSEM % static % solution(i,j,k,4,iEl)*T/(Tbar + T)

            myDGSEM % state % solution(i,j,k,1,iEl) = ( myDGSEM % state % solution(i,j,k,4,iEl) + &
                                                        myDGSEM % static % solution(i,j,k,4,iEl) )*myDGSEM % params % v0*&
                                                      tanh( (z-zc)/Hz )

          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateDevice( )
#endif

    !$OMP PARALLEL
    CALL myDGSEM % EquationOfState( )
    !$OMP END PARALLEL
    CALL myDGSEM % state % Calculate_Solution_At_Boundaries( myDGSEM % dgStorage )
    CALL myDGSEM % static % Calculate_Solution_At_Boundaries( myDGSEM % dgStorage )

#ifdef HAVE_CUDA
    CALL myDGSEM % state % UpdateHost( )
    CALL myDGSEM % static % UpdateHost( )
#endif

  END SUBROUTINE InitialCondition
!
END PROGRAM Fluid_InitialConditions
