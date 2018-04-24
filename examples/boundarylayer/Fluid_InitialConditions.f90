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
    INTEGER :: IFace, IFace2, e1, e2, s1, p2


    DO IFace = 1, myDGSEM % extComm % nBoundaries

      IFace2 = myDGSEM % extComm % boundaryIDs( IFace )
      e1    = myDGSEM % mesh % Faces % elementIDs(1,iFace2)
      s1    = myDGSEM % mesh % Faces % elementSides(1,iFace2)
      e2    = myDGSEM % mesh % Faces % elementIDs(2,iFace2)
      p2    = myDGSEM % extComm % extProcIDs( iFace )

      IF( e2 < 0 .AND. p2 == myeu % extComm % myRank )THEN

        IF( ABS(myDGSEM % mesh % elements % nHat(3, 0, 0, s1,e1) ) >&
            SQRT( myDGSEM % mesh % elements % nHat(1, 0, 0, s1,e1)**2 + &
                  myDGSEM % mesh % elements % nHat(2, 0, 0, s1,e1)**2 ) ) THEN ! Top or bottom surface

           IF( myDGSEM % mesh % elements % nHat(3, 0, 0, s1,e1) > 0.0_prec ) THEN ! Top

              myDGSEM % mesh % faces % elementIDs(2,iFace2) = PRESCRIBED
              myDGSEM % state % prescribedState(:,:,1,iFace) = myDGSEM % static % boundarySolution(:,:,4,s1,e1)*myDGSEM % params % v0

           ELSE ! Bottom

              myDGSEM % mesh % faces % elementIDs(2,iFace2) = NO_NORMAL_FLOW

           ENDIF

        ENDIF
        myDGSEM % mesh % faces % elementIDs(2,iFace2) = NO_NORMAL_FLOW

      ENDIF

    ENDDO


  END SUBROUTINE ResetBoundaryConditions
!
  SUBROUTINE InitialCondition( myDGSEM )

    TYPE( Fluid ), INTENT(inout) :: myDGSEM
    ! Local
    INTEGER    :: i, j, k, iEl, IFace
    REAL(prec) :: x, y, z, r, hScale, Lx, Ly, H, T, Tbar


    Lx = myDGSEM % params % xScale
    Ly = myDGSEM % params % yScale
    H  = myDGSEM % params % zScale

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

            myDGSEM % state % solution(i,j,k,1,iEl) = myDGSEM % static % solution(i,j,k,4,iEl)*myDGSEM % params % v0
            ! Drag profile
            myDGSEM % sourceTerms % drag(i,j,k,iEl) = myDGSEM % params % Cd*( 1.0_prec - tanh( z/myDGSEM % params % dragScale ) )
#ifdef PASSIVE_TRACERS
            myDGSEM % static % solution(i,j,k,6,iEl) = myDGSEM % static % solution(i,j,k,4,iEl)*( 1.0_prec - tanh( z/myDGSEM % params % dragScale ) )
#endif

          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA
    myDGSEM % state % solution_dev = myDGSEM % state % solution
#endif

    !$OMP PARALLEL
    CALL myDGSEM % EquationOfState( )
    !$OMP END PARALLEL

#ifdef HAVE_CUDA
    myDGSEM % state % solution = myDGSEM % state % solution_dev
#endif

  END SUBROUTINE InitialCondition
!
END PROGRAM Fluid_InitialConditions
