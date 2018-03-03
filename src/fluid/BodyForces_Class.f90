MODULE BodyForces_CLASS

  USE ModelPrecision

  IMPLICIT NONE



  TYPE BodyForces

    REAL(prec), ALLOCATABLE         :: drag(:,:,:,:)

#ifdef HAVE_CUDA
    REAL(prec), DEVICE, ALLOCATABLE :: drag_dev(:,:,:,:)
#endif

  CONTAINS

    PROCEDURE :: Build => Build_BodyForces
    PROCEDURE :: Trash => Trash_BodyForces

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateDevice => UpdateDevice_BodyForces
    PROCEDURE :: UpdateHost   => UpdateHost_BodyForces
#endif

  END TYPE BodyForces


CONTAINS

  SUBROUTINE Build_BodyForces( myForces, N, nEquations, nElements )
    IMPLICIT NONE
    CLASS( BodyForces ), INTENT(out) :: myForces
    INTEGER, INTENT(in)              :: N, nEquations, nElements


    ALLOCATE( myForces % drag(0:N,0:N,0:N,1:nElements ) )

    myForces % drag = 0.0_prec

#ifdef HAVE_CUDA

    ALLOCATE( myForces % drag_dev(0:N,0:N,0:N,1:nElements ) )
    myForces % drag_dev = 0.0_prec

#endif


  END SUBROUTINE Build_BodyForces

  SUBROUTINE Trash_BodyForces( myForces )
    IMPLICIT NONE
    CLASS( BodyForces ), INTENT(inout) :: myForces


    DEALLOCATE( myForces % drag )

#ifdef HAVE_CUDA

    DEALLOCATE( myForces % drag_dev )

#endif
  END SUBROUTINE Trash_BodyForces

#ifdef HAVE_CUDA

  SUBROUTINE UpdateDevice_BodyForces( myForces )
    IMPLICIT NONE
    CLASS( BodyForces ), INTENT(inout) :: myForces


    myForces % drag_dev = myForces % drag


  END SUBROUTINE UpdateDevice_BodyForces

  SUBROUTINE UpdateHost_BodyForces( myForces )
    IMPLICIT NONE
    CLASS( BodyForces ), INTENT(inout) :: myForces


    myForces % drag = myForces % drag_dev


  END SUBROUTINE UpdateHost_BodyForces

#endif

END MODULE BodyForces_CLASS
