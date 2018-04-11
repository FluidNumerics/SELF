! Voxels_Class.f90
!
! Copyright 2018 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE Voxels_Class

USE ModelPrecision
USE Lagrange_Class
USE HexMesh_Class
USE HexElements_Class

IMPLICIT NONE

  TYPE CharWrapper
    CHARACTER(20) :: text
  END TYPE CharWrapper

  TYPE Voxels
    INTEGER    :: nX, nY, nZ, nVoxels, nVars
    REAL(prec), ALLOCATABLE :: x(:), y(:), z(:) 
    INTEGER, ALLOCATABLE    :: elementIDs(:) 
    REAL(prec), ALLOCATABLE :: compCoords(:,:) 

    REAL(prec), ALLOCATABLE        :: variable(:,:,:,:)
    TYPE(CharWrapper), ALLOCATABLE :: variableName(:)


    CONTAINS

    PROCEDURE :: Build => Build_Voxels  
    PROCEDURE :: Trash => Trash_Voxels

    PROCEDURE :: MapMeshToVoxels
 
  END TYPE Voxels
 

CONTAINS

 SUBROUTINE Build_Voxels( myVoxels, nVars, nX, nY, nZ, x0, y0, z0, Lx, Ly, Lz )

   IMPLICIT NONE
   CLASS( Voxels ), INTENT(out) :: myVoxels
   INTEGER, INTENT(in)          :: nVars
   INTEGER, INTENT(in)          :: nX, nY, nZ
   REAL(prec), INTENT(in)       :: x0, y0, z0
   REAL(prec), INTENT(in)       :: Lx, Ly, Lz
   ! Local
   REAL(prec) :: xb(0:nX-1), yb(0:nY-1), zb(0:nZ-1)
   INTEGER    :: nEl

     myVoxels % nX = nX
     myVoxels % nY = nY
     myVoxels % nZ = nZ
     myVoxels % nVars = nVars

     nEl = nX*nY*nZ
     myVoxels % nVoxels = nEl

     ALLOCATE( myVoxels % x(0:nX-1), &
               myVoxels % y(0:nY-1), &
               myVoxels % z(0:nZ-1), &
               myVoxels % elementIDs(1:nEl), &
               myVoxels % compCoords(1:3,1:nEl), &
               myVoxels % variable(0:nX-1,0:nY-1,0:nZ-1,1:nVars), &
               myVoxels % variableName(1:nVars) )

     myVoxels % x = UniformPoints( x0, x0+Lx, nX )
     myVoxels % y = UniformPoints( y0, y0+Ly, nY )
     myVoxels % z = UniformPoints( z0, z0+Lz, nZ )

     myVoxels % elementIDs = 0
     myVoxels % compCoords = 0.0_prec 

 END SUBROUTINE Build_Voxels
!
 SUBROUTINE Trash_Voxels( myVoxels )
   IMPLICIT NONE
   CLASS( Voxels ), INTENT(inout) :: myVoxels

     DEALLOCATE( myVoxels % x, &
                 myVoxels % y, &
                 myVoxels % z, &
                 myVoxels % elementIDs, &
                 myVoxels % compCoords, &
                 myVoxels % variable, &
                 myVoxels % variableName )

 END SUBROUTINE Trash_Voxels
!
 SUBROUTINE MapMeshToVoxels( myVoxels, mesh, interp )
   IMPLICIT NONE
   CLASS( Voxels ), INTENT(inout) :: myVoxels
   TYPE( HexMesh ), INTENT(in)    :: mesh
   TYPE( Lagrange ), INTENT(in)   :: interp
   ! Local
   INTEGER :: iX, iY, iZ, iVoxel
   REAL(prec) :: coordinates(1:3,1:myVoxels % nVoxels)

     DO iZ = 0, myVoxels % nX-1
       DO iY = 0, myVoxels % nY-1
         DO iX = 0, myVoxels % nX-1
 
            iVoxel = iX + 1 + myVoxels % nX*( iY + iZ*myVoxels % nY )
            coordinates(1,iVoxel) = myVoxels % x(iX) 
            coordinates(2,iVoxel) = myVoxels % y(iY) 
            coordinates(3,iVoxel) = myVoxels % z(iZ) 

          ENDDO
        ENDDO
      ENDDO

      CALL mesh % elements % CalculateComputationalCoordinates( interp, coordinates, myVoxels % compCoords, myVoxels % elementIDs, myVoxels % nVoxels )

 END SUBROUTINE MapMeshToVoxels

END MODULE Voxels_Class
