MODULE DataViz_Class

USE ModelPrecision
USE HexMesh_Class
USE HexElements_Class


IMPLICIT NONE


  TYPE DataViz
    INTEGER :: nVizNodes, nVizElements
    INTEGER :: N, nElements           ! Number of points per spectral element
    REAL(prec), ALLOCATABLE :: x(:,:)
    INTEGER, ALLOCATABLE    :: specToVizNode(:,:,:,:)
    INTEGER, ALLOCATABLE    :: (:,:,:,:)
  END TYPE DataViz

CONTAINS

END MODULE DataViz_Class
