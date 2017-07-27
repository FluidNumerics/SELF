

MODULE RayCasting_Class


 ! src/common
 USE ModelPrecision
 ! src/geom/
 USE MappedGeometry_3D_Class 
 USE HexMesh_Class


 IMPLICIT NONE
 
    TYPE RayCasting 
       TYPE( Ray ), ALLOCATABLE :: 
       
    END TYPE RayCasting

END MODULE RayCasting_Class
