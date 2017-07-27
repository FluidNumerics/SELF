

MODULE ImagePlane_Class


 ! src/common
 USE ModelPrecision


 IMPLICIT NONE
 
    TYPE ImagePlane 
       INTEGER    :: nX, nY
       REAL(prec) :: x(:,:,:) ! (1:3,nX,nY)
       REAL(SP)   :: rgba(:,:,:) !(1:4,nX,nY) ! 32-bit, single-precision, values for bitmap
       
       CONTAINS
       
    END TYPE ImagePlane


 CONTAINS
 
 SUBROUTINE Build_ImagePlane( myImage, nX, nY )
   ! Allocates space for the image plane and sets default values
   IMPLICIT NONE
   CLASS( ImagePlane ), INTENT(inout) :: myImage
   INTEGER, INTENT(in)                :: nX, nY
   
   
      myImage % nX = nX
      myImage % nY = nY
      
      ALLOCATE( myImage % x(1:3, nX, nY), myImage % rgba(1:4,nX,nY) )
      myImage % x    = 0.0_prec
      myImage % rgba = 0.0_SP
      
      
 END SUBROUTINE Build_ImagePlane
!
END MODULE ImagePlane_Class
