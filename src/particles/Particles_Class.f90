! :::::::::::::::::::::::::::
! Author : Joseph Schoonover
! Copyright (2017)
! :::::::::::::::::::::::::::

MODULE Particles_Class

 ! src/common/
 USE ModelPrecision
 ! src/spectralops/
 USE NodalStorage
 ! src/geom/
 USE HexMesh_Class
 ! src/particles/
 USE ParticleParams
 
 IMPLICIT NONE
 
 ! The "Particles" data type is defined by arrays of particle positions
 ! and there tendencies. In the expectation that the particles may be 
 ! used on an eulerian mesh, they can also be associated with an element
 ! on an unstructured mesh. By keeping track of the element that it resides
 ! in, one can interpolate the eulerian velocity field onto the particle 
 ! position to aid in the computation of the particle tendency.
 !
 ! For this demonstration module, a hexahedral mesh is included in the 
 ! data structure. In production applications, this is not a requirement.
 !
 ! Currently, a fixed array size is used for the particles attributes. If 
 ! a particle is not found within the domain, its entry in an array of 
 ! logicals ("activeStatus") is set to false and the particulate is no 
 ! longer tracked.
 !
   TYPE Particles
      INTEGER                 :: nActive
      LOGICAL, ALLOCATABLE    :: activeStatus(:)
      INTEGER, ALLOCATABLE    :: element(:)
      REAL(prec), ALLOCATABLE :: mass(:)
      REAL(prec), ALLOCATABLE :: x(:,:)
      REAL(prec), ALLOCATABLE :: s(:,:)
      REAL(prec), ALLOCATABLE :: dxdt(:,:)
      REAL(prec), ALLOCATABLE :: meshVelocity(:,:,:,:)
      TYPE( HexMesh )         :: mesh
      TYPE( NodalStorage )    :: nodal
      TYPE( ParticleParams )  :: params
      
      CONTAINS
      
      PROCEDURE, PUBLIC  :: Build        => Build_Particles
      PROCEDURE, PRIVATE :: BuildHexMesh => BuildHexMesh_Particles
      PROCEDURE, PUBLIC  :: Trash        => Trash_Particles
      PROCEDURE, PRIVATE :: InitializeParticlePositions
      
      !PROCEDURE, PUBLIC  :: ForwardStepRK3 => ForwardStepRK3_Particles
      
      !PROCEDURE, PUBLIC  :: EulerianToParticles    => EulerianToParticles_Particles 
      PROCEDURE, PUBLIC  :: FindAssociatedElements => FindAssociatedElements_Particles
      
      PROCEDURE, PUBLIC  :: UpdateActiveParticleCount => UpdateActiveParticleCount_Particles
      
   END TYPE Particles

 CONTAINS 

 SUBROUTINE Build_Particles( myParticles )
 
   IMPLICIT NONE
   TYPE( Particles ), INTENT(inout) :: myParticles
   
      CALL myParticles % params % Build( )
      
      !Initially, all particles are active
      myParticles % nActive    = myParticles % params % nParticles
      
      ! For now, use only trilinear elements
      CALL myParticles % nodal % Build( 1, 2, GAUSS_LOBATTO, CG )
      
      ! Load in the hexahedral mesh
      CALL myParticles % BuildHexMesh( )
      
      ALLOCATE( myParticles % activeStatus(1:myParticles % params % nParticles), &
                myParticles % element(1:myParticles % params % nParticles), &
                myParticles % mass(1:myParticles % params % nParticles), &
                myParticles % x(1:3,1:myParticles % params % nParticles), &
                myParticles % s(1:3,1:myParticles % params % nParticles), &
                myParticles % dxdt(1:3,1:myParticles % params % nParticles),& 
                myParticles % meshVelocity(0:1,0:1,0:1,1:myParticles % mesh % nElems) )
                
      DO i = 1, myParticles % params % nParticles
         myParticles % activeStatus(i) = .TRUE.
         myParticles % element(i)      = 0
         myParticles % mass(i)         = 1.0_prec
         myParticles % x(1:3,i)        = 0.0_prec
         myParticles % s(1:3,i)        = 0.0_prec
         myParticles % dxdt(1:3,i)     = 0.0_prec
      ENDDO
      
      ! Set up initial particle positions
      CALL myParticles % InitializeParticlePositions( )
 
 END SUBROUTINE Build_Particles
!
 SUBROUTINE BuildHexMesh_Particles( myParticles )

   IMPLICIT NONE
   CLASS( Particles ), INTENT(inout) :: myParticles


      PRINT*,'Module Particles_Class.f90 : S/R BuildHexMesh :'

      PRINT*, 'Reading mesh from '//trim(myParticles % params % PeaceMeshFile)//'.'
      CALL myParticles % mesh % ReadPeaceMeshFile( myParticles % params % PeaceMeshFile )
      
      CALL myParticles % mesh % ScaleTheMesh( myParticles % nodal % interp, &
                                          myParticles % params % xScale, &
                                          myParticles % params % yScale, &
                                          myParticles % params % zScale )


 END SUBROUTINE BuildHexMesh_Particles
!
 SUBROUTINE Trash_Particles( myParticles )
 
   IMPLICIT NONE
   CLASS( Particles ), INTENT(inout) :: myParticles
   
      DEALLOCATE( myParticles % activeStatus, &
                  myParticles % element, &
                  myParticles % mass, &
                  myParticles % x, &
                  myParticles % s, &
                  myParticles % dxdt, &
                  myParticles % meshVelocity )
                  
      CALL myParticles % nodal % Trash( )
      CALL myParticles % mesh % Trash( )
 
 END SUBROUTINE Trash_Particles
!
 SUBROUTINE InitializeParticlePositions( myParticles )

   IMPLICIT NONE
   CLASS( Particles ), INTENT(inout) :: myParticles
   
   
 END SUBROUTINE InitializeParticlePositions
!
 SUBROUTINE FindAssociatedElements_Particles( myParticles )
 
   IMPLICIT NONE
   CLASS( Particles ), INTENT(inout) :: myParticles
   ! Local
   INTEGER    :: i, j, iEl
   LOGICAL    :: inElement
   
      DO i = 1, myParticles % nParticles
      
         isInElement = .FALSE.
         
         IF( myParticles % element(i) == 0 )THEN
            DO iEl = 1, myParticles % mesh % nElems
                      
               CALL myParticles % mesh % geom(iEl) % &
                 CalculateComputationalCoordinates( myParticles % nodal % interp, &
                                                    myParticles % x(1:3,i), &
                                                    myParticles % s(1:3,i), &
                                                    inElement )
               IF( inElement )THEN
                  myParticles % element(i) = iEl
                  BREAK
               ENDIF
            
            ENDDO
         
            IF( .NOT.(inElement) )THEN
               myParticles % activeStatus(i) = .FALSE.
            ENDIF
            
         ELSE
            
            ! Check if the particle is in its current element
            iEl = myParticles % element(i)
            CALL myParticles % mesh % geom(iEl) % &
                 CalculateComputationalCoordinates( myParticles % nodal % interp, &
                                                    myParticles % x(1:3,i), &
                                                    myParticles % s(1:3,i), &
                                                    inElement )
            IF( .NOT.(inElement) )THEN
               ! Search the element neighbors
            ENDIF
            
         ENDIF
      
      ENDDO
      
      
 END SUBROUTINE FindAssociatedElements_Particles
 
END MODULE Particles_Class
