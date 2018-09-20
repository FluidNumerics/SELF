PROGRAM Generate_Cubed_Sphere

USE ModelPrecision
USE ConstantsDictionary
USE Lagrange_Class
USE NodalDG_Class
USE HexMesh_Class
USE BoundaryCommunicator_Class

IMPLICIT NONE

  LOGICAL                 :: setupSuccess
  INTEGER, ALLOCATABLE    :: quads(:,:)
  INTEGER, ALLOCATABLE    :: hexes(:,:)
  INTEGER, ALLOCATABLE    :: uniqueNodeIDs(:)
  INTEGER                 :: polyDeg, N, nLayers, l1, l2, iEl
  INTEGER                 :: ii, jj, i, j, face, k, node, nBoundaryFaces
  INTEGER                 :: actualID, nNodes, nUnique, nQuads
  REAL(prec), ALLOCATABLE :: shell_radius(:)
  REAL(prec), ALLOCATABLE :: vertices(:,:)
  REAL(prec), ALLOCATABLE :: shell_vertices(:,:)
  REAL(prec), ALLOCATABLE :: unique_vertices(:,:)
  REAL(prec), ALLOCATABLE :: unique_shell_vertices(:,:)
  REAL(prec), ALLOCATABLE :: xGeom(:,:,:,:)
  REAL(prec), ALLOCATABLE :: xGeom_linear(:,:,:,:)
  REAL(prec), ALLOCATABLE :: xGeom_3D(:,:,:,:,:)
  REAL(prec), ALLOCATABLE :: vertices_3D(:,:)
  REAL(prec)              :: innerRadius, outerRadius
  REAL(prec)              :: cube_face_bounds(1:3,1:4,1:6)
  REAL(prec)              :: s, p, r
  REAL(prec)              :: r1(1:3), r2(1:3), r3(1:3), r4(1:3)
  REAL(prec)              :: x1(1:3), x2(1:3), x3(1:3), x4(1:3)
  REAL(prec)              :: x(1:3), rtp(1:3)
  CHARACTER(50)           :: meshName
  CHARACTER(7)            :: zoneID
  LOGICAL                 :: foundSameNode
  TYPE( NodalDG )         :: dgStorage
  TYPE( HexMesh )         :: mesh
  TYPE( BoundaryCommunicator ) :: bCom


    CALL InitializeFromCommandLine( N, polyDeg, innerRadius, outerRadius, nLayers, meshName )

    nNodes = 6*(N+1)*(N+1)
    nUnique = (N+1)*(N+1) + 2*N*(N+1) + (N-1)*(N+1) + 2*(N-1)*(N-1)
    nQuads  = 6*N*N

    ALLOCATE( quads(1:4,1:6*N*N), &
              vertices(1:3,1:6*(N+1)*(N+1)), &
              uniqueNodeIDs(1:6*(N+1)*(N+1)), &
              shell_vertices(1:3,1:6*(N+1)*(N+1)), &
              unique_shell_vertices(1:3,1:6*(N+1)*(N+1)), &
              unique_vertices(1:3,1:nUnique), &
              xGeom(1:3,0:polyDeg,0:polyDeg,1:6*N*N), &
              xGeom_linear(1:3,0:polyDeg,0:polyDeg,1:6*N*N), &
              shell_radius(0:nLayers), &
              vertices_3D(1:3,1:nUnique*(nLayers+1)), &
              hexes(1:8,1:6*N*N*nLayers), &
              xGeom_3D(1:3,0:polyDeg,0:polyDeg,0:polyDeg,1:6*N*N*nLayers) )


    CALL dgStorage % Build( UniformPoints(-1.0_prec, 1.0_prec, 12), &
                            polyDeg, 12, GAUSS_LOBATTO )

    CALL Cube_Faces_Setup( cube_face_bounds )

    CALL Setup_Shells( )    

    CALL Generate_3D_Mesh( )
    
    CALL mesh % WriteSELFMeshFile( meshName )

    CALL mesh % WriteTecplot( meshName )


    DEALLOCATE( quads, &
                vertices, &
                uniqueNodeIDs, &
                shell_vertices, &
                unique_shell_vertices, &
                unique_vertices, &
                xGeom, &
                xGeom_linear, &
                shell_radius, &
                vertices_3D, &
                hexes, &
                xGeom_3D )

    CALL dgStorage % Trash( )
    CALL mesh % Trash( )
    CALL bCom % Trash( )


CONTAINS

 SUBROUTINE InitializeFromCommandLine( N, polyDeg, innerRadius, outerRadius, nLayers, meshName )
   IMPLICIT NONE
   INTEGER, INTENT(out)      :: N, nLayers, polyDeg
   REAL(prec), INTENT(out)   :: innerRadius, outerRadius
   CHARACTER(*), INTENT(out) :: meshName
   ! Local
   INTEGER        :: nArg, argID
   CHARACTER(500) :: argname
   LOGICAL        :: resolutionGiven, polyGiven, iRadiusGiven, oRadiusGiven, radialGiven, meshGiven


     resolutionGiven = .FALSE.
     iRadiusGiven    = .FALSE.
     oRadiusGiven    = .FALSE.
     radialGiven     = .FALSE.
     meshGiven       = .FALSE.

     setupSuccess = .FALSE.

     N           = 5
     polyDeg     = 2
     innerRadius = 1
     outerRadius = 10
     nLayers     = 2
     meshName    = 'cubed_sphere'

     nArg = command_argument_count( )

     IF( nArg > 0 )THEN
 
        DO argID = 1, nArg
  
          CALL get_command_argument( argID, argName )

          SELECT CASE( TRIM(argName) )

             CASE("--resolution")
                resolutionGiven = .TRUE.
             CASE("--polynomial-degree")
                polyGiven = .TRUE.
             CASE("--inner-radius")
                iRadiusGiven = .TRUE.
             CASE("--outer-radius")
                oRadiusGiven = .TRUE.
             CASE("--n-radial-layers")
                radialGiven = .TRUE.
             CASE("--mesh-name")
                meshGiven = .TRUE.
             CASE DEFAULT

               IF( resolutionGiven )THEN

                  READ(argName,*) N 
                  resolutionGiven = .FALSE.

               ELSEIF( polyGiven )THEN

                  READ(argName,*) polyDeg
                  polyGiven = .FALSE.

               ELSEIF( iRadiusGiven )THEN

                  READ(argName,*) innerRadius
                  iRadiusGiven = .FALSE.

               ELSEIF( oRadiusGiven )THEN

                  READ(argName,*) outerRadius
                  oRadiusGiven = .FALSE.

               ELSEIF( radialGiven )THEN

                  READ(argName,*) nLayers 
                  radialGiven = .FALSE.

               ELSEIF( meshGiven )THEN

                  meshName  = TRIM(argName)
                  meshGiven = .FALSE.

               ENDIF

          END SELECT 
        ENDDO

     ELSE

        ! List possible options
        PRINT*, '  cubed_sphere '
        PRINT*, '    Quick spectral-element cubed-sphere mesh generation utility.'
        PRINT*, '--------------------------------------------------------------'
        PRINT*, '  Usage : cubed_sphere <inputs>                                      '
        PRINT*, ''
        PRINT*, ' --resolution <number> '
        PRINT*, ''
        PRINT*, '  The resolution number refers to the number of elements in each'
        PRINT*, '  direction of each face. For example --resolution 5 results in '
        PRINT*, '  5x5 elements on each cube face.'
        PRINT*, '  If resolution is not provided, the resolution defaults to 5.'
        PRINT*, ''
        PRINT*, ' --polynomial-degree <number> '
        PRINT*, ''
        PRINT*, '  Polynomial degree for the lagrange interpolating polynomial.'
        PRINT*, '  Defaults to 2.'
        PRINT*, ''
        PRINT*, ' --inner-radius <number> '
        PRINT*, ''
        PRINT*, '  Inner sphere radius defaults to 1 if not provided'
        PRINT*, ''
        PRINT*, ' --outer-radius <number> '
        PRINT*, ''
        PRINT*, '  Outer sphere radius defaults to 10 if not provided'
        PRINT*, ''
        PRINT*, ' --n-radial-layers <number> '
        PRINT*, ''
        PRINT*, '  The number of spherical shell layers between the inner and '
        PRINT*, '  outer radii. Number of layers defaults to 2'
        PRINT*, ''
        PRINT*, ' --mesh-name <Name of the mesh> '
        PRINT*, ''
        PRINT*, ' If the mesh-name is not provided, cubed_sphere.tec, '
        PRINT*, ' cubed_sphere.geom, and cubed_sphere.mesh, are created. '
        PRINT*, ''
        PRINT*, '--------------------------------------------------------------'
        
     ENDIF   


 END SUBROUTINE InitializeFromCommandLine

 FUNCTION Cartesian_to_Spherical( x ) RESULT( sph )
   IMPLICIT NONE
   REAL(prec) :: x(1:3)
   REAL(prec) :: sph(1:3)
   REAL(prec) :: r1

     sph(1) = sqrt( x(1)**2 + x(2)**2 + x(3)**2 )

     sph(2) = ATAN2( x(2), x(1) )

     r1 = sqrt( x(1)**2 + x(2)**2 )
     sph(3) = ATAN2( x(3), r1 )
     
 END FUNCTION Cartesian_to_Spherical
 
 SUBROUTINE Cube_Faces_Setup( cube_corners )
   IMPLICIT NONE
   REAL(prec), INTENT(out) :: cube_corners(1:3,1:4,1:6)

    ! South Face
    cube_corners(1:3,1,1) = (/ -1.0_prec, -1.0_prec, -1.0_prec /)    
    cube_corners(1:3,2,1) = (/ 1.0_prec, -1.0_prec, -1.0_prec /)    
    cube_corners(1:3,3,1) = (/ -1.0_prec, -1.0_prec, 1.0_prec /)    
    cube_corners(1:3,4,1) = (/ 1.0_prec, -1.0_prec, 1.0_prec /)    

    ! East Face
    cube_corners(1:3,1,2) = (/ 1.0_prec, -1.0_prec, -1.0_prec /)    
    cube_corners(1:3,2,2) = (/ 1.0_prec, 1.0_prec, -1.0_prec /)    
    cube_corners(1:3,3,2) = (/ 1.0_prec, -1.0_prec, 1.0_prec /)    
    cube_corners(1:3,4,2) = (/ 1.0_prec, 1.0_prec, 1.0_prec /)    

    ! North Face
    cube_corners(1:3,1,3) = (/ -1.0_prec, 1.0_prec, -1.0_prec /)    
    cube_corners(1:3,2,3) = (/ 1.0_prec, 1.0_prec, -1.0_prec /)    
    cube_corners(1:3,3,3) = (/ -1.0_prec, 1.0_prec, 1.0_prec /)    
    cube_corners(1:3,4,3) = (/ 1.0_prec, 1.0_prec, 1.0_prec /)    

    ! West Face
    cube_corners(1:3,1,4) = (/ -1.0_prec, -1.0_prec, -1.0_prec /)    
    cube_corners(1:3,2,4) = (/ -1.0_prec, 1.0_prec, -1.0_prec /)    
    cube_corners(1:3,3,4) = (/ -1.0_prec, -1.0_prec, 1.0_prec /)    
    cube_corners(1:3,4,4) = (/ -1.0_prec, 1.0_prec, 1.0_prec /)    

    ! Bottom Face
    cube_corners(1:3,1,5) = (/ -1.0_prec, -1.0_prec, -1.0_prec /)    
    cube_corners(1:3,2,5) = (/ 1.0_prec, -1.0_prec, -1.0_prec /)    
    cube_corners(1:3,3,5) = (/ -1.0_prec, 1.0_prec, -1.0_prec /)    
    cube_corners(1:3,4,5) = (/ 1.0_prec, 1.0_prec, -1.0_prec /)    

    ! Top Face
    cube_corners(1:3,1,6) = (/ -1.0_prec, -1.0_prec, 1.0_prec /)    
    cube_corners(1:3,2,6) = (/ 1.0_prec, -1.0_prec, 1.0_prec /)    
    cube_corners(1:3,3,6) = (/ -1.0_prec, 1.0_prec, 1.0_prec /)    
    cube_corners(1:3,4,6) = (/ 1.0_prec, 1.0_prec, 1.0_prec /)    


 END SUBROUTINE Cube_Faces_Setup
!
 SUBROUTINE Setup_Shells( )
   IMPLICIT NONE

    k = 0
    DO face = 1, 6
      DO j = 0, N

        p  = 2.0_prec*REAL(j,prec)/REAL(N,prec) -1.0_prec
        r2 = ( cube_face_bounds(1:3,4,face) - &
               cube_face_bounds(1:3,2,face) )*&
             ( p + 1.0_prec )/2.0_prec + &
             cube_face_bounds(1:3,2,face) 

        r4 = ( cube_face_bounds(1:3,3,face) - &
               cube_face_bounds(1:3,1,face) )*&
             ( p + 1.0_prec )/2.0_prec + &
             cube_face_bounds(1:3,1,face) 

        DO i= 0, N

          s  = 2.0_prec*REAL(i,prec)/REAL(N,prec) -1.0_prec

          r1 = ( cube_face_bounds(1:3,2,face) - &
                 cube_face_bounds(1:3,1,face) )*&
               ( s + 1.0_prec )/2.0_prec + &
               cube_face_bounds(1:3,1,face) 

          r3 = ( cube_face_bounds(1:3,4,face) - &
                 cube_face_bounds(1:3,3,face) )*&
               ( s + 1.0_prec )/2.0_prec + &
               cube_face_bounds(1:3,3,face) 


          x = 0.5_prec*( r3*(p+1.0_prec) - r1*(p-1.0_prec) )+&
              0.5_prec*( r2*(s+1.0_prec) - r4*(s-1.0_prec) )

          rtp = Cartesian_to_Spherical( x )

          k = k + 1
          vertices(1:3,k) = x
          shell_vertices(1:3,k) = x/rtp(1)


        ENDDO
      ENDDO

    ENDDO
   
    uniqueNodeIDs = 0
    k = 0
    DO j = 1, nNodes
     
      foundSameNode = .FALSE.
      DO i = 1, j-1

        foundSameNode = Nodes_are_same( vertices(1:3,i), vertices(1:3,j) )
        IF( foundSameNode )THEN
          uniqueNodeIDs(j) = uniqueNodeIDs(i)
          EXIT
        ENDIF

      ENDDO
      
      IF( .NOT. foundSameNode )THEN
        k = k + 1
        uniqueNodeIDs(j) = k
        unique_vertices(1:3,k) = vertices(1:3,j)
        unique_shell_vertices(1:3,k) = shell_vertices(1:3,j)
      ENDIF

    ENDDO
 
    k = 0
    DO face = 0, 5
      DO j = 0, N-1
        DO i = 0, N-1

          k = k+1
          quads(1,k) = uniqueNodeIDs( (face*(N+1) + j)*(N+1) + i+1 )
          quads(2,k) = uniqueNodeIDs( (face*(N+1) + j)*(N+1) + i+2 )
          quads(3,k) = uniqueNodeIDs( (face*(N+1) + j+1)*(N+1) + i+2 )
          quads(4,k) = uniqueNodeIDs( (face*(N+1) + j+1)*(N+1) + i+1 )

          x1 = unique_vertices(1:3,quads(1,k))
          x2 = unique_vertices(1:3,quads(2,k))
          x3 = unique_vertices(1:3,quads(3,k))
          x4 = unique_vertices(1:3,quads(4,k))

          DO jj = 0, polyDeg

            p  = dgStorage % interp % interpolationPoints(jj)
            r2 = ( x3 - x2 )*( p + 1.0_prec )/2.0_prec + x2
            r4 = ( x4 - x1 )*( p + 1.0_prec )/2.0_prec + x1

            DO ii = 0, polyDeg

              s  = dgStorage % interp % interpolationPoints(ii)
              r1 = ( x2 - x1 )*( s + 1.0_prec )/2.0_prec + x1 
              r3 = ( x3 - x4 )*( s + 1.0_prec )/2.0_prec + x4
            
              x = 0.5_prec*( r3*(p+1.0_prec) - r1*(p-1.0_prec) )+&
                  0.5_prec*( r2*(s+1.0_prec) - r4*(s-1.0_prec) )

              rtp = Cartesian_to_Spherical( x )

              xGeom(1:3,ii,jj,k) = x/rtp(1)

            ENDDO
          ENDDO

          x1 = unique_shell_vertices(1:3,quads(1,k))
          x2 = unique_shell_vertices(1:3,quads(2,k))
          x3 = unique_shell_vertices(1:3,quads(3,k))
          x4 = unique_shell_vertices(1:3,quads(4,k))

          DO jj = 0, polyDeg

            p  = dgStorage % interp % interpolationPoints(jj)
            r2 = ( x3 - x2 )*( p + 1.0_prec )/2.0_prec + x2
            r4 = ( x4 - x1 )*( p + 1.0_prec )/2.0_prec + x1

            DO ii = 0, polyDeg

              s  = dgStorage % interp % interpolationPoints(ii)
              r1 = ( x2 - x1 )*( s + 1.0_prec )/2.0_prec + x1 
              r3 = ( x3 - x4 )*( s + 1.0_prec )/2.0_prec + x4
            
              x = 0.25_prec*( ( r3*(p+1.0_prec) - r1*(p-1.0_prec) )+&
                              ( r2*(s+1.0_prec) - r4*(s-1.0_prec) ) )

              xGeom_linear(1:3,ii,jj,k) = x

            ENDDO
          ENDDO


        ENDDO
      ENDDO
    ENDDO

 END SUBROUTINE Setup_Shells

 SUBROUTINE Generate_3D_Mesh( )
   IMPLICIT NONE
   
    shell_radius = UniformPoints(innerRadius, outerRadius, nLayers)


    DO k = 1, nLayers

      l1 = 1+6*N*N*(k-1)
      l2 = 6*N*N*k

      IF( k == 1 )THEN

        DO j = 0, polyDeg
          s  = dgStorage % interp % interpolationPoints(j)
          xGeom_3D(1:3,0:polyDeg,0:polyDeg,j,l1:l2) = ( 0.5_prec*(xGeom_linear*shell_radius(k) - xGeom*shell_radius(k-1))*( s + 1.0_prec ) + xGeom*shell_radius(k-1) )
        ENDDO

      ELSEIF( k == nLayers )THEN

        DO j = 0, polyDeg
          s  = dgStorage % interp % interpolationPoints(j)
          xGeom_3D(1:3,0:polyDeg,0:polyDeg,j,l1:l2) = ( 0.5_prec*(xGeom*shell_radius(k) - xGeom_linear*shell_radius(k-1))*( s + 1.0_prec ) + xGeom_linear*shell_radius(k-1) )
        ENDDO

      ELSE

        DO j = 0, polyDeg
          s  = dgStorage % interp % interpolationPoints(j)
          xGeom_3D(1:3,0:polyDeg,0:polyDeg,j,l1:l2) = xGeom_linear*( 0.5_prec*(shell_radius(k) - shell_radius(k-1))*( s + 1.0_prec ) + shell_radius(k-1) )
        ENDDO

      ENDIF

      hexes(1:4,l1:l2) = quads(1:4,1:nQuads) + (k-1)*nUnique
      hexes(5:8,l1:l2) = quads(1:4,1:nQuads) + (k)*nUnique

    ENDDO

    DO k = 1, nLayers+1

      l1 = 1 + nUnique*(k-1)
      l2 = nUnique*k
      vertices_3D(1:3,l1:l2) =  unique_shell_vertices(1:3,1:nUnique)*shell_radius(k-1)

    ENDDO
   
    CALL mesh % Build( nUnique*(nLayers+1), 6*N*N*nLayers, 1, polyDeg )

    DO k = 1, nUnique*(nLayers+1)

      mesh % nodes % nodeID(k) = k
      mesh % nodes % x(1:3,k)  = vertices_3d(1:3,k)

    ENDDO

    DO iEl = 1, 6*N*N*nLayers

      mesh % elements % elementID(iEl)   = iEl
      mesh % elements % nodeIDs(1:8,iEl) = hexes(1:8,iEl)

      DO k = 0, polyDeg
        DO j = 0, polyDeg
          DO i = 0, polyDeg

             mesh % elements % x(i,j,k,1:3,iEl) = xGeom_3D(1:3,i,j,k,iEl)

          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CALL mesh % elements % GenerateMetrics( dgStorage % interp )
    CALL mesh % ConstructFaces( )
             

    ! Build the external/boundary communicator
    nBoundaryFaces = mesh % NumberOfBoundaryFaces(  )
    CALL bcom % Build( nBoundaryFaces )

    k = 0
    DO i = 1, mesh % faces % nFaces

      IF( mesh % faces % elementIDs(2,i) < 0 )THEN
        k = k + 1
        mesh % faces % boundaryID(i) = -k
        bcom % boundaryIDs(k)        = i
        bcom % boundaryGlobalIDs(k)  = i
        bcom % extProcIDs(k)         = 0
      ENDIF

    ENDDO

    CALL bcom % WritePickup( 'ExtComm.0000' )
    

 END SUBROUTINE Generate_3D_Mesh 

 FUNCTION Nodes_are_same( x1, x2 ) RESULT( sameNode )
   IMPLICIT NONE
   REAL(prec) :: x1(1:3), x2(1:3)
   LOGICAL    :: sameNode
   REAL(prec) :: r(1:3)

    sameNode = .FALSE.
    r = x2-x1
    IF( SQRT( r(1)**2 + r(2)**2 + r(3)**2 ) <= 10.0_prec**(-5) ) THEN
      sameNode = .TRUE.
    ENDIF

 END FUNCTION Nodes_are_same

END PROGRAM Generate_Cubed_Sphere

