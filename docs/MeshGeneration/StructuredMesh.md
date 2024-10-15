# Structured Mesh Generation in SELF


Although SELF uses unstructured mesh data structures, we have provided methods to create structured meshes and store them as an unstructured mesh. This can be quite useful as you are getting started with SELF or if the geometry for your problem can be defined using a structured mesh layout.


## One Dimension (1-D)
In one dimension, the only mesh we use is a structured mesh. To generate a structured mesh in one dimension, use the `StructuredMesh` generic in the [`Mesh1D`](../ford/type/mesh1d.html) class. 

At the moment, only uniformly spaced structured meshes of elements can be generated. This means that all of the elements are of the same width; keep in mind that within each element, there is a quadrature grid. The points in the quadrature grid are spaced so that spectral accuracy is guaranteed.

To generate a structured grid in 1-D, you need to provide the number of elements and the left and right endpoints of the mesh. 

In the example below, we create a 1-D mesh with 20 elements on the domain $x ∈ [0,1]$. The geometry fields are created from the mesh information and a $7^{th}$ degree interpolant through the Legendre-Gauss points. 


```fortran
  type(Mesh1D),target :: mesh
  type(Lagrange),target :: interp
  type(Geometry1D),target :: geometry

  call mesh % StructuredMesh(nElem=20, &
                             x=(/0.0_prec,1.0_prec/))

  call interp % Init(N=7, controlNodeType=GAUSS, &
                     M=10, targetNodeType=UNIFORM)

  call geometry % Init(interp,mesh%nElem)
  call geometry % GenerateFromMesh(mesh)

```

Notice that initializing the geometry requires an interpolant and the number of elements as input. 

!!! note
    Under the hood, the interpolant for the geometry (`geometry % interp` ) is associated with a pointer to `interp`, ie `geometry % interp => interp`.

Once the geometry is initialized, the physical positions and metric terms can be calculated and stored using the `GenerateFromMesh` method.

## Two Dimensions (2-D)
To generate a structured mesh in two dimensions, use the `StructuredMesh` generic in the [`Mesh2D`](../ford/type/mesh2d.html) class. 

At the moment, only uniformly spaced structured meshes of elements can be generated. This means that all of the elements are of the same width; keep in mind that within each element, there is a quadrature grid. The points in the quadrature grid are spaced so that spectral accuracy is guaranteed.

SELF uses a tiled structured grid. Tiled grids divide the 2-D grid into `nTilex`x`nTiley` tiles of size `nxPerTile`x`nyPerTile` . The width and height of the elements are defined as `dx` and `dy`. With these parameters,

* `nx = nTilex*nxPerTile` is the total number of elements in the x-direction
* `ny = nTiley*nyPerTile` is the total number of elements in the y-direction
* `nelem = nx*ny` is the total number of elements
* `Lx = dx*nx` is the domain length in the x-direction
* `Ly = dy*ny` is the domain length in the y-direction

You can set boundary conditions for each of the four sides of the structured mesh using a 1-D array of integers of length 4. The boundary conditions must be provided in counter-clockwise order, starting with the "south" boundary (south, east, north, west). The following built-in flags are available for setting boundary conditions

* `SELF_BC_NONORMALFLOW`
* `SELF_BC_PRESCRIBED`
* `SELF_BC_RADIATION`

The tiled layout is convenient for domain decomposition, when you are wanting to scale up your application for distributed memory platforms. You can further enable domain decomposition by setting the optional `enableDomainDecompisition` input to `.true.` . In this case, when you launch your application with `mpirun`, the domain will be automatically divided as evenly as possible across all MPI ranks.

!!! note
    It's good practice to set the total number of tiles equal to the number of MPI ranks that you are running with. Alternatively, you can use fairly small tiles when working with a large number of MPI ranks to increase the chance of minimizing point-to-point communications .

In the example below, we create a 2-D mesh with the following attributes

* $2 × 2$ tiles for the domain
* $10 × 10$ elements per tile
* Each element is has dimensions of $0.05 × 0.05$. The domain dimensions are then $L_x × L_y = 1 × 1$
* Domain decomposition is enabled

The geometry fields are created from the mesh information and a $7^{th}$ degree interpolant through the Legendre-Gauss points. 


```fortran
  type(Mesh2D),target :: mesh
  type(Lagrange),target :: interp
  type(SEMQuad),target :: geometry
  integer :: bcids(1:4)

  bcids(1:4) = [SELF_BC_PRESCRIBED,& ! south boundary
                  SELF_BC_PRESCRIBED,& ! east boundary
                  SELF_BC_PRESCRIBED,& ! north boundary
                  SELF_BC_PRESCRIBED]  ! west boundary

  call mesh % StructuredMesh( nxPerTile=10, nyPerTile=10,&
                              nTileX=2, nTileY=2,&
                              dx=0.05_prec, dy=0.05_prec, &
                              bcids)

  call interp % Init(N=7, controlNodeType=GAUSS, &
                     M=10, targetNodeType=UNIFORM)

  call geometry % Init(interp,mesh%nElem)
  call geometry % GenerateFromMesh(mesh)

```

Notice that initializing the geometry requires an interpolant and the number of elements as input. 

!!! note
    Under the hood, the interpolant for the geometry (`geometry % interp` ) is associated with a pointer to `interp`, ie `geometry % interp => interp`.

Once the geometry is initialized, the physical positions and metric terms can be calculated and stored using the `GenerateFromMesh` method.


## Three Dimensions (3-D)
To generate a structured mesh in three dimensions, use the `StructuredMesh` generic in the [`Mesh3D`](../ford/type/mesh3d.html) class. 

At the moment, only uniformly spaced structured meshes of elements can be generated. This means that all of the elements are of the same length, width, and height; though, the length, width, and height can each be their own value. Keep in mind that within each element, there is a quadrature grid. The points in the quadrature grid are spaced so that spectral accuracy is guaranteed.

SELF uses a tiled structured grid. Tiled grids divide the 3-D grid into `nTilex`x`nTiley`x`nTilez` tiles of size `nxPerTile`x`nyPerTile`x`nzPerTile` . The length, width, and height of the elements are defined as `dx`, `dy`, and `dz` respectively. With these parameters,

* `nx = nTilex*nxPerTile` is the total number of elements in the x-direction
* `ny = nTiley*nyPerTile` is the total number of elements in the y-direction
* `nz = nTilez*nzPerTile` is the total number of elements in the z-direction
* `nelem = nx*ny*nz` is the total number of elements
* `Lx = dx*nx` is the domain length in the x-direction
* `Ly = dy*ny` is the domain length in the y-direction
* `Lz = dz*nz` is the domain length in the z-direction

You can set boundary conditions for each of the four sides of the structured mesh using a 1-D array of integers of length 6. The boundary conditions must be provided in CGNS ordering (bottom,south, east, north, west,top). The following built-in flags are available for setting boundary conditions

* `SELF_BC_NONORMALFLOW`
* `SELF_BC_PRESCRIBED`
* `SELF_BC_RADIATION`

The tiled layout is convenient for domain decomposition, when you are wanting to scale up your application for distributed memory platforms. You can further enable domain decomposition by setting the optional `enableDomainDecompisition` input to `.true.` . In this case, when you launch your application with `mpirun`, the domain will be automatically divided as evenly as possible across all MPI ranks.

!!! note
    It's good practice to set the total number of tiles equal to the number of MPI ranks that you are running with. Alternatively, you can use fairly small tiles when working with a large number of MPI ranks to increase the chance of minimizing point-to-point communications .

In the example below, we create a 3-D mesh with the following attributes

* $2 × 2 × 2$ tiles for the domain
* $10 × 10 × 10$ elements per tile
* Each element is has dimensions of $0.05 × 0.05 × 0.05$. The domain dimensions are then $L_x × L_y × L_z = 1 × 1 × 1$
* Domain decomposition is enabled

The geometry fields are created from the mesh information and a $7^{th}$ degree interpolant through the Legendre-Gauss points. 


```fortran
  type(Mesh3D),target :: mesh
  type(Lagrange),target :: interp
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)

  bcids(1:6) = [SELF_BC_PRESCRIBED,&   ! bottom boundary
                  SELF_BC_PRESCRIBED,& ! south boundary
                  SELF_BC_PRESCRIBED,& ! east boundary
                  SELF_BC_PRESCRIBED,& ! north boundary
                  SELF_BC_PRESCRIBED,& ! west boundary
                  SELF_BC_PRESCRIBED]  ! yop boundary

  call mesh % StructuredMesh( nxPerTile=10, nyPerTile=10, nzPerTile=10&
                              nTileX=2, nTileY=2, nTileZ=2, &
                              dx=0.05_prec, dy=0.05_prec, dz=0.05_prec, &
                              bcids)

  call interp % Init(N=7, controlNodeType=GAUSS, &
                     M=10, targetNodeType=UNIFORM)

  call geometry % Init(interp,mesh%nElem)
  call geometry % GenerateFromMesh(mesh)

```

Notice that initializing the geometry requires an interpolant and the number of elements as input. 

!!! note
    Under the hood, the interpolant for the geometry (`geometry % interp` ) is associated with a pointer to `interp`, ie `geometry % interp => interp`.

Once the geometry is initialized, the physical positions and metric terms can be calculated and stored using the `GenerateFromMesh` method.
