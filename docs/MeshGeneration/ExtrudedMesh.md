# Periodic and Extruded Meshes

This page covers two mesh constructors that are used together to build layered volume meshes, such as those of a hydrostatic ocean model.
The first builds a 2-D structured mesh that is periodic in either horizontal direction.
The second extrudes any 2-D mesh downward through a stack of level surfaces into a 3-D mesh of hexahedra, whose columns follow the bathymetry.

## Per-axis periodic 2-D structured meshes

`Mesh2D % PeriodicStructuredMesh` builds the same tiled, uniformly spaced mesh as `Mesh2D % StructuredMesh` and differs only in its connectivity.
Each direction has its own periodic flag.
When `periodicX` is true, the east edge of the last column of elements and the west edge of the first column are wired as one interior edge, so the side exchange delivers the neighbor trace across $x = 0 \equiv x = L_x$.
The `periodicY` flag does the same for the south and north edges.
The `bcids` array is still given in the order south, east, north, west, and only the entries of a non-periodic direction are used.

```fortran
  type(Mesh2D),target :: mesh
  integer :: bcids(1:4)

  bcids(1:4) = [SELF_BC_NONORMALFLOW, & ! south wall
                SELF_BC_NONORMALFLOW, & ! east (unused, periodic)
                SELF_BC_NONORMALFLOW, & ! north wall
                SELF_BC_NONORMALFLOW]   ! west (unused, periodic)

  ! A channel of 16 x 8 elements, periodic in x and walled in y
  call mesh % PeriodicStructuredMesh(nxPerTile=8, nyPerTile=8, &
                                     nTileX=2, nTileY=1, &
                                     dx=0.0625_prec, dy=0.125_prec, &
                                     periodicX=.true., periodicY=.false., &
                                     bcids=bcids)
```

The connectivity follows the triply periodic `Mesh3D % PeriodicStructuredMesh`.
Every interior edge, wrapped or not, joins two edges whose trace coordinates run in the same direction, so every flip is zero.
The global side ids are assigned from the grid position of each edge and are dense in $1,\dots,$ `nUniqueSides`, with the two sides of a periodic edge sharing one id.
With $n_{Lx} = n_x$ in a periodic $x$ direction and $n_x+1$ otherwise, and likewise $n_{Ly}$, there are $n_{Lx} n_y + n_{Ly} n_x$ unique edges.
With neither direction periodic, the mesh is identical to the one `StructuredMesh` builds.

## Extruded 3-D meshes

`Mesh3D % ExtrudedMesh` turns a 2-D surface mesh into a 3-D mesh of `nz` layers.
The surface mesh can be any conforming `Mesh2D`: structured, periodic structured, read from HOPr, or read from HOHQMesh, with straight or curved isoparametric elements of any geometry degree $N_{geo}$.
The result is an ordinary `Mesh3D`, so the geometry, the mapped data types and the DG models use it unchanged.

### The map of an extruded element

The column above a surface element $Q$ is cut into $n_z$ layers by the level surfaces

$$
-h = z_0(x,y) < z_1(x,y) < \dots < z_{n_z}(x,y) = 0 .
$$

Let $(X(\xi,\eta), Y(\xi,\eta))$ be the isoparametric map of $Q$.
Element $k$ of the column has the map

$$
\vec{x}(\xi,\eta,\zeta) = \big(X(\xi,\eta),\ Y(\xi,\eta),\ Z(\xi,\eta,\zeta)\big), \qquad
Z = \tfrac{1-\zeta}{2}\, z_{k-1} + \tfrac{1+\zeta}{2}\, z_k ,
$$

where $Z$ is evaluated at the $N_{geo}+1$ geometry nodes of $Q$ in each horizontal direction.
The horizontal coordinates do not depend on $\zeta$, so lines of constant $(\xi,\eta)$ are vertical.
The 3-D mesh keeps the geometry degree and the geometry node type of the surface mesh.
Its node coordinates are sampled on the same $N_{geo}+1$ points in all three directions, which represents the linear vertical map exactly.
At the end nodes of a Lobatto-type or uniform node set the two weights are exactly 0 and 1, so the top face of layer $k$ and the bottom face of layer $k+1$ are sampled from the same $z_k$ and agree bitwise.

You give the levels at the surface geometry nodes, as `zLevels(1:nGeo+1,1:nGeo+1,0:nz,1:nElem2D)` for the elements of the surface mesh that this rank holds, in the length unit of the mesh coordinates and positive upward.
The constructor stops with an error unless `nz` is at least 1, the levels increase strictly with the level index at every node, and the top level is exactly zero.
The flat top is the reference sea surface the linear free surface condition is applied on.

The levels must also be continuous across lateral faces: two surface elements that share a geometry node must be given the same $z_k$ there.
This holds whenever the levels are computed pointwise from a bathymetry function sampled at the node coordinates of the surface mesh, as the two helpers below do.
The constructor cannot check it, because the neighbor across a lateral face may be on another rank.

### Building the levels from a bathymetry

Two helpers fill `zLevels` from the water depth `h(1:nGeo+1,1:nGeo+1,1:nElem2D)` at the surface geometry nodes, with $h > 0$.

`ExtrudedSigmaLevels(h, sigma, zLevels)` builds terrain-following levels $z_k = \sigma_k h$ from a stretching `sigma(0:nz)` with $\sigma_0 = -1$, $\sigma_{n_z} = 0$, increasing strictly.

`ExtrudedHybridLevels(h, depths, sigma, zLevels)` places $n_F$ flat interfaces at the depths $0 < d_1 < \dots < d_{n_F}$ below the surface and $n_S$ terrain-following layers between $z = -d_{n_F}$ and the sea floor below them.
With $d_0 = 0$ and $d_F = d_{n_F}$, the $n_z = n_S + n_F$ levels are

$$
z_k = -\big((1+\sigma_k)\, d_F - \sigma_k h\big), \quad k = 0,\dots,n_S, \qquad
z_{n_S+m} = -d_{n_F-m}, \quad m = 1,\dots,n_F .
$$

The deepest flat interface must lie above the sea floor everywhere, $d_F < \min h$, and the helper stops with an error otherwise.
With no flat interfaces this reduces to the sigma levels, and with a constant $h$ it gives flat $z$-levels throughout.

```fortran
  type(Mesh2D),target :: surface
  type(Mesh3D),target :: mesh
  real(prec),allocatable :: h(:,:,:),zLevels(:,:,:,:)
  integer :: nGeo,e,i,j

  call surface % Read_HOPr("basin_mesh.h5")
  nGeo = surface % nGeo

  ! Sample the bathymetry at the surface geometry nodes
  allocate(h(1:nGeo+1,1:nGeo+1,1:surface % nElem))
  do e = 1,surface % nElem
    do j = 1,nGeo+1
      do i = 1,nGeo+1
        h(i,j,e) = depth(surface % nodeCoords(1,i,j,e),surface % nodeCoords(2,i,j,e))
      enddo
    enddo
  enddo

  ! Two flat layers in the upper 100 m over four terrain-following layers
  allocate(zLevels(1:nGeo+1,1:nGeo+1,0:6,1:surface % nElem))
  call ExtrudedHybridLevels(h,[50.0_prec,100.0_prec], &
                            [-1.0_prec,-0.6_prec,-0.3_prec,-0.1_prec,0.0_prec],zLevels)

  call mesh % ExtrudedMesh(surface,zLevels, &
                           bottomBC=SELF_BC_NONORMALFLOW, &
                           topBC=SELF_BC_PRESCRIBED)
```

### Element ordering and domain decomposition

Elements are stored column by column.
Local element $e_3 = (e_2-1)\, n_z + k$ is layer $k$ of the column above local surface element $e_2$, with $k = 1$ at the bottom, and global element ids follow the same rule.
The mesh records the layer count in `mesh % nLayers`, and the surface element and layer of an element follow from its id alone, $e_2 = \lfloor (e_3-1)/n_z \rfloor + 1$ and $k = e_3 - (e_2-1)\, n_z$.

The 3-D decomposition is the surface decomposition expanded by $n_z$: rank $r$ owns the columns above the surface elements it owns, and `offsetElem` is $n_z$ times the surface offsets.
A column is therefore never split across ranks, and the local surface element $e_2$ of a column is element $e_2$ of the surface mesh you passed in.
The surface mesh is not copied, so keep it if your column operators need its geometry (see [Column operators](#column-operators)).

The HOHQMesh reader keeps every element on every rank rather than decomposing the mesh, and `ExtrudedMesh` stops with an error when the surface mesh does not hold exactly the elements its decomposition assigns to the rank.
For a parallel run, read the surface from HOPr or build it with one of the structured constructors.

### Connectivity, flips and global side ids

The 2-D sides south, east, north and west (1 to 4) become the lateral faces 2 to 5 of the hexahedron, and faces 1 and 6 are the bottom and the top.
A lateral face takes the neighbor, the neighbor side and the boundary condition id of its surface edge, in the same layer.
The bottom face of layer 1 takes `bottomBC`, the top face of layer $n_z$ takes `topBC`, and the faces between layers are interior.

The flip of a lateral face follows from its trace coordinates.
`BoundaryInterp` orders the nodes of the south and north faces by $(\xi, \zeta)$ and those of the east and west faces by $(\eta, \zeta)$, so on both sides of a lateral face the first trace coordinate is the coordinate along the surface edge and the second is $\zeta$, which points upward in both elements.
A surface edge flip of 0 therefore becomes face flip 0, and a surface edge flip of 1, an edge whose two sides run in opposite directions, becomes face flip 1, which reverses the first index and keeps the second.
The faces between layers have flip 0.

Global side ids are the MPI message tags of the side exchange, so they must be unique and agree across ranks.
The lateral face of layer $k$ on the surface edge with id $g$ has id $\mathrm{sign}(g)\big((|g|-1)\, n_z + k\big)$, which keeps the sign HOPr uses to mark the two sides of an edge.
The $n_z + 1$ horizontal faces of the column above global surface element $E$ are numbered after all lateral faces, as $n_z\, n_{s,2D} + (E-1)(n_z+1) + l$ for level $l-1$.
Here $n_{s,2D}$ is the `nUniqueSides` of the surface mesh, and the extruded mesh has

$$
n_{s,3D} = n_z\, n_{s,2D} + n_{E,2D}\,(n_z + 1)
$$

unique side ids, where $n_{E,2D}$ is the global number of surface elements.

### Metric identities

For the extruded map the exact metric terms are

$$
J\vec{a}^\xi = Z_\zeta\,(Y_\eta,\,-X_\eta,\,0), \qquad
J\vec{a}^\eta = Z_\zeta\,(-Y_\xi,\,X_\xi,\,0), \qquad
J\vec{a}^\zeta = \big(Y_\xi Z_\eta - Z_\xi Y_\eta,\ Z_\xi X_\eta - X_\xi Z_\eta,\ J_2\big),
$$

with $J_2 = X_\xi Y_\eta - X_\eta Y_\xi$ the Jacobian of the surface element.
`SEMHex` computes $J\vec{a}^i$ in Kopriva's curl-invariant form and $J$ as the determinant of the covariant basis.
The curl form differentiates products of two interpolated factors, and the products are of degree $2N_{geo}$ in each horizontal reference coordinate.
For instance, the horizontal components of $J\vec{a}^\xi$ contain $\partial_\eta(Y Z_\zeta)$, which the derivative matrix reproduces only when $Y Z_\zeta$ is resolved.
We therefore have the following, which the test `geometry_3d_extruded_identities` checks.

- With $N \ge 2N_{geo}$, the discrete metric terms satisfy the three identities to round-off: $J\vec{a}^\xi$ and $J\vec{a}^\eta$ have no vertical component and equal $Z_\zeta$ times the surface contravariant vectors, $J\vec{a}^\zeta_z = J_2$, and $\sum_i J\vec{a}^i D_i Z = J\hat{z}$.
- On an element with flat levels, $Z_\zeta$ is constant over the element, so the statements about $J\vec{a}^\xi$ and $J\vec{a}^\eta$ and the horizontal components of $\sum_i J\vec{a}^i D_i Z$ hold to round-off for every $N$.
- The vertical statements, $J\vec{a}^\zeta_z = J_2$ and the vertical component of $\sum_i J\vec{a}^i D_i Z = J\hat{z}$, compare the curl form of $J_2$ with the determinant $J_2$, and they hold to round-off only for $N \ge 2N_{geo}$ or an affine surface element, also on flat levels.
- Three weaker properties hold to round-off for every $N$, because each involves only products that are linear in $\zeta$ or constant along a vertical grid line: the vertical components of $J\vec{a}^\xi$ and $J\vec{a}^\eta$ vanish; the stored $J\vec{a}^\zeta_z$ is the same at every node of a vertical grid line, inside an element and from layer to layer; and the slope components $J\vec{a}^\zeta_x$, $J\vec{a}^\zeta_y$ are continuous across the faces between layers.
- Otherwise the identities hold to truncation error. On the curved Circle surface ($N_{geo} = 7$) with sigma levels at $N = 7$, the measured relative error of the horizontal statements is about $2\times 10^{-9}$ and that of the vertical ones about $3\times 10^{-11}$.

### Column operators

The module `SELF_ColumnOperators_3D` provides the operators that act along the vertical grid lines of an extruded mesh, the lines of constant $(\xi,\eta)$, which are the third index of SELF's `interior` arrays.
A column is never split across ranks, so none of them communicates.
Volume fields are `Scalar3D` or `MappedScalar3D` on the extruded mesh, and surface fields are `Scalar2D` or `MappedScalar2D` with one element per column, the same degree and the same number of variables.
The operators read only nodal values, so neither field needs a geometry of its own; the metric factor they need comes from the `SEMHex` passed to `Init`.

Write $\omega_l$ for the quadrature weights, $M = \mathrm{diag}(\omega)$, $D$ for the derivative matrix, $e_0$ and $e_N$ for the unit vectors of the bottom and top node, and $Z_\zeta = \partial z/\partial\zeta$, the entry `dxds % interior(i,j,l,iel,1,3,3)` of the covariant basis.
The column integral and its adjoint, the broadcast, are

$$
(\mathcal{C} q)_{ij} = \sum_{k=1}^{n_z} \sum_{l=0}^{N} \omega_l\, (Z_\zeta)^{(k)}_{ijl}\, q^{(k)}_{ijl}, \qquad
(\mathcal{B} Q)^{(k)}_{ijl} = Q_{ij} .
$$

The column integral of a field in units of $[q]$ is in units of $[q]$ metres, and $\mathcal{C}\mathbf{1}$ is the local depth.
Since the Jacobian of an extruded element is $J = J_2 Z_\zeta$ at every node, $\langle \mathcal{B}A, q\rangle_3 = \langle A, \mathcal{C}q\rangle_2$, where the inner products are weighted by $J\omega_i\omega_j\omega_l$ and $J_2\omega_i\omega_j$.
This holds to round-off both with the curl-form Jacobian $J^c = (J\vec{a}^\zeta)_z Z_\zeta$ and with the determinant that `SEMHex` stores, since the covariant basis gives $\det J = Z_\zeta \det J_2$ exactly.

The vertical sweeps solve the one dimensional upwind DG problem $\partial_\zeta q = r$ along every vertical line, layer by layer.
Upward, from a bottom trace, and downward, from a top trace,

$$
(D + M^{-1}e_0e_0^T)\, q^{(k)} = r^{(k)} + M^{-1}e_0\, q^{(k)}_{in}, \quad q^{(k)}_{in} = q^{(k-1)}_N ; \qquad
(D - M^{-1}e_Ne_N^T)\, q^{(k)} = r^{(k)} - M^{-1}e_N\, q^{(k)}_{in}, \quad q^{(k)}_{in} = q^{(k+1)}_0 .
$$

Both matrices map the constant vector to the inflow term, so the solution is the inflow value plus an element-local particular solution with a matrix inverse that `Init` computes once.
The right hand side $r$ is in reference scaling and the sweeps apply no metric factor: to solve $J^c_2\, \partial_\zeta w = f$, pass $r = f/J^c_2$.
The sweeps need nodes at $\zeta = \pm 1$, so they require Gauss-Lobatto nodes and stop with an error otherwise; the column integral and the broadcast work with any quadrature.
By the summation-by-parts property of Gauss-Lobatto nodes the two sweeps are adjoint: with $w$ from the upward sweep and $p$ from the downward sweep, on every vertical line,

$$
\sum_{k=1}^{n_z} \big( p^{(k)T} M r_w^{(k)} + w^{(k)T} M r_p^{(k)} \big) = p^{(n_z)}_{in}\, w^{(n_z)}_N - p^{(1)}_0\, w^{(1)}_{in} .
$$

```fortran
  type(ColumnOperators3D) :: columns
  type(MappedScalar3D) :: q,r,w
  type(Scalar2D) :: qbar,wBottom

  ! interp must use GAUSS_LOBATTO nodes for the sweeps
  call columns % Init(interp,mesh,geometry)
  call qbar % Init(interp,q % nVar,surface % nElem)

  call columns % ColumnIntegral(q,qbar)      ! qbar = C q
  call columns % Broadcast(qbar,q)           ! q = B qbar
  call columns % VerticalSweepUp(r,wBottom,w)
  call columns % Free()
```

On GPU builds the operators read and write the device copies of the fields, so upload the inputs with `UpdateDevice` and download the results with `UpdateHost` when you need them on the host.
The tests `column_integral_3d_extruded`, `column_broadcast_adjoint_3d` and `vertical_sweep_3d` check the exact integrals, the adjoint identities, and the exactness and inverse accuracy of the sweeps, on one and on two ranks.

### Relation to HOHQMesh extrusion

HOHQMesh can also sweep a quad mesh upward over a bottom topography `f(x,y)` to a flat top at height $H$ (the `SIMPLE_EXTRUSION` and `TOPOGRAPHY` blocks).
It spreads the topography over the whole column, placing level $j$ at $z_j = f\,(1 - j/n_z) + H\, j/n_z$, and it writes the elements layer by layer.
With $z' = z - H$ and $h = H - f$ these are the sigma levels $\sigma_j = j/n_z - 1$, so `ExtrudedSigmaLevels` with those $\sigma_j$ reproduces HOHQMesh's corner nodes and bottom topography, after reordering its elements column by column.
The interior nodes differ: the ISM file carries curved face data only for faces on the topography or on a curved boundary, so the faces of the upper layers inside the domain are rebuilt from their corners when the file is read.
The `share/mesh/Seamount` fixtures and the test `mesh3d_extruded_hohqmesh` show both.
`ExtrudedMesh` builds the mesh in memory instead, which is what an adaptive run needs when it rebuilds the columns above a refined surface.

### Current limitations

A surface mesh with nonconforming (mortar) edges is rejected.
A lateral face between columns at different refinement levels needs an anisotropic mortar, the 2-D mortar in the horizontal tensored with the identity in $\zeta$, which SELF does not provide yet.
