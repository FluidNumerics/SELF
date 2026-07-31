#!/usr/bin/env python3
"""Generate small synthetic 3-D ISM / ISM-MM fixture meshes for the SELF
`Read_HOHQMesh` (Mesh3D) reader tests.

Unlike ``generate_multimaterial3d.py`` (which mimics HOHQMesh's extrusion
output), these fixtures deliberately exercise reader paths that HOHQMesh's
own writer rarely or never produces:

1. ``ReaderFixture_ManyMaterials.mesh`` (ISM-MM, polyOrder 1)

   A 3x3x1 block of unit cubes where every element carries a unique
   material name (9 materials) and every physical boundary face carries a
   unique boundary name (30 names). This overflows the reader's initial
   material-name (8) and boundary-name (16) tables, exercising the table
   growth paths.

2. ``ReaderFixture_RotatedStack.mesh`` (ISM, polyOrder 2)

   A 1x1x3 stack of unit cubes where the middle element's corner numbering
   is rotated 90 degrees about z and the top element's by 180 degrees, so
   the shared horizontal faces pair with transposed side "flips" — the
   orientation combinations an extruded HOHQMesh mesh can never produce.
   Selected faces are flagged with (bilinear) face-point data in patterns
   HOHQMesh's sweep never emits (sidewalls flagged with bottom/top
   unflagged), exercising the curved-face-preferred edge extraction in the
   transfinite interpolation. All boundary names are "---" so the mesh has
   no named boundary conditions.

3. ``ReaderFixture_BadConnectivity.mesh`` (ISM, polyOrder 1)

   Two stacked cubes whose shared face is listed with corner orderings
   that share the same four corner nodes but are not related by any of
   the eight face orientations (the second element's face is "twisted").
   The reader must detect this and abort; used by a WILL_FAIL guard test.

See ``generate_multimaterial3d.py`` for the format description.
"""

import math
import os

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def cgl01(n):
    return [0.5 * (1.0 - math.cos(math.pi * i / n)) for i in range(n + 1)]


def fmt(x):
    return f"{x: .16E}"


def trilinear(corners, s, t, w):
    """Trilinear map of a hex from its 8 corners (CGNS order: bottom quad
    counter-clockwise, then the top quad above it)."""
    c = corners
    return tuple(
        (1 - s) * (1 - t) * (1 - w) * c[0][d] + s * (1 - t) * (1 - w) * c[1][d]
        + s * t * (1 - w) * c[2][d] + (1 - s) * t * (1 - w) * c[3][d]
        + (1 - s) * (1 - t) * w * c[4][d] + s * (1 - t) * w * c[5][d]
        + s * t * w * c[6][d] + (1 - s) * t * w * c[7][d]
        for d in range(3)
    )


def hohq_face_points(corners, hohq_face, u):
    """Sample a hex's trilinear map on one HOHQMesh face in the order of
    HOHQMesh's FaceFromVolume: faces 1=south(eta=0), 2=north(eta=1),
    3=bottom(zeta=0), 4=east(xi=1), 5=top(zeta=1), 6=west(xi=0); the first
    face index varies fastest."""
    pts = []
    for jj in u:
        for ii in u:
            if hohq_face == 1:
                pts.append(trilinear(corners, ii, 0.0, jj))
            elif hohq_face == 2:
                pts.append(trilinear(corners, ii, 1.0, jj))
            elif hohq_face == 3:
                pts.append(trilinear(corners, ii, jj, 0.0))
            elif hohq_face == 4:
                pts.append(trilinear(corners, 1.0, ii, jj))
            elif hohq_face == 5:
                pts.append(trilinear(corners, ii, jj, 1.0))
            elif hohq_face == 6:
                pts.append(trilinear(corners, 0.0, ii, jj))
    return pts


def write_fixture(path, order, nodes, elements):
    """elements: list of (nodeids[8], material_or_None, flags[6], names[6]),
    flags/names in HOHQMesh face order."""
    lines = [f" {len(nodes)} {len(elements)} {order}"]
    for (x, y, z) in nodes:
        lines.append(f" {fmt(x)} {fmt(y)} {fmt(z)}")
    u = cgl01(order)
    for (ids, material, flags, names) in elements:
        corner_line = " " + " ".join(str(i) for i in ids)
        if material is not None:
            corner_line += f" {material}"
        lines.append(corner_line)
        lines.append(" " + " ".join(str(f) for f in flags))
        corners = [nodes[i - 1] for i in ids]
        for f in range(1, 7):
            if flags[f - 1] == 1:
                for (x, y, z) in hohq_face_points(corners, f, u):
                    lines.append(f" {fmt(x)} {fmt(y)} {fmt(z)}")
        lines.append(" " + " ".join(names))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {path}: {len(nodes)} nodes, {len(elements)} elements, "
          f"order {order}")


def grid_nodes(nx, ny, zlevels):
    """(nx+1) x (ny+1) x len(zlevels) unit-spaced nodes, x fastest."""
    nodes = []
    for z in zlevels:
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                nodes.append((float(ix), float(iy), z))
    return nodes


def many_materials():
    nx = ny = 3
    nodes = grid_nodes(nx, ny, [0.0, 1.0])
    nid = lambda ix, iy, lv: 1 + ix + (nx + 1) * iy + (nx + 1) * (ny + 1) * lv
    elements = []
    e = 0
    for iy in range(ny):
        for ix in range(nx):
            e += 1
            ids = [nid(ix, iy, 0), nid(ix + 1, iy, 0),
                   nid(ix + 1, iy + 1, 0), nid(ix, iy + 1, 0),
                   nid(ix, iy, 1), nid(ix + 1, iy, 1),
                   nid(ix + 1, iy + 1, 1), nid(ix, iy + 1, 1)]
            # HOHQMesh face order: 1=south, 2=north, 3=bottom, 4=east,
            # 5=top, 6=west
            names = ["---"] * 6
            names[2] = f"bottom{e}"
            names[4] = f"top{e}"
            if iy == 0:
                names[0] = f"south{ix+1}"
            if iy == ny - 1:
                names[1] = f"north{ix+1}"
            if ix == nx - 1:
                names[3] = f"east{iy+1}"
            if ix == 0:
                names[5] = f"west{iy+1}"
            elements.append((ids, f"mat{e}", [0] * 6, names))
    write_fixture(os.path.join(OUTDIR, "ReaderFixture_ManyMaterials.mesh"),
                  1, nodes, elements)


def rotated_stack():
    # 4 nodes per level (unit square, counter-clockwise SW,SE,NE,NW),
    # 4 levels (z = 0..3)
    nodes = []
    for lv in range(4):
        nodes += [(0.0, 0.0, float(lv)), (1.0, 0.0, float(lv)),
                  (1.0, 1.0, float(lv)), (0.0, 1.0, float(lv))]
    n = lambda lv, c: 4 * lv + c  # c = 1(SW), 2(SE), 3(NE), 4(NW)

    # Element 1: unrotated
    e1 = [n(0, 1), n(0, 2), n(0, 3), n(0, 4),
          n(1, 1), n(1, 2), n(1, 3), n(1, 4)]
    # Element 2: corner numbering rotated 90 degrees (corner 1 at SE)
    e2 = [n(1, 2), n(1, 3), n(1, 4), n(1, 1),
          n(2, 2), n(2, 3), n(2, 4), n(2, 1)]
    # Element 3: rotated 180 degrees (corner 1 at NE)
    e3 = [n(2, 3), n(2, 4), n(2, 1), n(2, 2),
          n(3, 3), n(3, 4), n(3, 1), n(3, 2)]

    # Flags in HOHQMesh face order [S, N, B, E, T, W]: element 1 flags its
    # south/north/east sidewalls with bottom and top unflagged; element 2
    # flags only its (local) west face — combinations HOHQMesh's sweep
    # never writes, exercising the flagged-face-preferred edge extraction.
    elements = [
        (e1, None, [1, 1, 0, 1, 0, 0], ["---"] * 6),
        (e2, None, [0, 0, 0, 0, 0, 1], ["---"] * 6),
        (e3, None, [0] * 6, ["---"] * 6),
    ]
    write_fixture(os.path.join(OUTDIR, "ReaderFixture_RotatedStack.mesh"),
                  2, nodes, elements)


def bad_connectivity():
    nodes = []
    for lv in range(3):
        nodes += [(0.0, 0.0, float(lv)), (1.0, 0.0, float(lv)),
                  (1.0, 1.0, float(lv)), (0.0, 1.0, float(lv))]
    # Element 2 lists the shared face nodes (5,6,7,8) in a "twisted" order
    # (5,7,6,8) that no face orientation can produce.
    elements = [
        ([1, 2, 3, 4, 5, 6, 7, 8], None, [0] * 6, ["---"] * 6),
        ([5, 7, 6, 8, 9, 11, 10, 12], None, [0] * 6, ["---"] * 6),
    ]
    write_fixture(os.path.join(OUTDIR, "ReaderFixture_BadConnectivity.mesh"),
                  1, nodes, elements)


if __name__ == "__main__":
    many_materials()
    rotated_stack()
    bad_connectivity()
