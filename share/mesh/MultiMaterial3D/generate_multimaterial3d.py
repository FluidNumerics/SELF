#!/usr/bin/env python3
"""Generate 3-D HOHQMesh ISM / ISM-MM hexahedral test meshes for SELF.

This script reproduces the output of HOHQMesh's simple-extrusion pipeline
(`SimpleSweep.f90` -> `WriteISMHexMeshFile` in
`Source/3DSource/Mesh3DOutputMethods.f90`) for two small meshes:

1. ``InsulatedWire3D.mesh`` (ISM-MM, written next to this script)

   A two-material "insulated wire": the cross section is a classic O-grid
   with a 2x2 square core (material ``Copper``) surrounded by an 8-element
   ring (material ``Insulator``) whose outer boundary is the unit circle.
   The cross section is extruded in z over ``NZ`` layers. The outer
   cylinder faces are curved (polynomial order ``N``), so the mesh
   exercises SELF's curved-face transfinite interpolation as well as the
   ISM-MM material table.

2. ``Block3D.mesh`` (plain ISM, written to ../Block3D/Block3D.mesh)

   A single-material 2x2x2 unit cube with linear geometry (polyOrder 1),
   used to test the plain-ISM (no material string) path of the reader.

Format notes (matching the HOHQMesh writer exactly):

* 3-D ISM files have NO header line; the first line is the count line
  ``nNodes nElems polyOrder``. ISM-MM appends the element's material name
  to its corner-node line -- that is the only difference from plain ISM.
* Nodes are written layer-by-layer (all nodes of the bottom layer first).
* Elements are written layer-by-layer, quad-element-major within a layer.
* Per element: corner-node line (nodes 1-4 = bottom quad corners CCW,
  5-8 = the same corners on the layer above), a line of six boundary-face
  flags, one (N+1)^2 block of face points (x y z per line, first face
  index fastest) for every flagged face, and a line of six boundary names
  (``---`` marks an interior face).
* HOHQMesh hex face numbering: 1=south (eta=0), 2=north (eta=1),
  3=bottom (zeta=0), 4=east (xi=1), 5=top (zeta=1), 6=west (xi=0).
  Quad-edge k maps to hex face flagMap[k] = [1, 4, 2, 6] (k = 1..4).
* Face points are sampled at Chebyshev-Gauss-Lobatto points of degree N.
* Flag semantics from `sweepElements`: every element of the first layer
  flags faces 1,2,3,4,6; an element whose quad has a curved edge k flags
  faces flagMap[k], 3 and 5.
"""

import math
import os

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# Chebyshev-Gauss-Lobatto points mapped to [0,1]: u_i = (1 - cos(pi*i/N))/2


def cgl01(n):
    return [0.5 * (1.0 - math.cos(math.pi * i / n)) for i in range(n + 1)]


def lerp(p, q, t):
    return tuple((1.0 - t) * a + t * b for a, b in zip(p, q))


def fmt(x):
    return f"{x: .16E}"


class QuadTFI:
    """2-D transfinite (Coons) map of a quad with one optionally curved
    edge-2 (east) arc; edges are HOHQMesh-ordered: 1:(c1,c2), 2:(c2,c3),
    3:(c4,c3), 4:(c1,c4). Corners c1..c4 are CCW."""

    def __init__(self, corners, arc=None):
        self.c = corners  # list of 4 (x, y)
        self.arc = arc  # None or callable t->(x,y), arc(0)=c2, arc(1)=c3

    def edge1(self, s):
        return lerp(self.c[0], self.c[1], s)

    def edge2(self, t):
        if self.arc is not None:
            return self.arc(t)
        return lerp(self.c[1], self.c[2], t)

    def edge3(self, s):
        return lerp(self.c[3], self.c[2], s)

    def edge4(self, t):
        return lerp(self.c[0], self.c[3], t)

    def __call__(self, s, t):
        c1, c2, c3, c4 = self.c
        e1 = self.edge1(s)
        e3 = self.edge3(s)
        e4 = self.edge4(t)
        e2 = self.edge2(t)
        return tuple(
            (1.0 - t) * e1[d] + t * e3[d] + (1.0 - s) * e4[d] + s * e2[d]
            - ((1.0 - s) * (1.0 - t) * c1[d] + s * (1.0 - t) * c2[d]
               + s * t * c3[d] + (1.0 - s) * t * c4[d])
            for d in range(2)
        )


class QuadElement:
    def __init__(self, nodeids, tfi, material, edge_curved, edge_names):
        self.nodeids = nodeids  # 4 corner node ids (1-based, 2-D table)
        self.tfi = tfi
        self.material = material
        self.edge_curved = edge_curved  # [e1, e2, e3, e4] booleans
        self.edge_names = edge_names  # 4 boundary names or "---"


def face_points(quad, z0, z1, hohq_face, u):
    """Sample the volume mapping X(xi,eta,zeta) = (Q(xi,eta), z0+(z1-z0)*w)
    on one HOHQMesh face, in the order written by FaceFromVolume:
    xFace(:,ii,jj), ii fastest."""
    n1 = len(u)
    pts = []
    for jj in range(n1):
        for ii in range(n1):
            if hohq_face == 1:  # south: (xi, zeta)
                x, y = quad.tfi(u[ii], 0.0)
                z = z0 + (z1 - z0) * u[jj]
            elif hohq_face == 2:  # north: (xi, zeta)
                x, y = quad.tfi(u[ii], 1.0)
                z = z0 + (z1 - z0) * u[jj]
            elif hohq_face == 3:  # bottom: (xi, eta)
                x, y = quad.tfi(u[ii], u[jj])
                z = z0
            elif hohq_face == 4:  # east: (eta, zeta)
                x, y = quad.tfi(1.0, u[ii])
                z = z0 + (z1 - z0) * u[jj]
            elif hohq_face == 5:  # top: (xi, eta)
                x, y = quad.tfi(u[ii], u[jj])
                z = z1
            elif hohq_face == 6:  # west: (eta, zeta)
                x, y = quad.tfi(0.0, u[ii])
                z = z0 + (z1 - z0) * u[jj]
            pts.append((x, y, z))
    return pts


FLAGMAP = {1: 1, 2: 4, 3: 2, 4: 6}  # quad edge -> hex face


def write_mesh(path, nodes2d, quads, zlevels, order, ism_mm,
               bottom_name="bottom", top_name="top"):
    n1 = order + 1
    u = cgl01(order)
    nlayers = len(zlevels) - 1
    n2d = len(nodes2d)
    lines = []
    lines.append(f" {n2d*len(zlevels)} {len(quads)*nlayers} {order}")
    # Nodes, layer by layer
    for z in zlevels:
        for (x, y) in nodes2d:
            lines.append(f" {fmt(x)} {fmt(y)} {fmt(z)}")
    # Elements, layer by layer
    for j in range(1, nlayers + 1):
        z0 = zlevels[j - 1]
        z1 = zlevels[j]
        for q in quads:
            bottom_ids = [nid + n2d * (j - 1) for nid in q.nodeids]
            top_ids = [nid + n2d * j for nid in q.nodeids]
            ids = bottom_ids + top_ids
            corner_line = " " + " ".join(str(i) for i in ids)
            if ism_mm:
                corner_line += f" {q.material}"
            lines.append(corner_line)

            flags = [0] * 7  # 1-based
            names = ["---"] * 7  # 1-based
            if j == 1:
                for f in (1, 2, 4, 6):
                    flags[f] = 1
                flags[3] = 1
                names[3] = bottom_name
            if j == nlayers:
                names[5] = top_name
            for k in range(1, 5):
                if q.edge_curved[k - 1]:
                    flags[FLAGMAP[k]] = 1
                    flags[3] = 1
                    flags[5] = 1
                names[FLAGMAP[k]] = q.edge_names[k - 1]

            lines.append(" " + " ".join(str(flags[f]) for f in range(1, 7)))
            for f in range(1, 7):
                if flags[f] == 1:
                    for (x, y, z) in face_points(q, z0, z1, f, u):
                        lines.append(f" {fmt(x)} {fmt(y)} {fmt(z)}")
            lines.append(" " + " ".join(names[f] for f in range(1, 7)))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {path}: {n2d*len(zlevels)} nodes, "
          f"{len(quads)*nlayers} elements, order {order}")


def insulated_wire():
    """O-grid cross section: 2x2 square core [-0.5,0.5]^2 (Copper) plus an
    8-element ring out to the unit circle (Insulator)."""
    order = 4
    radius = 1.0
    # Square nodes: 3x3 grid, id = 1 + ix + 3*iy
    sq = {}
    nodes = []
    for iy in range(3):
        for ix in range(3):
            x = -0.5 + 0.5 * ix
            y = -0.5 + 0.5 * iy
            nodes.append((x, y))
            sq[(ix, iy)] = len(nodes)
    # Square boundary nodes CCW from angle 0 (mid-east)
    bd = [(2, 1), (2, 2), (1, 2), (0, 2), (0, 1), (0, 0), (1, 0), (2, 0)]
    bd_ids = [sq[p] for p in bd]
    angles = [math.pi * b / 4.0 for b in range(8)]
    circ_ids = []
    for th in angles:
        nodes.append((radius * math.cos(th), radius * math.sin(th)))
        circ_ids.append(len(nodes))

    quads = []
    # Core quads (Copper), straight edges, interior everywhere
    for iy in range(2):
        for ix in range(2):
            ids = [sq[(ix, iy)], sq[(ix + 1, iy)],
                   sq[(ix + 1, iy + 1)], sq[(ix, iy + 1)]]
            corners = [nodes[i - 1] for i in ids]
            quads.append(QuadElement(ids, QuadTFI(corners), "Copper",
                                     [False] * 4, ["---"] * 4))
    # Ring quads (Insulator): corner1 = inner at angle b, corner2 = outer at
    # angle b, corner3 = outer at angle b+1, corner4 = inner at angle b+1.
    # Edge 2 (corner2 -> corner3) is the outer circular arc.
    for b in range(8):
        bp = (b + 1) % 8
        ids = [bd_ids[b], circ_ids[b], circ_ids[bp], bd_ids[bp]]
        corners = [nodes[i - 1] for i in ids]
        th0 = angles[b]
        th1 = angles[b] + math.pi / 4.0

        def arc(t, th0=th0, th1=th1):
            th = th0 + t * (th1 - th0)
            return (radius * math.cos(th), radius * math.sin(th))

        quads.append(QuadElement(ids, QuadTFI(corners, arc), "Insulator",
                                 [False, True, False, False],
                                 ["---", "cylinder", "---", "---"]))

    zlevels = [0.0, 0.75, 1.5]
    write_mesh(os.path.join(OUTDIR, "InsulatedWire3D.mesh"),
               nodes, quads, zlevels, order, ism_mm=True)


def block3d():
    """Plain ISM 2x2x2 unit cube, linear geometry, single material."""
    order = 1
    sq = {}
    nodes = []
    for iy in range(3):
        for ix in range(3):
            nodes.append((0.5 * ix, 0.5 * iy))
            sq[(ix, iy)] = len(nodes)
    edge_bcs = {  # quad edge -> name for boundary edges of the unit square
        1: "south", 2: "east", 3: "north", 4: "west"}
    quads = []
    for iy in range(2):
        for ix in range(2):
            ids = [sq[(ix, iy)], sq[(ix + 1, iy)],
                   sq[(ix + 1, iy + 1)], sq[(ix, iy + 1)]]
            corners = [nodes[i - 1] for i in ids]
            names = ["---"] * 4
            if iy == 0:
                names[0] = edge_bcs[1]
            if ix == 1:
                names[1] = edge_bcs[2]
            if iy == 1:
                names[2] = edge_bcs[3]
            if ix == 0:
                names[3] = edge_bcs[4]
            quads.append(QuadElement(ids, QuadTFI(corners), "default",
                                     [False] * 4, names))
    zlevels = [0.0, 0.5, 1.0]
    write_mesh(os.path.join(OUTDIR, "..", "Block3D", "Block3D.mesh"),
               nodes, quads, zlevels, order, ism_mm=False)


if __name__ == "__main__":
    insulated_wire()
    block3d()
