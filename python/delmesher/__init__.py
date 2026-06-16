"""delmesher - robust tetrahedral meshing with surface chamfering.

This package wraps the C++ ``delmesher`` pipeline (Diazzi et al., *Surface
Chamfering for Robust Tetrahedral Meshing*, SIGGRAPH 2026). It turns a
triangulated input surface into a high-quality tetrahedral mesh of its interior
and is robust to sharp/acute angles and near-degenerate triangulations.

The main entry point is :func:`tetrahedralize`, which mirrors the command-line
tool ``delmesher`` (one keyword argument per CLI flag) but takes and returns
NumPy arrays instead of OFF / ``.tet`` files.

Example
-------
>>> import numpy as np, delmesher
>>> v, f = delmesher.read_off("input_models/boeing_part.off")
>>> result = delmesher.tetrahedralize(v, f, max_vertices=2000)
>>> result.vertices.shape, result.tetrahedra.shape
((..., 3), (..., 4))
"""

from __future__ import annotations

import numpy as np

from ._delmesher import Result, _tetrahedralize

__all__ = [
    "tetrahedralize",
    "Result",
    "read_off",
    "write_off",
    "write_tet",
]

try:  # populated by the build backend; optional at runtime
    from ._version import __version__
except Exception:  # pragma: no cover
    __version__ = "0.0.0"


# Expose the two "primary" meshes on Result with friendly aliases:
# by default (full pipeline) the enriched CDT is the final mesh; with
# enriched_cdt=False the angle-bounded Delaunay-refined mesh is the output.
def _result_vertices(self: Result) -> np.ndarray:
    """Vertices of the primary output mesh (enriched CDT if computed, else DR)."""
    return self.cdt_vertices if self.enriched_cdt_computed else self.dr_vertices


def _result_tetrahedra(self: Result) -> np.ndarray:
    """Tetrahedra of the primary output mesh (enriched CDT if computed, else DR)."""
    return self.cdt_tetrahedra if self.enriched_cdt_computed else self.dr_tetrahedra


Result.vertices = property(_result_vertices)
Result.tetrahedra = property(_result_tetrahedra)


def tetrahedralize(
    vertices,
    triangles,
    *,
    enriched_cdt: bool = True,
    sliver_removal: bool = False,
    lfs_exponent: int | None = None,
    max_vertices: int | None = None,
    compute_lfs: bool = False,
    verbose: bool = False,
) -> Result:
    """Tetrahedralize a triangulated surface.

    The input surface may be open or closed; it should be free of degenerate or
    self-intersecting triangles (duplicate vertices and triangles are removed
    automatically, matching the CLI). The pipeline runs
    *chamfering → Delaunay refinement → enriched CDT*.

    Parameters
    ----------
    vertices : array_like, shape (n, 3)
        Vertex coordinates. Converted to ``float64``.
    triangles : array_like, shape (m, 3)
        Triangle vertex indices into ``vertices``. Converted to ``uint32``.
    enriched_cdt : bool, default True
        Run the final enriched-CDT phase, making the output exactly conformal to
        the input surface. Set ``False`` to stop after Delaunay refinement (the
        angle-bounded mesh with min face angle > 14 degrees). CLI flag: ``-a``
        corresponds to ``enriched_cdt=False``.
    sliver_removal : bool, default False
        Enable sliver removal during Delaunay refinement. CLI: ``-d``.
    lfs_exponent : int, optional
        Ignore input constraints whose local feature size is below
        ``10**-lfs_exponent`` times the bounding-box diagonal (larger values
        push the output towards a plain CDT of the input). CLI: ``-e``.
    max_vertices : int, optional
        Stop Delaunay refinement after this many inserted vertices. ``None``
        runs to convergence (can take minutes). CLI: ``-m``.
    compute_lfs : bool, default False
        Compute the minimum local feature size of the input, reported as
        :attr:`Result.min_lfs`. CLI: ``-c``.
    verbose : bool, default False
        Print progress to stdout. CLI: ``-v``.

    Returns
    -------
    Result
        Holds the output meshes as NumPy arrays. See :attr:`Result.vertices`,
        :attr:`Result.tetrahedra` for the primary output, plus the per-phase
        meshes (``dr_*``, ``cdt_*``, ``chamfered_*``, ``dr_interface_*``).
    """
    v = np.ascontiguousarray(vertices, dtype=np.float64)
    f = np.ascontiguousarray(triangles, dtype=np.uint32)

    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"vertices must have shape (n, 3), got {v.shape}")
    if f.ndim != 2 or f.shape[1] != 3:
        raise ValueError(f"triangles must have shape (m, 3), got {f.shape}")
    if v.shape[0] == 0:
        raise ValueError("input has no vertices")
    if f.shape[0] == 0:
        raise ValueError("input has no triangles")
    if f.size and f.max() >= v.shape[0]:
        raise ValueError(
            f"triangle index {int(f.max())} out of range for {v.shape[0]} vertices"
        )
    if lfs_exponent is not None and lfs_exponent < 0:
        raise ValueError("lfs_exponent must be a non-negative integer")
    if max_vertices is not None and max_vertices < 0:
        raise ValueError("max_vertices must be a non-negative integer")

    return _tetrahedralize(
        v,
        f,
        enriched_cdt=enriched_cdt,
        sliver_removal=sliver_removal,
        lfs_exponent=-1 if lfs_exponent is None else int(lfs_exponent),
        max_vertices=-1 if max_vertices is None else int(max_vertices),
        compute_lfs=compute_lfs,
        verbose=verbose,
    )


# ---------------------------------------------------------------------------
# Lightweight OFF / .tet I/O. The pipeline re-cleans whatever vertices it is
# given (dedup + reorder), so read_off only needs to parse the raw file exactly
# as the C++ read_OFF_file does (no comment support, triangles only).
# ---------------------------------------------------------------------------
def read_off(path) -> tuple[np.ndarray, np.ndarray]:
    """Read an ASCII OFF triangle surface.

    Returns ``(vertices, triangles)`` as ``float64`` (n, 3) and ``uint32``
    (m, 3) arrays. Only triangular faces are supported, matching delmesher.
    """
    with open(path, "r") as fh:
        text = fh.read()
    tokens = _off_tokens(text)
    if next(tokens) != "OFF":
        raise ValueError(f"{path}: not an OFF file (missing 'OFF' header)")
    n_verts = int(next(tokens))
    n_faces = int(next(tokens))
    int(next(tokens))  # edge count, ignored (as in read_OFF_file)

    verts = np.empty((n_verts, 3), dtype=np.float64)
    for i in range(n_verts):
        verts[i, 0] = float(next(tokens))
        verts[i, 1] = float(next(tokens))
        verts[i, 2] = float(next(tokens))

    tris = np.empty((n_faces, 3), dtype=np.uint32)
    for i in range(n_faces):
        nv = int(next(tokens))
        if nv != 3:
            raise ValueError(f"{path}: non-triangular face not supported")
        tris[i, 0] = int(next(tokens))
        tris[i, 1] = int(next(tokens))
        tris[i, 2] = int(next(tokens))
    return verts, tris


def _off_tokens(text):
    # The 'OFF' tag may sit on its own line or be glued to the counts; split on
    # whitespace and also peel a leading 'OFF' off the first token.
    first = True
    for tok in text.split():
        if first and tok.startswith("OFF") and tok != "OFF":
            yield "OFF"
            tok = tok[3:]
            first = False
            if not tok:
                continue
        first = False
        yield tok


def write_off(path, vertices, faces) -> None:
    """Write a triangle/polygon surface to an ASCII OFF file."""
    v = np.ascontiguousarray(vertices, dtype=np.float64).reshape(-1, 3)
    f = np.ascontiguousarray(faces, dtype=np.int64)
    with open(path, "w") as fh:
        fh.write("OFF\n")
        fh.write(f"{len(v)} {len(f)} 0\n")
        for x, y, z in v:
            fh.write(f"{x:f} {y:f} {z:f}\n")
        for face in f:
            fh.write(f"{len(face)} " + " ".join(str(int(i)) for i in face) + "\n")


def write_tet(path, vertices, tetrahedra) -> None:
    """Write a tetrahedral mesh to a ``.tet`` file (delmesher's format).

    The format matches ``Tetrahedrization::saveTET``: a ``"<n> vertices"`` line,
    a ``"<m> tets"`` line, the vertex coordinates, then one ``"4 i j k l"`` line
    per tetrahedron. Coordinates use ``%f`` (six decimals), as in the C++ writer.
    """
    v = np.ascontiguousarray(vertices, dtype=np.float64).reshape(-1, 3)
    t = np.ascontiguousarray(tetrahedra, dtype=np.int64).reshape(-1, 4)
    with open(path, "w") as fh:
        fh.write(f"{len(v)} vertices\n")
        fh.write(f"{len(t)} tets\n")
        for x, y, z in v:
            fh.write(f"{x:f} {y:f} {z:f}\n")
        for a, b, c, d in t:
            fh.write(f"4 {int(a)} {int(b)} {int(c)} {int(d)}\n")
