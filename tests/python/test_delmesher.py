"""Python test suite for the delmesher bindings.

This mirrors the Catch2 CLI suite (tests/test_delmesher_cli.cpp): it runs the
meshing pipeline on every model in input_models/, switching on each accepted
option one at a time, and checks the run succeeds (the binding equivalent of
"the binary exits 0"). It additionally verifies that the bindings produce the
*same* mesh as the command-line tool, by comparing the returned NumPy arrays to
the .tet files the CLI writes for the same input.

Refinement is capped (DELMESHER_TEST_MAX_VERTICES, default 1000) so the suite
stays fast while still exercising every phase, exactly like the C++ suite.
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import tempfile

import numpy as np
import pytest

import delmesher

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
INPUT_DIR = REPO_ROOT / "input_models"

# Vertex cap (-m) applied to every run; "" means unbounded (very slow).
_cap_env = os.environ.get("DELMESHER_TEST_MAX_VERTICES", "1000")
CAP = int(_cap_env) if _cap_env.strip() else None


def discover_models():
    return sorted(INPUT_DIR.glob("*.off"))


def find_cli():
    """Locate the delmesher CLI for the parity tests, or None to skip them."""
    candidates = []
    env = os.environ.get("DELMESHER_BINARY")
    if env:
        candidates.append(pathlib.Path(env))
    candidates += [
        REPO_ROOT / "build" / "delmesher",
        REPO_ROOT / "build" / "Release" / "delmesher.exe",
        REPO_ROOT / "build" / "delmesher.exe",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


CLI = find_cli()
MODELS = discover_models()

# Option cases mirroring the C++ flag_cases(), restricted to flags that change
# what is computed (the file-output / logging flags -b -u -w -x -y -z -h -l do
# not alter the returned meshes and are covered by inspecting the result).
OPTION_CASES = {
    "no-flags": {},
    "a-skip-enriched": {"enriched_cdt": False},
    "c-min-lfs": {"compute_lfs": True},
    "d-sliver-removal": {"sliver_removal": True},
    "e-lfs-bound": {"lfs_exponent": 8},
    "m-max-vertices": {"max_vertices": 2000},
    "v-verbose": {"verbose": True},
}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def assert_valid_tet_mesh(vertices, tetrahedra):
    assert vertices.ndim == 2 and vertices.shape[1] == 3
    assert vertices.dtype == np.float64
    assert tetrahedra.ndim == 2 and tetrahedra.shape[1] == 4
    assert tetrahedra.dtype == np.uint32
    assert len(vertices) > 0
    assert len(tetrahedra) > 0
    assert np.isfinite(vertices).all()
    assert int(tetrahedra.max()) < len(vertices)
    # no degenerate tets (four distinct corners)
    assert (tetrahedra[:, 0] != tetrahedra[:, 1]).all()


def assert_valid_tri_surface(vertices, faces):
    assert vertices.ndim == 2 and vertices.shape[1] == 3
    assert faces.ndim == 2 and faces.shape[1] == 3
    if len(faces):
        assert int(faces.max()) < len(vertices)


def parse_tet_file(path):
    """Parse a delmesher .tet file into (vertices Nx3 float64, tets Mx4 int).

    Handles both header variants: Tetrahedrization::saveTET writes
    "<n> vertices" / "<m> tets"; TetMesh::saveTET writes "<n> vertices" /
    "<i> inner tets" / "<o> outer tets".
    """
    lines = [ln.strip() for ln in pathlib.Path(path).read_text().splitlines()
             if ln.strip()]
    n_verts = int(lines[0].split()[0])
    i = 1
    n_tets = 0
    while "tets" in lines[i]:
        n_tets += int(lines[i].split()[0])
        i += 1
    verts = np.array([[float(x) for x in lines[i + k].split()[:3]]
                      for k in range(n_verts)], dtype=np.float64)
    tets = np.array([[int(x) for x in lines[i + n_verts + k].split()[1:5]]
                     for k in range(n_tets)], dtype=np.int64)
    return verts, tets


def run_cli(model_path, args, workdir):
    """Run the CLI in workdir on a copy of the model; return its exit code."""
    local = pathlib.Path(workdir) / model_path.name
    local.write_bytes(model_path.read_bytes())
    cmd = [str(CLI)]
    if CAP is not None and "-m" not in args:
        cmd += ["-m", str(CAP)]
    cmd += list(args) + [str(local)]
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    return proc.returncode, proc


# --------------------------------------------------------------------------- #
# basic sanity / IO
# --------------------------------------------------------------------------- #
def test_models_present():
    assert MODELS, "no input models found in input_models/"


@pytest.mark.parametrize("model", MODELS, ids=lambda p: p.name)
def test_off_roundtrip(model, tmp_path):
    v, f = delmesher.read_off(model)
    assert v.shape[1] == 3 and f.shape[1] == 3
    assert v.dtype == np.float64 and f.dtype == np.uint32
    out = tmp_path / "rt.off"
    delmesher.write_off(out, v, f)
    v2, f2 = delmesher.read_off(out)
    assert np.allclose(v, v2, rtol=0, atol=1e-6)
    assert np.array_equal(f, f2)


def test_tet_roundtrip(tmp_path):
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    t = np.array([[0, 1, 2, 3]], dtype=np.uint32)
    out = tmp_path / "x.tet"
    delmesher.write_tet(out, v, t)
    v2, t2 = parse_tet_file(out)
    assert np.allclose(v, v2)
    assert np.array_equal(t, t2)


# --------------------------------------------------------------------------- #
# input validation
# --------------------------------------------------------------------------- #
def test_invalid_inputs():
    v = np.zeros((4, 3))
    with pytest.raises(ValueError):
        delmesher.tetrahedralize(np.zeros((4, 2)), np.zeros((1, 3), np.uint32))
    with pytest.raises(ValueError):
        delmesher.tetrahedralize(v, np.zeros((1, 4), np.uint32))
    with pytest.raises(ValueError):
        delmesher.tetrahedralize(np.zeros((0, 3)), np.zeros((0, 3), np.uint32))
    with pytest.raises(ValueError):  # index out of range
        delmesher.tetrahedralize(v, np.array([[0, 1, 99]], np.uint32))
    with pytest.raises(ValueError):
        delmesher.tetrahedralize(v, np.zeros((1, 3), np.uint32), max_vertices=-5)


# --------------------------------------------------------------------------- #
# the flag sweep (mirrors "delmesher exits 0 for each CLI flag")
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("model", MODELS, ids=lambda p: p.name)
@pytest.mark.parametrize("case_name", list(OPTION_CASES))
def test_option_sweep(model, case_name):
    opts = dict(OPTION_CASES[case_name])
    if CAP is not None and "max_vertices" not in opts:
        opts["max_vertices"] = CAP

    v, f = delmesher.read_off(model)
    result = delmesher.tetrahedralize(v, f, **opts)

    # The Delaunay-refined mesh is always produced and must be well-formed.
    assert_valid_tet_mesh(result.dr_vertices, result.dr_tetrahedra)
    assert_valid_tri_surface(result.dr_interface_vertices,
                             result.dr_interface_faces)
    assert_valid_tri_surface(result.chamfered_vertices, result.chamfered_faces)

    if opts.get("enriched_cdt", True):
        assert result.enriched_cdt_computed
        assert_valid_tet_mesh(result.cdt_vertices, result.cdt_tetrahedra)
        # primary output aliases the enriched CDT
        assert np.array_equal(result.vertices, result.cdt_vertices)
    else:
        assert not result.enriched_cdt_computed
        assert np.array_equal(result.vertices, result.dr_vertices)

    if case_name == "c-min-lfs":
        assert result.min_lfs > 0.0


# --------------------------------------------------------------------------- #
# parity with the command-line tool (the "same result" requirement)
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(CLI is None, reason="delmesher CLI not built; set DELMESHER_BINARY")
@pytest.mark.parametrize("model", MODELS, ids=lambda p: p.name)
def test_matches_cli(model, tmp_path):
    # One CLI run writes both meshes; one binding run produces both.
    rc, proc = run_cli(model, ["-x", "-z"], tmp_path)
    assert rc == 0, proc.stdout + proc.stderr
    cli_dr_v, cli_dr_t = parse_tet_file(tmp_path / "DR_mesh.tet")
    cli_cdt_v, cli_cdt_t = parse_tet_file(tmp_path / "enrichedCDT_mesh.tet")

    v, f = delmesher.read_off(model)
    r = delmesher.tetrahedralize(v, f, max_vertices=CAP, enriched_cdt=True)

    # Delaunay-refined mesh: written with %f (6 decimals).
    assert r.dr_vertices.shape[0] == cli_dr_v.shape[0]
    assert np.array_equal(r.dr_tetrahedra.astype(np.int64), cli_dr_t)
    assert np.allclose(r.dr_vertices, cli_dr_v, rtol=0, atol=2e-6)

    # Enriched CDT: written with operator<< (6 significant digits).
    assert r.cdt_vertices.shape[0] == cli_cdt_v.shape[0]
    assert np.array_equal(r.cdt_tetrahedra.astype(np.int64), cli_cdt_t)
    assert np.allclose(r.cdt_vertices, cli_cdt_v, rtol=1e-4, atol=1e-5)


@pytest.mark.parametrize("model", MODELS, ids=lambda p: p.name)
def test_determinism(model):
    """Repeated in-process calls must reproduce the same mesh bit-for-bit."""
    v, f = delmesher.read_off(model)
    r1 = delmesher.tetrahedralize(v, f, max_vertices=CAP)
    r2 = delmesher.tetrahedralize(v, f, max_vertices=CAP)
    assert np.array_equal(r1.cdt_tetrahedra, r2.cdt_tetrahedra)
    assert np.array_equal(r1.dr_tetrahedra, r2.dr_tetrahedra)
    assert np.array_equal(r1.cdt_vertices, r2.cdt_vertices)
