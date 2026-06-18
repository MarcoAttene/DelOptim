# Surface Chamfering for Robust Tetrahedral Meshing

[![CI](https://github.com/MarcoAttene/DelOptim/actions/workflows/ci.yml/badge.svg)](https://github.com/MarcoAttene/DelOptim/actions/workflows/ci.yml)

`delmesher` turns a triangulated input surface into a high-quality tetrahedral
mesh of its interior, and is robust to the inputs that defeat classical
meshers — in particular faces that meet at sharp/acute angles and faces with
arbitrarily low-quality (even near-degenerate) triangulations.

<p align="center"><img src="teaser_img.png" width="100%"></p>

> An input polyhedral surface whose faces meet at acute angles and carry a
> low-quality triangulation (left); the chamfered surface with the acute angles
> removed (mid-left); the volume meshed into tetrahedra with all face angles
> above 14° (mid-right); and the final mesh, made exactly conformal to the input
> by allowing a few tetrahedra with smaller angles (right).

This is the reference implementation of the paper:

> **Surface chamfering for robust tetrahedral meshing.**
> Lorenzo Diazzi, Jiacheng Dai, Daniele Panozzo, and Marco Attene.
> *ACM Transactions on Graphics* 45, 4 (SIGGRAPH 2026), Article 146.
> [doi:10.1145/3811395](https://doi.org/10.1145/3811395) ·
> [paper (PDF)](https://cims.nyu.edu/gcl/papers/2026-chamfering.pdf)

## Overview

Robust tetrahedral meshing usually forces a trade-off. *Boundary-conforming*
methods reproduce the input surface exactly but tend to choke on sharp features
and bad triangles; *boundary-approximating* methods are robust but no longer
match the input. This work takes a hybrid route that keeps robustness while
modifying the input only in a tightly controlled way. The pipeline has three
stages:

1. **Chamfering** — the input surface is locally and minimally modified to
   remove the acute angles between faces that are the main obstacle to robust,
   well-shaped meshing.
2. **Delaunay refinement** — the volume bounded by the chamfered surface is
   tetrahedralized into well-shaped elements, with a guaranteed lower bound on
   the face angles (> 14°).
3. **Enriched CDT** — exact conformance to the *original* input surface is
   restored, introducing a few lower-quality tetrahedra only where strictly
   necessary.

Every geometric decision is taken with exact predicates, so the pipeline is
numerically robust regardless of how degenerate the input is.

## Building

### Requirements

- [CMake](https://cmake.org/) ≥ 3.10
- A C++20 compiler — tested with GCC, Clang/AppleClang and MSVC 2019+
- Network access on the first configure: dependencies are fetched automatically
  ([SIMDe](https://github.com/simd-everywhere/simde) on ARM,
  [Catch2](https://github.com/catchorg/Catch2) for the tests)

### Build

```sh
git clone https://github.com/MarcoAttene/DelOptim
cd DelOptim
mkdir build && cd build
cmake ..
cmake --build . --config Release
```

The build defaults to an optimized `Release` configuration, so no
`-DCMAKE_BUILD_TYPE` is needed (the `--config Release` above only matters for
multi-config generators such as Visual Studio). This produces the `delmesher`
executable (`build/delmesher` on Linux/macOS, `build/Release/delmesher.exe` on
Windows).

> **Avoid Debug for real runs.** A `Debug` build (`-DCMAKE_BUILD_TYPE=Debug`)
> leaves `NDEBUG` unset, so the internal assertions stay enabled and the binary
> is several times slower.

### CMake options

| Option | Default | Description |
| --- | --- | --- |
| `CMAKE_BUILD_TYPE` | `Release` | Defaults to an optimized `Release` build; set `Debug` for assertions (much slower). |
| `DELMESHER_BUILD_TESTS` | `ON` | Build the Catch2 test suite (fetches Catch2). Set `OFF` to build offline / skip tests. |
| `LGPL` | `ON` | Build the variant covered solely by the LGPL (default). Set `OFF` to enable the `USE_MAROTS_METHOD` code path. |
| `DELMESHER_TEST_MAX_VERTICES` | `2000` | Refinement cap (`-m`) for the `boeing_part.off` flag sweep (and all models in smoke mode); empty string for an unbounded sweep. |

On x86, AVX2 is enabled by default; for older CPUs flip `ENABLE_AVX2` /
`ENABLE_SSE2` at the top of `CMakeLists.txt`. On ARM (e.g. Apple Silicon) those
x86 SIMD intrinsics are transparently emulated on NEON through SIMDe.

## Running

Run with no arguments to print the full list of options:

```sh
./build/delmesher
```

Process the bundled example model and write the angle-bounded tetrahedral mesh:

```sh
# chamfer -> Delaunay refinement; write the well-shaped output (out_mesh.tet)
./build/delmesher -a -z -v input_models/boeing_part.off
```

A full run refines until convergence and can take a few minutes; pass `-m N` to
stop after `N` inserted vertices (the test suite uses `-m 2000`). The program
prints `Execution correctly COMPLETED.` and returns `0` on success.

### Command-line options

| Flag | Argument | Effect |
| --- | --- | --- |
| *(none)* | | run the full pipeline: chamfer → Delaunay refinement → enriched CDT |
| `-a` | | skip the final *enriched CDT* phase; the Delaunay-refined output is angle-bounded (min face angle > 14°) but its constrained faces are not conformal to the input |
| `-c` | | compute the minimum local feature size (LFS) of the input triangulation |
| `-d` | | enable sliver removal during Delaunay refinement |
| `-e` | *int* | ignore input constraints whose LFS is below 10⁻ᵉ × the bounding-box diagonal (larger `e` ⇒ output closer to a plain CDT of the input) |
| `-m` | *int* | stop Delaunay refinement after this many inserted vertices |
| `-h` | | print angle histograms |
| `-l` | | logging mode: append a row to `delOpt_log.csv` |
| `-v` | | verbose output |

Output flags (each writes into the current directory):

| Flag | Needs | Output file | Type |
| --- | --- | --- | --- |
| `-b` | | `<input>_rebuilt.off` | surface |
| `-u` | | `chamfered_plc.off` | surface |
| `-w` | | `DR_interface.off` | surface |
| `-x` | | `DR_mesh.tet` | volume |
| `-y` | `-a` | `constrainedFaces.off` | surface |
| `-z` | `-a` | `out_mesh.tet` | volume |

## Python bindings

`delmesher` ships [nanobind](https://github.com/wjakob/nanobind) Python bindings
that drive the exact same pipeline as the CLI, but take and return
[NumPy](https://numpy.org/) arrays instead of OFF / `.tet` files.

### Install

```sh
pip install delmesher
```

Prebuilt binary wheels are published to
[PyPI](https://pypi.org/project/delmesher/) for Linux, macOS (Apple Silicon) and
Windows on CPython 3.9–3.13, so no compiler is required. To build from a source
checkout instead, run `pip install .`; this needs CMake and a C++20 compiler
(nanobind and, on ARM, SIMDe are fetched automatically).

### Quick start

```python
import delmesher

# vertices: (n, 3) float64;  triangles: (m, 3) integer indices
vertices, triangles = delmesher.read_off("input_models/boeing_part.off")

result = delmesher.tetrahedralize(vertices, triangles, max_vertices=2000)

result.vertices    # (V, 3) float64  -- final output mesh
result.tetrahedra  # (T, 4) uint32   -- indices into result.vertices
delmesher.write_tet("out_mesh.tet", result.vertices, result.tetrahedra)
```

The notebook [`examples/delmesher_demo.ipynb`](examples/delmesher_demo.ipynb)
walks through the whole pipeline and visualizes every output with
[PyVista](https://pyvista.org/) (`pip install ".[notebook]"`).

### Options

Every keyword argument of `tetrahedralize` maps to a CLI flag:

| Keyword | Default | CLI flag | Effect |
| --- | --- | --- | --- |
| `enriched_cdt` | `True` | `-a` ⇒ `False` | run the final enriched-CDT phase (exact conformance); `False` stops after Delaunay refinement |
| `sliver_removal` | `False` | `-d` | enable sliver removal during refinement |
| `lfs_exponent` | `None` | `-e` | ignore constraints with LFS below `10⁻ᵉ × bbox diagonal` |
| `max_vertices` | `None` | `-m` | stop refinement after this many inserted vertices |
| `compute_lfs` | `False` | `-c` | compute the input's minimum LFS (see `result.min_lfs`) |
| `verbose` | `False` | `-v` | print progress to stdout |

The returned `Result` exposes every mesh the CLI can write, as NumPy arrays:
`vertices`/`tetrahedra` (the primary output — enriched CDT, or the
Delaunay-refined mesh when `enriched_cdt=False`), `dr_vertices`/`dr_tetrahedra`
(`-x`), `cdt_vertices`/`cdt_tetrahedra` (`-z`), `chamfered_vertices`/
`chamfered_faces` (`-u`) and `dr_interface_vertices`/`dr_interface_faces` (`-w`),
plus the scalars `input_is_manifold`, `input_has_interior` and `min_lfs`.

Each call runs the pipeline in an isolated child process, so results are
reproducible regardless of call history and invalid input raises a Python
exception instead of terminating the interpreter.

### Python tests

`tests/python/` mirrors the C++ suite: it sweeps every option on the reference
model `boeing_part.off` (capped), runs every other model in `input_models/` once
uncapped, and additionally checks that the bindings reproduce the CLI's `.tet`
output bit-for-bit. The helper script builds everything into a virtual
environment under `build/venv` and runs the suite:

```sh
scripts/build_python_venv.sh
```

## Testing

The [Catch2](https://github.com/catchorg/Catch2) suite in `tests/` checks that
`delmesher` exits with code `0`, in two complementary ways:

- **Reference model** (`boeing_part.off`): every accepted command-line flag is
  switched on one at a time, with refinement capped
  (`DELMESHER_TEST_MAX_VERTICES`, default `2000`) so the per-flag sweep stays
  fast while still exercising every downstream phase.
- **Every other model** in `input_models/`: a single run, no flags and no
  vertex cap — a full-pipeline check that the mesher handles a variety of real
  inputs end to end.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Configure with `-DDELMESHER_TEST_MAX_VERTICES=""` for an unbounded boeing sweep,
`-DDELMESHER_TEST_SMOKE=ON` to reduce the suite to the single capped boeing run
(used for the slow unoptimized Debug builds, which skip the uncapped other-model
runs entirely), or `-DDELMESHER_BUILD_TESTS=OFF` to skip the tests. Continuous
integration builds and runs the suite on Linux (x86-64), macOS (Apple Silicon)
and Windows (x86-64) in both Debug and Release on every push; the Debug builds
run only the capped boeing sweep.

### Test runtimes

Wall-clock for one uncapped full run of each non-reference model in an optimized
Release build (single core, Apple Silicon M-series):

| model | input verts | Release (uncapped) |
| --- | ---: | ---: |
| `48416.off` | 268 | 0.3s |
| `423818.off` | 125 | 0.6s |
| `343542.off` | 433 | 1.7s |
| `49529.off` | 840 | 2.5s |
| `85110.off` | 1130 | 3.8s |
| `56099.off` | 776 | 4.8s |
| `227616.off` | 1498 | 7.4s |
| `1558739.off` | 1026 | 8.0s |
| `39641.off` | 540 | 9.1s |
| `350273.off` | 2080 | 12.7s |
| `110907.off` | 1702 | 14.3s |
| `472100.off` | 152 | 15.3s |
| `80363.off` | 827 | 17.5s |
| `67924.off` | 13698 | 19.7s |
| `1130082.off` | 2947 | 21.0s |
| `72669.off` | 1131 | 27.9s |
| `62572.off` | 1694 | 28.2s |
| `58007.off` | 1136 | 30.8s |
| `61242.off` | 3718 | 42.5s |
| `1207651.off` | 360 | 152.8s |

Runtime is driven by the refinement geometry, not the input size. These uncapped
runs are unbearably slow in an unoptimized Debug build (minutes each, several
exceeding five), which is why the Debug builds skip them and run only the capped
`boeing_part.off` sweep — that stays well under a minute per flag even in Debug.

## Citation

If you use this code in your research, please cite:

```bibtex
@article{diazzi2026chamfering,
  author    = {Diazzi, Lorenzo and Dai, Jiacheng and Panozzo, Daniele and Attene, Marco},
  title     = {Surface Chamfering for Robust Tetrahedral Meshing},
  journal   = {ACM Transactions on Graphics},
  year      = {2026},
  volume    = {45},
  number    = {4},
  articleno = {146},
  month     = jul,
  publisher = {Association for Computing Machinery},
  doi       = {10.1145/3811395},
  url       = {https://doi.org/10.1145/3811395}
}
```

## License

This program is free software, released under the **GNU Lesser General Public
License (LGPL), version 3 or (at your option) any later version**. It is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See [`lgpl.txt`](lgpl.txt) and <https://www.gnu.org/licenses/lgpl.txt>
for the full terms.

Copyright © 2019–2026 IMATI-GE / CNR and the authors. By default the build is
covered solely by the LGPL; configure with `-DLGPL=OFF` to enable the
`USE_MAROTS_METHOD` code path.
