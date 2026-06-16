#!/usr/bin/env python3
"""Minimal post-build smoke test for a freshly built delmesher wheel.

cibuildwheel runs this against every wheel it produces (in an isolated venv with
only the wheel + numpy installed), so it must not depend on the source tree
beyond the input model passed on the command line. Cross-platform on purpose:
invoked as `python smoke_test.py <model.off>` so the same line works in bash and
cmd, unlike a multi-statement `python -c "..."`.
"""
import sys

import numpy as np  # noqa: F401  (ensures the runtime dependency is importable)
import delmesher

model = sys.argv[1] if len(sys.argv) > 1 else "input_models/boeing_part.off"

v, f = delmesher.read_off(model)
result = delmesher.tetrahedralize(v, f, max_vertices=500)

assert result.tetrahedra.ndim == 2, result.tetrahedra.shape
assert result.tetrahedra.shape[1] == 4, result.tetrahedra.shape
assert len(result.tetrahedra) > 0, "no tetrahedra produced"

print("smoke ok", result.vertices.shape, result.tetrahedra.shape)
