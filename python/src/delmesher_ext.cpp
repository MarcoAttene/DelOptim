// Python bindings for delmesher (nanobind).
//
// This file exposes the exact same meshing pipeline that src/main.cpp drives
// from the command line, but operating on in-memory NumPy arrays instead of
// reading/writing OFF and .tet files. The options map one-to-one to the CLI
// flags (see the mapping in tetrahedralize() below), and the meshes returned
// are extracted from the very same C++ data structures that the CLI writers
// (TetMesh::saveTET, Tetrahedrization::saveTET, saveOFFInterface, ...) read, so
// a binding run and a CLI run on the same input produce identical results.

#ifdef _MSC_VER  // mirror src/main.cpp
#define _HAS_STD_BYTE 0
#endif

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#ifndef _WIN32
#include <sys/wait.h>
#include <unistd.h>
#include <cerrno>
#include <iostream>
#endif

#include "interfaces.h"

namespace nb = nanobind;

namespace {

// ---------------------------------------------------------------------------
// Plain, NumPy-friendly mesh container: flat coordinate / connectivity buffers.
// ---------------------------------------------------------------------------
struct PyMesh {
    std::vector<double> coords;     // x0,y0,z0, x1,y1,z1, ...
    std::vector<uint32_t> elems;    // flat connectivity (cols entries per element)
    uint32_t cols = 4;              // 3 for triangles, 4 for tetrahedra

    size_t num_vertices() const { return coords.size() / 3; }
    size_t num_elements() const { return cols ? elems.size() / cols : 0; }
};

// The full result of one pipeline run. Mirrors what each CLI output flag writes:
//   dr_mesh      <-> -x  DR_mesh.tet              (Delaunay-refined tet mesh)
//   cdt_mesh     <-> -z  enrichedCDT_mesh.tet     (enriched CDT tet mesh)
//   chamfered    <-> -u  chamfered_plc.off        (chamfered surface)
//   dr_interface <-> -w  DR_interface.off         (DR mesh outer interface)
struct PipelineResult {
    PyMesh dr_mesh;        // tetrahedra
    PyMesh cdt_mesh;       // tetrahedra (only if enriched CDT was computed)
    PyMesh chamfered;      // triangles
    PyMesh dr_interface;   // triangles

    bool input_is_manifold = false;
    bool input_has_interior = false;
    bool enriched_cdt_computed = false;
    double min_lfs = -1.0;  // input PLC LFS / bbox diagonal, or -1 if not computed
};

// ---------------------------------------------------------------------------
// Extraction helpers - each mirrors the corresponding C++ file writer so the
// vertex/element ordering (and hence the integer connectivity) is identical.
// ---------------------------------------------------------------------------

// Mirrors Tetrahedrization::saveTET (include/tetmesh.h): vertices in V order,
// tets in T order, indices into V.
void extract_tetrahedrization(const Tetrahedrization& m, PyMesh& out) {
    const TetVertices& V = m.vrts();
    const Tetrahedra& T = m.tets();

    out.cols = 4;
    out.coords.resize(3 * V.size());
    std::unordered_map<const TetVertex*, uint32_t> idx;
    idx.reserve(V.size() * 2);

    uint32_t k = 0;
    for (const TetVertex* v : V) {
        double x, y, z;
        v->getPoint()->getApproxXYZCoordinates(x, y, z, true);
        out.coords[3 * k] = x;
        out.coords[3 * k + 1] = y;
        out.coords[3 * k + 2] = z;
        idx[v] = k++;
    }

    out.elems.clear();
    out.elems.reserve(4 * T.size());
    for (const Tetrahedron* t : T) {
        out.elems.push_back(idx[t->v0()]);
        out.elems.push_back(idx[t->v1()]);
        out.elems.push_back(idx[t->v2()]);
        out.elems.push_back(idx[t->v3()]);
    }
}

// Mirrors TetMesh::saveTET with inner_only == false (src/delaunay.hpp): all
// vertices, then inner (DT_IN) tets followed by outer non-ghost tets.
void extract_tetmesh(const TetMesh& m, PyMesh& out) {
    out.cols = 4;
    const uint32_t nv = m.numVertices();
    out.coords.resize(3 * (size_t)nv);
    for (uint32_t i = 0; i < nv; i++) {
        double x, y, z;
        m.vertices[i]->getApproxXYZCoordinates(x, y, z, true);
        out.coords[3 * i] = x;
        out.coords[3 * i + 1] = y;
        out.coords[3 * i + 2] = z;
    }

    const uint32_t nt = m.numTets();
    out.elems.clear();
    auto push_tet = [&](uint32_t i) {
        out.elems.push_back(m.tet_node[i * 4]);
        out.elems.push_back(m.tet_node[i * 4 + 1]);
        out.elems.push_back(m.tet_node[i * 4 + 2]);
        out.elems.push_back(m.tet_node[i * 4 + 3]);
    };
    for (uint32_t i = 0; i < nt; i++)
        if (m.mark_tetrahedra[i] == DT_IN) push_tet(i);
    for (uint32_t i = 0; i < nt; i++)
        if (!m.isGhost(i) && m.mark_tetrahedra[i] != DT_IN) push_tet(i);
}

// Chamfered surface: the chamfering interface already exposes the vertices and
// triangulated faces it feeds into Delaunay refinement (the -u output).
void extract_chamfered(chamfering_interface& cham, PyMesh& out) {
    out.cols = 3;
    const auto& verts = cham.out_vertices;
    out.coords.resize(3 * verts.size());
    for (size_t i = 0; i < verts.size(); i++) {
        double x, y, z;
        verts[i]->getApproxXYZCoordinates(x, y, z, true);
        out.coords[3 * i] = x;
        out.coords[3 * i + 1] = y;
        out.coords[3 * i + 2] = z;
    }
    out.elems.clear();
    for (const std::vector<uint32_t>& f : cham.out_faces)
        if (f.size() == 3)
            out.elems.insert(out.elems.end(), {f[0], f[1], f[2]});
}

// Mirrors Tetrahedrization::saveOFFInterface (include/tetmesh.h): the interface
// triangles between internal and external tets, compacted to used vertices,
// oriented outward (the -w output).
void extract_dr_interface(const Tetrahedrization& m, PyMesh& out) {
    out.cols = 3;
    std::unordered_map<const TetVertex*, uint32_t> idx;
    std::vector<const TetVertex*> used;

    auto vid = [&](const TetVertex* v) -> uint32_t {
        auto it = idx.find(v);
        if (it != idx.end()) return it->second;
        uint32_t n = (uint32_t)used.size();
        idx[v] = n;
        used.push_back(v);
        return n;
    };

    out.elems.clear();
    for (const TetFace* f : m.faces()) {
        if (!f->isInterface()) continue;
        uint32_t a = vid(f->v0());
        uint32_t b = vid(f->v1());
        uint32_t c = vid(f->v2());
        if (f->t1()->is_internal)
            out.elems.insert(out.elems.end(), {a, b, c});
        else
            out.elems.insert(out.elems.end(), {c, b, a});
    }

    out.coords.resize(3 * used.size());
    for (size_t i = 0; i < used.size(); i++) {
        double x, y, z;
        used[i]->getPoint()->getApproxXYZCoordinates(x, y, z, true);
        out.coords[3 * i] = x;
        out.coords[3 * i + 1] = y;
        out.coords[3 * i + 2] = z;
    }
}

// ---------------------------------------------------------------------------
// The pipeline. A faithful, file-free transcription of src/main.cpp's body
// (non-tetgen path). Option arguments correspond to the CLI flags one-to-one.
// ---------------------------------------------------------------------------
PipelineResult run_pipeline(std::vector<double> coords, uint32_t n_pts,
                            std::vector<uint32_t> tris, uint32_t n_tris,
                            bool enriched_cdt,    // CLI: NOT -a
                            bool comp_inLFS,      // CLI: -c
                            bool remove_slivers,  // CLI: -d
                            uint32_t min_dist_exp,// CLI: -e   (UINT32_MAX = unset)
                            uint32_t max_vrts,    // CLI: -m   (UINT32_MAX = unset)
                            bool verbose) {       // CLI: -v
    initFPU();
    // A fresh CLI process always starts the segment-recovery PRNG from h == 1;
    // reset it here so every binding call reproduces a fresh run bit-for-bit,
    // independent of how many times the pipeline has already run in-process.
    reset_myrand();

    PipelineResult R;

    // ---- Load PLC (main.cpp: initFromFile + setBoundingBox) ---------------
    inputPLC plc;
    plc.initFromVectors(coords.data(), n_pts, tris.data(), n_tris, verbose);
    plc.setBoundingBox(1.0);  // bbox_factor = 1.0

    // ---- Internal parameters (identical defaults to main.cpp) -------------
    const double bbox_factor = 1.0;  (void)bbox_factor;
    double epsilon = DBL_MAX;
    const bool cham_simpl = true;
    const bool cham_safe = false;
    const double optim_ratio = 2.0;
    double toll = DBL_MAX;

    // ---- Optional Steiner CDT for LFS / small-feature handling ------------
    cdt_interface input_cdt;
    bool empty_plc = false;
    if (comp_inLFS || min_dist_exp != UINT32_MAX)
        input_cdt.createSteinerCDT(plc, comp_inLFS, min_dist_exp != UINT32_MAX);

    if (min_dist_exp != UINT32_MAX) {
        toll = std::pow(10.0, -1.0 * (double)min_dist_exp);
        input_cdt.mark_small_tris(toll, plc);
        empty_plc = (plc.count_ignored_tris() == plc.numTriangles());
    }

    // ---- Chamfering -------------------------------------------------------
    chamfering_interface cham(cham_safe, cham_simpl, /*save_output=*/false,
                              empty_plc);
    cham.perform_chamfering(plc, epsilon, toll, verbose);

    // ---- Delaunay refinement ---------------------------------------------
    delRef_interface DR(plc, cham, input_cdt, verbose, /*log=*/false);
    DR.init(comp_inLFS);
    DR.recover_segments();
    DR.recover_faces();
    DR.optimize_tetrahedra(optim_ratio, remove_slivers, max_vrts);
    DR.intExt_classification();
    DR.get_mesh_statistics(/*log=*/false, /*histo=*/false, /*DR_skin=*/false,
                           /*DR_outmesh=*/false);

    R.input_is_manifold = cham.input_is_manifold;
    R.input_has_interior = cham.input_has_interior;
    if (comp_inLFS && input_cdt.is_defined())
        R.min_lfs = input_cdt.inputPLC_LFS() / plc.get_BBox_diag();

    extract_tetrahedrization(DR.mesh, R.dr_mesh);
    extract_chamfered(cham, R.chamfered);
    extract_dr_interface(DR.mesh, R.dr_interface);

    // ---- Enriched CDT (inlined from get_enriched_CDT, file writes removed) -
    if (enriched_cdt) {
        recyled_surface_mesh cdt_input(DR.mesh, plc, verbose);
        cdt_input.set_ref_vrts(cham.out_conn_vertices);
        cdt_input.rebuild_surface(toll);

        cdt_interface e_cdt;
        e_cdt.createSteinerCDT(cdt_input.cdt_vrts, cdt_input.cdt_tris, false);
        extract_tetmesh(*e_cdt.mesh, R.cdt_mesh);
        R.enriched_cdt_computed = true;
    }

    return R;
}

// ---------------------------------------------------------------------------
// Process isolation.
//
// The underlying C++ pipeline is written as a one-shot program: several data
// structures share ownership of the same geometric points, and on teardown the
// destructors would double-free them. That is harmless for the CLI (it just
// exits), but in a long-lived Python process it corrupts the heap and breaks
// the next call. The pipeline also calls exit() on invalid input, which would
// take the interpreter down with it.
//
// To get exactly the CLI's behaviour (a clean, independent run every time) we
// run the pipeline in a forked child and ship the resulting meshes back over a
// pipe. The child _exit()s before any destructor runs, so there is no
// double-free and no cross-call state; a child that exit()s or crashes on bad
// input becomes a Python exception in the parent. On Windows (no fork) we fall
// back to an in-process call, which is reliable for one call per process.
// ---------------------------------------------------------------------------

namespace {

void put_bytes(std::vector<char>& b, const void* p, size_t n) {
    const char* c = (const char*)p;
    b.insert(b.end(), c, c + n);
}
template <typename T>
void put_pod(std::vector<char>& b, T v) { put_bytes(b, &v, sizeof(T)); }

void put_mesh(std::vector<char>& b, const PyMesh& m) {
    put_pod<uint32_t>(b, m.cols);
    put_pod<uint64_t>(b, m.coords.size());
    put_bytes(b, m.coords.data(), m.coords.size() * sizeof(double));
    put_pod<uint64_t>(b, m.elems.size());
    put_bytes(b, m.elems.data(), m.elems.size() * sizeof(uint32_t));
}

std::vector<char> serialize(const PipelineResult& R) {
    std::vector<char> b;
    put_pod<uint8_t>(b, R.input_is_manifold);
    put_pod<uint8_t>(b, R.input_has_interior);
    put_pod<uint8_t>(b, R.enriched_cdt_computed);
    put_pod<double>(b, R.min_lfs);
    put_mesh(b, R.dr_mesh);
    put_mesh(b, R.cdt_mesh);
    put_mesh(b, R.chamfered);
    put_mesh(b, R.dr_interface);
    return b;
}

struct Cursor {
    const char* p;
    const char* end;
    void take(void* dst, size_t n) {
        if (p + n > end) throw std::runtime_error("truncated pipeline result");
        std::memcpy(dst, p, n);
        p += n;
    }
    template <typename T> T pod() { T v; take(&v, sizeof(T)); return v; }
    void mesh(PyMesh& m) {
        m.cols = pod<uint32_t>();
        m.coords.resize(pod<uint64_t>());
        take(m.coords.data(), m.coords.size() * sizeof(double));
        m.elems.resize(pod<uint64_t>());
        take(m.elems.data(), m.elems.size() * sizeof(uint32_t));
    }
};

PipelineResult deserialize(const std::vector<char>& b) {
    Cursor c{b.data(), b.data() + b.size()};
    PipelineResult R;
    R.input_is_manifold = c.pod<uint8_t>() != 0;
    R.input_has_interior = c.pod<uint8_t>() != 0;
    R.enriched_cdt_computed = c.pod<uint8_t>() != 0;
    R.min_lfs = c.pod<double>();
    c.mesh(R.dr_mesh);
    c.mesh(R.cdt_mesh);
    c.mesh(R.chamfered);
    c.mesh(R.dr_interface);
    return R;
}

#ifndef _WIN32
bool write_all(int fd, const void* buf, size_t n) {
    const char* p = (const char*)buf;
    while (n) {
        ssize_t w = ::write(fd, p, n);
        if (w < 0) { if (errno == EINTR) continue; return false; }
        p += w;
        n -= (size_t)w;
    }
    return true;
}
bool read_all(int fd, void* buf, size_t n) {
    char* p = (char*)buf;
    while (n) {
        ssize_t r = ::read(fd, p, n);
        if (r < 0) { if (errno == EINTR) continue; return false; }
        if (r == 0) return false;  // EOF before all bytes (child died early)
        p += r;
        n -= (size_t)r;
    }
    return true;
}
#endif

// Forwarding wrapper: argument list mirrors the CLI options exactly.
template <typename... Args>
PipelineResult run_pipeline_isolated(Args&&... args) {
#ifdef _WIN32
    return run_pipeline(std::forward<Args>(args)...);
#else
    int fds[2];
    if (::pipe(fds) != 0)
        throw std::runtime_error("delmesher: pipe() failed");

    pid_t pid = ::fork();
    if (pid < 0) {
        ::close(fds[0]);
        ::close(fds[1]);
        throw std::runtime_error("delmesher: fork() failed");
    }

    if (pid == 0) {
        // Child: compute, ship the result, and _exit() before any destructor
        // (or Python finalizer) runs. Never returns to Python.
        ::close(fds[0]);
        std::vector<char> buf;
        try {
            PipelineResult R = run_pipeline(std::forward<Args>(args)...);
            std::cout.flush();
            buf = serialize(R);
        } catch (...) {
            ::_exit(2);  // unexpected C++ exception in the pipeline
        }
        uint64_t len = buf.size();
        bool ok = write_all(fds[1], &len, sizeof(len)) &&
                  write_all(fds[1], buf.data(), buf.size());
        ::close(fds[1]);
        ::_exit(ok ? 0 : 3);
    }

    // Parent: drain the pipe (keeps the child unblocked), then reap it.
    ::close(fds[1]);
    uint64_t len = 0;
    bool ok = read_all(fds[0], &len, sizeof(len));
    std::vector<char> buf;
    if (ok) {
        buf.resize(len);
        ok = read_all(fds[0], buf.data(), len);
    }
    ::close(fds[0]);

    int status = 0;
    while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) {}

    if (!ok || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        if (WIFSIGNALED(status))
            throw std::runtime_error(
                "delmesher: meshing process crashed (signal " +
                std::to_string(WTERMSIG(status)) +
                "); the input is likely invalid (non-manifold, self-intersecting "
                "or degenerate).");
        throw std::runtime_error(
            "delmesher: meshing failed; the input is likely invalid "
            "(non-manifold, self-intersecting or degenerate).");
    }
    return deserialize(buf);
#endif
}

}  // namespace

// ---------------------------------------------------------------------------
// nanobind glue
// ---------------------------------------------------------------------------

using F64Array = nb::ndarray<const double, nb::ndim<2>, nb::c_contig,
                             nb::device::cpu>;
using U32Array = nb::ndarray<const uint32_t, nb::ndim<2>, nb::c_contig,
                             nb::device::cpu>;

// Build a NumPy array (rows x cols) that owns a fresh copy of `src`.
template <typename T>
nb::ndarray<nb::numpy, T, nb::ndim<2>> make_2d(const std::vector<T>& src,
                                               size_t cols) {
    size_t rows = cols ? src.size() / cols : 0;
    T* data = new T[src.size() ? src.size() : 1];
    if (!src.empty()) std::memcpy(data, src.data(), src.size() * sizeof(T));
    nb::capsule owner(data, [](void* p) noexcept { delete[] (T*)p; });
    return nb::ndarray<nb::numpy, T, nb::ndim<2>>(data, {rows, cols}, owner);
}

nb::object coords_array(const PyMesh& m) { return nb::cast(make_2d(m.coords, 3)); }
nb::object elems_array(const PyMesh& m) {
    return nb::cast(make_2d(m.elems, m.cols));
}

}  // namespace

NB_MODULE(_delmesher, m) {
    m.doc() = "Robust tetrahedral meshing (delmesher) - native bindings";

    nb::class_<PipelineResult>(m, "Result",
        "Output of delmesher.tetrahedralize(). Meshes are NumPy arrays: "
        "vertex coordinates are (n, 3) float64; element connectivity is "
        "(m, k) uint32 (k=4 tetrahedra, k=3 triangles).")
        .def_prop_ro("dr_vertices", [](const PipelineResult& r) { return coords_array(r.dr_mesh); },
                     "Vertices of the Delaunay-refined tet mesh (CLI -x).")
        .def_prop_ro("dr_tetrahedra", [](const PipelineResult& r) { return elems_array(r.dr_mesh); },
                     "Tetrahedra of the Delaunay-refined mesh (CLI -x).")
        .def_prop_ro("cdt_vertices", [](const PipelineResult& r) { return coords_array(r.cdt_mesh); },
                     "Vertices of the enriched-CDT tet mesh (CLI -z); empty if not computed.")
        .def_prop_ro("cdt_tetrahedra", [](const PipelineResult& r) { return elems_array(r.cdt_mesh); },
                     "Tetrahedra of the enriched-CDT mesh (CLI -z); empty if not computed.")
        .def_prop_ro("chamfered_vertices", [](const PipelineResult& r) { return coords_array(r.chamfered); },
                     "Vertices of the chamfered input surface (CLI -u).")
        .def_prop_ro("chamfered_faces", [](const PipelineResult& r) { return elems_array(r.chamfered); },
                     "Triangles of the chamfered input surface (CLI -u).")
        .def_prop_ro("dr_interface_vertices", [](const PipelineResult& r) { return coords_array(r.dr_interface); },
                     "Vertices of the DR mesh outer interface surface (CLI -w).")
        .def_prop_ro("dr_interface_faces", [](const PipelineResult& r) { return elems_array(r.dr_interface); },
                     "Triangles of the DR mesh outer interface surface (CLI -w).")
        .def_ro("input_is_manifold", &PipelineResult::input_is_manifold)
        .def_ro("input_has_interior", &PipelineResult::input_has_interior)
        .def_ro("enriched_cdt_computed", &PipelineResult::enriched_cdt_computed)
        .def_ro("min_lfs", &PipelineResult::min_lfs,
                "Input PLC local feature size / bbox diagonal (CLI -c), or -1 if not computed.");

    m.def(
        "_tetrahedralize",
        [](F64Array verts, U32Array tris, bool enriched_cdt, bool sliver_removal,
           long lfs_exponent, long max_vertices, bool compute_lfs,
           bool verbose) {
            const uint32_t n_pts = (uint32_t)verts.shape(0);
            const uint32_t n_tris = (uint32_t)tris.shape(0);

            std::vector<double> coords(verts.data(), verts.data() + 3 * (size_t)n_pts);
            std::vector<uint32_t> t(tris.data(), tris.data() + 3 * (size_t)n_tris);

            const uint32_t min_dist_exp =
                lfs_exponent < 0 ? UINT32_MAX : (uint32_t)lfs_exponent;
            const uint32_t max_vrts =
                max_vertices < 0 ? UINT32_MAX : (uint32_t)max_vertices;

            PipelineResult r;
            {
                // The pipeline is single-threaded but heavy; release the GIL
                // while the child computes and we wait on it.
                nb::gil_scoped_release release;
                r = run_pipeline_isolated(std::move(coords), n_pts, std::move(t),
                                          n_tris, enriched_cdt, compute_lfs,
                                          sliver_removal, min_dist_exp, max_vrts,
                                          verbose);
            }
            return r;
        },
        nb::arg("vertices"), nb::arg("triangles"), nb::arg("enriched_cdt") = true,
        nb::arg("sliver_removal") = false, nb::arg("lfs_exponent") = -1,
        nb::arg("max_vertices") = -1, nb::arg("compute_lfs") = false,
        nb::arg("verbose") = false,
        "Low-level pipeline entry point; prefer delmesher.tetrahedralize().");
}
