// Integration tests for the delmesher command-line interface.
//
// Two complementary test cases drive the compiled `delmesher` binary. The
// binary documents that it "returns 0 when the whole execution terminates
// correctly", so each invocation passes iff it exits with code 0.
//
//   1. The reference model (boeing_part.off) gets the exhaustive treatment:
//      every accepted command-line flag is switched on one at a time, with a
//      vertex cap (-m) so the sweep stays fast while still driving every
//      downstream phase.
//   2. Every other model in input_models/ is run exactly once, with no flags
//      and -- in optimized builds -- no vertex cap, as a full-pipeline check
//      that the mesher handles a variety of real inputs end to end.
//
// The binary path (DELMESHER_BINARY) and the input-model directory
// (DELMESHER_INPUT_DIR) are injected by CMake as compile definitions.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#ifndef _WIN32
#include <sys/wait.h>
#endif

namespace fs = std::filesystem;

namespace {

// The reference model that gets the exhaustive per-flag, capped sweep. Every
// other model is run once and uncapped (see the second test case below).
constexpr const char* PRIMARY_MODEL = "boeing_part.off";

// One CLI invocation: a human-readable name plus the argument list inserted
// between the binary and the input filename.
struct FlagCase {
    std::string name;
    std::vector<std::string> args;
};

// Every command-line flag accepted by delmesher, switched on one at a time.
// Flags that take an integer operand (-e, -m) are given a valid value; the two
// flags documented as "needs [-a]" (-y, -z) are paired with -a, which is the
// documented way to enable them.
const std::vector<FlagCase>& flag_cases() {
    static const std::vector<FlagCase> cases = {
        {"no-flags",         {}},
        {"a-skip-enriched",  {"-a"}},
        {"b-save-cdt-input", {"-b"}},
        {"c-min-lfs",        {"-c"}},
        {"d-sliver-removal", {"-d"}},
        {"e-lfs-bound",      {"-e", "8"}},
        {"m-max-vertices",   {"-m", "2000"}},
        {"h-histograms",     {"-h"}},
        {"l-logging",        {"-l"}},
        {"v-verbose",        {"-v"}},
        {"u-save-chamfer",   {"-u"}},
        {"w-save-dr-skin",   {"-w"}},
        {"x-save-dr-mesh",   {"-x"}},
        {"y-save-out-skin",  {"-a", "-y"}},
        {"z-save-out-mesh",  {"-a", "-z"}},
    };
#if defined(DELMESHER_TEST_SMOKE)
    // Smoke mode: a single capped run is enough to confirm the binary works end
    // to end. Used for the slow unoptimized Debug CI builds; the Release builds
    // run the full sweep above.
    static const std::vector<FlagCase> smoke(cases.begin(), cases.begin() + 1);
    return smoke;
#else
    return cases;
#endif
}

// Discover every .off model in the input directory, so dropping new models into
// input_models/ extends the test automatically (no reconfigure needed).
std::vector<fs::path> discover_input_models() {
    std::vector<fs::path> models;
    for (const auto& entry : fs::directory_iterator(DELMESHER_INPUT_DIR)) {
        if (entry.is_regular_file() && entry.path().extension() == ".off")
            models.push_back(entry.path());
    }
    std::sort(models.begin(), models.end());
    return models;
}

// Every model except the reference model -- the ones run once and uncapped.
std::vector<fs::path> other_models() {
    std::vector<fs::path> models = discover_input_models();
    models.erase(std::remove_if(models.begin(), models.end(),
                                [](const fs::path& p) {
                                    return p.filename() == PRIMARY_MODEL;
                                }),
                 models.end());
    return models;
}

std::string quote(const std::string& s) { return "\"" + s + "\""; }

// Run a shell command and return the process exit code (-1 if it never ran).
int run(const std::string& cmd) {
#ifdef _WIN32
    // cmd.exe strips the outermost pair of quotes from its argument, so wrap the
    // whole command to keep quoted paths (which may contain spaces) intact.
    const std::string line = "\"" + cmd + "\"";
#else
    const std::string& line = cmd;
#endif
    const int raw = std::system(line.c_str());
    if (raw == -1) return -1;
#ifdef _WIN32
    return raw;
#else
    return WIFEXITED(raw) ? WEXITSTATUS(raw) : -1;
#endif
}

// Copy a model into the (throwaway) working directory and return the local
// path. Running against a copy keeps any artifact the binary derives from the
// input name -- e.g. <model>_rebuilt.off from -b -- out of the source tree,
// where, ending in .off, it would be re-discovered as an input.
fs::path stage_model(const fs::path& model) {
    const fs::path local = fs::current_path() / model.filename();
    fs::copy_file(model, local, fs::copy_options::overwrite_existing);
    return local;
}

}  // namespace

TEST_CASE("delmesher exits 0 for each CLI flag (boeing_part, capped)",
          "[cli][integration]") {
    const fs::path primary = fs::path(DELMESHER_INPUT_DIR) / PRIMARY_MODEL;
    REQUIRE(fs::exists(primary));  // the reference model must be present

    const auto& flags = GENERATE_REF(from_range(flag_cases()));
    const fs::path local = stage_model(primary);

    // Build:  "<binary>" [-m <cap>] <flags...> "<model>"
    // The optional cap bounds Delaunay refinement so the suite stays fast; it
    // still drives every downstream phase the flags touch. See CMakeLists.txt.
    // Skip it for the case that exercises -m itself, so the two don't collide.
    const bool sets_max_vertices =
        std::find(flags.args.begin(), flags.args.end(), "-m") != flags.args.end();
    std::string cmd = quote(DELMESHER_BINARY);
    const std::string cap = DELMESHER_TEST_MAX_VERTICES;
    if (!cap.empty() && !sets_max_vertices) cmd += " -m " + cap;
    for (const auto& arg : flags.args) cmd += " " + arg;
    cmd += " " + quote(local.string());

    INFO("model : " << primary.filename().string());
    INFO("flag  : " << flags.name);
    INFO("command: " << cmd);

    CHECK(run(cmd) == 0);
}

TEST_CASE("delmesher exits 0 for each other model (single uncapped run)",
          "[cli][integration]") {
    // An uncapped full run of these models is expensive (minutes each), and
    // unbearably so in an unoptimized Debug build. Smoke mode -- used for the
    // slow Debug CI builds -- therefore skips them entirely: only the capped
    // boeing_part sweep runs there. They run uncapped in the Release builds.
#if defined(DELMESHER_TEST_SMOKE)
    SKIP("smoke mode: only the reference model is exercised");
#else
    const std::vector<fs::path> models = other_models();
    REQUIRE_FALSE(models.empty());  // there must be at least one other model

    const auto& model = GENERATE_REF(from_range(models));
    const fs::path local = stage_model(model);

    // No flags and no vertex cap: a full meshing run of the real input.
    const std::string cmd =
        quote(DELMESHER_BINARY) + " " + quote(local.string());

    INFO("model : " << model.filename().string());
    INFO("command: " << cmd);

    CHECK(run(cmd) == 0);
#endif
}
