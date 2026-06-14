// Integration tests for the delmesher command-line interface.
//
// A single generic test case runs the compiled `delmesher` binary on every
// input model found in the input_models/ directory, switching on each accepted
// command-line flag one at a time. The binary documents that it "returns 0 when
// the whole execution terminates correctly", so each invocation passes iff it
// exits with code 0.
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
    return cases;
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

std::string quote(const std::string& s) { return "\"" + s + "\""; }

// Run a shell command and return the process exit code (-1 if it never ran).
int run(const std::string& cmd) {
    const int raw = std::system(cmd.c_str());
    if (raw == -1) return -1;
#ifdef _WIN32
    return raw;
#else
    return WIFEXITED(raw) ? WEXITSTATUS(raw) : -1;
#endif
}

}  // namespace

TEST_CASE("delmesher exits 0 for each CLI flag", "[cli][integration]") {
    const std::vector<fs::path> models = discover_input_models();
    REQUIRE_FALSE(models.empty());  // there must be at least one input model

    // Cartesian product: every model crossed with every flag case.
    const auto& model = GENERATE_REF(from_range(models));
    const auto& flags = GENERATE_REF(from_range(flag_cases()));

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
    cmd += " " + quote(model.string());

    INFO("model : " << model.filename().string());
    INFO("flag  : " << flags.name);
    INFO("command: " << cmd);

    CHECK(run(cmd) == 0);
}
