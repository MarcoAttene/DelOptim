# Unit/integration tests for delmesher, built on Catch2 (fetched on demand).

include(FetchContent)
FetchContent_Declare(Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG v3.7.1
)
FetchContent_MakeAvailable(Catch2)

add_executable(delmesher_tests test_delmesher_cli.cpp)
target_compile_features(delmesher_tests PRIVATE cxx_std_20)
target_link_libraries(delmesher_tests PRIVATE Catch2::Catch2WithMain)

# Run each test in a dedicated working directory so the meshes the binary writes
# (out_mesh.tet, delOpt_log.csv, ...) do not pollute the build tree.
set(DELMESHER_TEST_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/workdir")
file(MAKE_DIRECTORY "${DELMESHER_TEST_WORKDIR}")

# A full meshing run of the shipped model is expensive (tens of minutes), and
# the test runs the binary once per accepted flag. By default we cap Delaunay
# refinement with `-m <N>` so the suite stays in the minutes range while still
# exercising every flag and every downstream phase (segment/face recovery,
# enriched CDT, and all the output writers run regardless of the cap). Set to an
# empty string for unbounded full runs:
#   cmake -DDELMESHER_TEST_MAX_VERTICES="" ..
set(DELMESHER_TEST_MAX_VERTICES "2000" CACHE STRING
	"Vertex cap (-m) applied to every CLI test run; empty for unbounded full runs")

target_compile_definitions(delmesher_tests PRIVATE
	DELMESHER_BINARY="$<TARGET_FILE:delmesher>"
	DELMESHER_INPUT_DIR="${CMAKE_SOURCE_DIR}/input_models"
	DELMESHER_TEST_MAX_VERTICES="${DELMESHER_TEST_MAX_VERTICES}"
)

# The tests shell out to the binary, so make sure it is built first.
add_dependencies(delmesher_tests delmesher)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(delmesher_tests WORKING_DIRECTORY "${DELMESHER_TEST_WORKDIR}")
