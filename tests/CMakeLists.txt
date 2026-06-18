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

# The reference model (boeing_part.off) is run once per accepted flag, capped
# with `-m <N>` so that sweep stays in the minutes range while still exercising
# every flag and every downstream phase (segment/face recovery, enriched CDT,
# and all the output writers run regardless of the cap). Every other model is
# run once, uncapped, as a full-pipeline check. This cap applies only to the
# boeing sweep; set it to an empty string for an unbounded sweep:
#   cmake -DDELMESHER_TEST_MAX_VERTICES="" ..
set(DELMESHER_TEST_MAX_VERTICES "2000" CACHE STRING
	"Vertex cap (-m) for the boeing_part flag sweep; empty for an unbounded sweep")

# Smoke mode reduces the boeing flag sweep to a single capped invocation and
# skips the (uncapped, minutes-each) other-model runs entirely, so the slow
# unoptimized Debug builds in CI stay fast. Release builds run the full per-flag
# sweep and every other model uncapped.
option(DELMESHER_TEST_SMOKE "Reduce the suite to the single capped boeing run (for slow Debug builds)" OFF)

target_compile_definitions(delmesher_tests PRIVATE
	DELMESHER_BINARY="$<TARGET_FILE:delmesher>"
	DELMESHER_INPUT_DIR="${CMAKE_SOURCE_DIR}/input_models"
	DELMESHER_TEST_MAX_VERTICES="${DELMESHER_TEST_MAX_VERTICES}"
)
if(DELMESHER_TEST_SMOKE)
	target_compile_definitions(delmesher_tests PRIVATE DELMESHER_TEST_SMOKE=1)
endif()

# The tests shell out to the binary, so make sure it is built first.
add_dependencies(delmesher_tests delmesher)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(delmesher_tests WORKING_DIRECTORY "${DELMESHER_TEST_WORKDIR}")
