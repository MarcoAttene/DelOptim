cmake_minimum_required(VERSION 3.10)

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64
# processors, do not support any of the two.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)

#####################

# use gcc and g++ instead of clang
#set(CMAKE_C_COMPILER "gcc-12")
#set(CMAKE_CXX_COMPILER "g++-12")

# Build the LGPL-only variant by default (disables the USE_MAROTS_METHOD code
# path). Use CMake -DLGPL=OFF to opt back into USE_MAROTS_METHOD.
option(LGPL "LGPL" ON)

# Build the command-line executable (delmesher).
option(DELMESHER_BUILD_CLI "Build the delmesher command-line executable" ON)

# Build the Catch2-based test suite (fetches Catch2 on first configure)
option(DELMESHER_BUILD_TESTS "Build the delmesher test suite" ON)

# Build the Python extension module (driven by scikit-build-core / pyproject.toml).
option(DELMESHER_BUILD_PYTHON "Build the delmesher Python bindings" OFF)

# set the project name
project(delmesher)

# Default to a Release (optimized) build when the user hasn't asked for a
# specific configuration. Only applies to single-config generators; multi-config
# generators (Visual Studio, Ninja Multi-Config) select the config at build time.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build (Debug Release RelWithDebInfo MinSizeRel)" FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
		Debug Release RelWithDebInfo MinSizeRel)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# On ARM (e.g. Apple Silicon) the exact predicates rely on SIMDe to emulate the
# x86 AVX2/FMA intrinsics (<x86/avx2.h>, <x86/fma.h>) on top of NEON. Fetched
# once here and added to every target by delmesher_configure_target().
set(DELMESHER_IS_ARM FALSE)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
	set(DELMESHER_IS_ARM TRUE)
	include(FetchContent)
	FetchContent_Declare(simde
		GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
		GIT_TAG v0.8.2
	)
	FetchContent_MakeAvailable(simde)
endif()

# Apply the compile flags, include paths and definitions shared by the CLI
# executable and the Python module to the given target.
function(delmesher_configure_target tgt)
	target_include_directories(${tgt} PRIVATE src include)

	if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# optimize for speed
		target_compile_options(${tgt} PRIVATE /Ot)
	# grant IEEE 754 compliance
		target_compile_options(${tgt} PRIVATE "/fp:strict")
	# use intrinsic functions
		target_compile_options(${tgt} PRIVATE "/Oi")
	# turn off annoying warnings
		target_compile_options(${tgt} PRIVATE "/D _CRT_SECURE_NO_WARNINGS")
	# stop <windows.h> (pulled in by logger.h / getRSS.c) from defining min/max
	# macros, which otherwise clash with std::min/std::max (e.g. in cham.cpp)
		target_compile_options(${tgt} PRIVATE "/D NOMINMAX")
	# set target architecture
		if(ENABLE_AVX2)
			target_compile_options(${tgt} PRIVATE "/arch:AVX2")
		elseif(ENABLE_SSE2)
			target_compile_options(${tgt} PRIVATE "/arch:SSE2")
		endif()
		if(NOT LGPL)
			target_compile_options(${tgt} PRIVATE "/D USE_MAROTS_METHOD")
		endif()
	else()
	# set standard optimization level
		target_compile_options(${tgt} PRIVATE -O2)
	# grant IEEE 754 compliance
		target_compile_options(${tgt} PRIVATE -frounding-math)
	# the exact predicates rely on error-free transformations (e.g. Dekker
	# two-product): forbid the compiler from contracting a*b+c into an FMA, which
	# rounds once instead of twice and silently breaks those identities. This
	# mirrors the /fp:strict used on the MSVC path.
		target_compile_options(${tgt} PRIVATE -ffp-contract=off)
	# set target architecture
		if(DELMESHER_IS_ARM)
			# ARM (e.g. Apple Silicon): the x86 -mavx2/-mfma/-msse2 flags are
			# invalid here. The code emulates the x86 AVX2/FMA intrinsics on top
			# of NEON via SIMDe (fetched above).
		else()
			if(ENABLE_AVX2)
				target_compile_options(${tgt} PRIVATE "-mavx2")
				target_compile_options(${tgt} PRIVATE "-mfma")
			elseif(ENABLE_SSE2)
				target_compile_options(${tgt} PRIVATE "-msse2")
			endif()
		endif()
		if(NOT LGPL)
			target_compile_options(${tgt} PRIVATE "-DUSE_MAROTS_METHOD")
		endif()
	endif()

	if(DELMESHER_IS_ARM)
		target_include_directories(${tgt} PRIVATE ${simde_SOURCE_DIR}/simde)
	endif()
endfunction()

# ---- Command-line executable ------------------------------------------------
if(DELMESHER_BUILD_CLI)
	add_executable(${PROJECT_NAME}
		src/main.cpp
		src/getRSS.c
	)
	delmesher_configure_target(${PROJECT_NAME})

	if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# reserve enough stack size
		target_link_options(${PROJECT_NAME} PRIVATE "/STACK:8421376")
	else()
	# reserve enough stack size (GNU ld only; the macOS linker uses the 8 MB default)
		if(NOT APPLE)
			target_compile_options(${PROJECT_NAME} PRIVATE -Wl,-z,stacksize=8421376)
		endif()
	endif()
endif()

# ---- Python extension module ------------------------------------------------
if(DELMESHER_BUILD_PYTHON)
	find_package(Python 3.8
		COMPONENTS Interpreter Development.Module REQUIRED)

	# Locate nanobind's CMake config via the active interpreter (works whether
	# nanobind is installed in a venv or provided by scikit-build-core). The
	# project targets CMake 3.10, which predates CMP0074, so nanobind_ROOT alone
	# is not used as a search hint -- pass the directory explicitly via PATHS.
	execute_process(
		COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
		OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NANOBIND_CMAKE_DIR)
	find_package(nanobind CONFIG REQUIRED PATHS "${NANOBIND_CMAKE_DIR}")

	nanobind_add_module(_delmesher NB_STATIC
		python/src/delmesher_ext.cpp
		src/getRSS.c
	)
	delmesher_configure_target(_delmesher)

	# Install next to the pure-Python package (delmesher/_delmesher*.so).
	install(TARGETS _delmesher LIBRARY DESTINATION delmesher)
endif()

# ---- Test suite -------------------------------------------------------------
if(DELMESHER_BUILD_CLI AND DELMESHER_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
