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

# Use CMake -DLGPL=ON if you opt for the LGPL license
option(LGPL "LGPL" OFF)

# set the project name
project(delmesher)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(${PROJECT_NAME} 
    src/main.cpp
    src/getRSS.c
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
# optimize for speed
	target_compile_options(${PROJECT_NAME} PUBLIC /Ot)
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
# use intrinsic functions
	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
# turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
	endif()
# reserve enough stack size
	target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")
	if(NOT LGPL)
		target_compile_options(${PROJECT_NAME} PUBLIC "/D USE_MAROTS_METHOD")
	endif()
else()
# set standard optimization level
	target_compile_options(${PROJECT_NAME} PUBLIC -O2)
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
	# set target architecture
	if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
		# ARM (e.g. Apple Silicon): the x86 -mavx2/-mfma/-msse2 flags are
		# invalid here. The code emulates the x86 AVX2/FMA intrinsics on top
		# of NEON via SIMDe, which is fetched below.
	else()
		if(ENABLE_AVX2)
			target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
			target_compile_options(${PROJECT_NAME} PUBLIC "-mfma")
		elseif(ENABLE_SSE2)
			target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
		endif()
	endif()
# reserve enough stack size (GNU ld only; the macOS linker uses the 8 MB default)
	if(NOT APPLE)
		target_compile_options(${PROJECT_NAME} PUBLIC -Wl,-z,stacksize=8421376)
	endif()
	if(NOT LGPL)
		target_compile_options(${PROJECT_NAME} PUBLIC "-DUSE_MAROTS_METHOD")
	endif()
endif()

# On ARM (e.g. Apple Silicon) the exact predicates rely on SIMDe to emulate the
# x86 AVX2/FMA intrinsics (<x86/avx2.h>, <x86/fma.h>) on top of NEON.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
	include(FetchContent)
	FetchContent_Declare(simde
		GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
		GIT_TAG v0.8.2
	)
	FetchContent_MakeAvailable(simde)
	target_include_directories(${PROJECT_NAME} PUBLIC ${simde_SOURCE_DIR}/simde)
endif()

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC
	src
	include
)

