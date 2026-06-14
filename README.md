# Surface chamfering for robust tetrahedral meshing

This code implements our quality meshing algorithm described in "**Surface chamfering for robust tetrahedral meshing**" by L. Diazzi, J. Dai, D. Panozzo and <a href="http://saturno.ge.imati.cnr.it/ima/personal-old/attene/PersonalPage/attene.html">M. Attene</a> (ACM Trans Graphics Vol 45, N. 4, SIGGRAPH 2026).
Please refer to our paper for details here: https://doi.org/10.1145/3811395

<p align="center"><img src="teaser_img.png"></p>

## Usage
Clone this repository, including submodules, with:
```
git clone https://github.com/MarcoAttene/DelOptim
```
You may build the executable as follows:
```
mkdir build
cd build
cmake ..
```
This will produce an appropriate building configuration for your system.
On Windows MSVC, this will produce a delmesher.sln file.
On Linux/MacOS, this will produce a Makefile. 
Use it as usual to compile delmesher.

When compiled, the code generates an executable called ``delmesher``.
Launch it with no command line parameters to have a list of supported options.

We tested our code on MacOS (GCC-10) and Windows (MSVC 2019).
It should work on Linux-GCC and MacOS-Clang too, but we have not tested it on these configurations.

On Apple Silicon (arm64) and other ARM targets the x86 SIMD intrinsics used by
the exact predicates are emulated on top of NEON via
[SIMDe](https://github.com/simd-everywhere/simde), which CMake fetches
automatically; no extra setup is required.

| --- |

## Testing
A small [Catch2](https://github.com/catchorg/Catch2)-based test suite lives in
`tests/`. It is built by default (CMake fetches Catch2 automatically) and runs
``delmesher`` on every model in `input_models/`, switching on each accepted
command-line flag one at a time and checking that the binary exits with code 0.
Build and run it with:
```
cd build
cmake ..
cmake --build .
ctest
```
A full meshing run is expensive, so by default each test run caps Delaunay
refinement (`-m 2000`) to keep the suite to a few minutes; the cap still
exercises every downstream phase and output writer. Use
``-DDELMESHER_TEST_MAX_VERTICES=""`` for unbounded full runs, or
``-DDELMESHER_BUILD_TESTS=OFF`` to skip the tests entirely (e.g. when building
offline).

| --- |

## Citing us
If you use our code in your academic projects, please cite our paper using the following BibTeX entry:
```
@article{delmesher2026,
  title   = {Surface chamfering for robust tetrahedral meshing},
  author  = {Diazzi, Lorenzo and Dai, Jiacheng and Panozzo, Daniele and Attene, Marco},
  journal = {ACM Transactions on Graphics (SIGGRAPH 2026)},
  year    = {2026},
  volume  = {45},
  number  = {4}
}
```

## License
This program is free software; you can redistribute it and/or modify 
it under the terms of the GNU Lesser General Public License (LGPL) as 
published by the Free Software Foundation; either version 3 of the 
License, or (at your option) any later version.                          
                
This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
GNU LGPL (http://www.gnu.org/licenses/lgpl.txt) for more details.
