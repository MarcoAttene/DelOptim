# Surface chamfering for robust tetrahedral meshing

This code implements our quality meshing algorithm described in "**Surface chamfering for robust tetrahedral meshing**" by L. Diazzi, J. Dai, D. Panozzo and M. Attene (ACM Trans Graphics Vol 45, N. 4, SIGGRAPH 2026).
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
