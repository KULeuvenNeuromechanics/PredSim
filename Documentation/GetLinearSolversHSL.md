
Get a library with HSL linear solvers for use with ipopt
========================================================

Before getting the source code make sure to read the terms of the licence, especially regarding sharing of the software.

1. Get the source code for the HSL solvers from https://licences.stfc.ac.uk/product/coin-hsl. 
You will receive `coinhsl-2021.05.05.zip` (or maybe a newer version?).

2. Install [MSYS2](https://www.msys2.org/)

3. Run MinGW64 (came with MSYS2), this will open a command window.
    > Use right-click or shift+insert to paste text into this window, ctrl+v does not work.
    > Make sure all file paths use `/`, not `\`.

4. Install packages
- In the MinGW64 command window, run:
    - General utility `pacman -S git make pkg-config`
    - GNU compilers (c and fortran) `pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran`
    - Dependencies for solvers (lapack for math, metis for matrix ordering) `pacman -S mingw-w64-x86_64-lapack mingw-w64-x86_64-metis`

5. Organise the source code
- In the MinGW64 command window, run:
    - Go to a base folder `cd c:/documents/coin-or/coinhsl`
    - `mkdir source`
    - `cd source`
    - `git clone https://github.com/coin-or-tools/ThirdParty-HSL.git`
- In file explorer:
    - extract `coinhsl-2021.05.05.zip` into `c:/documents/coin-or/coinhsl/source/ThirdParty-HSL`
    - rename `coinhsl-2021.05.05` to `coinhsl`

6. Compile the code
- In the MinGW64 command window, run:
    - Go to the folder with source code `cd c:/documents/coin-or/coinhsl/source/ThirdParty-HSL`
    - Configure the compiler `./configure --prefix="/c/documents/coin-or/coinhsl/build"--enable-openmp`
    - `make`
    - `make install`
- In `c:/documents/coin-or/coinhsl/build/bin/`, rename `libcoinhsl-2.dll` to `libhsl.dll`. 

7. Use the new solver
- In `/PredSim/main.m`, add:
```
S.solver.ipopt_options.hsllib = 'c:/documents/coin-or/coinhsl/build/bin/libhsl.dll';
S.solver.linear_solver = 'ma97';
```


