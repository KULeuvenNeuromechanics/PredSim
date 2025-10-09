
Get a library with HSL linear solvers for use with ipopt
========================================================

1. Get the source code for the HSL solvers from https://licences.stfc.ac.uk/product/coin-hsl. 
You will receive `coinhsl-2021.05.05.zip` (or maybe a newer version?).

2. Install [MSYS2](https://www.msys2.org/)

3. Run MinGW64 (came with MSYS2), this will open a command window.
    > Use right-click or shift+insert to paste text into this window, ctrl+v does not work.
    > Make sure all file paths use `/`, not `\`.

4. Install packages
- General utility `pacman -S git make`
- GNU compilers (c and fortran) `pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran`
- Dependencies for solvers (lapack for math, metis for matrix ordering) `pacman -S mingw-w64-x86_64-lapack mingw-w64-x86_64-metis`



cd c:/GBW_MyPrograms/coin-or/coinhsl/source

git clone https://github.com/coin-or-tools/ThirdParty-HSL.git

put coinhsl source code in c:/GBW_MyPrograms/coin-or/coinhsl/source/ThirdParty-HSL


