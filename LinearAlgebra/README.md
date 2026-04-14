During an attempt to port PredSim to Linux (in the context of trying to run it
on the VSC supercomputer), dependencies such as CasADi, OpenSim, IPopt were
compiled using OpenBLAS and ScaLAPACK. It turned out that this clashes with
the way MATLAB handles linear algebra; MATLAB ships with Intel MKL libraries
which contain BLAS and LAPACK wrappers, so in principle they should be
compatible with OpenBLAS/ScaLAPACK. In practice, that is not the case; a
hypothesis is that this happens because MATLAB has its own additional layer
in between MATLAB code and BLAS/LAPACK calls, but it is hard to debug a
closed-source code. The solution to get the dependencies such as
CasADi/OpenSim/IPopt to work from within MATLAB, was to set the following
variables:

```bash
export BLAS_VERSION=${EBROOTOPENBLAS}/lib64/libopenblas.so
export LAPACK_VERSION=${EBROOTOPENBLAS}/lib64/libopenblas.so
```

(Note that the LAPACK library from OpenBLAS is used instead of the one from
ScaLAPACK, not sure why, but ScaLAPACK did not seem to work). With these
variables set, the external dependencies work correctly from MATLAB. However,
they do mess up MATLAB internal linear algebra calls. For example, solving a
linear system of equations with `mldivide` gives a segmentation violation
(which is again hard to debug, since the stack trace points to MATLAB libraries
which we cannot compile with debug symbols). The solution proposed here is to
replace MATLAB linear algebra calls with functions compiled using `mex` that
directly call the correct BLAS library. The command to compile such a function
is, after loading MATLAB and OpenBLAS modules:

```bash
mex mldivide_lapack.c -lopenblas
```

Note that when solving an overdetermined system of equations, (hopefully)
small deviations between the output of `mldivide` and `mldivide_lapack` are
observed.

What is described above should be considered a workaround and in the future
a better solution could be found. A first possibility is to try and rebuild
the external dependencies with Intel MKL instead of OpenBLAS/ScaLAPACK. If the
Intel MKL version is close enough to the one shipped with MATLAB, that might
solve incompatibilities. A second option is that the problematic part of
PredSim in the context of this problem (which is generating the polynomial
fit) is done separately. This would require some changes to the MATLAB scripts
to do (part of) the preprocessing without the need to include any external
dependencies, and then there should be no interference from BLAS/LAPACK
libraries outside of the MATLAB installation.

## TODO

- Check if the `mex` output is compatible across architectures and different
  OpenBLAS versions
- Try to compile OpenSim and OpenSimAD with `intel` toolchain instead of `foss`
