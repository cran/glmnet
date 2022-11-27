# glmnetpp C++ Core Engine

## Overview

The `glmnetpp` C++ core engine implements the core components of the
system.

## Dependencies

If you have set up a conda environment following the instructions in the main repo README, you should already have a C++ toolchain installed along with Bazel. If not, we require a C++ toolchain that supports C++-14 and [OpenMP](https://www.openmp.org/) and an installation of [Bazel](https://bazel.build/)

Suggested compilers:
- [GCC >= 9.3.0](https://gcc.gnu.org/)
- [Clang >= 10.0.0](https://clang.llvm.org/)

## Build

```
mamba update -y conda
mamba env create
conda activate glmnetpp
```

__Note__: On Linux, it's best to specify whether you want to use `clang` or `gcc`.
Add the appropriate flag to each `bazel` call below:
```
# For gcc
# For clang
bazel ... --config=gcc
bazel ... --config=clang
```

__To build `glmnetpp`__:
```
bazel build //:glmnetpp 
```
Note that `glmnetpp` is a header-only library,
so this will simply collect all the headers and register its dependencies.
For release mode, add the flag `-c opt` after `build`.
For debug mode, add the flag `-c dbg` after `build`.

__To run all tests__:
```
bazel test -c dbg //test/... 
```

__To run a particular test__:
```
bazel test -c dbg //test:name-of-test
```

__To run the benchmarks__:
```
bazel run -c opt //benchmark:name-of-benchmark
```