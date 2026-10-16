# Using SELF to Build your own applications

When you install SELF, it installs a CMake package alongside the library and its Fortran module files.
A downstream project finds that package with `find_package` and links a single imported target, `SELF::self`.
This page explains what the package provides and shows a minimal project that uses it.

## What the package installs

Running `make install` (or `cmake --install`) on a SELF build writes the following under the install prefix.

* `lib/libself.so`, the SELF library.
* `include/`, the Fortran module files (`self_constants.mod`, `self_mesh_3d.mod`, and the rest) and `SELF_Macros.h`.
* `lib/cmake/SELF/`, the package files `SELFConfig.cmake`, `SELFConfigVersion.cmake`, and `SELFTargets.cmake`.

The imported target `SELF::self` carries everything a consumer needs to compile and link against SELF.
Its include directories point at the installed `include/` directory, so the SELF modules are found without any further flags.
Its link interface names MPI (`MPI::MPI_Fortran`), HDF5 (`HDF5::HDF5`), and FEQ-Parse, and `SELFConfig.cmake` finds those dependencies again on your machine when you call `find_package(SELF)`.
When SELF was built with GPU support, the package also finds HIP or the CUDA toolkit, whichever backend SELF was built with.

## A minimal project

The project below builds one program against an installed SELF.

```cmake
cmake_minimum_required(VERSION 3.21)
project(myapp LANGUAGES Fortran C)

find_package(SELF CONFIG REQUIRED)

add_executable(myapp myapp.f90)
target_link_libraries(myapp PRIVATE SELF::self)
```

The program can use any of the SELF modules directly.

```fortran
program myapp
  use SELF_Constants
  use SELF_Mesh_3D
  implicit none
  type(Mesh3D) :: mesh
  integer :: bcids(1:6)

  bcids(1:6) = 0
  call mesh%StructuredMesh(2,2,2,2,2,2,0.25_prec,0.25_prec,0.25_prec,bcids)
  print*,"Number of elements : ",mesh%nElem
  call mesh%Free()
endprogram myapp
```

Point CMake at the SELF install prefix when you configure the project.
If SELF's dependencies were installed with Spack, activate the same Spack environment first so that MPI, HDF5, and FEQ-Parse are found where SELF found them.

```shell
cmake -S . -B build -DCMAKE_PREFIX_PATH=${HOME}/opt/self
cmake --build build
```

The `test/external_consumer` directory of the SELF repository is a complete example of this kind of project, and continuous integration builds it against a fresh SELF install.

## Matching precision and backend

SELF is compiled either in double precision, the default, or in single precision when it is configured with `-DSELF_ENABLE_DOUBLE_PRECISION=OFF`.
The real kind `prec` in `SELF_Constants` is fixed when the library is compiled, and every SELF module that you `use` is compiled with it.
Your own code must use the same precision, and a program that hands single precision reals to a double precision SELF will not compile.
Code that selects a precision with `#ifdef DOUBLE_PRECISION`, for instance to choose an MPI datatype, is more dangerous, because it compiles either way and silently disagrees with the library when the definition is missing.

To make this automatic, `SELF::self` passes the preprocessor definitions that SELF itself was compiled with on to every target that links it.
These are `DOUBLE_PRECISION` for a double precision build, `ENABLE_GPU` or `ENABLE_APU` for a device build, `HAVE_HIP` or `HAVE_CUDA` for the device backend, and `MULTITHREADING` when CPU multithreading is enabled.
A source file only sees them when it is preprocessed, so name it with a `.F90` suffix or set the `Fortran_PREPROCESS` source property.
Declaring reals as `real(prec)` rather than with a literal kind is the simplest way to stay consistent, because `prec` always follows the library.

The package also records the configuration in variables that you can test after `find_package(SELF)`.

* `SELF_DOUBLE_PRECISION` is `ON` when SELF was built in double precision.
* `SELF_ENABLE_GPU` and `SELF_ENABLE_APU` say whether a device backend was built.
* `SELF_GPU_BACKEND` is `HIP`, `CUDA`, or empty for a CPU build.
* `SELF_COMPILE_DEFINITIONS` lists the definitions carried by `SELF::self`.
* `SELF_Fortran_COMPILER_ID` and `SELF_Fortran_COMPILER_VERSION` name the compiler that wrote the module files.

A project with its own precision option can compare it against `SELF_DOUBLE_PRECISION` and stop with an error at configure time rather than fail later.

Fortran module files can only be read by the compiler, and usually only by the compiler version, that wrote them.
Build your project with the same Fortran compiler that built SELF; the package warns at configure time when the compiler vendor differs.
