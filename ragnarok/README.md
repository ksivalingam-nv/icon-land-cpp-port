# ICON RAGNAROK

Ragnarok is the C++ implementation of the AES physics based on the [Kokkos](https://github.com/kokkos/kokkos) programming model.


## Project Structure

```
└── 📁ragnarok
    └── 📁cmake
        └── 📁caches
            ├── levante_cpu_gcc.cmake
            ├── levante_cpu_intel.cmake
            ├── levante_gpu_nvhpc.cmake
            ├── lumi_cpu_cray.cmake
            ├── lumi_gpu_cray.cmake
    └── 📁common
        ├── CMakeLists.txt
        ├── types_test.cpp
        ├── types.hpp
    └── 📁support
        └── 📁icon_bridge
        └── 📁testing
        └── 📁timer
        ├── CMakeLists.txt
    └── 📁aes_component
        ├── CMakeLists.txt
    ├── CMakeLists.txt
    ├── CONTRIBUTING.md
    └── README.md
```

Shared functionality that is used across physics components belongs in the `common/` directory.

The `support/` directory provides infrastructure used throughout Ragnarok:

* `icon_bridge` provides access to the data and functionality available in the main ICON code. Additionally, it also provides initialize and finalize routines for Ragnarok. Its support will be extended as needed. Design and workflow notes live in [`support/icon_bridge/docs/`](support/icon_bridge/docs/): [patch_recon.md](support/icon_bridge/docs/patch_recon.md) establishes the general framework for representing ICON state in C++; [intp_recon.md](support/icon_bridge/docs/intp_recon.md), [config_recon.md](support/icon_bridge/docs/config_recon.md) and [nh_metrics_recon.md](support/icon_bridge/docs/nh_metrics_recon.md) build on it for interpolation, configuration and nonhydrostatic-metrics state.
* `kernels` provides the switchable-parallelism loop layer for writing ICON block/index loop nests once, with the build deciding which loop runs parallel; usage is described in [icon_loops.md](support/kernels/docs/icon_loops.md).
* `memory` provides temporary work arrays for physics kernels from a reusable buffer; usage is described in [scratch_memory.md](support/memory/docs/scratch_memory.md).
* `timer` mirrors the ICON timer infrastructure with one implementation for all builds: measurements start on a ragnarok-internal wall-clock backend; in the ICON build the bridge setup switches to the ICON Fortran timers, so Ragnarok timers appear in the experiment logfile. Usage is described in [timer_usage.md](support/timer/docs/timer_usage.md).
* `testing` provides the entry point for standalone unit testing and is also used in the GitLab CI pipelines.

## Usage

Ragnarok is integrated into ICON and is built automatically when ICON is configured with `--enable-ragnarok`. The ICON build system also ensures that the correct compilation flags and the appropriate Kokkos backend (from the `externals/kokkos` submodule) are used.

## Available Compile Options

> **Note:** default options are marked in __bold__

* _Unit-test_ - compile unit-tests
  * BUILD_TESTING=__ON__/OFF
* Standalone - compile standalone components
  * RGK_ENABLE_STANDALONE=__ON__/OFF
* Bindings - interoperability with other programming languages; Fortran is supported by default.
  * RGK_ENABLE_PYTHON_BINDINGS=ON/**OFF**
* Precision - either single or double precision
  * RGK_ENABLE_SINGLE_PRECISION=ON/**OFF**

### Build as part of ICON

The ICON configuration scripts automatically configure both Ragnarok and Kokkos.

### Ragnarok standalone

The standalone has a dependency of [NetCDF for CXX](https://github.com/Unidata/netcdf-cxx4). For Levante: `spack load netcdf-cxx4@4.3.1`.

To build ragnarok standalone, one can use the cache files provided for levante under `ragnarok/cmake/caches/`, e.g.:
  ```
  cmake -B build \
        -S ragnarok \
        -C ragnarok/cmake/caches/levante_cpu.cmake \
        && cmake --build build --parallel
  ```

The aes_microphysics executable can be run in single or double precision, e.g.:
  ```
    ./build/aes_microphysics/standalone/aes_graupel <input.nc> double  
  ```
The results are saved to NetCDF file called `output.nc`.

### Interoperability with python

The aes_microphysics provides python bindings for connecting the library to python environments. To ensure that, python3 (&pybind11) and netcdf-cxx4 to be accessible. For Levante:
  ```
  module load python3
  python3 -m venv venv
  source venv/bin/activate
  pip install -r ragnarok/aes_microphysics/bindings/python/requirements.txt
  cmake -B build_py \
        -S ragnarok \
        -C ragnarok/cmake/caches/levante_cpu.cmake \
        -DRGK_ENABLE_PYTHON_BINDINGS=ON \
        -Dpybind11_DIR=<path-to-venv>/venv/lib/python3.10/site-packages/pybind11/share/cmake/pybind11/
  cmake --build build_py --parallel
  ```

The following shared library will be generated `build_py/aes_microphysics/aes_muphys_py.so`, which can be later imported by a python source code.

## Testing

Ragnarok uses GoogleTest for unit testing. To run the tests:

  ```
  cd build && ctest
  ```
