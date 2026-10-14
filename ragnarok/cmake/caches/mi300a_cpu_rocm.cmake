# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

set(ROCM_ROOT "/opt/rocm")

# cmake options
set(CMAKE_CXX_COMPILER
    "${ROCM_ROOT}/bin/amdclang++"
    CACHE STRING "C++ compiler" FORCE
)
set(CMAKE_C_COMPILER
    "${ROCM_ROOT}/bin/amdclang"
    CACHE STRING "C compiler" FORCE
)
set(CMAKE_Fortran_COMPILER
    "gfortran"
    CACHE STRING "Fortran compiler" FORCE
)
# -Wno-deprecated-attributes: see the same flag in mi300a_gpu_rocm.cmake.
# -ffp-contract=off is required, not a preference: ragnarok's unit tests compare
# against hardcoded reference values with a 4-ULP tolerance, and this compiler
# contracts a*b+c into an FMA by default, which does not match the arithmetic
# the references were recorded with.
set(CMAKE_CXX_FLAGS
    "-Wno-deprecated-attributes -ffp-contract=off"
    CACHE STRING "CXX flags" FORCE
)
set(CMAKE_PREFIX_PATH
    "${ROCM_ROOT}"
    CACHE STRING "ROCm install prefix" FORCE
)
set(CMAKE_BUILD_TYPE
    "RelWithDebInfo"
    CACHE STRING "Build type" FORCE
)

# ragnarok options
set(BUILD_TESTING
    ON
    CACHE BOOL "Enable testing" FORCE
)
# Avoid using NetCDF: this host has none installed
set(RGK_ENABLE_STANDALONE
    OFF
    CACHE BOOL "Enable standalone" FORCE
)

# kokkos options
set(Kokkos_ENABLE_SERIAL
    ON
    CACHE BOOL "Build Kokkos Serial backend" FORCE
)
set(Kokkos_ENABLE_OPENMP
    ON
    CACHE BOOL "Build Kokkos OpenMP backend" FORCE
)
set(Kokkos_ENABLE_HIP
    OFF
    CACHE BOOL "Build Kokkos HIP backend" FORCE
)
set(Kokkos_ARCH_ZEN4
    ON
    CACHE BOOL "Build for AMD ZEN4 arch" FORCE
)
