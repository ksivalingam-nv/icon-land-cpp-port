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

# system paths on levante
set(GTest_ROOT
    "/sw/spack-levante/googletest-1.10.0-opzgcq"
    CACHE PATH "Path to Googletest on Levante" FORCE
)

set(GCC_ROOT
    "/sw/spack-levante/gcc-11.2.0-bcn7mb"
    CACHE PATH "Path to GCC on Levante" FORCE
)
set(INTEL_ROOT
    "/sw/spack-levante/intel-oneapi-compilers-2025.3.2-hrqlez/compiler/2025.3"
    CACHE PATH "Path to Intel oneAPI Compilers on Levante" FORCE
)

# cmake options
set(CMAKE_CXX_COMPILER
    "${INTEL_ROOT}/bin/icpx"
    CACHE STRING "C++ compiler" FORCE
)
set(CMAKE_C_COMPILER
    "${INTEL_ROOT}/bin/icx"
    CACHE STRING "C compiler" FORCE
)
set(CMAKE_Fortran_COMPILER
    "${INTEL_ROOT}/bin/ifx"
    CACHE STRING "Fortran compiler" FORCE
)
set(CMAKE_BUILD_TYPE
    "RelWithDebInfo"
    CACHE STRING "Build type" FORCE
)
set(CMAKE_BUILD_RPATH
    "${CMAKE_BUILD_RPATH}:${GCC_ROOT}/lib64"
    CACHE INTERNAL "Extending -rpath"
)
set(CMAKE_CXX_FLAGS
    "-g -fp-model=strict --gcc-toolchain=${GCC_ROOT} -Wall -Werror"
    CACHE STRING "C++ compiler flags" FORCE
)

# ragnarok options
set(BUILD_TESTING
    ON
    CACHE BOOL "Enable testing" FORCE
)
set(RGK_ENABLE_STANDALONE
    ON
    CACHE BOOL "Enable standalone" FORCE
)

# Kokkos options
set(Kokkos_ENABLE_SERIAL
    ON
    CACHE BOOL "Build Kokkos Serial backend" FORCE
)
set(Kokkos_ENABLE_OPENMP
    ON
    CACHE BOOL "Build Kokkos OpenMP backend" FORCE
)
