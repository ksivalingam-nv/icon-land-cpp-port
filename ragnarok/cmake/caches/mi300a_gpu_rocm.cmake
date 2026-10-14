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

# system paths on the MI300A node of DKRZ's `gh` cluster (mi300a.d.gh.dkrz.de).
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
set(CMAKE_HIP_COMPILER
    "${ROCM_ROOT}/bin/amdclang++"
    CACHE STRING "HIP compiler" FORCE
)
# Kokkos is compiled as the HIP language
# (Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE with the HIP backend), so
# CMAKE_CXX_FLAGS alone would not reach those objects.
set(CMAKE_HIP_FLAGS
    "-Wno-deprecated-attributes"
    CACHE STRING "HIP flags" FORCE
)
set(CMAKE_CXX_FLAGS
    "-Wno-deprecated-attributes"
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
set(RGK_ENABLE_STANDALONE
    OFF
    CACHE BOOL "Enable standalone" FORCE
)

# kokkos options
set(Kokkos_ENABLE_SERIAL
    ON
    CACHE BOOL "Build Kokkos Serial backend" FORCE
)
set(Kokkos_ENABLE_HIP
    ON
    CACHE BOOL "Build Kokkos HIP backend" FORCE
)
set(Kokkos_ARCH_AMD_GFX942_APU
    ON
    CACHE BOOL "Build Kokkos for AMD MI300A APU arch" FORCE
)
