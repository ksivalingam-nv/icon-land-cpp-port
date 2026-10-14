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

# system paths on lumi

set(CRAY_ROOT "/opt/cray/pe/craype/2.7.34")

# cmake options
set(CMAKE_CXX_COMPILER
    "${CRAY_ROOT}/bin/CC"
    CACHE STRING "C++ compiler" FORCE
)
set(CMAKE_C_COMPILER
    "${CRAY_ROOT}/bin/cc"
    CACHE STRING "C compiler" FORCE
)
set(CMAKE_Fortran_COMPILER
    "${CRAY_ROOT}/bin/ftn"
    CACHE STRING "Fortran compiler" FORCE
)
set(CMAKE_CXX_FLAGS
    "-O0 -ffp-model=strict"
    CACHE STRING "CXX flags" FORCE
)
set(CMAKE_HIP_ARCHITECTURES
    gfx90a
    CACHE STRING "Hip Arch" FORCE
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
    ON
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
set(Kokkos_ARCH_AMD_GFX90A
    ON
    CACHE BOOL "Build Kokkos for AMD  arch" FORCE
)
