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

if(BUILD_TESTING)
  set(kokkos_dependents tests_runner.cpp)

  set_source_files_properties(
    ${kokkos_dependents} TARGET_DIRECTORY ragnarok_test
    PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
  )

  target_sources(ragnarok_test PRIVATE ${kokkos_dependents})

  if(CMAKE_CXX_COMPILER_ID STREQUAL "CrayClang")
    set_target_properties(ragnarok_test PROPERTIES LINKER_LANGUAGE Fortran)
  endif()
endif()
