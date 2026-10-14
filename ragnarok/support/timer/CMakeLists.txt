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

set(kokkos_dependents timer_sync.cpp)

set_source_files_properties(
  ${kokkos_dependents} TARGET_DIRECTORY ragnarok
  PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
)

target_sources(
  ragnarok
  PRIVATE timer.cpp timer.hpp timer_ops.hpp timer_sync.hpp ${kokkos_dependents}
)

if(BUILD_TESTING)
  set(kokkos_dependents timer_test.cpp)

  set_source_files_properties(
    ${kokkos_dependents} TARGET_DIRECTORY ragnarok_test
    PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE}
  )

  target_sources(ragnarok_test PRIVATE ${kokkos_dependents})
endif()
