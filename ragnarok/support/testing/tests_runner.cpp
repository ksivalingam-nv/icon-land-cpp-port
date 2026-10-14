// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------
#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <cstdlib>
#include <string>

#include "support/timer/timer.hpp"

#ifndef __STANDALONE
#include "support/icon_bridge/ragnarok.hpp"
#endif

#include "support/common/execution_space.hpp"

int main(int argc, char** argv) {
  // Initialize Kokkos before any tests run.
  Kokkos::initialize(argc, argv);

#ifndef __STANDALONE
  ragnarok::initialize_serial_backend();
#endif

  // test builds have no namelist; the runner states the timers intent -
  // RAGNAROK_TIMERS_LEVEL, or the minimal active level
  const char* level = std::getenv("RAGNAROK_TIMERS_LEVEL");
  timer::configure(level ? std::stoi(level) : 1);

  // for the purpose of the unit-tests with no ICON connected to ragnarok,
  // we use the Kokkos::DefaultExecutionSpace as the default execution space instance
  icon::create_default_execution_space();

  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

#ifndef __STANDALONE
  ragnarok::finalize_serial_backend();
#endif

  icon::destroy_default_execution_space();

  // Finalize Kokkos after all tests have completed.
  Kokkos::finalize();
  return result;
}
