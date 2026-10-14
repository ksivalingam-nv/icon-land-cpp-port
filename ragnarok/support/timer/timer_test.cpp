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
#include "timer.hpp"

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>

namespace {
using namespace timer;

// constructed before main: covers id issuance at static-initialization time
Timer early_timer("early_timer");

// the test runner adopts RAGNAROK_TIMERS_LEVEL before any test runs; the
// tests written against the default timers_level skip under an override
bool level_overridden() { return std::getenv("RAGNAROK_TIMERS_LEVEL") != nullptr; }

TEST(TimerTests, CheckConfiguredLevel) {
  if (not level_overridden()) GTEST_SKIP() << "needs RAGNAROK_TIMERS_LEVEL";
  // a fresh timer's activity reflects the adopted level:
  const int level = std::stoi(std::getenv("RAGNAROK_TIMERS_LEVEL"));
  Timer threshold1_probe("threshold1_probe");
  EXPECT_EQ(threshold1_probe.active(), level >= 1);
  if (level < 1) return;  // thresholds below 1 are rejected
  Timer at_level_probe("at_level_probe", level);
  Timer above_level_probe("above_level_probe", level + 1);
  EXPECT_TRUE(at_level_probe.active());
  EXPECT_FALSE(above_level_probe.active());
}

TEST(TimerTests, CheckSecondConfigurationAborts) {
  // the runner has configured already; timers_level is a runtime constant
  EXPECT_DEATH(configure(1), "configured more than once");
}

TEST(TimerTests, CheckNestingEnforcement) {
  if (level_overridden()) GTEST_SKIP() << "assumes the default timers_level";
  static Timer outer_timer("outer_timer");
  static Timer inner_timer("inner_timer");

  // properly nested measurements pass:
  outer_timer.start();
  inner_timer.start();
  inner_timer.stop();
  outer_timer.stop();

  // an ill-nested stop aborts naming the innermost running timer:
  EXPECT_DEATH(
      {
        outer_timer.start();
        inner_timer.start();
        outer_timer.stop();
      },
      "innermost running timer");
}

TEST(TimerTests, CheckDestructionWhileRunning) {
  if (level_overridden()) GTEST_SKIP() << "assumes the default timers_level";
  // a timer leaving scope with its measurement open aborts:
  EXPECT_DEATH(
      {
        Timer doomed_timer("doomed_timer");
        doomed_timer.start();
      },
      "during its measurement");
}

TEST(TimerTests, CheckTimerMeasurement) {
  if (level_overridden()) GTEST_SKIP() << "assumes the default timers_level";
  Timer timer1("timer1");
  timer1.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(42));
  timer1.stop();

  // test rough precision of timer measurement:
  double value1 = timer1.value();
  EXPECT_GE(value1, 0.042);
  EXPECT_LT(value1, 1.0);
}

TEST(TimerTests, CheckEarlyConstructedTimer) {
  early_timer.start();
  early_timer.stop();
  EXPECT_TRUE(early_timer.value() >= 0.0);
}

TEST(TimerTests, CheckLevelFilter) {
  if (level_overridden()) GTEST_SKIP() << "assumes the default timers_level";
  // with the default timers_level 1, a timer activating at level 1
  // records and one activating at level 2 is filtered
  Timer level1_timer("level1_timer", 1);
  Timer level2_timer("level2_timer", 2);
  EXPECT_TRUE(level1_timer.active());
  EXPECT_FALSE(level2_timer.active());
  level1_timer.start();
  level1_timer.stop();
  level2_timer.start();
  level2_timer.stop();

  std::stringstream report;
  std::streambuf* cout_buf = std::cout.rdbuf(report.rdbuf());
  report_all_timers();
  std::cout.rdbuf(cout_buf);
  EXPECT_TRUE(report.str().find("level1_timer: ncalls=1,") != std::string::npos);
  EXPECT_TRUE(report.str().find("level2_timer") == std::string::npos);
}

TEST(TimerTests, CheckParallelKernelMeasurement) {
  if (level_overridden()) GTEST_SKIP() << "assumes the default timers_level";
  // measure inside a kernel on the host backend (multithreaded when OpenMP):
  Timer kernel_timer("kernel_timer");
  constexpr int n = 1000;
  Kokkos::parallel_for("timed_kernel", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
                       [&kernel_timer](int) {
                         kernel_timer.start();
                         kernel_timer.stop();
                       });

  // all threads' measurements appear in the report:
  std::stringstream report;
  std::streambuf* cout_buf = std::cout.rdbuf(report.rdbuf());
  report_all_timers();
  std::cout.rdbuf(cout_buf);

  // every pool thread gets iterations (static schedule, n far above the
  // thread count), so the report names exactly the configured team:
  const std::string expected = "kernel_timer: ncalls=" + std::to_string(n) + ", ";
  const int threads          = Kokkos::DefaultHostExecutionSpace().concurrency();
  if (threads > 1) {
    EXPECT_TRUE(report.str().find(expected + "threads=" + std::to_string(threads) + ",") != std::string::npos);
  } else {
    EXPECT_TRUE(report.str().find(expected + "total=") != std::string::npos);
  }
}

}  // namespace
