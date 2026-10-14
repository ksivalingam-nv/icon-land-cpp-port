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
#include "execution_space.hpp"

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <type_traits>

namespace {

constexpr int stream_id = 0;

class ExecutionSpaceTest : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {
    // destroy the default execution space instance if it was created by the tests runner
    icon::destroy_default_execution_space();
  }

  static void TearDownTestSuite() {
    // recreate the default execution space instance for the tests runner
    icon::create_default_execution_space();
  }
};

TEST_F(ExecutionSpaceTest, CreateExecSpace) {
  // The singleton instance must be of the correct type for the build variant.
  icon::create_default_execution_space(stream_id);

  EXPECT_NE(icon::detail::default_execution_space_instance, nullptr);
  EXPECT_EQ(&icon::default_execution_space_instance(), icon::detail::default_execution_space_instance.get());

#if defined(KOKKOS_ENABLE_CUDA)
  EXPECT_EQ(icon::default_execution_space_instance().cuda_stream(), reinterpret_cast<cudaStream_t>(stream_id));
#elif defined(KOKKOS_ENABLE_HIP)
  EXPECT_EQ(icon::default_execution_space_instance().hip_stream(), reinterpret_cast<hipStream_t>(stream_id));
#elif defined(KOKKOS_ENABLE_OPENMP)
  EXPECT_EQ(icon::default_execution_space_instance(), Kokkos::OpenMP::execution_space());
#else
  EXPECT_EQ(icon::default_execution_space_instance(), Kokkos::Serial::execution_space());
#endif
  icon::destroy_default_execution_space();
}

TEST_F(ExecutionSpaceTest, NoDoubleInstantiation) {
  icon::create_default_execution_space(stream_id);
  // Attempting to create the singleton instance again should abort via finish().
  EXPECT_DEATH((icon::create_default_execution_space(stream_id)),
               "create_default_execution_space\\(\\) called more than once");
  icon::destroy_default_execution_space();
}

TEST_F(ExecutionSpaceTest, DestroyExecSpace) {
  icon::create_default_execution_space(stream_id);
  icon::destroy_default_execution_space();
  EXPECT_EQ(icon::detail::default_execution_space_instance, nullptr);
}

}  // namespace
