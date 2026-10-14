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
#ifndef RAGNAROK_SUPPORT_COMMON_EXECUTION_SPACE_H_
#define RAGNAROK_SUPPORT_COMMON_EXECUTION_SPACE_H_

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <memory>

#include "support/common/finish.hpp"

/// @file
/// @brief The execution space instance ragnarok kernels dispatch on by default.
///
/// icon::default_execution_space_instance() is bound, on CUDA/HIP builds, to the CUDA/HIP
/// stream behind OpenACC async queue 1 - the queue ICON's physics kernels run on. Kokkos
/// kernels launched through this instance are therefore ordered against the surrounding
/// !$ACC ... ASYNC(1) region by the stream itself, without an explicit !$ACC WAIT(1) before
/// a HOST_DATA USE_DEVICE handoff.
namespace icon {

namespace detail {
inline std::unique_ptr<Kokkos::DefaultExecutionSpace> default_execution_space_instance;
}  // namespace detail

inline Kokkos::DefaultExecutionSpace& default_execution_space_instance() {
  icon::require(detail::default_execution_space_instance.get() != nullptr,
                "default_execution_space_instance() called before create_default_execution_space()");
  return *detail::default_execution_space_instance;
}

inline void create_default_execution_space() {
  icon::require(detail::default_execution_space_instance.get() == nullptr,
                "create_default_execution_space() called more than once");
  detail::default_execution_space_instance = std::make_unique<Kokkos::DefaultExecutionSpace>();
}

inline void create_default_execution_space(std::intptr_t stream_id) {
  icon::require(detail::default_execution_space_instance.get() == nullptr,
                "create_default_execution_space() called more than once");
#if defined(KOKKOS_ENABLE_CUDA)
  detail::default_execution_space_instance =
      std::make_unique<Kokkos::DefaultExecutionSpace>(reinterpret_cast<cudaStream_t>(stream_id));
#elif defined(KOKKOS_ENABLE_HIP)
  detail::default_execution_space_instance =
      std::make_unique<Kokkos::DefaultExecutionSpace>(reinterpret_cast<hipStream_t>(stream_id));
#else
  create_default_execution_space();
#endif
}

/// Drops the held instance. Must run before Kokkos::finalize() - destroying a stream-bound
/// execution space instance after Kokkos has torn down its runtime is undefined behaviour.
inline void destroy_default_execution_space() {
  icon::require(detail::default_execution_space_instance.get() != nullptr,
                "destroy_default_execution_space() called before create_default_execution_space()");
  detail::default_execution_space_instance.reset();
}

}  // namespace icon

#endif  // RAGNAROK_SUPPORT_COMMON_EXECUTION_SPACE_H_
