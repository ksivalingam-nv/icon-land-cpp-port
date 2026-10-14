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

#include "ragnarok.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <cassert>
#include <string>

#include "finish.hpp"
#include "logging.hpp"
#include "support/common/execution_space.hpp"
#include "support/memory/memory_pool.hpp"

namespace ragnarok {

enum class InitState { pre_init = 0, initialized, finalized, problem };

InitState init_state = InitState::pre_init;

bool is_initialized() { return init_state == InitState::initialized; }

namespace {

std::string construct_kokkos_version() {
  return std::to_string(KOKKOS_VERSION_MAJOR) + "." + std::to_string(KOKKOS_VERSION_MINOR) + "." +
         std::to_string(KOKKOS_VERSION_PATCH);
}

}  // namespace

void init(int device_id, std::intptr_t stream_id) {
  icon::require(init_state == InitState::pre_init);
  // Kokkos must not pick its own GPU: its local-rank heuristic can disagree with
  // the per-rank device ICON already bound for the OpenACC runtime, and a
  // mid-init cudaSetDevice/hipSetDevice to another device breaks the OpenACC async queues.
  if (device_id >= 0) {
    Kokkos::initialize(Kokkos::InitializationSettings().set_device_id(device_id));
  } else {
    Kokkos::initialize();
  }

  icon::create_default_execution_space(stream_id);
  if (ragnarok::isCPU()) {
    initialize_serial_backend();
  }
  init_state = InitState::initialized;
}

void finalize() {
  assert(init_state == InitState::initialized);
  if (ragnarok::isCPU()) {
    finalize_serial_backend();
  }
  icon::destroy_default_execution_space();
  Kokkos::finalize();
  // enter final state - we are not even allowed to call Kokkos::initialize again:
  init_state = InitState::finalized;
}

const std::string& retrieve_kokkos_version() {
  const static std::string version = construct_kokkos_version();
  return version;
}

}  // namespace ragnarok

extern "C" {

void init_ragnarok_c(int is_stdio_proc, int device_id, std::intptr_t stream_id) {
  ragnarok::init(device_id, stream_id);
  icon::logging::set_stdio_proc(is_stdio_proc != 0);
}

void finalize_ragnarok() { ragnarok::finalize(); }

const char* retrieve_kokkos_version_c(int& length) {
  const std::string& version = ragnarok::retrieve_kokkos_version();
  length                     = static_cast<int>(version.length());
  return version.c_str();
}

void ragnarok_ensure_empty_memory_pool() {
  ragnarok::MemoryPool<Kokkos::DefaultExecutionSpace>::ensure_empty_memory_pool();
}

}  // extern "C"
