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

/// @file
/// @brief C++ and C interface for initializing, finalizing, and querying the ragnarok backend, including Kokkos/OpenMP
/// integration and Fortran interoperability.
///
/// This header provides:
/// - Initialization/finalization routines for the ragnarok backend
/// - Management of serial execution spaces (Kokkos/OpenMP)
/// - Functions for querying backend state and Kokkos version
/// - C bindings for Fortran interoperability

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_RAGNAROK_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_RAGNAROK_H_

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <type_traits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ragnarok {

/// @brief: Returns true if the code was compiled with ICON for CPU target
bool constexpr isCPU() {
#ifndef __STANDALONE
  using MemorySpace            = Kokkos::DefaultExecutionSpace::memory_space;
  constexpr bool is_host_space = std::is_same<MemorySpace, Kokkos::HostSpace>::value;
  if constexpr (is_host_space) {
    return true;
  } else {
    return false;
  }
#else
  return false;
#endif
}

inline std::vector<Kokkos::Serial> serial_spaces;
inline void initialize_serial_backend() {
#ifdef _OPENMP
  std::vector<int> weights(omp_get_max_threads(), 1);
  serial_spaces = Kokkos::Experimental::partition_space(Kokkos::Serial(), weights);
#else
  Kokkos::Serial s;
  serial_spaces.push_back(s);
#endif
}

inline void finalize_serial_backend() { serial_spaces.clear(); }

inline auto& get_serial_exec_space() {
#ifdef _OPENMP
  int tid = omp_get_thread_num();
  return serial_spaces[tid];
#else
  return serial_spaces[0];
#endif
}

bool is_initialized();

// device_id >= 0 pins Kokkos to that GPU; -1 lets Kokkos choose (host builds)
// stream_id is used when a GPU backend is used;
void init(int device_id, std::intptr_t stream_id);

const std::string& retrieve_kokkos_version();

}  // namespace ragnarok

extern "C" {
void init_ragnarok_c(int, int, std::intptr_t);
const char* retrieve_kokkos_version_c(int&);
void ragnarok_ensure_empty_memory_pool();
}

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_RAGNAROK_H_
