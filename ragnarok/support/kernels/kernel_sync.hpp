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

#ifndef RAGNAROK_SUPPORT_KERNELS_KERNEL_SYNC_H_
#define RAGNAROK_SUPPORT_KERNELS_KERNEL_SYNC_H_

#include <Kokkos_Core.hpp>

#include "icon_loops.hpp"
#include "support/common/execution_space.hpp"
#include "support/common/finish.hpp"

// Completion and readback primitives that are legal inside a block-loop body
// under every role assignment (see docs/icon_loops.md):
// - block-parallel builds: the body may run inside a parallel region where no
//   Kokkos call is allowed - here the primitives compile to plain reads or
//   nothing (index loops ran in place, mirrors alias their originals).
// - index-parallel builds: the block loop is a serial host loop and
//   index_loops launched asynchronous kernels - here the primitives emit the
//   fence, and a copy where the data is a distinct allocation, on the default
//   instance, the instance the kernels run on.
// The prefers_coalesced_access() branching below is infrastructure; physics
// code uses these verbs and never branches itself.

namespace icon {

// completion point for all previously launched kernels - bridge entry
// points call this before returning to Fortran
inline void kernel_fence() {
  if constexpr (prefers_coalesced_access()) default_execution_space_instance().fence();
}

// read back a kernel-computed scalar inside a block body
template <typename T, typename... P>
[[nodiscard]] inline T fetch_scalar(const Kokkos::View<T, P...>& v) {
  if constexpr (prefers_coalesced_access()) {
    T value;
    Kokkos::deep_copy(default_execution_space_instance(), value, v);
    default_execution_space_instance().fence();
    return value;
  } else {
    return v();
  }
}

// staging transfer of the same logical data (mirror pattern); a copy between
// *distinct* views is compute, not transfer - use an index_loop with an
// assignment body instead
template <class DstView, class SrcView>
inline void fetch(const DstView& dst, const SrcView& src) {
  if constexpr (prefers_coalesced_access()) {
    // on unified-memory systems the mirror pair may alias - then the fence
    // alone completes the kernels that produced src
    if (static_cast<const void*>(dst.data()) != static_cast<const void*>(src.data()))
      Kokkos::deep_copy(default_execution_space_instance(), dst, src);
    default_execution_space_instance().fence();
  } else {
    require(dst.data() == src.data());
  }
}

}  // end of namespace icon
#endif
