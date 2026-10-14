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

#ifndef RAGNAROK_SUPPORT_KERNELS_ICON_LOOPS_H_
#define RAGNAROK_SUPPORT_KERNELS_ICON_LOOPS_H_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <type_traits>

#include "grid_range.hpp"
#include "support/common/execution_space.hpp"

// capture markers for the switchable loops below; both capture by value, so
// enclosing scalars are read-only copies in every build - order-dependent
// state belongs in serial_loop. index_lambda carries the device annotation a
// kernel needs on device builds (and is plain [=] otherwise). Deliberately
// lowercase and unprefixed: each marker pairs with its loop so a nest reads
// as one construct, and the declaration-in-argument syntax exposes the macro
// at every use site.
#define block_lambda [=]
#define index_lambda KOKKOS_LAMBDA

namespace icon {

// the one build-mode predicate of the kernels layer: execution spaces whose
// performance lives on coalesced access across the parallel index
// parallelize the stride-1 (index) dimension of (..., nlev, nproma) fields;
// host spaces parallelize blocks. Equivalently: true iff index_loop launches
// kernels (which kernel_sync.hpp completes). Explicit list: extend when
// enabling a new device backend.
constexpr bool prefers_coalesced_access() {
  bool coalescing = false;
#if defined(KOKKOS_ENABLE_CUDA)
  coalescing |= std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::Cuda>;
#endif
#if defined(KOKKOS_ENABLE_HIP)
  coalescing |= std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::HIP>;
#endif
  return coalescing;
}

// a device backend outside the list above must not fall through to the host
// roles - fail the build instead of silently mis-parallelizing
static_assert(
    prefers_coalesced_access() ||
        Kokkos::SpaceAccessibility<Kokkos::HostSpace, Kokkos::DefaultExecutionSpace::memory_space>::accessible,
    "unsupported device execution space - extend prefers_coalesced_access()");

// block_loop/index_loop swap parallelization roles with the build: host
// builds parallelize the block loop (index body serial), device builds the
// index loop (block body serial, one kernel launch per block - good for few
// large blocks, poor for many small ones). Loop bodies must be valid under
// both roles; a label reaches the profiler only on the side that is parallel
// in the current build.

// run body(jb) over the blocks of grid_range. Ordering, not completion:
// kernels launched by nested index_loop calls execute in submission order on
// the default instance but may still be in flight on return - host access to
// their results goes through the kernel_sync verbs
template <class Body>
inline void block_loop(const char* label, const GridRange& grid_range, const Body& body) {
  if constexpr (prefers_coalesced_access()) {
    const auto brange = grid_range.block_range();
    for (int jb = brange.start; jb < brange.end; ++jb) body(jb);
  } else {
    Kokkos::parallel_for(label, grid_range.block_range().range_policy(icon::default_execution_space_instance()), body);
  }
}

template <class Body>
inline void block_loop(const GridRange& grid_range, const Body& body) {
  block_loop("icon::block_loop", grid_range, body);
}

// run body(jc) over indices; asynchronous on device builds - complete via
// the kernel_sync verbs before host access. The IndexRange form covers call
// sites where the block loop lives elsewhere (e.g. on the Fortran side).
template <class Body>
inline void index_loop(const char* label, const IndexRange indices, const Body& body) {
  if constexpr (prefers_coalesced_access()) {
    Kokkos::parallel_for(label, indices.range_policy(icon::default_execution_space_instance()), body);
  } else {
    for (int jc = indices.start; jc < indices.end; ++jc) body(jc);
  }
}

template <class Body>
inline void index_loop(const IndexRange indices, const Body& body) {
  index_loop("icon::index_loop", indices, body);
}

// collapsed level+index loop (OpenACC COLLAPSE(2) analog): run body(jk, jc)
// over levels x indices, jc fastest (stride 1 under LayoutRight); asynchronous
// on device builds like the single-index index_loop. One flat RangePolicy, jc
// recovered by div/mod: coalesced access without MDRangePolicy tiling
// heuristics.
template <class Body>
inline void index_loop(const char* label, const LevelRange levels, const IndexRange indices, const Body& body) {
  const int istart = indices.start;
  const int kstart = levels.start;
  const int nidx   = indices.size();
  const int nlev   = levels.size();
  if constexpr (prefers_coalesced_access()) {
    Kokkos::parallel_for(
        label, Kokkos::RangePolicy(icon::default_execution_space_instance(), 0, nidx * nlev),
        KOKKOS_LAMBDA(const int f) {
          const int jk = kstart + f / nidx;
          const int jc = istart + f % nidx;
          body(jk, jc);
        });
  } else {
    for (int jk = levels.start; jk < levels.end; ++jk)
      for (int jc = indices.start; jc < indices.end; ++jc) body(jk, jc);
  }
}

template <class Body>
inline void index_loop(const LevelRange levels, const IndexRange indices, const Body& body) {
  index_loop("icon::index_loop", levels, indices, body);
}

// collapsed extra x index loop: like the level x index form above, but the
// slow dimension carries no grid meaning (e.g. prognostics); it must be
// spelled as a plain Range - the exact-type constraint keeps the tagged
// ranges out of the untagged slot
template <std::same_as<Range> ExtraRange, class Body>
inline void index_loop(const char* label, const ExtraRange extra, const IndexRange indices, const Body& body) {
  const int istart = indices.start;
  const int estart = extra.start;
  const int nidx   = indices.size();
  if constexpr (prefers_coalesced_access()) {
    Kokkos::parallel_for(
        label, Kokkos::RangePolicy(icon::default_execution_space_instance(), 0, extra.size() * nidx),
        KOKKOS_LAMBDA(const int f) {
          const int jc = istart + f % nidx;
          const int i  = estart + f / nidx;
          body(i, jc);
        });
  } else {
    for (int i = extra.start; i < extra.end; ++i)
      for (int jc = indices.start; jc < indices.end; ++jc) body(i, jc);
  }
}

template <std::same_as<Range> ExtraRange, class Body>
inline void index_loop(const ExtraRange extra, const IndexRange indices, const Body& body) {
  index_loop("icon::index_loop", extra, indices, body);
}

// collapsed extra x level x index loop (COLLAPSE(3) analog): run
// body(i, jk, jc) over extra x levels x indices, jc fastest as above. The
// extra dimension is the slowest (leftmost view extent) and carries no grid
// meaning (e.g. tracers); it must be spelled as a plain Range - the
// exact-type constraint keeps the tagged ranges out of the untagged slot.
template <std::same_as<Range> ExtraRange, class Body>
inline void index_loop(const char* label, const ExtraRange extra, const LevelRange levels, const IndexRange indices,
                       const Body& body) {
  const int istart = indices.start;
  const int kstart = levels.start;
  const int estart = extra.start;
  const int nidx   = indices.size();
  const int nlev   = levels.size();
  if constexpr (prefers_coalesced_access()) {
    Kokkos::parallel_for(
        label, Kokkos::RangePolicy(icon::default_execution_space_instance(), 0, extra.size() * nlev * nidx),
        KOKKOS_LAMBDA(const int f) {
          const int jc = istart + f % nidx;
          const int t  = f / nidx;
          const int jk = kstart + t % nlev;
          const int i  = estart + t / nlev;
          body(i, jk, jc);
        });
  } else {
    for (int i = extra.start; i < extra.end; ++i)
      for (int jk = levels.start; jk < levels.end; ++jk)
        for (int jc = indices.start; jc < indices.end; ++jc) body(i, jk, jc);
  }
}

template <std::same_as<Range> ExtraRange, class Body>
inline void index_loop(const ExtraRange extra, const LevelRange levels, const IndexRange indices, const Body& body) {
  index_loop("icon::index_loop", extra, levels, indices, body);
}

// run body(jc) over the indices of block jb of grid_range
template <class Body>
inline void index_loop(const char* label, const GridRange& grid_range, const int jb, const Body& body) {
  index_loop(label, grid_range.index_range(jb), body);
}

template <class Body>
inline void index_loop(const GridRange& grid_range, const int jb, const Body& body) {
  index_loop("icon::index_loop", grid_range.index_range(jb), body);
}

// collapsed loop over levels x the indices of block jb of grid_range
template <class Body>
inline void index_loop(const char* label, const GridRange& grid_range, const int jb, const LevelRange levels,
                       const Body& body) {
  index_loop(label, levels, grid_range.index_range(jb), body);
}

template <class Body>
inline void index_loop(const GridRange& grid_range, const int jb, const LevelRange levels, const Body& body) {
  index_loop("icon::index_loop", levels, grid_range.index_range(jb), body);
}

// serial_loop: no parallelism at all, in contrast to the switchable loops
// above - a plain host loop over the same iteration spaces, ascending order
// with the rightmost body argument fastest, single-threaded in every build;
// the device is never involved. The body is any callable
// ([&], [=], ...): there are no build-dependent capture rules to hide, so no
// marker and no label. The body must not launch kernels. For CPU-first
// development stages: order-dependent state and host-only calls are fine
// here. Porting to the GPU = rename to index_loop, spell the lambda
// index_lambda, give the data a device home - the compiler then flags each
// scalar carry and host-only call left in the body.

// run body(j) over indices
template <class Body>
inline void serial_loop(const IndexRange indices, const Body& body) {
  for (int j = indices.start; j < indices.end; ++j) body(j);
}

// run body(jk, jc) over levels x indices
template <class Body>
inline void serial_loop(const LevelRange levels, const IndexRange indices, const Body& body) {
  for (int jk = levels.start; jk < levels.end; ++jk)
    for (int jc = indices.start; jc < indices.end; ++jc) body(jk, jc);
}

// run body(jc) over the indices of block jb of grid_range
template <class Body>
inline void serial_loop(const GridRange& grid_range, const int jb, const Body& body) {
  serial_loop(grid_range.index_range(jb), body);
}

// run body(jk, jc) over levels x the indices of block jb of grid_range
template <class Body>
inline void serial_loop(const GridRange& grid_range, const int jb, const LevelRange levels, const Body& body) {
  serial_loop(levels, grid_range.index_range(jb), body);
}

// run body(jb, jc) over the whole grid_range
template <class Body>
inline void serial_loop(const GridRange& grid_range, const Body& body) {
  const auto brange = grid_range.block_range();
  for (int jb = brange.start; jb < brange.end; ++jb) {
    const auto irange = grid_range.index_range(jb);
    for (int jc = irange.start; jc < irange.end; ++jc) body(jb, jc);
  }
}

// run body(jb, jk, jc) over the whole grid_range x levels
template <class Body>
inline void serial_loop(const GridRange& grid_range, const LevelRange levels, const Body& body) {
  const auto brange = grid_range.block_range();
  for (int jb = brange.start; jb < brange.end; ++jb) {
    const auto irange = grid_range.index_range(jb);
    for (int jk = levels.start; jk < levels.end; ++jk)
      for (int jc = irange.start; jc < irange.end; ++jc) body(jb, jk, jc);
  }
}

}  // end of namespace icon
#endif
