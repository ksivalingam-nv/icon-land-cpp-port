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
#ifndef RAGNAROK_SUPPORT_MEMORY_SCRATCH_VIEW_H_
#define RAGNAROK_SUPPORT_MEMORY_SCRATCH_VIEW_H_

/// @file scratch_view.hpp
/// @brief Allocate unmanaged Kokkos views from a ScratchSpace (a MemoryPool frame, see memory_pool.hpp).
///
/// The concrete Kokkos::View type on the LHS determines element type, layout, memory space,
/// and memory traits.  No template argument is needed on the RHS:
///
///   Kokkos::View<int*, ES, Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::Restrict>> index_iv
///       = ragnarok::scratch_view(local_space, nvol);
///
/// Do not use `auto` on the LHS: `auto` suppresses the implicit conversion and the
/// variable would hold a ScratchViewProxy, not a Kokkos::View. The proxy converts only
/// as an rvalue and is non-copyable, so any later lvalue use of such a variable fails to
/// compile (the stray `auto` declaration itself does compile).

#include <concepts>
#include <cstddef>
#include <tuple>
#include <type_traits>

#include "memory_pool.hpp"

namespace ragnarok {

namespace detail {

// alias template collapsing a pack element to T; indirection through a struct works around NVCC
template <typename T, typename>
struct Always {
  using type = T;
};
template <typename T, typename U>
using always_t = typename Always<T, U>::type;

}  // namespace detail

// ─────────────────────────────────────────────────────────────────────────────
// ScratchViewProxy — return type of scratch_view()
// ─────────────────────────────────────────────────────────────────────────────

/// Proxy returned by scratch_view(local_space, dims...).
///
/// `operator ViewType()` fires on assignment to a concrete Kokkos::View.  ViewType is then
/// fully known: element type, layout, memory space, and traits are all read from it.
/// The frame is asked for sizeof(element) * product(dims) bytes, and the view is
/// constructed from (ptr, dims...).
///
/// Extents must be fully dynamic (DataType = T*, T**, …); the allocation size is the
/// product of the runtime arguments — enforced below.  To support mixed static/dynamic
/// extents (Kokkos DataType T*[4] or mdspan std::extents<…>), query ViewType::static_extent(i)
/// inside operator ViewType() and multiply those compile-time factors into the allocation;
/// the constructor call is unchanged because it only ever receives the dynamic dims.
template <class MemSpace, std::integral... Dims>
struct [[nodiscard]] ScratchViewProxy {
  ScratchSpace<MemSpace>& frame;
  std::tuple<detail::always_t<std::size_t, Dims>...> dims;

  ScratchViewProxy(ScratchSpace<MemSpace>& frame_, Dims... dims_) noexcept : frame(frame_), dims(dims_...) {}

  // single-use: each conversion allocates, so a stored (auto) proxy converting twice
  // would silently produce two distinct buffers
  ScratchViewProxy(const ScratchViewProxy&)            = delete;
  ScratchViewProxy& operator=(const ScratchViewProxy&) = delete;

  // constrained so the proxy never competes in overload resolution for non-View targets
  template <typename ViewType>
    requires(Kokkos::is_view<ViewType>::value)
  operator ViewType() && noexcept {
    // static extents would make the view span more than we allocate — silent out-of-bounds
    static_assert(std::size_t(ViewType::rank_dynamic) == std::size_t(ViewType::rank),
                  "scratch_view: static extents are not supported");
    static_assert(sizeof...(Dims) == std::size_t(ViewType::rank),
                  "scratch_view: number of dims must match the view rank");
    // a mismatched space would silently wrap e.g. a host pointer in a device view
    static_assert(std::is_same_v<typename ViewType::memory_space, MemSpace>,
                  "scratch_view: view memory space must match the frame's memory space");
    // scratch memory starts uninitialized, so a const view of it could only ever read garbage
    static_assert(!std::is_const_v<typename ViewType::value_type>, "scratch_view: const element type is not supported");
    static_assert(ViewType::memory_traits::is_unmanaged,
                  "scratch_view: the view type must carry the Unmanaged memory trait");
    using T = typename ViewType::value_type;
    // alloc hands out raw bytes; it never constructs or destroys elements
    static_assert(std::is_trivially_default_constructible_v<T> && std::is_trivially_destructible_v<T>,
                  "scratch_view: element type must be trivial");
    // saturating product: a plain fold could wrap past the limit and back under it
    constexpr std::size_t kMax = detail::kMaxScratchElements;
    std::size_t n              = 1;
    std::apply([&n](auto... d) { ((n = (d != 0 && n > kMax / d) ? kMax + 1 : n * d), ...); }, dims);
    if (n > kMax) Kokkos::abort("scratch_view: more than 2^31-1 elements — negative or overflowed dim?");
    T* ptr = frame.template alloc<T>(n);
    return std::apply([ptr](auto... d) { return ViewType(ptr, d...); }, dims);
  }
};

/// View type and element type are inferred from the assignment target — do not use `auto`.
template <class MemSpace, std::integral... Dims>
ScratchViewProxy<MemSpace, Dims...> scratch_view(ScratchSpace<MemSpace>& frame, Dims... dims) noexcept {
  return ScratchViewProxy<MemSpace, Dims...>(frame, dims...);
}

/// Alias: the ScratchSpace type associated with an execution space.
/// Allows the syntax `ragnarok::LocalScratchSpace<ExecSpace>` instead of
/// `ragnarok::ScratchSpace<typename ExecSpace::memory_space>`.
template <typename ExecSpace>
using LocalScratchSpace = ScratchSpace<typename ExecSpace::memory_space>;

/// Wrapper function for get_local_space()
template <typename ExecSpace>
[[nodiscard]] LocalScratchSpace<ExecSpace> get_local_space() noexcept {
  return MemoryPool<ExecSpace>::get_local_space();
}

}  // namespace ragnarok

#endif  // RAGNAROK_SUPPORT_MEMORY_SCRATCH_VIEW_H_
