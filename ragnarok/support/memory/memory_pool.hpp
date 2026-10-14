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
#ifndef RAGNAROK_SUPPORT_MEMORY_MEMORY_POOL_H_
#define RAGNAROK_SUPPORT_MEMORY_MEMORY_POOL_H_

static_assert(__cplusplus >= 202002L, "memory_pool.hpp requires C++20 or later");

/// @file memory_pool.hpp
/// @brief Scratch allocator for per-call physics kernel scratch arrays.
///
/// Terms: the pool owns one large allocation per memory space, the *slab*. A *frame*
/// (ScratchSpace) hands out slab memory by advancing the slab's *offset*; destroying the
/// frame rewinds the offset, releasing all its allocations at once.
///
/// Typical usage:
///   MemoryPool<ES>::ensure_empty_memory_pool();   // once per timestep; aborts if a frame leaked
///   auto local_space = MemoryPool<ES>::get_local_space();
///   See scratch_view.hpp for application-level usage.
/// Nesting works: each get_local_space() starts from the current offset.
/// A single request is limited to 2^31-1 elements, matching the model's int indexing.
/// Slab auto-grows past the previous cycle's peak on the first call of a new cycle.
/// Host: per-thread slab (thread-safe). Device: global static; must not be called inside OpenMP.
///
/// Stream ordering: a frame's destructor rewinds the offset WITHOUT fencing. With a single
/// stream this is safe — kernels serialize. With more than one stream, going out of scope
/// before all kernels finished is undefined behaviour: synchronize via fence first.

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "support/common/execution_space.hpp"

// set this in your configure wrapper if you want to debug the memory usage
// #define MEMORY_POOL_DEBUG

namespace ragnarok {

// ─────────────────────────────────────────────────────────────────────────────
// Utilities
// ─────────────────────────────────────────────────────────────────────────────

namespace detail {

/// Round `n` up to the next multiple of `align` (must be a power of two).
constexpr std::size_t align_up(std::size_t n, std::size_t align) noexcept { return (n + align - 1) & ~(align - 1); }

/// Alignment of the slab base and of every allocation carved from it: covers device
/// coalescing granularity, host cachelines, and any over-aligned element type.
inline constexpr std::size_t kAllocAlign         = 128;

/// Upper bound on the elements of a single request: ICON indexes fields with int, so a
/// larger request is always a bug — typically a negative or int-overflowed dimension.
inline constexpr std::size_t kMaxScratchElements = 0x7fffffff;

/// True only for plain HostSpace. Host-accessibility is deliberately NOT the criterion:
/// unified spaces (CudaUVMSpace, SharedSpace, ...) are host-accessible but must still be
/// allocated via kokkos_malloc and treated as device memory (global slab, OpenMP guard).
template <class MemSpace>
inline constexpr bool is_host_space_v = std::is_same_v<MemSpace, Kokkos::HostSpace>;

/// User override (0 = unset); written by set_scratch_buffer_initial_size().
/// Atomic: host worker threads read it lazily when creating their thread-local slabs.
inline std::atomic<std::size_t> scratch_user_size_{0};

/// Initial slab size: set_scratch_buffer_initial_size() > RAGNAROK_SCRATCH_INIT_SIZE env > 64 MiB.
inline std::size_t scratch_initial_size() noexcept {
  const std::size_t user_size = scratch_user_size_.load(std::memory_order_relaxed);
  if (user_size != 0) return user_size;
  static const std::size_t env_value = []() noexcept -> std::size_t {
    constexpr std::size_t kDefault = 64ULL << 20;  // 64 MiB
    const char* env                = std::getenv("RAGNAROK_SCRATCH_INIT_SIZE");
    if (!env || env[0] == '\0') return kDefault;
    // strict plain-bytes integer: a typo like "64M" or "-1" must abort, not silently
    // become a tiny or huge slab
    errno                      = 0;
    char* end                  = nullptr;
    const unsigned long long v = std::strtoull(env, &end, 10);
    if (std::strchr(env, '-') || end == env || *end != '\0' || errno == ERANGE) {
      std::fprintf(stderr, "RAGNAROK_SCRATCH_INIT_SIZE: invalid value '%s' (plain bytes expected)\n", env);
      Kokkos::abort("RAGNAROK_SCRATCH_INIT_SIZE: invalid value");
    }
    return v != 0 ? static_cast<std::size_t>(v) : kDefault;
  }();
  return env_value;
}

}  // namespace detail

template <class MemSpace>
struct Slab;

namespace detail {

/// Slab storage, keyed by memory space: exec spaces sharing a memory space (e.g. Serial and
/// OpenMP on HostSpace) share one pool. Per-thread on host, global otherwise.
template <class MemSpace>
Slab<MemSpace>& pool_slab() noexcept {
  if constexpr (is_host_space_v<MemSpace>) {
    static thread_local Slab<MemSpace> s;
    return s;
  } else {
    static Slab<MemSpace> s;
    return s;
  }
}

}  // namespace detail

/// Override the initial slab size (bytes). Call before the first get_local_space().
inline void set_scratch_buffer_initial_size(std::size_t bytes) noexcept {
  detail::scratch_user_size_.store(bytes, std::memory_order_relaxed);
}

// ─────────────────────────────────────────────────────────────────────────────
// MemBlock — RAII owner of a single raw allocation
// ─────────────────────────────────────────────────────────────────────────────

/// RAII owner of a single malloc/kokkos_malloc allocation. Move-only.
template <class MemSpace>
struct MemBlock {
  char* ptr           = nullptr;

  MemBlock() noexcept = default;

  explicit MemBlock(std::size_t bytes) noexcept {
    // aligned_alloc requires a size that is a multiple of the alignment
    const std::size_t asize = detail::align_up(bytes > 0 ? bytes : 1, detail::kAllocAlign);
    if (asize < bytes) Kokkos::abort("MemBlock: size overflow");
    if constexpr (detail::is_host_space_v<MemSpace>)
      ptr = static_cast<char*>(std::aligned_alloc(detail::kAllocAlign, asize));
    else
      ptr = static_cast<char*>(Kokkos::kokkos_malloc<MemSpace>(asize));
    if (!ptr) Kokkos::abort("MemBlock: allocation failed");
    if (reinterpret_cast<std::uintptr_t>(ptr) % detail::kAllocAlign != 0)
      Kokkos::abort("MemBlock: allocator returned insufficient alignment");
  }

  MemBlock(const MemBlock&)            = delete;
  MemBlock& operator=(const MemBlock&) = delete;

  MemBlock(MemBlock&& o) noexcept : ptr(o.ptr) { o.ptr = nullptr; }

  MemBlock& operator=(MemBlock&& o) noexcept {
    free_ptr();
    ptr   = o.ptr;
    o.ptr = nullptr;
    return *this;
  }

  ~MemBlock() noexcept { free_ptr(); }

 private:
  void free_ptr() noexcept {
    if (!ptr) return;
    if constexpr (detail::is_host_space_v<MemSpace>)
      std::free(ptr);
    else
      Kokkos::kokkos_free<MemSpace>(ptr);
  }
};

// ─────────────────────────────────────────────────────────────────────────────
// OverflowAllocator — fallback for allocations that exceed the slab
// ─────────────────────────────────────────────────────────────────────────────

/// Fallback allocator when the slab is exhausted. Holds up to kMaxBlocks MemBlocks.
template <class MemSpace>
class OverflowAllocator {
 public:
  static constexpr int kMaxBlocks                        = 16;  // increase as needed

  OverflowAllocator() noexcept                           = default;

  OverflowAllocator(const OverflowAllocator&)            = delete;
  OverflowAllocator(OverflowAllocator&&)                 = delete;
  OverflowAllocator& operator=(const OverflowAllocator&) = delete;
  OverflowAllocator& operator=(OverflowAllocator&&)      = delete;

  // unlike the rewind-only frame destructor, this actually frees memory a still-running
  // kernel may touch; not every backend's free synchronizes on its own (SYCL's does not).
  // Default-constructing MemSpace::execution_space here would fence the wrong stream on a
  // device backend: ragnarok kernels dispatch through icon::default_execution_space_instance()
  // (see support/common/execution_space.hpp), a distinct instance from any default-constructed one.
  ~OverflowAllocator() {
    if constexpr (!detail::is_host_space_v<MemSpace>) {
      if (n_ > 0) {
        static_assert(std::is_same_v<MemSpace, typename Kokkos::DefaultExecutionSpace::memory_space>,
                      "OverflowAllocator: MemSpace must be the default execution space's memory space");
        icon::default_execution_space_instance().fence("ragnarok::OverflowAllocator free");
      }
    }
  }

  char* alloc(std::size_t bytes) noexcept {
    if (n_ >= kMaxBlocks) Kokkos::abort("OverflowAllocator: too many overflow allocs — increase kMaxBlocks");
    blocks_[n_] = MemBlock<MemSpace>{bytes};
    return blocks_[n_++].ptr;
  }

 private:
  MemBlock<MemSpace> blocks_[kMaxBlocks];
  int n_ = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// Slab — plain data: pointer, capacity, offset, peak
// ─────────────────────────────────────────────────────────────────────────────

/// Slab: a MemBlock plus current offset and peak (largest offset seen).
template <class MemSpace>
struct Slab {
  MemBlock<MemSpace> block;
  std::size_t size   = 0;
  std::size_t offset = 0;
  std::size_t peak   = 0;
  int live_frames    = 0;  // ScratchSpaces currently alive on this slab

  Slab() noexcept    = default;

  explicit Slab(std::size_t bytes) noexcept : block(bytes), size(bytes) {}

  Slab(const Slab&)            = delete;
  Slab& operator=(const Slab&) = delete;

  // defaulted moves would leave the moved-from scalars stale behind a null block pointer
  Slab(Slab&& o) noexcept
      : block(std::move(o.block)), size(o.size), offset(o.offset), peak(o.peak), live_frames(o.live_frames) {
    o.reset_scalars();
  }

  Slab& operator=(Slab&& o) noexcept {
    block       = std::move(o.block);
    size        = o.size;
    offset      = o.offset;
    peak        = o.peak;
    live_frames = o.live_frames;
    o.reset_scalars();
    return *this;
  }

 private:
  void reset_scalars() noexcept {
    size        = 0;
    offset      = 0;
    peak        = 0;
    live_frames = 0;
  }
};

// ─────────────────────────────────────────────────────────────────────────────
// ScratchSpace — frame over a Slab
// ─────────────────────────────────────────────────────────────────────────────

/// Frame over a Slab: rewinds the slab offset on destruction.
/// Use as a local variable only, so frames die in reverse creation order and the offset
/// stays consistent. Copy/move are deleted to enforce this; live_frames catches the rest.
template <class MemSpace>
class ScratchSpace {
 public:
  explicit ScratchSpace(Slab<MemSpace>& buf) noexcept : frame_base_(buf.offset), buf_(buf) { ++buf.live_frames; }

  ScratchSpace(const ScratchSpace&)            = delete;
  ScratchSpace(ScratchSpace&&)                 = delete;
  ScratchSpace& operator=(const ScratchSpace&) = delete;
  ScratchSpace& operator=(ScratchSpace&&)      = delete;

  ~ScratchSpace() {
    buf_.offset = frame_base_;
    --buf_.live_frames;
    // overflow_ destructs here, freeing any overflow blocks
  }

  /// Advances the slab offset; falls back to a fresh overflow block when the slab is full.
  /// Every allocation starts kAllocAlign-aligned (slab base and overflow blocks are too).
  template <typename U>
  [[nodiscard]] U* alloc(std::size_t count) noexcept {
    static_assert(alignof(U) <= detail::kAllocAlign, "ScratchSpace: element type over-aligned beyond kAllocAlign");
    // with count bounded, sizeof(U)*count cannot wrap; nor can the offset, which only ever
    // advances in step with physically backed allocations (slab or overflow blocks)
    if (count > detail::kMaxScratchElements) Kokkos::abort("ScratchSpace: more than 2^31-1 elements requested");
    const std::size_t bytes = sizeof(U) * count;
    const std::size_t off   = detail::align_up(buf_.offset, detail::kAllocAlign);
    buf_.offset             = off + bytes;
#ifdef MEMORY_POOL_DEBUG
    if (buf_.offset > buf_.peak)
      std::fprintf(stderr, "[ragnarok::ScratchSpace] watermark: %zu -> %zu bytes\n", buf_.peak, buf_.offset);
#endif
    buf_.peak = std::max(buf_.peak, buf_.offset);
    if (buf_.offset <= buf_.size) return reinterpret_cast<U*>(buf_.block.ptr + off);
    return reinterpret_cast<U*>(overflow_.alloc(bytes));
  }

 private:
  std::size_t frame_base_;
  Slab<MemSpace>& buf_;
  OverflowAllocator<MemSpace> overflow_;
};

// ─────────────────────────────────────────────────────────────────────────────
// MemoryPool — factory that hands out ScratchSpaces
// ─────────────────────────────────────────────────────────────────────────────

/// Hands out frames. Host: per-thread slab (safe in OpenMP). Device: one global slab.
template <class ExecSpace>
struct MemoryPool {
  using MemSpace = typename ExecSpace::memory_space;

  /// Call once per timestep. Aborts if a frame was not released (offset != 0).
  /// Also applies any pending slab growth from the previous cycle: doing the (synchronizing)
  /// device realloc here keeps it at the controlled timestep boundary. get_local_space() has
  /// the same growth logic because host worker threads never pass through this function.
  static void ensure_empty_memory_pool() noexcept {
    Slab<MemSpace>& slab = get_slab();
    if (slab.offset != 0) Kokkos::abort("MemoryPool: offset != 0 at timestep start — a frame was not released");
    if (slab.live_frames != 0) Kokkos::abort("MemoryPool: a frame is still alive at timestep start");
    if (slab.size != 0 && slab.peak > slab.size) grow_slab(slab);
  }

  /// Returns a frame starting at the current offset.
  /// When no frame is open (offset == 0), grows the slab if the previous peak exceeded capacity.
  [[nodiscard]] static ScratchSpace<MemSpace> get_local_space() noexcept {
#ifdef _OPENMP
    if constexpr (!detail::is_host_space_v<MemSpace>)
      if (omp_in_parallel()) Kokkos::abort("MemoryPool<device>: called inside OpenMP parallel region");
#endif
    Slab<MemSpace>& slab = get_slab();
    // no frame alive but offset nonzero: a previous frame on this thread leaked.
    // Complements ensure_empty_memory_pool(), which cannot see other threads' slabs.
    if (slab.live_frames == 0 && slab.offset != 0)
      Kokkos::abort("MemoryPool: offset != 0 with no live ScratchSpace — a previous frame leaked");
    // grow/init only when no frame is alive: replacing the slab would free the block under
    // a live frame and clobber its live_frames count. Growth just waits for a clean entry.
    if (slab.offset == 0 && slab.live_frames == 0) {
      if (slab.size == 0) {
        slab = Slab<MemSpace>{detail::scratch_initial_size()};
#ifdef MEMORY_POOL_DEBUG
        std::fprintf(stderr, "[ragnarok::MemoryPool] initial alloc: %zu bytes\n", slab.size);
#endif
        if constexpr (!detail::is_host_space_v<MemSpace>) Kokkos::push_finalize_hook([]() { get_slab() = {}; });
      } else if (slab.peak > slab.size) {
        grow_slab(slab);
      }
    }
    return ScratchSpace<MemSpace>{slab};
  }

 private:
  /// Per-thread (host) or global (device) slab reference.
  static Slab<MemSpace>& get_slab() noexcept { return detail::pool_slab<MemSpace>(); }

  static void grow_slab(Slab<MemSpace>& slab) noexcept {
    // 25% headroom: an exact fit would trigger a realloc every cycle when the peak creeps up
    const std::size_t new_size = slab.peak + slab.peak / 4;
#ifdef MEMORY_POOL_DEBUG
    std::fprintf(stderr, "[ragnarok::MemoryPool] slab grew: %zu -> %zu bytes\n", slab.size, new_size);
#endif
    slab = Slab<MemSpace>{new_size};  // move-assigns: frees old block, allocs new one
  }
};

}  // namespace ragnarok

#endif  // RAGNAROK_SUPPORT_MEMORY_MEMORY_POOL_H_
