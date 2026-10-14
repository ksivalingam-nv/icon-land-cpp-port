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
#include <algorithm>
#include <vector>

#include "common/types.hpp"
#include "icon_loops.hpp"
#include "kernel_sync.hpp"

namespace {

using icon::GridRange;

// a small grid range with a trimmed first and last block and a full block in
// between; inclusive constructor bounds
constexpr int block_size = 8;
constexpr int nblocks    = 3;
GridRange make_test_range() { return GridRange(icon::Entity::cells, block_size, 0, nblocks - 1, 3, 5); }

// kernel dispatch lives in a free function: nvcc rejects extended lambdas
// inside the private TestBody member (see types_test.cpp)
void apply_block_index(const GridRange grid_range, const ConstView2D<double, MemorySpace> d_in,
                       View2D<double, MemorySpace> d_out) {
  icon::block_loop(
      "kernels_block_loop", grid_range, block_lambda(const int jb) {
        icon::index_loop(
            "kernels_index_loop", grid_range, jb,
            index_lambda(const int jc) { d_out(jb, jc) = 2.0 * d_in(jb, jc) + jc; });
      });
  icon::kernel_fence();
}

// the block/index role swap over a trimmed GridRange: the same elemental
// expression through block_loop/index_loop must reproduce the plain nested
// reference loops bitwise - and must not touch anything outside the range
// (checked via poisoned output)
TEST(KernelsTest, BlockIndexLoop) {
  const GridRange grid_range = make_test_range();
  constexpr double poison    = -1.0;

  ManagedHostView2D<double> h_in("h_in", nblocks, block_size);
  ManagedHostView2D<double> h_ref("h_ref", nblocks, block_size);
  for (int jb = 0; jb < nblocks; ++jb)
    for (int jc = 0; jc < block_size; ++jc) {
      h_in(jb, jc)  = jb * block_size + jc;
      h_ref(jb, jc) = poison;
    }

  // reference: plain nested loops over the range accessors
  const auto brange = grid_range.block_range();
  for (int jb = brange.start; jb < brange.end; ++jb) {
    const auto irange = grid_range.index_range(jb);
    for (int jc = irange.start; jc < irange.end; ++jc) h_ref(jb, jc) = 2.0 * h_in(jb, jc) + jc;
  }

  // candidate: the same expression through the loop abstraction, staged in
  // the default memory space
  ManagedView2D<double> d_in(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_in"), nblocks, block_size);
  ManagedView2D<double> d_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_out"), nblocks, block_size);
  Kokkos::deep_copy(d_in, h_in);
  Kokkos::deep_copy(d_out, poison);
  apply_block_index(grid_range, d_in, d_out);
  ManagedHostView2D<double> h_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "h_out"), nblocks, block_size);
  Kokkos::deep_copy(h_out, d_out);

  // full-array comparison: poison outside the range, exact values inside
  for (int jb = 0; jb < nblocks; ++jb)
    for (int jc = 0; jc < block_size; ++jc) EXPECT_EQ(h_out(jb, jc), h_ref(jb, jc)) << "jb=" << jb << " jc=" << jc;
}

void apply_collapsed(const GridRange grid_range, const icon::LevelRange levels,
                     const ConstView3D<double, MemorySpace> d_in, View3D<double, MemorySpace> d_out) {
  icon::block_loop(
      "kernels_collapsed_blocks", grid_range, block_lambda(const int jb) {
        icon::index_loop(
            "kernels_collapsed", grid_range, jb, levels,
            index_lambda(const int jk, const int jc) { d_out(jb, jk, jc) = d_in(jb, jk, jc) + 10.0 * jk + jc; });
      });
  icon::kernel_fence();
}

// collapsed level+index loop with a nonzero level start: a level-dependent
// elemental expression must match the reference exactly - a wrong (jk, jc) /
// stride mapping or level offset shows up in the values, an out-of-range
// write in the poison
TEST(KernelsTest, CollapsedBlockIndexLoop) {
  const GridRange grid_range = make_test_range();
  constexpr int nlev         = 5;
  const icon::LevelRange levels(1, nlev);  // skip level 0
  constexpr double poison = -1.0;

  ManagedHostView3D<double> h_in("h_in", nblocks, nlev, block_size);
  ManagedHostView3D<double> h_ref("h_ref", nblocks, nlev, block_size);
  for (int jb = 0; jb < nblocks; ++jb)
    for (int jk = 0; jk < nlev; ++jk)
      for (int jc = 0; jc < block_size; ++jc) {
        h_in(jb, jk, jc)  = (jb * nlev + jk) * block_size + jc;
        h_ref(jb, jk, jc) = poison;
      }

  const auto brange = grid_range.block_range();
  for (int jb = brange.start; jb < brange.end; ++jb) {
    const auto irange = grid_range.index_range(jb);
    for (int jk = levels.start; jk < levels.end; ++jk)
      for (int jc = irange.start; jc < irange.end; ++jc) h_ref(jb, jk, jc) = h_in(jb, jk, jc) + 10.0 * jk + jc;
  }

  ManagedView3D<double> d_in(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_in"), nblocks, nlev, block_size);
  ManagedView3D<double> d_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_out"), nblocks, nlev, block_size);
  Kokkos::deep_copy(d_in, h_in);
  Kokkos::deep_copy(d_out, poison);
  apply_collapsed(grid_range, levels, d_in, d_out);
  ManagedHostView3D<double> h_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "h_out"), nblocks, nlev, block_size);
  Kokkos::deep_copy(h_out, d_out);

  for (int jb = 0; jb < nblocks; ++jb)
    for (int jk = 0; jk < nlev; ++jk)
      for (int jc = 0; jc < block_size; ++jc)
        EXPECT_EQ(h_out(jb, jk, jc), h_ref(jb, jk, jc)) << "jb=" << jb << " jk=" << jk << " jc=" << jc;
}

void apply_collapsed_flat(const icon::LevelRange levels, const icon::IndexRange indices,
                          const ConstView2D<double, MemorySpace> d_in, View2D<double, MemorySpace> d_out) {
  icon::index_loop(
      "kernels_collapsed_flat", levels, indices,
      index_lambda(const int jk, const int jc) { d_out(jk, jc) = d_in(jk, jc) + 10.0 * jk + jc; });
  icon::kernel_fence();
}

// the block-free collapsed form (the saturation_adjustment_2d shape): no
// GridRange, no block loop - completion is the caller's job via
// kernel_fence
TEST(KernelsTest, CollapsedFlatLoop) {
  constexpr int n    = 10;
  constexpr int nlev = 4;
  const icon::LevelRange levels(1, nlev);
  const icon::IndexRange indices(2, 9);
  constexpr double poison = -1.0;

  ManagedHostView2D<double> h_in("h_in", nlev, n);
  ManagedHostView2D<double> h_ref("h_ref", nlev, n);
  for (int jk = 0; jk < nlev; ++jk)
    for (int jc = 0; jc < n; ++jc) {
      h_in(jk, jc)  = jk * n + jc;
      h_ref(jk, jc) = poison;
    }
  for (int jk = levels.start; jk < levels.end; ++jk)
    for (int jc = indices.start; jc < indices.end; ++jc) h_ref(jk, jc) = h_in(jk, jc) + 10.0 * jk + jc;

  ManagedView2D<double> d_in(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_in"), nlev, n);
  ManagedView2D<double> d_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_out"), nlev, n);
  Kokkos::deep_copy(d_in, h_in);
  Kokkos::deep_copy(d_out, poison);
  apply_collapsed_flat(levels, indices, d_in, d_out);
  ManagedHostView2D<double> h_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "h_out"), nlev, n);
  Kokkos::deep_copy(h_out, d_out);

  for (int jk = 0; jk < nlev; ++jk)
    for (int jc = 0; jc < n; ++jc) EXPECT_EQ(h_out(jk, jc), h_ref(jk, jc)) << "jk=" << jk << " jc=" << jc;
}

void apply_collapsed_flat2(const icon::Range extra, const icon::IndexRange indices,
                           const ConstView2D<double, MemorySpace> d_in, View2D<double, MemorySpace> d_out) {
  icon::index_loop(
      "kernels_collapsed_flat2", extra, indices,
      index_lambda(const int i, const int jc) { d_out(i, jc) = d_in(i, jc) + 10.0 * i + jc; });
  icon::kernel_fence();
}

// the grid-free 2d collapsed form (the prognostics-init shape): a plain
// Range as the slow dimension over the indices of a block; nonzero starts
// on both dimensions catch a wrong div/mod recovery, the poison an
// out-of-range write
TEST(KernelsTest, CollapsedFlatLoopGridFree) {
  constexpr int n_extra = 4;
  constexpr int n       = 10;
  const icon::Range extra(1, n_extra);
  const icon::IndexRange indices(2, 9);
  constexpr double poison = -1.0;

  ManagedHostView2D<double> h_in("h_in", n_extra, n);
  ManagedHostView2D<double> h_ref("h_ref", n_extra, n);
  for (int i = 0; i < n_extra; ++i)
    for (int jc = 0; jc < n; ++jc) {
      h_in(i, jc)  = i * n + jc;
      h_ref(i, jc) = poison;
    }
  for (int i = extra.start; i < extra.end; ++i)
    for (int jc = indices.start; jc < indices.end; ++jc) h_ref(i, jc) = h_in(i, jc) + 10.0 * i + jc;

  ManagedView2D<double> d_in(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_in"), n_extra, n);
  ManagedView2D<double> d_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_out"), n_extra, n);
  Kokkos::deep_copy(d_in, h_in);
  Kokkos::deep_copy(d_out, poison);
  apply_collapsed_flat2(extra, indices, d_in, d_out);
  ManagedHostView2D<double> h_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "h_out"), n_extra, n);
  Kokkos::deep_copy(h_out, d_out);

  for (int i = 0; i < n_extra; ++i)
    for (int jc = 0; jc < n; ++jc) EXPECT_EQ(h_out(i, jc), h_ref(i, jc)) << "i=" << i << " jc=" << jc;
}

void apply_collapsed_flat3(const icon::Range extra, const icon::LevelRange levels, const icon::IndexRange indices,
                           const ConstView3D<double, MemorySpace> d_in, View3D<double, MemorySpace> d_out) {
  icon::index_loop(
      "kernels_collapsed_flat3", extra, levels, indices, index_lambda(const int i, const int jk, const int jc) {
        d_out(i, jk, jc) = d_in(i, jk, jc) + 100.0 * i + 10.0 * jk + jc;
      });
  icon::kernel_fence();
}

// the 3d collapsed block-free form (the tracer-copy shape): a grid-free
// slowest dimension in front of levels x indices; nonzero starts on all
// three dimensions catch a wrong div/mod recovery, the poison an
// out-of-range write
TEST(KernelsTest, CollapsedFlatLoop3D) {
  constexpr int n_extra = 4;
  constexpr int n       = 10;
  constexpr int nlev    = 5;
  const icon::Range extra(1, n_extra);
  const icon::LevelRange levels(1, nlev);
  const icon::IndexRange indices(2, 9);
  constexpr double poison = -1.0;

  ManagedHostView3D<double> h_in("h_in", n_extra, nlev, n);
  ManagedHostView3D<double> h_ref("h_ref", n_extra, nlev, n);
  for (int i = 0; i < n_extra; ++i)
    for (int jk = 0; jk < nlev; ++jk)
      for (int jc = 0; jc < n; ++jc) {
        h_in(i, jk, jc)  = (i * nlev + jk) * n + jc;
        h_ref(i, jk, jc) = poison;
      }
  for (int i = extra.start; i < extra.end; ++i)
    for (int jk = levels.start; jk < levels.end; ++jk)
      for (int jc = indices.start; jc < indices.end; ++jc)
        h_ref(i, jk, jc) = h_in(i, jk, jc) + 100.0 * i + 10.0 * jk + jc;

  ManagedView3D<double> d_in(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_in"), n_extra, nlev, n);
  ManagedView3D<double> d_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_out"), n_extra, nlev, n);
  Kokkos::deep_copy(d_in, h_in);
  Kokkos::deep_copy(d_out, poison);
  apply_collapsed_flat3(extra, levels, indices, d_in, d_out);
  ManagedHostView3D<double> h_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "h_out"), n_extra, nlev, n);
  Kokkos::deep_copy(h_out, d_out);

  for (int i = 0; i < n_extra; ++i)
    for (int jk = 0; jk < nlev; ++jk)
      for (int jc = 0; jc < n; ++jc)
        EXPECT_EQ(h_out(i, jk, jc), h_ref(i, jk, jc)) << "i=" << i << " jk=" << jk << " jc=" << jc;
}

// per-block count-and-compact via an atomic counter, fetch_scalar readback,
// and a counter-sized follow-up kernel - the graupel jmx pattern; the block
// loop lives at the caller (a plain host loop, as if on the Fortran side)
int compact_evens(const GridRange grid_range, const int jb, const ConstView2D<double, MemorySpace> d_in,
                  View1D<double, MemorySpace> d_sel) {
  Kokkos::View<int> d_count("count");
  icon::index_loop(
      "kernels_count", grid_range, jb, index_lambda(const int jc) {
        if (static_cast<int>(d_in(jb, jc)) % 2 == 0) {
          const int pos = Kokkos::atomic_fetch_add(&d_count(), 1);
          d_sel(pos)    = d_in(jb, jc);
        }
      });
  const int n = icon::fetch_scalar(d_count);
  icon::index_loop(
      "kernels_transform", icon::IndexRange(0, n), index_lambda(const int j) { d_sel(j) = d_sel(j) + 1.0; });
  icon::kernel_fence();
  return n;
}

TEST(KernelsTest, FetchScalarReadback) {
  const GridRange grid_range = make_test_range();

  ManagedHostView2D<double> h_in("h_in", nblocks, block_size);
  for (int jb = 0; jb < nblocks; ++jb)
    for (int jc = 0; jc < block_size; ++jc) h_in(jb, jc) = jb * block_size + jc;
  ManagedView2D<double> d_in(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_in"), nblocks, block_size);
  Kokkos::deep_copy(d_in, h_in);

  const auto brange = grid_range.block_range();
  for (int jb = brange.start; jb < brange.end; ++jb) {
    // reference: the even values of block jb, transformed, in ascending order
    std::vector<double> expected;
    const auto irange = grid_range.index_range(jb);
    for (int jc = irange.start; jc < irange.end; ++jc) {
      if (static_cast<int>(h_in(jb, jc)) % 2 == 0) expected.push_back(h_in(jb, jc) + 1.0);
    }

    ManagedView1D<double> d_sel(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_sel"), block_size);
    const int n = compact_evens(grid_range, jb, d_in, d_sel);
    ASSERT_EQ(n, static_cast<int>(expected.size())) << "jb=" << jb;

    ManagedHostView1D<double> h_sel(Kokkos::view_alloc(Kokkos::WithoutInitializing, "h_sel"), block_size);
    Kokkos::deep_copy(h_sel, d_sel);
    // the atomic compaction order is nondeterministic - compare sorted
    std::vector<double> selected(h_sel.data(), h_sel.data() + n);
    std::sort(selected.begin(), selected.end());
    EXPECT_EQ(selected, expected) << "jb=" << jb;
  }
}

void apply_fill(const icon::IndexRange indices, View1D<double, MemorySpace> d_data) {
  icon::index_loop("kernels_fetch_fill", indices, index_lambda(const int j) { d_data(j) = 1.0 + j; });
}

// the mirror staging pattern fetch() serves: an asynchronous index_loop fills
// the view, fetch() through the host mirror completes it and makes the values
// host-readable - the aliasing no-op on host builds, the completing copy on
// device builds (poison proves the copy happened)
TEST(KernelsTest, FetchMirrorStaging) {
  constexpr int n         = 10;
  constexpr double poison = -1.0;

  ManagedView1D<double> d_data(Kokkos::view_alloc(Kokkos::WithoutInitializing, "d_data"), n);
  auto h_data = Kokkos::create_mirror_view(d_data);  // aliases d_data on host builds
  for (int j = 0; j < n; ++j) h_data(j) = poison;

  apply_fill(icon::IndexRange(0, n), d_data);
  icon::fetch(h_data, d_data);

  for (int j = 0; j < n; ++j) EXPECT_EQ(h_data(j), 1.0 + j) << "j=" << j;
}

// running sum over the whole grid_range - loop-carried across jc and jb; the
// ordered whole-range serial_loop takes a plain by-reference lambda, so an
// ordinary scalar carries the sum
void apply_scan(const GridRange grid_range, const HostConstView2D<double> h_in, HostView2D<double> h_out) {
  double sum = 0.0;
  icon::serial_loop(grid_range, [&](const int jb, const int jc) {
    sum += h_in(jb, jc);
    h_out(jb, jc) = sum;
  });
}

// the ordering guarantee of serial_loop: an order-dependent scan must
// reproduce the reference bitwise - any parallelization or reordering breaks
// the non-associative float sums. Host-only views: on device builds this also
// proves serial_loop stays on the host (the data has no device home).
TEST(KernelsTest, SerialScanOrdering) {
  const GridRange grid_range = make_test_range();
  constexpr double poison    = -1.0;

  ManagedHostView2D<double> h_in("h_in", nblocks, block_size);
  ManagedHostView2D<double> h_ref("h_ref", nblocks, block_size);
  ManagedHostView2D<double> h_out("h_out", nblocks, block_size);
  for (int jb = 0; jb < nblocks; ++jb)
    for (int jc = 0; jc < block_size; ++jc) {
      h_in(jb, jc)  = 1.0 / (1.0 + jb * block_size + jc);  // order-sensitive sums
      h_ref(jb, jc) = poison;
      h_out(jb, jc) = poison;
    }

  double sum        = 0.0;
  const auto brange = grid_range.block_range();
  for (int jb = brange.start; jb < brange.end; ++jb) {
    const auto irange = grid_range.index_range(jb);
    for (int jc = irange.start; jc < irange.end; ++jc) {
      sum += h_in(jb, jc);
      h_ref(jb, jc) = sum;
    }
  }

  apply_scan(grid_range, h_in, h_out);

  for (int jb = 0; jb < nblocks; ++jb)
    for (int jc = 0; jc < block_size; ++jc) EXPECT_EQ(h_out(jb, jc), h_ref(jb, jc)) << "jb=" << jb << " jc=" << jc;
}

}  // namespace
