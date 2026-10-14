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
#include "grid_range.hpp"

#include <gtest/gtest.h>

#include <type_traits>

namespace {

using icon::BlockRange;
using icon::Entity;
using icon::GridRange;
using icon::IndexRange;
using icon::LevelRange;
using icon::Range;
using icon::unchecked;

// --- Range: normalization, unchecked passthrough, size

constexpr Range empty_range(5, 2);  // arbitrary start > end is empty input
static_assert(empty_range.start == 5 && empty_range.end == 5 && empty_range.size() == 0);

constexpr Range plain_range(icon::unchecked, 3, 7);
static_assert(plain_range.start == 3 && plain_range.end == 7 && plain_range.size() == 4);

// --- dimension tags: LevelRange/IndexRange/BlockRange are Ranges, but no
// Range (or foreign flavor) can slip into a tagged slot

static_assert(std::is_convertible_v<LevelRange, Range>);
static_assert(std::is_convertible_v<IndexRange, Range>);
static_assert(std::is_convertible_v<BlockRange, Range>);
static_assert(!std::is_convertible_v<Range, IndexRange>);
static_assert(!std::is_convertible_v<Range, LevelRange>);
static_assert(!std::is_convertible_v<Range, BlockRange>);
static_assert(!std::is_convertible_v<LevelRange, IndexRange>);
static_assert(!std::is_convertible_v<IndexRange, LevelRange>);
static_assert(!std::is_convertible_v<BlockRange, IndexRange>);
static_assert(!std::is_convertible_v<IndexRange, BlockRange>);

// --- GridRange: inclusive 6-int constructor, per-block index trims

constexpr int bs = 8;
constexpr GridRange grid6(Entity::cells, bs, 0, 2, 3, 5);
static_assert(grid6.block_range().start == 0 && grid6.block_range().end == 3);
static_assert(grid6.index_range(0).start == 3 && grid6.index_range(0).end == bs);  // trimmed first block
static_assert(grid6.index_range(1).start == 0 && grid6.index_range(1).end == bs);  // full middle block
static_assert(grid6.index_range(2).start == 0 && grid6.index_range(2).end == 6);   // trimmed last block

constexpr GridRange single_block(Entity::cells, bs, 1, 1, 2, 4);
static_assert(single_block.block_range().start == 1 && single_block.block_range().end == 2);
static_assert(single_block.index_range(1).start == 2 && single_block.index_range(1).end == 5);
static_assert(single_block.get_flat_range().start == bs + 2 && single_block.get_flat_range().end == bs + 5);

// normalized empty (arbitrary junk index bounds must not leak out)
constexpr GridRange empty_grid(Entity::cells, bs, 5, 3, 2, 1000);
static_assert(empty_grid.block_range().size() == 0);
static_assert(empty_grid.get_flat_range().start == 0 && empty_grid.get_flat_range().end == 0);

// index trims are clamped to [0, block_size]
constexpr GridRange clamped(Entity::cells, bs, 0, 1, -7, 1000);
static_assert(clamped.index_range(0).start == 0 && clamped.index_range(1).end == bs);
constexpr GridRange clamped_single(Entity::cells, bs, 2, 2, 3, 999);
static_assert(clamped_single.index_range(2).start == 3 && clamped_single.index_range(2).end == bs);

// --- flat constructor <-> get_flat_range round-trip

constexpr int fbs = 48;

// non-empty: exact round-trip
constexpr GridRange from_flat(Entity::cells, fbs, Range(5, 100));
static_assert(from_flat.get_flat_range().start == 5 && from_flat.get_flat_range().end == 100);

// single element
constexpr GridRange from_flat_one(Entity::cells, fbs, Range(50, 51));
static_assert(from_flat_one.get_flat_range().start == 50 && from_flat_one.get_flat_range().end == 51);
static_assert(from_flat_one.index_range(1).start == 2 && from_flat_one.index_range(1).end == 3);

// in-block empty: exact round-trip
constexpr GridRange from_flat_empty_mid(Entity::cells, fbs, Range(50, 50));
static_assert(from_flat_empty_mid.get_flat_range().start == 50 && from_flat_empty_mid.get_flat_range().size() == 0);

// block-boundary empties: empty in, empty out
constexpr GridRange from_flat_empty_zero(Entity::cells, fbs, Range(0, 0));
static_assert(from_flat_empty_zero.get_flat_range().size() == 0);
constexpr GridRange from_flat_empty_edge(Entity::cells, fbs, Range(fbs, fbs));
static_assert(from_flat_empty_edge.get_flat_range().size() == 0);

// runtime visibility of the compile-time invariants (ctest lists a Kernels
// test; the static_asserts above are the actual coverage)
TEST(KernelsGridRange, FlatRoundtrip) {
  const GridRange g(Entity::cells, fbs, Range(5, 100));
  EXPECT_EQ(g.get_flat_range().start, 5);
  EXPECT_EQ(g.get_flat_range().end, 100);
  EXPECT_EQ(empty_grid.get_flat_range().size(), 0);
}

}  // namespace
