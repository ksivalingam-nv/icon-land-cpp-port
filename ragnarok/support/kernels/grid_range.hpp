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

#ifndef RAGNAROK_SUPPORT_KERNELS_GRID_RANGE_H_
#define RAGNAROK_SUPPORT_KERNELS_GRID_RANGE_H_

#include <Kokkos_Core.hpp>

#include "common/entity.hpp"

namespace icon {

// construction tag: skip bounds normalization
struct Unchecked {
  explicit Unchecked() = default;
};
inline constexpr Unchecked unchecked{};

class Range {
 public:
  const int start;  // inclusive start index
  const int end;    // exclusive end index

  // normalizes to a valid - possibly empty - [start, end) range;
  // accepts arbitrary start > end as empty input
  constexpr Range(int start, int end) : Range(unchecked, start, Kokkos::max(start, end)) {}

  // caller guarantees start <= end
  constexpr Range(Unchecked, int start, int end) : start(start), end(end) {}

  [[nodiscard]] constexpr int size() const noexcept { return end - start; }

  template <typename ExecSpace>
  [[nodiscard]] auto range_policy(const ExecSpace& space) const {
    return Kokkos::RangePolicy(space, start, end);
  }
};

// dimension-tagged Range flavors: plain Ranges distinguished at the type
// level, so interfaces can tell the loop dimensions apart and an argument
// cannot slip into the wrong slot silently

// vertical levels
struct LevelRange : Range {
  using Range::Range;
};

// horizontal indices (within a block)
struct IndexRange : Range {
  using Range::Range;
};

// blocks of the horizontal grid
struct BlockRange : Range {
  using Range::Range;
};

// iteration support on grids
// similar to mo_model_domain::t_patch::t_subset_range
class GridRange {
 private:
  const Entity ent;
  const int blk_size;
  // blocks are stored inclusive - they act as identities (jb == end_blk)
  const int start_blk;
  const int end_blk;
  const int start_blk_start_idx;
  // exclusive upper index bound of the end block
  const int end_blk_end_idx;

 public:
  // - constructor arguments are all inclusive bounds - the form in which
  //   bounds arrive at the bridge; the inclusive->exclusive
  //   conversion happens once, here
  // - bounds are normalized such that block_range() and index_range()
  //   always form valid - possibly empty - [start, end) ranges:
  //   arbitrary start > end is accepted as empty input, index trims are
  //   clamped to [0, block_size]
  constexpr GridRange(Entity ent, int block_size, int start_block, int end_block, int start_index, int end_index)
      : ent(ent),
        blk_size(block_size),
        start_blk(start_block),
        end_blk(Kokkos::max(end_block, start_block - 1)),
        start_blk_start_idx(Kokkos::min(blk_size, Kokkos::max(0, start_index))),
        end_blk_end_idx(
            Kokkos::min(blk_size, Kokkos::max(end_index + 1, (start_blk == end_blk) ? start_blk_start_idx : 0))) {}

  // construct from the flat range: the iteration space in flattened
  // block-major coordinates; requires flat.start >= 0
  constexpr GridRange(Entity ent, int block_size, Range flat)
      : GridRange(ent, block_size, flat.start / block_size, (flat.end - 1) / block_size, flat.start % block_size,
                  (flat.end - 1) % block_size) {}

  // the blocks of this range
  [[nodiscard]] constexpr BlockRange block_range() const noexcept {
    return BlockRange(unchecked, start_blk, end_blk + 1);
  }

  // indices within block jb
  [[nodiscard]] constexpr IndexRange index_range(int jb) const noexcept {
    const int start = (jb == start_blk) ? start_blk_start_idx : 0;
    const int end   = (jb == end_blk) ? end_blk_end_idx : blk_size;
    return IndexRange(unchecked, start, end);
  }

  // the whole iteration space as one flat range (see the flat constructor);
  // empty grid ranges yield [0, 0) - their stored index trims are unspecified
  [[nodiscard]] constexpr Range get_flat_range() const noexcept {
    if (start_blk > end_blk) return Range(unchecked, 0, 0);
    return Range(unchecked, start_blk * blk_size + start_blk_start_idx, end_blk * blk_size + end_blk_end_idx);
  }

  // state infos:

  [[nodiscard]] constexpr Entity get_entity() const noexcept { return ent; }

  [[nodiscard]] constexpr int get_block_size() const noexcept { return blk_size; }
};

}  // end of namespace icon
#endif
