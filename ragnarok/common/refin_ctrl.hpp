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

#ifndef RAGNAROK_COMMON_REFIN_CTRL_H_
#define RAGNAROK_COMMON_REFIN_CTRL_H_

#include "entity.hpp"

namespace icon {

namespace grid_constants {

// Maximum halo width
// see mo_impl_constants.f90.
inline constexpr int max_hw = 2;

// Refinement control (refin_ctrl) level bounds per grid entity.
// see mo_impl_constants.f90
template <Entity Ent>
struct GridRL;
template <>
struct GridRL<Entity::cells> {
  static constexpr int grf_bdywidth    = 4;
  static constexpr int min_rl_int      = -4;
  static constexpr int min_rl          = min_rl_int - 2 * max_hw;
  static constexpr int max_rl          = 5;
  static constexpr int start_prog      = grf_bdywidth + 1;
  static constexpr int end_prog        = min_rl_int;
  static constexpr int start_halo      = min_rl_int - 1;
  static constexpr int end_halo        = min_rl;
  static constexpr int start_halo_lev1 = min_rl_int - 1;
  static constexpr int end_halo_lev1   = min_rl_int - 1;
};

template <>
struct GridRL<Entity::edges> {
 private:
  using Cells = GridRL<Entity::cells>;

 public:
  static constexpr int grf_bdywidth = 9;
  static constexpr int min_rl_int   = 2 * Cells::min_rl_int;
  static constexpr int min_rl       = min_rl_int - (2 * max_hw + 1);
  static constexpr int max_rl       = 2 * Cells::max_rl;
  static constexpr int start_prog   = grf_bdywidth + 1;
  static constexpr int end_prog     = min_rl_int;
  static constexpr int start_halo   = min_rl_int - 1;
  static constexpr int end_halo     = min_rl;
};

template <>
struct GridRL<Entity::vertices> {
 private:
  using Cells = GridRL<Entity::cells>;

 public:
  static constexpr int min_rl_int = Cells::min_rl_int;
  static constexpr int min_rl     = min_rl_int - (max_hw + 1);
  static constexpr int max_rl     = Cells::max_rl;
  // - start_prog and grf_bdywidth_v are not defined in ICON
  // - vertex loops choose rl_start per stencil, typically 2
  //   (skipping the outermost vertex row on the lateral boundary)
  // - end_prog mirrors mo_impl_constants::end_prog_verts (separation from halos)
  static constexpr int end_prog   = min_rl_int;
  static constexpr int start_halo = min_rl_int - 1;
  static constexpr int end_halo   = min_rl;
};

// Lock the values to their Fortran counterparts. These mirror the documented
// results in src/shared/mo_impl_constants.f90 and must stay in lockstep.
static_assert(GridRL<Entity::cells>::min_rl == -8);
static_assert(GridRL<Entity::edges>::min_rl == -13);
static_assert(GridRL<Entity::edges>::min_rl_int == -8);
static_assert(GridRL<Entity::edges>::max_rl == 10);
static_assert(GridRL<Entity::vertices>::min_rl == -7);

// General ordering invariant: halo bound below interior bound below max.
static_assert(GridRL<Entity::cells>::min_rl < GridRL<Entity::cells>::min_rl_int &&
              GridRL<Entity::cells>::min_rl_int < GridRL<Entity::cells>::max_rl);
static_assert(GridRL<Entity::edges>::min_rl < GridRL<Entity::edges>::min_rl_int &&
              GridRL<Entity::edges>::min_rl_int < GridRL<Entity::edges>::max_rl);
static_assert(GridRL<Entity::vertices>::min_rl < GridRL<Entity::vertices>::min_rl_int &&
              GridRL<Entity::vertices>::min_rl_int < GridRL<Entity::vertices>::max_rl);

}  // namespace grid_constants

// Public alias for grid refin_ctrl levels:
template <Entity Ent>
using GridRL = grid_constants::GridRL<Ent>;

}  // namespace icon

#endif  // RAGNAROK_COMMON_REFIN_CTRL_H_
