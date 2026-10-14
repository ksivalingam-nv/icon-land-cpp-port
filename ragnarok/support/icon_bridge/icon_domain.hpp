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

///
/// @file
/// @brief top level t_patch initialization
///

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_DOMAIN_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_DOMAIN_H_

#include <Kokkos_Core.hpp>
#include <Kokkos_OffsetView.hpp>

#include "common/entity.hpp"
#include "common/refin_ctrl.hpp"
#include "domain_views.hpp"
#include "icon_f2c.hpp"
#include "support/kernels/grid_range.hpp"

namespace icon {

// see src/grid/mo_model_domain.f90 for documentation about Patch components

class CommPatch : public f2c::CommPatch {
 public:
  CommPatch() = delete;
  CommPatch(const f2c::PatchDescr f2c_patch_descr);
};

// The four refin_ctrl marker arrays of a grid entity, grouped so GridBase takes one
// named bundle instead of four same-typed views (which are easy to transpose positionally).
struct RefinCtrlBounds {
  ManagedHostConstOffsetView1D<int> start_block;
  ManagedHostConstOffsetView1D<int> end_block;
  ManagedHostConstOffsetView1D<int> start_index;
  ManagedHostConstOffsetView1D<int> end_index;
};

template <Entity Ent>
class GridBase {
 public:
  // grid-type specific constant:
  static constexpr Entity ent = Ent;
  // block size (nproma) of the blocked field views
  const int block_size;

  // start/end markers of index regions of index fields
  // support: refin_ctrl levels = {min_rl..max_rl of cells}
  // preserving the model semantics of prognostic/halo/... refin_ctrl levels
  // using OffsetViews to handle negative integers
  // example: Fortran:  patch%cells%start_block(-8) == jbs
  //          Ragnarok: patch.cells.start_block(-8) == jbs - 1;
  const ManagedHostConstOffsetView1D<int> start_block;
  const ManagedHostConstOffsetView1D<int> end_block;
  const ManagedHostConstOffsetView1D<int> start_index;
  const ManagedHostConstOffsetView1D<int> end_index;
  //
  GridBase(int block_size, const RefinCtrlBounds& bounds)
      : block_size(block_size),
        start_block(bounds.start_block),
        end_block(bounds.end_block),
        start_index(bounds.start_index),
        end_index(bounds.end_index) {}
  // declare grid as invariant and unique: built once, then never copied or moved.
  GridBase(const GridBase&)            = delete;
  GridBase(GridBase&&)                 = delete;
  GridBase& operator=(const GridBase&) = delete;
  GridBase& operator=(GridBase&&)      = delete;

  [[nodiscard]] GridRange make_grid_range(int rls, int rle = GridRL<ent>::min_rl) const noexcept {
    debug_require(GridRL<ent>::min_rl <= rls && rls <= GridRL<ent>::max_rl, "make_grid_range: rl_start out of range");
    debug_require(GridRL<ent>::min_rl <= rle && rle <= GridRL<ent>::max_rl, "make_grid_range: rl_end out of range");
    return GridRange(ent, block_size, start_block(rls), end_block(rle), start_index(rls), end_index(rle));
  }
};

// grid cells data for ragnarok framework
class GridCells : public GridBase<Entity::cells> {
 public:
  // predefined ranges:
  const GridRange prog_range = make_grid_range(GridRL<ent>::start_prog, GridRL<ent>::end_prog);
  // zero based data fields:
  // Fortran arrays have shape (nproma, nblks_c);
  // Ragnarok uses LayoutRight views with reversed extents:
  // view(jb-1, jc-1) is the same as Fortran_array(jc,jb)
  const HostConstView2D<int> num_edges;
  const HostConstView2D<Real> area;
  const HostConstView2D<GeographicCoords> center;
  const ManagedHostConstView2D<Real> center_lon;
  const ManagedHostConstView2D<Real> center_lat;
  //
  // zero-based index fields:
  // The extents are reversed with respect to Fortran (see data fields above).
  // example: Fortran:  patch%cells%edge_idx(ia, ib, m) == iedge
  //          Ragnarok: patch.cells.edge_idx(m-1, ib-1, ia-1) == iedge - 1;
  // => ragnarok_edge_value(ragnarok_edge_idx(m-1,ib-1,ia-1)) == fortran_edge_value(fortran_edge_idx(ia, ib, m))
  const ManagedHostConstView3D<int> edge_idx;
  const ManagedHostConstView3D<int> edge_blk;
  //
  GridCells() = delete;
  GridCells(const f2c::GridCells& f2c_cells);
  GridCells(const GridCells&)            = delete;
  GridCells(GridCells&&)                 = delete;
  GridCells& operator=(const GridCells&) = delete;
  GridCells& operator=(GridCells&&)      = delete;
};

// grid edges data for ragnarok framework
class GridEdges : public GridBase<Entity::edges> {
 public:
  // predefined ranges:
  const GridRange prog_range = make_grid_range(GridRL<ent>::start_prog, GridRL<ent>::end_prog);
  // zero-based data fields:
  const HostConstView3D<TangentVec> primal_normal_cell;
  const ManagedHostConstView3D<Real> primal_normal_cell_v1;
  const ManagedHostConstView3D<Real> primal_normal_cell_v2;
  const HostConstView3D<TangentVec> primal_normal_vert;
  const ManagedHostConstView3D<Real> primal_normal_vert_v1;
  const ManagedHostConstView3D<Real> primal_normal_vert_v2;
  const HostConstView3D<TangentVec> dual_normal_cell;
  const ManagedHostConstView3D<Real> dual_normal_cell_v1;
  const ManagedHostConstView3D<Real> dual_normal_cell_v2;
  const HostConstView3D<TangentVec> dual_normal_vert;
  const ManagedHostConstView3D<Real> dual_normal_vert_v1;
  const ManagedHostConstView3D<Real> dual_normal_vert_v2;
  const HostConstView2D<Real> tangent_orientation;
  const HostConstView2D<Real> inv_primal_edge_length;
  const HostConstView2D<Real> inv_dual_edge_length;
  const HostConstView2D<Real> inv_vert_vert_length;
  const HostConstView3D<Real> edge_cell_length;
  //
  // zero-based index fields:
  const ManagedHostConstView3D<int> cell_idx;
  const ManagedHostConstView3D<int> cell_blk;
  const ManagedHostConstView3D<int> vertex_idx;
  const ManagedHostConstView3D<int> vertex_blk;
  GridEdges() = delete;
  GridEdges(const f2c::GridEdges& f2c_edges);
  GridEdges(const GridEdges&)            = delete;
  GridEdges(GridEdges&&)                 = delete;
  GridEdges& operator=(const GridEdges&) = delete;
  GridEdges& operator=(GridEdges&&)      = delete;
};

// grid vertices data for ragnarok framework
class GridVertices : public GridBase<Entity::vertices> {
 public:
  // no prog_range for vertices
  // zero-based data fields:
  const HostConstView2D<GeographicCoords> vertex;
  const ManagedHostConstView2D<Real> vertex_lon;
  const ManagedHostConstView2D<Real> vertex_lat;
  //
  // zero-based index fields:
  const ManagedHostConstView3D<int> cell_idx;
  const ManagedHostConstView3D<int> cell_blk;
  GridVertices() = delete;
  GridVertices(const f2c::GridVertices& f2c_vertices);
  GridVertices(const GridVertices&)            = delete;
  GridVertices(GridVertices&&)                 = delete;
  GridVertices& operator=(const GridVertices&) = delete;
  GridVertices& operator=(GridVertices&&)      = delete;
};

class Patch : public f2c::PatchInfo {
 public:
  CommPatch comm_patch;
  const GridCells cells;
  const GridEdges edges;
  const GridVertices verts;
  Patch() = delete;
  Patch(const f2c::PatchDescr f2c_patch_descr);
  // A Patch is built once in domain::init() and then lives unchanged for the whole
  // execution. It is never copied or moved.
  Patch(const Patch&)            = delete;
  Patch(Patch&&)                 = delete;
  Patch& operator=(const Patch&) = delete;
  Patch& operator=(Patch&&)      = delete;
};

// an accidental copy/move ability is meant to be a compile error.
static_assert(not std::is_copy_constructible_v<Patch>, "Patch must NOT be copy constructible!");
static_assert(not std::is_move_constructible_v<Patch>, "Patch must NOT be move constructible!");

// Host-only orchestration: the registry mirrors ICON grid metadata in HostSpace.
// A device consumer must deep_copy into its own mirrors -- Kokkos flags a HostSpace
// view dereferenced in a kernel only at runtime (Kokkos::abort), not at compile time.
namespace domain {

void init();

}  // namespace domain

const Patch& get_patch(int dom_id);

}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_DOMAIN_H_
