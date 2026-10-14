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

#include "icon_domain.hpp"

#include <memory>
#include <vector>

#include "finish.hpp"
#include "icon_f2c.hpp"
namespace icon {
namespace {

// container for all domains
std::vector<std::unique_ptr<Patch>> p_patch;
int p_patch_id_offset = 0;

// nvhpc-24.7: nvc++ cannot parse explicit template arguments on names found only via
// ADL (C++20 P0846); make the f2c view factories visible to ordinary lookup.
using f2c::make_component_view2d;
using f2c::make_component_view3d;
using f2c::make_constview2d;
using f2c::make_constview3d;
using f2c::make_index_view3d;

// --- index conversion:

// we use a proxy object to avoid unwanted conversions
struct IndexOffsetViewProxy {
  ManagedHostOffsetView1D<int> view;
  operator ManagedHostOffsetView1D<const int>() const { return view; }
};

// The 1D offset view is a one-off and reads the repr directly; the 3D factory
// lives in domain_views.hpp.
IndexOffsetViewProxy make_index_offsetview1d(const f2c::FieldRepr& repr) {
  require(repr.rank == 1);
  // Allocate the OffsetView uninitialized directly; the index range is inclusive,
  // so [lb, lb + shape[0] - 1] gives extent == shape[0] with lower bound lb.
  const int lb = repr.lbounds[0];
  ManagedHostOffsetView1D<int> ov(Kokkos::view_alloc(Kokkos::WithoutInitializing, "index_offsetview1d"),
                                  {lb, lb + repr.shape[0] - 1});
  f2c::transform_fortran_to_c_indices(f2c::get_data_ptr<int>(repr), ov.data(), ov.size());
  return {ov};
}

// Bundle the four refin_ctrl marker arrays of an f2c grid struct into a RefinCtrlBounds.
// The three f2c grid structs share these field names, so one template covers all of them,
// and the designated initializers name each marker exactly once (no positional mix-ups).
template <typename F2cGrid>
RefinCtrlBounds make_refin_ctrl_bounds(const F2cGrid& g) {
  return {.start_block = make_index_offsetview1d(g.start_block_f),
          .end_block   = make_index_offsetview1d(g.end_block_f),
          .start_index = make_index_offsetview1d(g.start_index_f),
          .end_index   = make_index_offsetview1d(g.end_index_f)};
}

// --- f2c data fetch:

// retrieve grid cells from icon:
f2c::GridCells get_f2c_grid_cells(const f2c::PatchDescr f2c_patch_descr) {
  f2c::GridCells f2c_cells;
  f2c::get_fun_table().get_grid_cells(f2c_patch_descr, &f2c_cells);
  return f2c_cells;
}
f2c::GridEdges get_f2c_grid_edges(const f2c::PatchDescr f2c_patch_descr) {
  f2c::GridEdges f2c_edges;
  f2c::get_fun_table().get_grid_edges(f2c_patch_descr, &f2c_edges);
  return f2c_edges;
}
f2c::GridVertices get_f2c_grid_vertices(const f2c::PatchDescr f2c_patch_descr) {
  f2c::GridVertices f2c_vertices;
  f2c::get_fun_table().get_grid_vertices(f2c_patch_descr, &f2c_vertices);
  return f2c_vertices;
}

f2c::CommPatch gen_f2c_comm_patch(const f2c::PatchDescr f2c_patch_descr) {
  f2c::CommPatch f2c_comm_patch;
  f2c::get_fun_table().get_comm_patch(f2c_patch_descr, &f2c_comm_patch);
  return f2c_comm_patch;
}

f2c::PatchInfo gen_patch_info(const f2c::PatchDescr f2c_patch_descr) {
  f2c::PatchInfo patch_info;
  f2c::get_fun_table().get_patch_info(f2c_patch_descr, &patch_info);
  return patch_info;
}

}  // namespace

// --- constructors:

CommPatch::CommPatch(const f2c::PatchDescr f2c_patch_descr) : f2c::CommPatch(gen_f2c_comm_patch(f2c_patch_descr)) {}

GridCells::GridCells(const f2c::GridCells& f2c_cells)
    : GridBase(f2c_cells.edge_idx_f.shape[0], make_refin_ctrl_bounds(f2c_cells)),
      // data fields:
      num_edges(make_constview2d<int>(f2c_cells.num_edges)),
      area(make_constview2d<Real>(f2c_cells.area)),
      center(make_constview2d<GeographicCoords>(f2c_cells.center)),
      center_lon(make_component_view2d<GeographicCoords>(f2c_cells.center, &GeographicCoords::lon)),
      center_lat(make_component_view2d<GeographicCoords>(f2c_cells.center, &GeographicCoords::lat)),
      // index fields:
      edge_idx(make_index_view3d(f2c_cells.edge_idx_f)),
      edge_blk(make_index_view3d(f2c_cells.edge_blk_f)) {}

GridEdges::GridEdges(const f2c::GridEdges& f2c_edges)
    : GridBase(f2c_edges.cell_idx_f.shape[0], make_refin_ctrl_bounds(f2c_edges)),
      // data fields:
      primal_normal_cell(make_constview3d<TangentVec>(f2c_edges.primal_normal_cell)),
      primal_normal_cell_v1(make_component_view3d<TangentVec>(f2c_edges.primal_normal_cell, &TangentVec::v1)),
      primal_normal_cell_v2(make_component_view3d<TangentVec>(f2c_edges.primal_normal_cell, &TangentVec::v2)),
      primal_normal_vert(make_constview3d<TangentVec>(f2c_edges.primal_normal_vert)),
      primal_normal_vert_v1(make_component_view3d<TangentVec>(f2c_edges.primal_normal_vert, &TangentVec::v1)),
      primal_normal_vert_v2(make_component_view3d<TangentVec>(f2c_edges.primal_normal_vert, &TangentVec::v2)),
      dual_normal_cell(make_constview3d<TangentVec>(f2c_edges.dual_normal_cell)),
      dual_normal_cell_v1(make_component_view3d<TangentVec>(f2c_edges.dual_normal_cell, &TangentVec::v1)),
      dual_normal_cell_v2(make_component_view3d<TangentVec>(f2c_edges.dual_normal_cell, &TangentVec::v2)),
      dual_normal_vert(make_constview3d<TangentVec>(f2c_edges.dual_normal_vert)),
      dual_normal_vert_v1(make_component_view3d<TangentVec>(f2c_edges.dual_normal_vert, &TangentVec::v1)),
      dual_normal_vert_v2(make_component_view3d<TangentVec>(f2c_edges.dual_normal_vert, &TangentVec::v2)),
      tangent_orientation(make_constview2d<Real>(f2c_edges.tangent_orientation)),
      inv_primal_edge_length(make_constview2d<Real>(f2c_edges.inv_primal_edge_length)),
      inv_dual_edge_length(make_constview2d<Real>(f2c_edges.inv_dual_edge_length)),
      inv_vert_vert_length(make_constview2d<Real>(f2c_edges.inv_vert_vert_length)),
      edge_cell_length(make_constview3d<Real>(f2c_edges.edge_cell_length)),
      // index fields:
      cell_idx(make_index_view3d(f2c_edges.cell_idx_f)),
      cell_blk(make_index_view3d(f2c_edges.cell_blk_f)),
      vertex_idx(make_index_view3d(f2c_edges.vertex_idx_f)),
      vertex_blk(make_index_view3d(f2c_edges.vertex_blk_f)) {}

GridVertices::GridVertices(const f2c::GridVertices& f2c_vertices)
    : GridBase(f2c_vertices.cell_idx_f.shape[0], make_refin_ctrl_bounds(f2c_vertices)),
      // data fields:
      vertex(make_constview2d<GeographicCoords>(f2c_vertices.vertex)),
      vertex_lon(make_component_view2d<GeographicCoords>(f2c_vertices.vertex, &GeographicCoords::lon)),
      vertex_lat(make_component_view2d<GeographicCoords>(f2c_vertices.vertex, &GeographicCoords::lat)),
      // index fields:
      cell_idx(make_index_view3d(f2c_vertices.cell_idx_f)),
      cell_blk(make_index_view3d(f2c_vertices.cell_blk_f)) {}

Patch::Patch(const f2c::PatchDescr patch_descr)
    : f2c::PatchInfo(gen_patch_info(patch_descr)),
      comm_patch(patch_descr),
      cells(get_f2c_grid_cells(patch_descr)),
      edges(get_f2c_grid_edges(patch_descr)),
      verts(get_f2c_grid_vertices(patch_descr)) {}

const Patch& get_patch(int dom_id) {
  // dom_id follows model conventions: 1-based, with 0 reserved for the optional radiation grid.
  int pos = p_patch_id_offset + dom_id;
  require(pos >= 0 and pos < static_cast<int>(p_patch.size()));
  return *p_patch[pos];
}

namespace domain {

void init() {
  require(f2c::is_initialized());
  require(p_patch.size() == 0);
  const f2c::FunTable& fun_table  = f2c::get_fun_table();
  const f2c::DomainInfo& dom_info = f2c::get_domain_info();
  int nids                        = dom_info.id_max - dom_info.id_min + 1;
  require(nids > 0);
  p_patch.reserve(nids);
  require(dom_info.id_min == 0 or dom_info.id_min == 1);
  // id_min is 0 when a radiation grid is present, 1 otherwise; the offset maps dom_id to p_patch index.
  p_patch_id_offset = -dom_info.id_min;
  for (int id = dom_info.id_min; id <= dom_info.id_max; id++) {
    const f2c::PatchDescr f2c_patch_descr = fun_table.get_mo_model_domain_p_patch_descr(id);
    const auto& patch                     = *p_patch.emplace_back(std::make_unique<Patch>(f2c_patch_descr));
    require(patch.id == id);
  }
  require(static_cast<int>(p_patch.size()) == nids);

  // Views must be released before Kokkos::finalize().
  Kokkos::push_finalize_hook([]() { p_patch.clear(); });
}

}  // namespace domain

}  // namespace icon
