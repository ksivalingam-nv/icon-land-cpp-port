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
/// @brief support forFortran <=> C access.
///

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_H_

#include <cstddef>

#include "finish.hpp"
#include "icon_f2c_basic_types.hpp"

namespace icon {

namespace f2c {

// opaque pointer to a Fortran t_patch instance
struct PatchDescr {
  void* cptr = NULL;
};

// portable part of top level t_patch components
struct PatchInfo {
  int id       = -1;
  int nlev     = 0;
  int nblks_c  = 0;
  int npromz_c = 0;
  int nblks_e  = 0;
  int npromz_e = 0;
  int nblks_v  = 0;
  int npromz_v = 0;
  // to be extended
};

// domain info:
struct DomainInfo {
  int nproma = 0;
  int id_min = -1;
  int id_max = -2;
};

// opaque pointer to a Fortran ICON CommPattern instance
struct CommPatternDescr {
  void* cptr = NULL;
};

inline constexpr CommPatternDescr comm_pat_null = {NULL};

using CommPattern                               = CommPatternDescr;

inline constexpr int f2c_repr_max_ndims         = 5;

struct FieldRepr {
  void* data;
  int f2c_kind;
  int rank;
  int shape[f2c_repr_max_ndims];
  int lbounds[f2c_repr_max_ndims];
};

// Guard binary compatibility with Fortran t_f2c_repr5 (mo_ragnarok_f2c.f90).
// Both sides follow C struct-layout rules (BIND(C)); check that no unexpected
// internal padding has been inserted by verifying each field's offset.
static_assert(offsetof(FieldRepr, data) == 0);
static_assert(offsetof(FieldRepr, f2c_kind) == sizeof(void*));
static_assert(offsetof(FieldRepr, rank) == sizeof(void*) + sizeof(int));
static_assert(offsetof(FieldRepr, shape) == sizeof(void*) + 2 * sizeof(int));
static_assert(offsetof(FieldRepr, lbounds) == sizeof(void*) + 7 * sizeof(int));
// Also lock the total size (no trailing padding): FieldRepr members sit back-to-back
// in the Grid* structs below, so their offsets are multiples of sizeof(FieldRepr).
static_assert(sizeof(FieldRepr) == sizeof(void*) + 12 * sizeof(int));

// Returns a typed pointer to the data pointer in a FieldRepr structure.
// Evaluation of f2c_typekind<T> makes this a safe conversion if successful
// or a failed conversion with a fatal error.
template <typename T>
T* get_data_ptr(const FieldRepr& repr) {
  static_assert(f2c_typekind<T> != 0, "get_data_ptr: Unsupported type for f2c_typekind");
  require(repr.f2c_kind == f2c_typekind<T>, "bad type conversion");
  return static_cast<T*>(repr.data);
}

struct GridCells {
  FieldRepr num_edges;
  FieldRepr area;
  FieldRepr center;
  FieldRepr edge_idx_f;
  FieldRepr edge_blk_f;
  FieldRepr start_block_f;
  FieldRepr end_block_f;
  FieldRepr start_index_f;
  FieldRepr end_index_f;
};

struct GridEdges {
  FieldRepr primal_normal_cell;
  FieldRepr primal_normal_vert;
  FieldRepr dual_normal_cell;
  FieldRepr dual_normal_vert;
  FieldRepr tangent_orientation;
  FieldRepr inv_primal_edge_length;
  FieldRepr inv_dual_edge_length;
  FieldRepr inv_vert_vert_length;
  FieldRepr edge_cell_length;
  FieldRepr cell_idx_f;
  FieldRepr cell_blk_f;
  FieldRepr vertex_idx_f;
  FieldRepr vertex_blk_f;
  FieldRepr start_block_f;
  FieldRepr end_block_f;
  FieldRepr start_index_f;
  FieldRepr end_index_f;
};

struct GridVertices {
  FieldRepr vertex;
  FieldRepr cell_idx_f;
  FieldRepr cell_blk_f;
  FieldRepr start_block_f;
  FieldRepr end_block_f;
  FieldRepr start_index_f;
  FieldRepr end_index_f;
};

// communication patch
struct CommPatch {
  CommPattern comm_pat_c;
  CommPattern comm_pat_e;
  CommPattern comm_pat_v;
  CommPattern comm_pat_c1;
};

struct ProcessInfo {
  bool is_mpi_parallel;
};

struct FunTable {
  void (*get_domain_info)(DomainInfo* dom_info)                                               = NULL;
  void (*get_patch_info)(const PatchDescr patch_descr, PatchInfo* patch_info)                 = NULL;
  PatchDescr (*get_mo_model_domain_p_patch_descr)(int id)                                     = NULL;
  void (*get_comm_patch)(const PatchDescr patch_descr, CommPatch* comm_patch)                 = NULL;
  void (*get_grid_cells)(const PatchDescr patch_descr, GridCells* grid_cells)                 = NULL;
  void (*get_grid_edges)(const PatchDescr patch_descr, GridEdges* grid_edges)                 = NULL;
  void (*get_grid_vertices)(const PatchDescr patch_descr, GridVertices* grid_vertices)        = NULL;
  void (*message)(const char* name, int name_len, const char* text, int text_len)             = NULL;
  void (*finish)(const char* name, int name_len, const char* text, int text_len)              = NULL;
  void (*exchange_data_r3d)(CommPatternDescr descr, bool lacc, double* recv, int* recv_shape) = NULL;
  ProcessInfo (*get_process_info)()                                                           = NULL;
  int (*new_timer)(const char* name, int name_len)                                            = NULL;
  void (*timer_start)(int timer)                                                              = NULL;
  void (*timer_stop)(int timer)                                                               = NULL;
  double (*timer_value)(int timer)                                                            = NULL;
  int (*get_timers_level)()                                                                   = NULL;
};

// see docs/friend-gated-details.md for rationale of this construct
void init(const FunTable* funtab);
bool is_initialized();
const FunTable& get_fun_table();
const DomainInfo& get_domain_info();
bool have_radiation_grid();

class Internals {
  Internals() = delete;
  // access functions:
  friend void init(const FunTable* funtab) {
    require(not init_state, "init_internals called more than once");
    fun_table = *funtab;
    funtab->get_domain_info(&dom_info);
    require(dom_info.id_max >= dom_info.id_min, "no domain available");
    init_state = true;
  }
  friend bool is_initialized() { return init_state; }
  friend const FunTable& get_fun_table() {
    debug_require(init_state);
    return fun_table;
  }
  friend const DomainInfo& get_domain_info() {
    debug_require(init_state);
    return dom_info;
  }
  friend bool have_radiation_grid() {
    debug_require(init_state);
    return dom_info.id_min == 0;
  }
  // state:
  inline static bool init_state = false;
  inline static FunTable fun_table{};
  inline static DomainInfo dom_info{};
};

}  // namespace f2c
}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_H_
