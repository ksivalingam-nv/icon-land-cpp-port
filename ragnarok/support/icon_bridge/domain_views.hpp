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

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_DOMAIN_VIEWS_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_DOMAIN_VIEWS_H_

#include <Kokkos_OffsetView.hpp>

#include "common/types.hpp"
#include "finish.hpp"
#include "icon_f2c.hpp"

// offsetview types needed to capture the established domain classification scheme using negative indices:
// ManagedHostOffsetView1D<const int> and ManagedHostConstOffsetView1D<int> must stay the same type;
// icon_domain.cpp relies on it to reach a RefinCtrlBounds member with a single user-defined conversion.

template <typename T>
using ManagedHostConstOffsetView1D =
    Kokkos::Experimental::OffsetView<const T*, Kokkos::HostSpace, ManagedRestrictedTrait>;

template <typename T>
using ManagedHostOffsetView1D = Kokkos::Experimental::OffsetView<T*, Kokkos::HostSpace, ManagedRestrictedTrait>;

namespace icon {
namespace f2c {

// view related helper:

// true if every dimension of repr uses a 1-based lower bound (i.e. no Fortran index offset).
inline bool have_default_offsets(const f2c::FieldRepr& repr) {
  const auto& lb     = repr.lbounds;
  bool have_all_ones = true;
  for (int i = 0; i < repr.rank; ++i) {
    have_all_ones &= (lb[i] == 1);
  }
  return have_all_ones;
}

// Builds a mutable, unmanaged 2D view over the data represented by repr.
// Reversed extents -> LayoutRight
template <typename T>
[[nodiscard]] HostView2D<T> make_view2d(const f2c::FieldRepr& repr) {
  require(repr.rank == 2);
  require(have_default_offsets(repr));
  return HostView2D<T>(get_data_ptr<T>(repr), repr.shape[1], repr.shape[0]);
}

// Builds a mutable, unmanaged 3D view over the data represented by repr.
// Reversed extents -> LayoutRight
template <typename T>
[[nodiscard]] HostView3D<T> make_view3d(const f2c::FieldRepr& repr) {
  require(repr.rank == 3);
  require(have_default_offsets(repr));
  return HostView3D<T>(get_data_ptr<T>(repr), repr.shape[2], repr.shape[1], repr.shape[0]);
}

// Const-view counterparts, expressed in terms of the mutable builders above.

template <typename T>
[[nodiscard]] HostConstView2D<T> make_constview2d(const f2c::FieldRepr& repr) {
  return make_view2d<T>(repr);
}

template <typename T>
[[nodiscard]] HostConstView3D<T> make_constview3d(const f2c::FieldRepr& repr) {
  return make_view3d<T>(repr);
}

// Reversed extents -> LayoutRight (see make_view3d).
template <typename T>
[[nodiscard]] HostConstView4D<T> make_constview4d(const f2c::FieldRepr& repr) {
  require(repr.rank == 4);
  require(have_default_offsets(repr));
  return HostConstView4D<T>(get_data_ptr<T>(repr), repr.shape[3], repr.shape[2], repr.shape[1], repr.shape[0]);
}

// --- index conversion (Fortran 1-based -> C zero-based):

inline void transform_fortran_to_c_indices(const int* s_ptr, int* d_ptr, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    d_ptr[i] = s_ptr[i] - 1;
  }
}

// Index-converting view factory; goes through the zero-cost const view above,
// which owns the FieldRepr->LayoutRight convention.
[[nodiscard]] inline ManagedHostView3D<int> make_index_view3d(const f2c::FieldRepr& repr) {
  const HostConstView3D<int> src = make_constview3d<int>(repr);
  ManagedHostView3D<int> dst(Kokkos::view_alloc(Kokkos::WithoutInitializing, "index_view3d"), src.extent(0),
                             src.extent(1), src.extent(2));
  transform_fortran_to_c_indices(src.data(), dst.data(), src.size());
  return dst;
}

// --- struct component extraction (array-of-structs -> struct-of-arrays):

template <typename StructT, typename MemberT>
void extract_member(const StructT* s_ptr, MemberT* d_ptr, size_t n, MemberT StructT::* member) {
  for (size_t i = 0; i < n; ++i) d_ptr[i] = s_ptr[i].*member;
}

template <typename StructT, typename MemberT>
[[nodiscard]] ManagedHostView2D<MemberT> make_component_view2d(const f2c::FieldRepr& repr, MemberT StructT::* member) {
  const HostConstView2D<StructT> src = make_constview2d<StructT>(repr);
  ManagedHostView2D<MemberT> dst(Kokkos::view_alloc(Kokkos::WithoutInitializing, "component_view2d"), src.extent(0),
                                 src.extent(1));
  extract_member(src.data(), dst.data(), src.size(), member);
  return dst;
}

template <typename StructT, typename MemberT>
[[nodiscard]] ManagedHostView3D<MemberT> make_component_view3d(const f2c::FieldRepr& repr, MemberT StructT::* member) {
  const HostConstView3D<StructT> src = make_constview3d<StructT>(repr);
  ManagedHostView3D<MemberT> dst(Kokkos::view_alloc(Kokkos::WithoutInitializing, "component_view3d"), src.extent(0),
                                 src.extent(1), src.extent(2));
  extract_member(src.data(), dst.data(), src.size(), member);
  return dst;
}

}  // namespace f2c
}  // namespace icon

#endif
