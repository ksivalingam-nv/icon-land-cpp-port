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

/// @file
/// numerical element types used in ICON

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ELEMENT_TYPES_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ELEMENT_TYPES_H_

#include <concepts>
#include <type_traits>

#include "common/real.hpp"

namespace icon {

/// @brief Constrains template parameters to the two REAL kinds we use.
template <typename T>
concept RealTypes = std::same_as<T, float> || std::same_as<T, double>;

// icon math types: t_cartesian_coordinates, t_geographical_coordinates, t_tangent_vectors.
// All three structs are standard-layout by construction and therefore
// are interoperable with appropriate BIND(C) Fortran types.

template <RealTypes T = Real>
struct CartCoords {
  T x[3];
};

template <RealTypes T = Real>
struct GeoCoords {
  T lon;
  T lat;
};

template <RealTypes T = Real>
struct TangVec {
  T v1;
  T v2;
};

// working-precision aliases — equivalent to, e.g., CartCoords<>
using CartesianCoords  = CartCoords<Real>;
using GeographicCoords = GeoCoords<Real>;
using TangentVec       = TangVec<Real>;

}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_ELEMENT_TYPES_H_
