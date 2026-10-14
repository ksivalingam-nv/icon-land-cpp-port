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
/// Basic C++ types for ICON-Fortran C interoperability.
/// Maps basic f2c types to integer and string identifiers.

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_BASIC_TYPES_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_BASIC_TYPES_H_

#include "element_types.hpp"

namespace icon {
namespace f2c {

// C++ bool and C _Bool are the same underlying type per the Itanium ABI
// (GCC/Clang on Linux), which is what Fortran's ISO_C_BINDING targets.
using c_bool = bool;
static_assert(sizeof(c_bool) == 1, "c_bool must be 1 byte to match C _Bool ABI");
static_assert(alignof(c_bool) == 1, "c_bool must have 1-byte alignment to match C _Bool ABI");

/// @brief Maps C++ types to integer type IDs.
template <typename T>
inline constexpr int f2c_typekind = 0;

template <>
inline constexpr int f2c_typekind<c_bool> = 1;
template <>
inline constexpr int f2c_typekind<int> = 2;
template <>
inline constexpr int f2c_typekind<float> = 3;
template <>
inline constexpr int f2c_typekind<double> = 4;

template <RealTypes T>
inline constexpr int f2c_typekind<CartCoords<T>> = std::same_as<T, float> ? 5 : 8;
template <RealTypes T>
inline constexpr int f2c_typekind<GeoCoords<T>> = std::same_as<T, float> ? 6 : 9;
template <RealTypes T>
inline constexpr int f2c_typekind<TangVec<T>> = std::same_as<T, float> ? 7 : 10;

/// @brief String names for each F2C typekind integer ID.
constexpr const char* typekind_str[11]        = {
    "undef",         // 0
    "c_bool",        // 1
    "int",           // 2
    "float",         // 3
    "double",        // 4
    "CartCoordsSp",  // 5
    "GeoCoordsSp",   // 6
    "TangVecSp",     // 7
    "CartCoordsDp",  // 8
    "GeoCoordsDp",   // 9
    "TangVecDp",     // 10
};

}  // namespace f2c
}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_F2C_BASIC_TYPES_H_
