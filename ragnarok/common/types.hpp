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
#ifndef RAGNAROK_COMMON_TYPES_H_
#define RAGNAROK_COMMON_TYPES_H_

#include <Kokkos_Core.hpp>
#include <type_traits>

using UnmanagedRestrictedTrait = Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>;
using ManagedRestrictedTrait   = Kokkos::MemoryTraits<Kokkos::Restrict>;

// Note: Rank 1 deliberately carries no layout - LayoutRight and LayoutLeft coincide there,
// and spelling it out would silently make these types distinct.

// unmanaged mutable views
template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using View0D = Kokkos::View<T, ExecSpace, UnmanagedRestrictedTrait>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using View1D = Kokkos::View<T*, ExecSpace, UnmanagedRestrictedTrait>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using View2D = Kokkos::View<T**, Kokkos::LayoutRight, ExecSpace, UnmanagedRestrictedTrait>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using View3D = Kokkos::View<T***, Kokkos::LayoutRight, ExecSpace, UnmanagedRestrictedTrait>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using View4D = Kokkos::View<T****, Kokkos::LayoutRight, ExecSpace, UnmanagedRestrictedTrait>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using View5D = Kokkos::View<T*****, Kokkos::LayoutRight, ExecSpace, UnmanagedRestrictedTrait>;

// Rank-indexed selector for unmanaged mutable views
template <typename T, int Rank, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using RankedView = std::conditional_t<
    Rank == 0, View0D<T, ExecSpace>,
    std::conditional_t<
        Rank == 1, View1D<T, ExecSpace>,
        std::conditional_t<
            Rank == 2, View2D<T, ExecSpace>,
            std::conditional_t<Rank == 3, View3D<T, ExecSpace>,
                               std::conditional_t<Rank == 4, View4D<T, ExecSpace>,
                                                  std::conditional_t<Rank == 5, View5D<T, ExecSpace>, void>>>>>>;

// unmanaged immutable views
template <typename T, typename ExecSpace>
using ConstView0D = View0D<const T, ExecSpace>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ConstView1D = View1D<const T, ExecSpace>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ConstView2D = View2D<const T, ExecSpace>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ConstView3D = View3D<const T, ExecSpace>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ConstView4D = View4D<const T, ExecSpace>;

// managed mutable views
template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ManagedView0D = Kokkos::View<T, ExecSpace, ManagedRestrictedTrait>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ManagedView1D = Kokkos::View<T*, ExecSpace, ManagedRestrictedTrait>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ManagedView2D = Kokkos::View<T**, Kokkos::LayoutRight, ExecSpace, ManagedRestrictedTrait>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ManagedView3D = Kokkos::View<T***, Kokkos::LayoutRight, ExecSpace, ManagedRestrictedTrait>;

// managed immutable views
template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ManagedConstView1D = ManagedView1D<const T, ExecSpace>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ManagedConstView2D = ManagedView2D<const T, ExecSpace>;

template <typename T, typename ExecSpace = Kokkos::DefaultExecutionSpace>
using ManagedConstView3D = ManagedView3D<const T, ExecSpace>;

// Host shorthands. The space is always spelled out: relying on the default above would make the
// managed ones allocate in device memory as soon as the default execution space is a device.

// unmanaged mutable host views
template <typename T>
using HostView0D = View0D<T, Kokkos::HostSpace>;

template <typename T>
using HostView1D = View1D<T, Kokkos::HostSpace>;

template <typename T>
using HostView2D = View2D<T, Kokkos::HostSpace>;

template <typename T>
using HostView3D = View3D<T, Kokkos::HostSpace>;

// unmanaged immutable host views
template <typename T>
using HostConstView1D = ConstView1D<T, Kokkos::HostSpace>;

template <typename T>
using HostConstView2D = ConstView2D<T, Kokkos::HostSpace>;

template <typename T>
using HostConstView3D = ConstView3D<T, Kokkos::HostSpace>;

template <typename T>
using HostConstView4D = ConstView4D<T, Kokkos::HostSpace>;

// managed mutable host views
template <typename T>
using ManagedHostView1D = ManagedView1D<T, Kokkos::HostSpace>;

template <typename T>
using ManagedHostView2D = ManagedView2D<T, Kokkos::HostSpace>;

template <typename T>
using ManagedHostView3D = ManagedView3D<T, Kokkos::HostSpace>;

// managed immutable host views
template <typename T>
using ManagedHostConstView1D = ManagedConstView1D<T, Kokkos::HostSpace>;

template <typename T>
using ManagedHostConstView2D = ManagedConstView2D<T, Kokkos::HostSpace>;

template <typename T>
using ManagedHostConstView3D = ManagedConstView3D<T, Kokkos::HostSpace>;

// strided views
template <typename T>
using StridedView1D = Kokkos::View<T*, Kokkos::LayoutStride, UnmanagedRestrictedTrait>;

template <typename T>
using StridedView2D = Kokkos::View<T**, Kokkos::LayoutStride, UnmanagedRestrictedTrait>;

template <typename T>
using StridedView3D = Kokkos::View<T***, Kokkos::LayoutStride, UnmanagedRestrictedTrait>;

template <typename T>
using StridedView4D = Kokkos::View<T****, Kokkos::LayoutStride, UnmanagedRestrictedTrait>;

// Rank-indexed selector for strided views.
template <typename T, int Rank>
using StridedRankedView =
    std::conditional_t<Rank == 1, StridedView1D<T>,
                       std::conditional_t<Rank == 2, StridedView2D<T>,
                                          std::conditional_t<Rank == 3, StridedView3D<T>,
                                                             std::conditional_t<Rank == 4, StridedView4D<T>, void>>>>;

// Kokkos default memory space
using MemorySpace   = Kokkos::DefaultExecutionSpace::memory_space;

using RangePolicy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>;

#endif  // RAGNAROK_COMMON_TYPES_H_
