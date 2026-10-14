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

#ifndef RAGNAROK_SUPPORT_TIMER_TIMER_SYNC_HPP_
#define RAGNAROK_SUPPORT_TIMER_TIMER_SYNC_HPP_

namespace timer {

// the timers_level at which device sync before measuring begins;
// identical behavior to mo_real_timer, which syncs above level 10
constexpr int timers_level_sync_threshold = 11;

// device sync before each measurement in the synced regime; the definition
// (icon::kernel_fence) lives in timer_sync.cpp, keeping Kokkos out of the
// timer headers
void sync();

}  // namespace timer

#endif  // RAGNAROK_SUPPORT_TIMER_TIMER_SYNC_HPP_
