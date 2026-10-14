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

// Bridge wiring: switches the timer backend to the ICON Fortran timers.

#include "icon_f2c.hpp"
#include "support/timer/timer.hpp"

namespace timer {

void connect_icon_backend() {
  const icon::f2c::FunTable& ftab = icon::f2c::get_fun_table();
  icon::require(ftab.get_timers_level and ftab.new_timer and ftab.timer_start and ftab.timer_stop and ftab.timer_value,
                "timer::connect_icon_backend(): incomplete timer entries in the ftable");

  // switch to the bridged backend; the getter maps ltimer=false to
  // timers_level 0, deactivating every timer, and floors it at 1 otherwise:
  Timer::set_backend(TimerOps(ftab.new_timer, ftab.timer_start, ftab.timer_stop, ftab.timer_value),
                     ftab.get_timers_level());
}

}  // namespace timer

extern "C" void init_ragnarok_timer() { timer::connect_icon_backend(); }
