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

#ifndef RAGNAROK_SUPPORT_TIMER_TIMER_OPS_HPP_
#define RAGNAROK_SUPPORT_TIMER_TIMER_OPS_HPP_

namespace timer {

// The installed backend's operations; Timer dispatches through them.
struct TimerOps {
  using gen_id_t       = int (*)(const char* name, int name_len);
  using timer_switch_t = void (*)(int timer);
  using timer_value_t  = double (*)(int timer);

  constexpr TimerOps(gen_id_t gen_id, timer_switch_t start, timer_switch_t stop, timer_value_t value)
      : gen_id(gen_id), start(start), stop(stop), value(value) {}

  gen_id_t gen_id;
  timer_switch_t start;
  timer_switch_t stop;
  timer_value_t value;
};

}  // namespace timer

#endif  // RAGNAROK_SUPPORT_TIMER_TIMER_OPS_HPP_
