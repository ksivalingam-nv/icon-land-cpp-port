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

// Timer with switchable backend: timers construct from program start and
// record once a timers_level is configured. Standalone builds measure on the
// ragnarok-internal backend; in the ICON build the icon bridge's
// connect_icon_backend() switches to the ICON Fortran timers (mo_timer.f90,
// mo_real_timer.f90) before any measurement.
// Only start/stop are performance-relevant and stay inline; everything else
// lives in timer.cpp, together with the internal backend.
// See docs/timer_usage.md for usage and docs/timer_design.md for the model.
#ifndef RAGNAROK_SUPPORT_TIMER_TIMER_HPP_
#define RAGNAROK_SUPPORT_TIMER_TIMER_HPP_

#include <cstddef>
#include <string>

#include "timer_ops.hpp"
#include "timer_sync.hpp"

namespace timer {

// a single named, accumulating timer
class Timer {
  // the configuration entry points (defined in timer.cpp and
  // icon_bridge/icon_timer.cpp) drive the private statics:
  friend void connect_icon_backend();
  friend void configure(int timers_level);

 public:
  // threshold_level: the timers_level at which the timer becomes active,
  // fixed at declaration; the default reproduces a plain IF (ltimer) guard
  explicit Timer(const std::string& name, int threshold_level = 1);

  // destroyed timers drop out of re-registration; their completed internal
  // measurements remain in the report. Destruction during a running
  // measurement aborts (checked on the internal backend).
  ~Timer();

  // copies would alias the pooled registration:
  Timer(const Timer&)            = delete;
  Timer& operator=(const Timer&) = delete;

  void start() {
    if (not active()) return;
    if (synced) sync();
    backend.start(id);
  }

  void stop() {
    if (not active()) return;
    if (synced) sync();
    backend.stop(id);
  }

  // the configured timers_level has reached the timer's threshold_level
  bool active() const { return timers_level >= threshold_level; }

  // the calling thread's accumulated time
  double value() const { return backend.value(id); }

  // the key under which the current backend knows this timer; after the
  // switch it is the ICON timer id, accepted by mo_real_timer's queries
  int backend_id() const { return id; }

 private:
  int id = -1;  // issued at construction, re-issued at the switch
  std::string name;
  int threshold_level;

  // our slot in the pool, kept current by the O(1) deregistration:
  // pool()[t->pool_pos] == t
  std::size_t pool_pos = 0;

  // the current backend's operations; constant-initialized in timer.cpp,
  // so timers work from program start
  static TimerOps backend;

  // the configured level; unconfigured (0) equals timers off, so no
  // measurement precedes the configured intent
  static inline int timers_level = 0;

  // detailed runs sync before every measurement (mo_real_timer regime)
  static inline bool synced      = false;

  void issue_id();
  static void apply_config(int level);
  static void set_backend(TimerOps ops, int level);
};

// Standalone counterpart of connect_icon_backend(): adopts the timers_level
// for the pre-installed internal backend. Exactly one configuration per
// program - a second one, by either entry point, aborts.
void configure(int timers_level);

// internal-backend report; after the switch it defers to the ICON report
void report_all_timers();

// bridge socket: defined in icon_bridge/icon_timer.cpp, where all f2c
// interop lives
void connect_icon_backend();

}  // namespace timer

#endif  // RAGNAROK_SUPPORT_TIMER_TIMER_HPP_
