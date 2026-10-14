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

// Everything that is not performance-relevant: the internal backend, the
// configuration entry points, and Timer's cold paths.

#include "timer.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "support/common/finish.hpp"

namespace timer {

namespace {

// cacheline size of the targeted systems (x86-64); architectures with
// longer lines still keep most of the isolation below
constexpr std::size_t cacheline_size = 64;

// stats tables grow in chunks of eight timers: bounded slack instead of
// geometric growth, and late-constructed timers (function-local statics)
// rewrite the tables only at every eighth construction
constexpr int growth_granularity     = 8;
constexpr int rounded_capacity(int count) {
  return (count + growth_granularity - 1) / growth_granularity * growth_granularity;
}

// Self-contained wall-clock backend with per-thread measurement state.
// Threading model (OpenMP only; other threads are undefined behavior): one
// flat team, fixed at configure() - the contract mo_real_timer already
// imposes, whose new_timer initializes every thread's slot from the master.
// All per-thread state is keyed by team rank in rows the master lays out in
// serial context; each row is written only by its owning thread. Outside
// parallel regions exactly one thread executes, so all serial-context code
// below is mutually exclusive by contract, and the region join publishes
// the team's rows. No shared structure is ever written inside a parallel
// region.
class InternalBackend {
 public:
  // per-thread measurement state of one timer; one cacheline each, so
  // measurement writes never share a line across threads (adjacent tables
  // would otherwise touch at their boundaries):
  struct alignas(cacheline_size) Stats {
    std::chrono::steady_clock::time_point start_ts;  // start of the ongoing measurement
    double time_sum = 0.0;                           // sum over all completed measurements
    int call_num    = 0;                             // number of completed measurements
    bool is_on      = false;                         // true while a measurement is ongoing
  };

  // Per-thread stack of running measurements, enforcing LIFO nesting like
  // mo_real_timer's active_timers: the innermost running timer must stop
  // first. Fixed depth - deeper nesting aborts, mirroring the Fortran cap.
  struct ActiveStack {
    static constexpr int max_depth = 32;
    int ids[max_depth];
    int depth = 0;
  };

  // one team rank's state; cacheline-aligned so that a rank's stack writes
  // never share a line with its neighbor's row
  struct alignas(cacheline_size) ThreadState {
    std::vector<Stats> stats;  // always sized to count
    ActiveStack stack;
  };

  struct SharedState {
    std::vector<std::string> names;    // grown by gen_id only (serial context)
    int count = 0;                     // issued ids (serial context)
    std::vector<ThreadState> threads;  // one row per team rank; laid out by configure(), empty before
  };

  // immortal: stays usable during static destruction
  static SharedState& shared_state() {
    static SharedState& instance = *new SharedState;
    return instance;
  }

  // the calling thread's row index
  static int thread_rank() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
  }

  // the team configure() lays out rows for
  static int team_size() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
  }

  // calling thread's row; rows exist for the whole team from configure()
  // on, so the measurement path stays free of locks and size checks
  static ThreadState& my_state() {
    SharedState& s = shared_state();
    const int rank = thread_rank();
    icon::debug_require(static_cast<std::size_t>(rank) < s.threads.size(),
                        "timer::InternalBackend: thread outside the configured team");
    return s.threads[rank];
  }

  // chunk-rounded growth. The reserve is load-bearing: it fixes the
  // capacity at the chunk rounding, where resize alone would grow
  // geometrically.
  static void grow(std::vector<Stats>& stats, int count) {
    stats.reserve(rounded_capacity(count));
    stats.resize(count);
  }

  static int gen_id(const char* name, int name_len) {
    require_serial_context("timer::InternalBackend: gen_id() within a parallel region");
    SharedState& s = shared_state();
    s.names.emplace_back(name, name_len);
    s.count++;
    // grow every row now, in serial context:
    for (ThreadState& t : s.threads) grow(t.stats, s.count);
    return s.count - 1;
  }

  // one fresh, fully sized row per team rank (serial context). Precondition:
  // Kokkos is initialized - it fixes the team size, and the OpenMP runtime
  // must not be first queried before its first region (Kokkos documents
  // the Intel/hwloc process-mask pitfall).
  static void lay_out_team() {
    require_serial_context("timer::InternalBackend: team layout within a parallel region");
    SharedState& s = shared_state();
    s.threads.assign(team_size(), ThreadState{});
    for (ThreadState& t : s.threads) grow(t.stats, s.count);
  }

  static void start(int id) {
    ThreadState& t = my_state();
    Stats& s       = t.stats[id];
    icon::require(not s.is_on, "timer::InternalBackend: start() on a running timer");
    ActiveStack& a = t.stack;
    icon::require(a.depth < ActiveStack::max_depth, "timer::InternalBackend: measurements nested deeper than 32");
    a.ids[a.depth++] = id;
    s.is_on          = true;
    s.start_ts       = std::chrono::steady_clock::now();
  }

  static void stop(int id) {
    ThreadState& t = my_state();
    Stats& s       = t.stats[id];
    icon::require(s.is_on, "timer::InternalBackend: stop() without start()");
    // is_on implies this thread started id, so the stack is non-empty:
    ActiveStack& a = t.stack;
    if (a.ids[a.depth - 1] != id) {
      icon::finish("timer::InternalBackend", "stop('" + shared_state().names[id] + "') while '" +
                                                 shared_state().names[a.ids[a.depth - 1]] +
                                                 "' is the innermost running timer - stop that first");
    }
    a.depth--;
    s.time_sum += seconds_since(s.start_ts);
    s.is_on = false;
    s.call_num++;
  }

  // calling thread's time sum plus its ongoing measurement, without
  // stopping it; nothing before the configuration laid out the rows
  static double value(int id) {
    if (static_cast<std::size_t>(thread_rank()) >= shared_state().threads.size()) return 0.0;
    const Stats& s = my_state().stats[id];
    double val     = s.time_sum;
    if (s.is_on) val += seconds_since(s.start_ts);
    return val;
  }

  static double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  // abort when used within a parallel region (mirrors the Fortran backend):
  static void require_serial_context([[maybe_unused]] const char* what) {
#ifdef _OPENMP
    icon::require(not omp_in_parallel(), what);
#endif
  }

  // abort when id has a running measurement on any thread: a destroyed
  // timer must not leave an open measurement behind in the tables
  static void require_closed(int id, const std::string& name) {
    for (const ThreadState& t : shared_state().threads) {
      if (t.stats[id].is_on) icon::finish("timer::Timer", "destroying timer '" + name + "' during its measurement");
    }
  }

  // give back all internal-backend memory; only legal after the switch,
  // when the internal operations are unreachable
  static void release() {
    SharedState& s = shared_state();
    s.names        = {};
    s.threads      = {};
    s.count        = 0;
  }

  // completed measurements only; ncalls summed over threads, total as
  // min/max over the threads that used the timer
  static void report() {
    require_serial_context("timer::InternalBackend: report() within a parallel region");
    SharedState& s = shared_state();
    std::cout << "Timer report:\n";
    for (std::size_t id = 0; id < s.names.size(); id++) {
      int ncalls = 0, nthreads = 0;
      double tot_min = 0., tot_max = 0.;
      for (const ThreadState& t : s.threads) {
        const Stats& st = t.stats[id];
        if (st.call_num == 0) continue;
        ncalls += st.call_num;
        tot_min = nthreads ? std::min(tot_min, st.time_sum) : st.time_sum;
        tot_max = nthreads ? std::max(tot_max, st.time_sum) : st.time_sum;
        nthreads++;
      }
      if (ncalls == 0) continue;
      std::cout << s.names[id] << ": ncalls=" << ncalls;
      if (nthreads > 1) {
        std::cout << ", threads=" << nthreads << ", total=[" << tot_min << ".." << tot_max << "]\n";
      } else {
        std::cout << ", total=" << tot_max << "\n";
      }
    }
  }
};

// timers_level is a runtime constant: exactly one configuration
bool configured = false;

// true once the backend switch happened; the switch is one-way
bool switched   = false;

// pool of all live timers; only accessed single-threaded (construction,
// destruction, and the backend switch); immortal like shared_state()
std::vector<Timer*>& pool() {
  static std::vector<Timer*>& instance = *new std::vector<Timer*>;
  return instance;
}

}  // namespace

// pre-installed internal backend; constant-initialized, so namespace-scope
// timers constructed before main() dispatch safely
constinit TimerOps Timer::backend{InternalBackend::gen_id, InternalBackend::start, InternalBackend::stop,
                                  InternalBackend::value};

Timer::Timer(const std::string& name, int threshold_level) : name(name), threshold_level(threshold_level) {
  // level 0 would stay active even with timers switched off (ltimer=false
  // arrives as timers_level 0):
  icon::require(threshold_level >= 1, "timer::Timer: threshold_level below 1");
  issue_id();
  pool_pos = pool().size();
  pool().push_back(this);
}

// O(1) deregistration: swap-and-pop, updating the moved timer's slot. Only
// the internal tables are queryable for an open measurement; after the
// switch the ICON backend owns that state.
Timer::~Timer() {
  if (not switched) InternalBackend::require_closed(id, name);
  std::swap(pool()[pool_pos], pool().back());
  pool()[pool_pos]->pool_pos = pool_pos;
  pool().pop_back();
}

// id generation always succeeds (both backends abort internally on failure)
void Timer::issue_id() { id = backend.gen_id(name.c_str(), static_cast<int>(name.length())); }

// Adopt the timers_level, once. Until then the level is 0, so no
// measurement precedes the configuration. Both configuration entry points
// run after Kokkos initialization, so `synced` never predates it - sync()
// is safe.
void Timer::apply_config(int level) {
  InternalBackend::require_serial_context("timer::Timer: configuration within a parallel region");
  icon::require(not configured, "timer::Timer: configured more than once");
  configured   = true;
  timers_level = level;
  synced       = timers_level >= timers_level_sync_threshold;
}

// Switch from the internal to the bridged backend. Ids are
// backend-specific, so all pooled timers get new ids.
void Timer::set_backend(TimerOps ops, int level) {
  apply_config(level);
  backend  = ops;
  switched = true;
  for (Timer* t : pool()) t->issue_id();
  InternalBackend::release();
}

void configure(int timers_level) {
  Timer::apply_config(timers_level);
  // the team's rows exist from here on, so measurements inside parallel
  // regions stay allocation-free. Not done at the switch: bridged runs
  // never use the internal tables.
  InternalBackend::lay_out_team();
}

void report_all_timers() {
  if (switched) {
    std::cout << "report_all_timers() not applicable: see the ICON timer report\n";
  } else {
    InternalBackend::report();
  }
}

}  // namespace timer
