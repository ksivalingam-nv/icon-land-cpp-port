<!--
ICON

---------------------------------------------------------------
Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
Contact information: icon-model.org

See AUTHORS.TXT for a list of authors
See LICENSES/ for license information
SPDX-License-Identifier: BSD-3-Clause
---------------------------------------------------------------
-->

# Timer usage

Timers measure the wall-clock time of code sections. In the ICON build the
results appear in the ICON timer report at the end of the experiment
logfile, next to the Fortran timers; in standalone builds (e.g. unit tests)
`timer::report_all_timers()` prints them.

Include `support/timer/timer.hpp`.

## Defining a timer

Every timer object gets its own line in the single timer report — define
exactly one object per timer. The recommended place is namespace scope,
next to the code it measures: such a timer is constructed at program
start, in serial context, before any measurement can happen.

```cpp
timer::Timer mytimer("myname");     // namespace scope in a .cpp file
inline timer::Timer iotimer("io");  // in a header
```

A function-local `static` is an alternative where namespace scope is
inconvenient. Its first call constructs the timer, so that call must
happen outside parallel regions (see the rules below).

```cpp
void physics_step() {
  static timer::Timer t("physics_step");
  t.start();
  // ... work ...
  t.stop();
}
```

## Measuring

```cpp
mytimer.start();
// ... work ...
mytimer.stop();
```

Each start/stop pair adds one measurement; repeated pairs accumulate.
Running measurements on a thread must nest: stop the innermost running
timer first. Ill-formed sequences abort with an error message: double
start, stop without start, ill-nested stop, destroying a timer during its
measurement (standalone builds only).

`mytimer.value()` returns the calling thread's accumulated time in seconds,
including an ongoing measurement (without stopping it).

## Levels

A timer declares its activation threshold: the `timers_level` at which it
becomes active. The ICON build takes `timers_level` from the namelist. A
standalone program states it by calling `timer::configure(level)` once
after Kokkos is initialized; the test runner does so with the
`RAGNAROK_TIMERS_LEVEL` environment variable (default 1). Until then
`timers_level` is 0 and every timer is inactive. `timers_level` is a
runtime constant: a second configuration is a fatal error.

```cpp
timer::Timer coarse("phases");                       // default threshold_level 1
timer::Timer fine("detail", /*threshold_level=*/5);  // active from timers_level 5 on
```

An inactive timer's start/stop do nothing and it is absent from the
report. The default `threshold_level` 1 is the minimal valid choice and
mirrors a plain `IF (ltimer)` guard: in the ICON build such timers measure
whenever `ltimer` is on, and `ltimer = .false.` disables all timers.

## Timers around kernels

On GPU builds `index_loop` kernels run asynchronously to the host, so a
meaningful kernel measurement must synchronize. This requires
`timers_level >= timer::timers_level_sync_threshold`: from that level on,
every active timer synchronizes the device before recording. Define timers
for kernel measurements with a `threshold_level` of at least
`timer::timers_level_sync_threshold` — this ensures the measurement is
always synchronized when the timer is active.

```cpp
timer::Timer graupel_phase("graupel_phase", timer::timers_level_sync_threshold);
```

## Rules for parallel code

Measurements must stay portable across the loop layer's role swap
([icon_loops.md](../../kernels/docs/icon_loops.md)): a `block_loop` body
runs on the host in every build (OpenMP-parallel on CPU builds, serial on
GPU builds), whereas an `index_loop` body may become a GPU kernel.

- start/stop belong outside or inside a `block_loop` body, never inside an
  `index_loop` body — on GPU builds that is a compile error. Time an
  `index_loop` from the enclosing block body or around the whole loop nest.
- A measurement is per thread: the thread that starts it must stop it, and a
  measurement inside a `block_loop` body records on every block-loop thread
  (the report sums the counts and shows the totals per thread).
- Measurements inside a parallel region must run in the OpenMP team fixed
  at configuration: not in a nested region, and not in a region whose
  `num_threads` clause exceeds that team (undefined in `mo_real_timer` as
  well). Regions without measurements are unconstrained.
- OpenMP threads only: measuring from any other thread type
  (e.g. `std::thread`) is undefined behavior.
- Define timers and print reports outside parallel regions.

## Reports

- **ICON build**: nothing to do — the ICON timer report covers all timers.
- **Standalone / tests**: `timer::report_all_timers()` prints call counts
  and total times; timers used by several threads show the count summed and
  the total as `[min..max]` over the threads.
