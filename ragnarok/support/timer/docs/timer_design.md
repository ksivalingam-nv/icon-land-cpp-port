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

# Timer design

High-level view of the implementation. Usage is in
[timer_usage.md](timer_usage.md).

## One Timer, two backends

`Timer` dispatches through a small table of function pointers. An internal
wall-clock backend is in place from program start. In the ICON build the
bridge switches once to the ICON Fortran timers; test builds keep the
internal backend and set the level via `timer::configure()`. Exactly one
of these happens per program; a second configuration aborts. Until it
happens the level is 0 and all timers are off, so every measurement is
governed by the one configuration and belongs to the one backend it ran
on.

## The pool

Ids are backend keys, so the switch re-issues them for all live timers. That
is the pool's only purpose.

## Levels

Per timer: active iff `timers_level >= threshold_level`. Fortran puts this
condition at each call site, ragnarok at construction. Global: from
`timers_level_sync_threshold` on, every measurement syncs the device first,
like `mo_real_timer`'s `!$ACC WAIT` above level 10.

## Per-thread statistics

The internal backend mirrors `mo_real_timer`'s `threadprivate` state: each
thread accumulates its own values, nesting is enforced per thread, only the
report aggregates. Like `mo_real_timer` it assumes one fixed OpenMP team:
`configure()` creates one statistics table per team rank up front, so a
measuring thread only indexes its own table and needs no locks. Kokkos
OpenMP kernels are ordinary parallel regions and fit this model.

After the switch the same holds for the ICON backend, whose `threadprivate`
state `new_timer` initializes for the Fortran OpenMP team. Kokkos host
kernels must therefore run on that team, which is the case when Fortran and
C++ share one OpenMP runtime. The testbed checks this property: it measures
on every thread of a Kokkos host kernel and reads the counts back per
thread from a Fortran parallel region.
