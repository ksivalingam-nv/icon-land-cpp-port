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

# Scratch memory: temporary work arrays

A physics routine often needs temporary arrays that live only for one call.
Allocating memory on every call can be slow, so ragnarok hands out these work
arrays from one big pre-allocated buffer instead. You get ordinary Kokkos
views, and there is nothing to free — everything is released automatically at
the end of the routine.

## How to use it

Three steps, taken from `aes_microphysics/graupel.ipp`:

**1. Open a scratch space** at the top of your routine, as a plain local variable:

```c++
ragnarok::LocalScratchSpace<ExecSpace> local_space = ragnarok::get_local_space<ExecSpace>();
```

**2. Create views** from it. The element type and rank come from the left-hand
side; the dimensions are the arguments:

```c++
View0D<int, ExecSpace> d_jmx    = ragnarok::scratch_view(local_space);
View1D<int, ExecSpace> index_iv = ragnarok::scratch_view(local_space, nvol);
View2D<T, ExecSpace> d_vt       = ragnarok::scratch_view(local_space, idx::np, nvec);
```

Use the `View0D`/`View1D`/`View2D` aliases from `common/types.hpp` — they are
the right view types for scratch memory.

**3. Do nothing at the end.** When `local_space` goes out of scope, all the
memory it handed out is released at once.

## Rules

- **Don't write `auto`** on the left of `scratch_view(...)` — always spell out
  the view type. (An `auto` variable would not be a view; using it fails to
  compile.)
- **Scratch memory starts with undefined content.** Initialize before reading,
  e.g. `Kokkos::deep_copy(execSpace, d_vt, ZERO<T>);` in graupel.
- **Keep `local_space` local.** Don't store it in a class member, don't return
  it. The views die with it, so don't keep them past the routine either.
- **Nesting is fine.** A routine you call may open its own scratch space.
- Mismatches (wrong view type, wrong number of dimensions, wrong memory space,
  const or non-trivial element type) are caught by the compiler.

## Good to know

The buffer starts at 64 MiB and grows automatically to whatever your run
actually needs, so an "out of scratch memory" situation cannot occur. ICON
checks once per timestep that everything was released; a physics routine never
manages the buffer itself. The initial size can be overridden with the
environment variable `RAGNAROK_SCRATCH_INIT_SIZE` (plain bytes) or with
`ragnarok::set_scratch_buffer_initial_size()`. To watch the memory usage,
compile with `MEMORY_POOL_DEBUG` (see `memory_pool.hpp`).
