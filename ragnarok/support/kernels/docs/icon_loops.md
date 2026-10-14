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

# icon_loops: writing ICON loop nests in C++

How to write the block/index loop nests of the AES physics with
`icon_loops.hpp`. Working examples: the tests next door
(`kernels_test.cpp`) and the testbed (`ragnarok_testbed.cpp`).

## One loop nest, two parallelization modes

ICON runs every block/index loop nest in one of two modes:

- **CPU:** parallel over blocks `jb` (OpenMP), serial over indices `jc`
- **GPU:** serial over blocks (ideally one huge block), parallel over `jc`

In Fortran both modes are spelled out with directives (`!$OMP` outside,
`!$ACC` inside). In Ragnarok you write the nest once, and the **build**
decides which loop is parallel:

```cpp
block_loop(grid_range, block_lambda(const int jb) {
  index_loop(grid_range, jb, index_lambda(const int jc) {
    d_out(jb, jc) = 2.0 * d_in(jb, jc);
  });
});
```

|              | CPU build         | GPU build              |
| ------------ | ----------------- | ---------------------- |
| `block_loop` | parallel (OpenMP) | serial (host)          |
| `index_loop` | serial (host)     | parallel (GPU kernel)  |

What this means for the body you write:

- **The body must be correct in both modes.** Whatever runs inside may not
  assume a particular loop is the parallel one. Order-dependent bodies
  (scans, running sums) belong in a `serial_loop`, below.
- **There is no nested parallelism** — ICON parallelizes over blocks *or*
  over the indices of one block, never both. The API cannot express the
  illegal combination, so you cannot get it wrong.
- **The loops promise ordering, not completion.** Everything runs in
  submission order on one queue, so a later loop sees the results of an
  earlier one without any waiting on your part. But on GPU builds kernels
  may still be running when a loop returns — *host* code may only touch
  their results through the verbs of `kernel_sync.hpp`, below. This is the
  Fortran `ASYNC(1)`/`WAIT` scheme: launch freely, wait only at consumption
  points.

Every loop takes an optional label as first argument
(`block_loop("my_label", grid_range, …)`); the label shows up in profiler
output for the loop that is parallel in the current build. Without a label
a generic default is used.

## Collapsed level×index loop

Most of the physics runs the index loop *collapsed* with the vertical level
loop. Pass a `LevelRange` and take two indices, with `jc` — the stride-1
dimension of the views — in the last slot:

```cpp
block_loop(grid_range, block_lambda(const int jb) {
  index_loop(grid_range, jb, LevelRange(kstart, kend),
             index_lambda(const int jk, const int jc) {
    d_out(jb, jk, jc) = d_in(jb, jk, jc) + fac(jk);
  });
});
```

`LevelRange` is a half-open range `[kstart, kend)`: level loops often start
below the model top (nonzero `kstart`). On GPU builds the collapsed loop
runs as one kernel per call with `jc` as the fastest index — matching the
memory layout, so access stays coalesced.

There is also a block-free form for callers whose block loop lives
elsewhere: `index_loop(levels, indices, body)` with an `IndexRange` for the
horizontal indices.

A third, grid-free dimension can be collapsed in front (`COLLAPSE(3)`
analog) — the typical shape of a tracer loop:

```cpp
index_loop(Range(0, ntracer), levels, indices,
           index_lambda(const int i, const int jk, const int jc) {
  zqx(i, jk, jc) = qx(i, jk, jc);
});
```

The extra dimension is the slowest one (the leftmost view extent) and has
no grid meaning, so it is spelled as a plain `Range` — deliberately: the
tagged ranges (`LevelRange`, `IndexRange`) do not convert into this slot,
so loop dimensions cannot trade places silently.

The same works one rank lower, when the slow dimension of a 2D collapse is
not levels (e.g. initializing per-column prognostics):

```cpp
index_loop(Range(0, prog_idx::count), indices,
           index_lambda(const int p, const int jc) {
  d_qp(p, jc) = 0.0;
});
```

`LevelRange` in a collapsed signature therefore always means levels; a slow
dimension that is not levels is always a plain `Range`.

## Reading results back (`kernel_sync.hpp`)

On GPU builds `index_loop` kernels run asynchronously to the host (among
themselves they execute in launch order); on CPU builds a block body may
execute inside an OpenMP parallel region where no Kokkos synchronization is
allowed. `kernel_sync.hpp` provides verbs that do the
right thing in both modes — use them instead of any hand-written fence or
`deep_copy`:

- `kernel_fence()` — completion point for all previously launched kernels.
  Call it before host code reads their results; bridge entry points end
  with it before returning to Fortran.
- `fetch_scalar(counter)` — read back a single kernel-computed value, held
  in a rank-0 view (e.g. an atomic compaction counter), inside a block
  body — e.g. to size a follow-up `index_loop`.
- `fetch(dst, src)` — bring the *same logical data* across the mirror pair
  `dst`/`src` (staging pattern). A copy between *distinct* fields is
  compute, not transfer: spell it as an `index_loop` with an assignment
  body.

The rule behind all of this: **never write `prefers_coalesced_access()` (or
any build-mode branch) in physics code.** If a body seems to need one, a
verb is missing — it belongs in `kernel_sync.hpp` or `icon_loops.hpp`.

## `serial_loop`: CPU-first porting

Freshly ported code often works on host data only and may depend on
iteration order (scans, running sums). Such code must still build and run
correctly in a GPU build — merely without using the GPU. `serial_loop` is
the loop for this stage: no parallelism at all, in any build — a plain host
loop over the same iteration spaces as the switchable loops:

```cpp
double sum = 0.0;
serial_loop(grid_range, [&](const int jb, const int jc) {
  sum += h_in(jb, jc);
  h_out(jb, jc) = sum;
});
```

- Runs on the host, single-threaded, in every build, iterating in ascending
  order with the rightmost body argument fastest (`jb` outer, `jc` inner in
  the example above). The device stays uninvolved; host-only data is the
  point, not a limitation.
- Takes any C++ lambda — `[&]`, `[=]`, whatever the body needs. There is no
  capture marker (nothing build-dependent to hide) and no label (the loop is
  never parallel, never profiled). Scalar carries like `sum` above are
  ordinary C++.
- Forms: per-block `serial_loop(grid_range, jb, body(jc))` — the usual shape
  inside a `block_loop`, so blocks still run OpenMP-parallel on CPU builds
  while the ported body stays serial; whole-range
  `serial_loop(grid_range, body(jb, jc))` for scans and fully serial nests;
  a `LevelRange` argument adds the collapsed level dimension to either; and
  block-free forms over plain ranges.
- The body must not launch kernels — a serial loop *around* kernels is a
  `block_loop` in its device role, not a `serial_loop`.
- `grep -rn serial_loop` lists every not-yet-ported loop.

Porting a `serial_loop` to the GPU is a mechanical recipe: rename it to
`index_loop`, spell the lambda `index_lambda`, and give the data a device
home — the compiler then flags every leftover scalar carry and host-only
call in the body, so nothing ships silently broken.

## Capture markers

The switchable loops take their lambdas through markers:

| loop         | marker         |
| ------------ | -------------- |
| `block_loop` | `block_lambda` |
| `index_loop` | `index_lambda` |

Both capture by value (`index_lambda` additionally carries the device
annotation a GPU kernel needs), so enclosing scalars are read-only copies
inside these bodies: `sum += x` does not compile. That is deliberate — in a
loop that may run in parallel, a scalar update would race or be lost.
Order-dependent state belongs in a `serial_loop`, whose plain lambdas are
free to capture by reference. Inside member functions, don't touch data
members from a `block_lambda`/`index_lambda` body — copy them into locals
first (a by-value lambda captures the object *pointer*, which a GPU cannot
follow).

## Naming rationale

The switchable loops are named for *what they iterate* — blocks, indices —
not for *where they sit* in the nest. Positional names (`outer_loop`/
`inner_loop`) were considered and rejected; the choice was deliberate, so
here is the reasoning.

- **Block and index are ICON's own words.** Fortran ICON speaks in blocks
  and indices: loop nests run over `jb` and `jc`, `get_indices_c` yields
  the index bounds, `mo_omp_block_loop` parallelizes block loops. Porting
  therefore maps one to one — an `!$OMP` loop over `jb` becomes a
  `block_loop`, a `jc` loop under `!$ACC` becomes an `index_loop` — with
  no new vocabulary to learn.
- **`index_loop` is named for its defining dimension, not its arity.** The
  collapsed forms iterate more than indices, but every one of them keeps
  `jc` as the fastest dimension — the one that is device-parallel and
  coalesced. That invariant *is* the concept; collapsed level or extra
  dimensions ride along and do not change what the loop is.
- **Positional names do not survive the block-free forms.** "Inner" only
  means something while an outer loop is in view, but
  `index_loop(levels, indices, body)` is legal with the block loop on the
  Fortran side — there is no outer loop in the same code to be inner to.
  `index_loop` describes itself wherever it appears.
- **Positional names suggest nested parallelism.** In the GPU world,
  outer/inner is the vocabulary of hierarchical parallelism — teams and
  vectors, both parallel at once. ICON's design point is the opposite:
  exactly one loop is parallel per build, the other is plain serial.
  Naming the pair by their nesting relation would emphasize the one
  property that is deliberately not exploited.

## GridRange

`GridRange` (`grid_range.hpp`) describes the iteration space the loops run
over: the blocks of a range of cells, edges or vertices, and the indices
within the trimmed first and last block. Physics code receives it ready-made
from the bridge and passes it through to the loops. For hand-written
reference loops (tests), `block_range()` and `index_range(jb)` provide the
same bounds as half-open ranges with `start` and `end` members.
