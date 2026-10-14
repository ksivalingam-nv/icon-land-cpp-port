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

# Patch Reconstruction

How to extend the Fortran→C++ patch bridge in ICON-Ragnarok.

## Architecture

```
mo_model_domain.f90  →[gen_repr5]→  mo_ragnarok_f2c.f90 / icon_f2c.hpp  →  icon_domain.hpp / .cpp
```

`icon::Patch` exposes zero-based C++ views (`cells`, `edges`, `verts`).

Struct-typed fields (`TangentVec`, `GeographicCoords`) get per-member SoA companions
(`_v1`/`_v2`, `_lon`/`_lat`) extracted on the C++ side.

## View layout

Fortran is column-major (first index contiguous); Kokkos `LayoutRight` requires the last
index contiguous. `make_constview2d` / `make_constview3d` reverse the extent order to
alias the Fortran buffer without copying:

| Rank | Fortran shape | C++ extents | C++ access |
|------|--------------|-------------|------------|
| 2-D  | `(nproma, nblks)` | `(nblks, nproma)` | `view(jb, jc)` |
| 3-D  | `(nproma, nblks, nconn)` | `(nconn, nblks, nproma)` | `view(m, jb, jc)` |

`GridBase::block_size` is `extent(1)` (= `nproma`). The connectivity dimension is
leading for 3-D fields — the only layout-valid choice for the unmanaged `_f` views.
Managed and unmanaged views share the reversed layout and are indexed identically.

The `start/end_block/index` arrays keep their Fortran lower bounds (negative refin_ctrl
levels), transported via `lbounds` and represented as `OffsetView`s
(`make_index_offsetview1d`), with the stored values shifted to 0-based.

## Adding a field

1. `mo_ragnarok_f2c.f90` — add `TYPE(t_f2c_repr5) :: <name>` to `t_f2c_grid_*`
2. `icon_f2c.hpp` — add `FieldRepr <name>;` to `f2c::Grid*`
3. `mo_ragnarok_support.f90` — call `gen_repr5(patch%…%<field>)` in `get_f2c_grid_*`
4. `icon_domain.hpp` — add the view field(s) to `Grid*`
5. `icon_domain.cpp` — wire up the constructor with `make_constview*` / `make_index_*` / `make_component_view*`

For a new Fortran type: add a `gen_repr5` overload in `mo_ragnarok_f2c.f90` and a
`f2c_typekind<T>` specialisation in `icon_f2c_basic_types.hpp` / `element_types.hpp`.
