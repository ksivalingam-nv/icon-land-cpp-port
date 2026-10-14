# AGENTS.md

Guidance for AI coding agents working in the ICON-Land (JSBACH) repository. The [ICON-Land documentation](https://jsbach.gitlab-pages.dkrz.de/jsbach) (FORD-generated from inline comments plus `doc/guide/`) and the [GitLab wiki](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/home) hold the full technical/scientific docs. The wiki is the richest source; where it and the source tree disagree, prefer the source tree.

This directory is part of the reduced ICON tree it sits in (`externals/jsbach`); it is an ordinary directory there, not a submodule. It was vendored from the jsbach repository (`gitlab.dkrz.de:jsbach/jsbach`) at commit `03f30edbb` (`icon-land-2026.04-25-g03f30edbb`, i.e. release `icon-land-2026.04` plus 25 commits) and then reduced. Changes here are committed to the ICON tree; nothing is synced with the jsbach repository. See the parent `AGENTS.md` for ICON-wide guidance.

## What this is

ICON-Land — historically and still in repo/namelist naming "JSBACH" — is the land component of ICON: a Fortran framework for land-surface processes. **This working copy is a reduced version**: it contains only what the `jsbach_lite` usecase needs in *standalone* mode (forcing from files), i.e. the bio-geophysical processes energy/water balance, soil, snow, hydrology, radiation and climatological phenology. Its purpose is to serve as a small, bit-reproducible reference for a C++ port (see `ragnarok/` in the ICON tree).

It is still a **framework plus processes**, but the process set is fixed (see "Processes"). Removed relative to the full ICON-Land (do not look for it, do not re-add it without being asked):

- all QUINCY/IQJ code (`VEG_`, `SB_`, `Q_ASSIMI_`, `Q_PHENO_`, `Q_SYL_`, `Q_AGR_`, `MODEL_QUINCY`, BGCM classes/store, isotope utilities, QUINCY model config/output classes, `__QUINCY_STANDALONE__` blocks),
- the JSBACH processes `ASSIMI_`, `CARBON_`, `DISTURB_`/`FUEL_`/`FLCC_`/`WLCC_`, `CLIMBUF_`, `PPLCC_`, `FAGE_`, `ALCC_`, `NLCC_`,
- the land-cover-change (lcc) infrastructure (`mo_jsb_lcc*`, `mo_jsb_cqt_class`, conserved quantities),
- PFT tiles and the land-cover-type library (`mo_jsb_lctlib_class`, `data/lctlib_*.def`, `model%lctlib`),
- all usecases except `jsbach_lite`, the `_tmx` usecase variants and all `use_tmx` code (explicit coupling to the AES atmosphere),
- phenology schemes other than `climatology`, the soil-albedo scheme based on carbon pools, PAR/canopy-layer radiation,
- the hydrological discharge process `HD_` (`src/hd/`, river routing and its ocean coupling) and the `__NO_JSBACH_HD__` guards.

ICON-Land **cannot be compiled or run on its own**. It always needs a *host model* that provides I/O, domain decomposition, parallelization and time control — currently ICON, the only supported host, via the adapter file `src/shared/mo_util_jsbach.f90` on the ICON side. See "Adapters" below.

Building and testing happens through the ICON tree this directory is part of.

## Core concepts

You cannot usefully edit this code without these. They are all wiki-documented; the authoritative code lives in `src/base/`.

### Model → tiles → process memory

A **model** (`t_jsb_model`, `src/base/mo_jsb_model_class.f90`) owns a **hierarchical tile tree**, a list of **processes**, a **task queue**, and the **task schedule** derived from them. Tiles are nodes of the tile tree (`t_jsb_tile_tree_node`, `src/base/mo_jsb_tile_tree_class.f90`: `parent`, `first_child`, `next_sibling`, `Has_children`, ...). The model keeps all tiles in precomputed lists, `model%tiles(:)` (pre-order, parents first) and `model%tiles_postorder(:)` (children first), built once by `Build_tile_lists` in `jsbach_setup_tiles`; to loop over tiles, iterate over one of these (`DO itile=1,SIZE(model%tiles); tile => model%tiles(itile)%p`).

**Tiles** (`t_jsb_tile_abstract`, `src/base/mo_jsb_tile_class.f90`, `src/base/mo_jsb_tile.f90`) subdivide each grid cell by surface type: `box → {land, lake}`, `land → {veg, glac}` (there are no PFT tiles in this reduced version). Each tile carries a fraction `tile%fract`; children sum to their parent, and the `box` fraction is the land(+lake) fraction of the ICON grid cell. Each tile holds **its own memory per process**: `tile%mem(PROC_)%p`.

**Process memory** (`t_jsb_memory`, `src/base/mo_jsb_memory_class.f90`) is a list of variables of type `t_jsb_var_real1d` (1D - scalar/monitoring), `t_jsb_var_real2d` (2D - horizontal layer) and `t_jsb_var_real3d` (3D - horizontal and vertical layer), registered with `mem%Add_var(...)`. Note: Output variable names are composed as `PROC_variable_tile` so grepping for an output field name will fail; grep for the bare `variable` part inside the owning process directory.

### Processes

Each scientific process lives in its own `src/<name>/` directory and is usually identified by an enumerator ending in `_`, whereby some processes have subdivisions represented by additional enumerators. The enumerators are defined in `src/base/mo_jsb_process_class.f90`.

- Processes in this reduced version: `srf_energy_bal/` (`SEB_`), `turbulence/` (`TURB_`), `soil_snow_energy/` (`SSE_`), `hydrology/` (`HYDRO_`), `radiation/` (`RAD_`), `phenology/` (`PHENO_`, scheme `climatology` only), `atmland/` (`A2L_`, `L2A_`).

Conventional module set per process `PROC` (see `templates/` and `templates/README.md`. Note: the templates might not be up-to-date and copying a comparable existing process is the recommended route):

| module | required | role |
| --- | --- | --- |
| `mo_PROC_interface.f90` | yes | defines the process **tasks**: their `Integrate` (`update_*`) and `Aggregate` (`aggregate_*`) procedures |
| `mo_PROC_config_class.f90` | yes | namelist parsing (`jsb_PROC_nml`), the `active` switch, and creation of process-specific vertical grids |
| `mo_PROC_memory_class.f90` | if it has state | variable declarations + `Init_PROC_memory` with the `Add_var` calls |
| `mo_PROC_init.f90` | often | reads ic/bc files and sets initial/boundary conditions; called from `jsbach_init` |
| `mo_PROC_process.f90` | optional | the actual formulas, infrastructure-free |
| `mo_PROC_constants.f90` | optional | process parameters or constants, `PARAMETER` or `SAVE` + `init_PROC_constants()`;  |
| `mo_PROC_parameters.f90` | optional | process parameters — in addition to or instead of `mo_PROC_constants.f90`, `SAVE` + `init_PROC_parameters()` |
| `mo_PROC_util.f90` | optional | helper routines |

Some specifics:
- `A2L_` and `L2A_` share `src/atmland/mo_atmland_interface.f90`, have no config class, and their update routines are called directly from `mo_jsb_interface.f90` (`jsbach_interface`, the entry point for a coupled host), not via the task queue. The standalone driver does not call `jsbach_interface`: it fills the `A2L_` memory from the forcing and runs the task queue itself.
- `PHENO_` has no `mo_pheno_process.f90` any more; `mo_pheno_parameters.f90` only keeps `pheno_param_jsbach%EG_SG%tau_pseudo_soil` (used by `SEB_` init).

**Registering a new process** touches several infrastructure modules — miss one and it either crashes or silently does nothing:

1. add the ID + name mapping in `src/base/mo_jsb_process_class.f90` (`Get_process_name`, `Get_process_id`),
2. associate it with tiles in `src/base/mo_jsb_model_usecases.f90`,
3. queue its tasks in `jsbach_init` (`src/base/mo_jsb_model_init.f90`) — order matters, see below,
4. add interface + config class (+ memory class) to `src/base/mo_jsb_process_factory.f90` and bump `max_no_of_processes`; set `has_memory = .FALSE.` in `Create_process` if it has none,
5. call its `PROC_init` from `jsbach_init`,
6. call `init_PROC_constants()` and / or `init_PROC_parameters()` from `jsbach_setup_models` (`src/base/mo_jsb_base.f90`) if applicable.

### Process actions

Where a process runs on the tile tree (see associations in `src/base/mo_jsb_model_usecases.f90`) is an *action* per (process, tile), enumerated in `src/base/mo_jsb_process_class.f90`: `SKIP_`, `INHERIT_`, `AGGREGATE_`, `ON_LEAFS_`, `ON_TILE_`, `ON_SUBTREE_`. `ON_TILE_` on a tile gives descendants `INHERIT_` (their memory points at the parent's); `ON_LEAFS_` gives leaves `ON_TILE_` and intermediate tiles `AGGREGATE_`.

Query helpers (see FAQ): `model%Is_process_enabled` (configured at all), `tile%Is_process_active` (integrated *or* aggregated here), `tile%Is_process_calculated` (integrated here), `tile%Has_process_memory`. In addition, the config of a process `PROC` is often directly used to test if a process is active — using `dsl4jsb_Config(PROC_)%active` or (mostly in the infrastructure) `model%processes(PROC_)%p%config%active`.

### Task queue and schedule

The task queue (`src/base/mo_jsb_task_class.f90`, built in `jsbach_init` in `src/base/mo_jsb_model_init.f90`) fixes the order of `update_*` calls within a timestep. At the end of `jsbach_init`, `model%Build_task_schedule` turns it into the **task schedule** `model%schedule(:)`: per task, the ordered steps (tile, `TASK_INTEGRATE`/`TASK_AGGREGATE`) obtained by a depth-first walk of the tile tree — a tile with `SKIP_`/`INHERIT_` is skipped with its sub-tree, `ON_TILE_`/`ON_SUBTREE_` or a leaf is integrated, any other tile has its children processed first and is then aggregated (**leaves → top**). Each timestep, `model%Run_tasks` just iterates over these steps and calls `task%Do_it` (with OpenACC, consecutive integrate steps of a task run on separate async queues that are joined before the next aggregate step). Set `debug_level` ≥ 1 to print the schedule at initialization, ≥ 3 to trace every step.

There is a single queue (`model%config%model_scheme` is always `MODEL_JSBACH`), ordered by scientific dependency: `RAD_ surface_radiation`, `PHENO_ phenology`, `PHENO_ foliage_projected_cover`, `HYDRO_ snow_and_wet_fraction`, `HYDRO_ soil_properties`, `SEB_ surface_energy`, `HYDRO_ evaporation`, `SEB_ surface_fluxes`, `HYDRO_ surface_hydrology`, `HYDRO_ snow_and_wet_fraction`, `SSE_ soil_and_snow_properties`, `SSE_ soil_and_snow_temperature`, `HYDRO_ soil_hydrology`, `HYDRO_ canopy_cond_unstressed`, `HYDRO_ water_stress`, `HYDRO_ canopy_cond_stressed`, `TURB_ humidity_scaling`, `TURB_ surface_roughness`, `RAD_ albedo`, `HYDRO_ water_balance`, `SEB_ asselin_filter`, `SEB_ energy_balance`.

Tasks are queued regardless of the process's `active` namelist flag; the schedule drops tasks of inactive processes and tiles where `tile%Is_process_active` is false. Still, **every `update_`/`aggregate_` routine should start with `IF (.NOT. tile%Is_process_active(XYZ_)) RETURN`**.

### Aggregation

`Aggregate_weighted_by_fract_2d`/`_3d` (`src/base/mo_jsb_tile_class.f90`) average child values weighted by tile fraction. The `l_aggregate_all` flag on `Add_var` decides whether children *lacking* the variable still contribute their fraction to the denominator (`.TRUE.`, correct for fluxes that must be per whole grid box, e.g. what the atmosphere receives) or are ignored entirely (`.FALSE.`, default, correct for intensive quantities like temperature).

A variable only reaches `_box` output if it is actually listed in the process's `aggregate_` routine — otherwise `_box` output will stay the initial value `initval_r` as specified in the `Add_var` call in the process memory module — which is often but not always zero. Note: technically it is possible to only add variables for certain tiles — aggregation routines can of course only be called for variables added along the tile tree.

### Usecases

`src/base/mo_jsb_model_usecases.f90` — selected via `usecase` in `jsb_model_nml`. Only `jsbach_lite` exists (minimal physics-only setup); any other value stops the model. It defines the tile hierarchy `box → {land, lake}`, `land → {veg, glac}` and the process/action association.

Tile construction passes name, `lct` (land cover type, `src/base/mo_jsb_lct_class.f90`), parent, `fract_varname` (which field in the fractions file, `bc_land_frac.nc`, supplies the initial fraction; `'equal'` means "split evenly, read nothing"), and the process/action arrays. Processes can also be attached conditionally afterwards via `tile%Set_process_action(...)`. `t_jsb_lct%lib_id` (formerly the index into the lctlib) and the `lib_id` argument of the `Init_PROC_memory` routines are still there but always `0`.

### Decomposition and grids

Horizontal: the global cell vector is split across MPI domains, each domain into blocks of `nproma` cells. 2D variables are `var(nproma, nblks)`, 3D are `var(nproma, nlev, nblks)`. Inside a task you work on one chunk, addressed by `options%iblk`, `options%ics`, `options%ice` — i.e. `var(options%ics:options%ice, options%iblk)`.

Vertical: `t_jsb_vgrid` (`src/base/mo_jsb_grid_class.f90`), created with `new_vgrid(...)` in a process config class, published with `Register_vgrid`, retrieved anywhere with `Get_vgrid('name')`. Convention: **`ubounds > lbounds` always**, and layer 1 is always the one at the ground with `lbounds(1) = 0.0_wp`, independent of whether the axis points down (soil) or up (canopy).

### Adapters (host-model abstraction)

ICON-Land should not `USE` host-model modules directly. Instead the host provides a fixed set of generic interface modules — `mo_jsb_parallel_iface`, `mo_jsb_domain_iface`, `mo_jsb_namelist_iface`, `mo_read_namelists_iface`, `mo_jsb_time_iface`, `mo_jsb_io_iface`, `mo_jsb_io_netcdf_iface`, `mo_jsb_grid_iface`, `mo_jsb_utils_iface`, `mo_jsb_varlist_iface`, `mo_jsb_vertical_axes_iface`, `mo_physical_constants_iface` — all implemented in a single host-side file, `src/shared/mo_util_jsbach.f90` in ICON. ICON-Land *requests* a feature through the generic interface; the host decides how to provide it.

**Rule: only files in `src/base/` may `USE` the adapter modules.** Everything outside `base` must go through a `base` module. If you need new host functionality you must define the host-agnostic interface here *and* implement it in the host adapter file (same commit — see "Repository layout").

## The DSL (`dsl4jsb`)

Because everything is an object behind pointers, process code is written against a small macro language, expanded by `scripts/dsl4jsb/dsl4jsb.py` (regex substitution) **before compilation**.

```fortran
dsl4jsb_Use_processes HYDRO_
dsl4jsb_Use_memory(HYDRO_)
dsl4jsb_Def_memory(HYDRO_)
dsl4jsb_Real2D_onChunk :: evapotrans
dsl4jsb_Get_memory(HYDRO_)
dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)
```

expands to `USE` statements, a `TYPE(t_hydro_memory), POINTER :: hydro__mem`, a `SELECT TYPE` on `tile%mem(hydro_)%p`, and `evapotrans => hydro__mem%evapotrans%ptr(ics:ice,iblk)`.

The macros assume that the names `tile`, `options`, `iblk`, `ics`, `ice` exist and are assigned, which can lead to segmentation faults at runtime if this is not the case. In a task routine they therefore need to be declared and set from `options` first. `_tile` variants (`dsl4jsb_Get_memory_tile`, `dsl4jsb_Get_var2D_onChunk_tile`, ...) let you name a different tile. `_onDomain` variants point at the whole array and are used in `mo_PROC_init` modules, which run before decomposition on the full domain.

Full macro inventory: grep `dsl4jsb_` in `scripts/dsl4jsb/dsl4jsb.py`. Families are `Use_*`/`Def_*`/`Get_*` for processes, memory, config and plain variables (`var2D/var3D` × `onChunk/onDomain`). The script is unchanged and still knows the macros for lctlib parameters, lcc relocations and BGCM matrices, but nothing they expand to exists any more: do not use them.

**Consequences you will hit:**

- Preprocessed files are `<name>.pp-jsb.f90` under `<build-dir>/externals/jsbach/src/...`. **Never edit those** — edits are lost on the next build.
- Compiler and runtime **line numbers are those of the original `.f90`**: `dsl4jsb.py` writes a `#line 1 "<original file>"` directive after its 5-line header, so the offset is compensated.
- If a rebuild seems not to pick up your change, check that the corresponding `.pp-jsb.f90` was regenerated. `make clean` in the ICON root removes the preprocessed files and fixes stale-DSL problems.

## Source layout (`src/`)

Infrastructure:

- `base/` — the framework: model/tile/process/task/memory/variable classes, tile tree (`mo_jsb_tile_tree_class.f90`), usecases, task-queue and schedule construction, grids, time, I/O, parallel, `mo_jsb_version.f90`
- `interfaces/` — `mo_jsb_interface.f90`, the entry points the host calls: `jsbach_interface` (coupled host only, unused in standalone runs), `jsbach_start_timestep`, `jsbach_finish_timestep`, `jsbach_get_var`, `invalidate_cuda_graphs`
- `drivers/` — `mo_jsbach_model.f90` (standalone driver, invoked from ICON's `src/drivers/icon.f90`; also initializes/finalizes ragnarok inside `#ifndef __NO_RAGNAROK__`) and `mo_jsb4_forcing.f90` (offline forcing)
- `shared/` — cross-process helpers: physical/math constants (still containing constants of removed code), `mo_phy_schemes.f90`, `mo_orbit_solar.f90`, OpenACC helpers

Process directories: see the inventory under "Processes" above.

Other top-level directories: `doc/` (FORD guide pages, media, CSS), `scripts/` (DSL, CI checks, preprocessing, restart helpers), `templates/` (new-process skeletons), `LICENSES/`. These were not reduced and still describe the full ICON-Land.

## Configuration (namelists)

- `jsb_control_nml` and one `jsb_model_nml` **per model instance** live in the *master* namelist file: technical settings (model id/name/namelist filename, debug and timer levels, `is_standalone`, `restart_jsbach`), parsed in `mo_jsb_control.f90`.
- A second, differently-scoped `jsb_model_nml` in the *model* namelist file carries the scientific configuration (`model_scheme` (only `jsbach`), `usecase` (only `jsbach_lite`), `npft` (without effect, 0), `use_lakes`, `use_glacier`, ...), as does `jsb_grid_nml`; grep `mo_jsb_config_class.f90` for the full list. `use_tmx`, `quincy_model` and the QUINCY element/spinup switches are gone: a namelist that still sets them fails on read.
- Per-process namelists: **`jsb_<proc>_nml`** (`jsb_seb_nml`, `jsb_turb_nml`, `jsb_sse_nml`, `jsb_hydro_nml`, `jsb_rad_nml`, `jsb_pheno_nml`). Each carries at least `active`, plus `bc_filename`/`ic_filename`. `jsb_pheno_nml` must set `scheme = 'climatology'` (the code default is still `'logrop'`, which no longer exists). Namelist groups of removed processes are simply not read.
- In coupled runs ICON-Land output namelists belong to `NAMELIST_atm`, not `NAMELIST_lnd` as in standalone runs.

`__NO_JSBACH__` guards the whole component; other cpp keys in use are `_OPENACC`, `_OPENMP`, `HAVE_CDI_PIO`, `HAVE_F2003`, `HAVE_FC_ATTRIBUTE_CONTIGUOUS`. Keep preprocessor use to technical concerns only, and mind capitalization.

## Build and run (via the host model; here for levante)

From an ICON(-MPIM) clone — ICON-Land is built as part of ICON, there is no build system in ICON-Land. The `config/dkrz/levante.*` wrappers pin their whole toolchain by absolute path, so **no module needs to be loaded to build**:

```bash
mkdir build && cd build
../config/dkrz/levante.intel     # or another config/dkrz/levante.* wrapper
make -j32
```

`--enable-jsbach` is on by default; `--enable-jsbach-hd` still exists in ICON's `configure` but has no effect, since the HD process is removed. ICON's `icon.mk.in` runs `scripts/dsl4jsb/dsl4jsb.py` over every `externals/jsbach/src/**/*.f90` before compiling. If running ICON-Land standalone, i.e. without ICON atmosphere and with forcing from a file, more code can be disabled (`--disable-les --disable-upatmo --disable-ocean --disable-coupling --disable-yaxt --disable-atmo --disable-comin --disable-art --disable-waves`).

Running:

- **mkexp** (preferred for real experiments): work in ICON's `run/` directory, copy an example config, then run mkexp *by path* setting `ACCOUNT=<account-id>` where `account-id` is the account to be used for accounting compute time, if required, e.g. on levante. In many cases, `account-id` will be your "effective group ID" (the first listed in the terminal with command `groups`). Then submit the generated `run_start`:

  ```bash
  cd run
  cp examples/JSB-standalone-R2B4.config myexp.config
  ../utils/mkexp/mkexp myexp.config ACCOUNT=<account-id>
  cd ../experiments/myexp/scripts && sbatch myexp.run_start
  ```

  Nothing needs to be on `PATH` and no module is needed. For an out-of-source build, `export MKEXP_PATH=.:<build-dir>/run` so mkexp finds the build's `run/SETUP.config`.

  The ICON-Land example is `run/examples/JSB-standalone-R2B4.config`, with the short test variant `run/checksuite.icon-land/test_jsb_standalone.config`; both select `LAND_USECASE = jsbach_lite`. Their `EXP_TYPE` resolves under `run/mkexp/types/` (`land/icon-land-standalone-R02B04.config`).

- The `exp.*` runscript templates of the full ICON-Land (`exp.land_jsbach_R2B4_test`, `exp.land_quincy_canopy_R2B4_test`, ...) do not apply to this reduced version.

## Documentation

The [ICON-Land documentation](https://jsbach.gitlab-pages.dkrz.de/jsbach) is generated by FORD from the `!>`/`!<` inline comments plus the guide pages in `doc/guide/`, and published by the scheduled `pages` job — which runs only on the default branch, from a schedule, never on a merge request.

Build it locally with `./build-doc.sh`: `-f` skips the search index and is much faster, `-c` removes the output again. The result is `doc/html/index.html` (`output_dir` in `project.md`).

**The script needs outbound network access** — it bootstraps `uv` from `astral.sh` into `3rdparty/` and then installs a FORD fork from GitHub, so it cannot run offline or in a sandbox. `graphviz` is needed in addition if `project.md` sets `graph: true` (currently `false`).

Nothing verifies the rendered result automatically: the `Test build doc` CI job is `manual`. So when you can run the build, do, and check that the annotations you added come out as intended — FORD markup mistakes stay invisible until someone renders the pages.

## Conventions

- **License header**: every source file starts with a one-line module summary, then the ICON-Land block using `!>` comment markers — `ICON-Land`, `Copyright (C) 2013-2026, MPI-M, MPI-BGC`, `Contact: icon-model.org`, `Authors: AUTHORS.md`, `See LICENSES/ for license information`, `SPDX-License-Identifier: BSD-3-Clause`. **This differs from ICON's own header** (different years, institutions and comment marker) — copy from a neighbouring file in this repository, not from ICON. See `templates/module_header_with_BSD-3-license.template`; enforced by `scripts/ci_scripts/check_license.sh`.
- **Model name in prose**: always write "ICON-Land" — never "ICON-L" or "ICON-LAND" — in documentation, wiki pages, commit messages and merge-request text. (Code identifiers, paths and namelist names keep their historical "jsbach" spelling.)
- **File hygiene**: ASCII only, no non-UTF-8 bytes, no trailing blanks (`scripts/ci_scripts/check_file_format.sh`). Fortran line limit 132 chars, but keep to ≤110 for readable FORD pages. Two-space indentation, no tabs.
- **Style**: spaces around operators; real literals as `_wp` (`2._wp`); cast integers explicitly (`REAL(iexp, wp)`) rather than mixing types; `==`/`<`/`<=` rather than `.EQ.`/`.LT.`/`.LE.`.
- **Naming**: variables `small_caps`; 3D variables carry a suffix naming their third dimension (`_cl` canopy layer, `_sl` soil layer, `_snl` snow layer); enumerators `ALL_CAPS` (`ASSIMI_`, `MODEL_QUINCY`); type-bound procedures `Capitalised_small_caps`, usually verb-first (`Init_vars`, `Configure_model_processes`).
- **Argument order** for scientific routines: dimensions → lctlib parameters → config variables → `in` → `inout` → `out`, with bgc_material structures after the plain 2D/3D variables in each group.
- **Inline documentation is part of the deliverable.** FORD uses `!>` before a code line and `!<` after it. Every memory variable and every constant needs a long name and `[unit]` at declaration; every constant needs a literature reference (or an explicit "tuned against X") at its definition. See "Documentation" for checking how it renders.
- **Only `update_*` (task Integrate) routines should touch model infrastructure**; where possible keep formulas in infrastructure-free routines so they stay testable and portable.

### GPU (OpenACC)

GPU support is OpenACC. `Add_var` already emits the `!$ACC ENTER DATA COPYIN` / `!$ACC UPDATE DEVICE` for registered variables — anything not created via `Add_var` must be placed on the device explicitly. No array/`WHERE` notation inside accelerated regions: write loops out with `!$ACC PARALLEL DEFAULT(PRESENT)` + `!$ACC LOOP GANG VECTOR COLLAPSE(n)`, cell index innermost. Never issue `!$ACC UPDATE HOST/DEVICE` inside the time loop — it is very expensive. Copy an existing accelerated routine rather than inventing a pattern.

## CI and testing

The `.gitlab-ci.yml` in this directory is the upstream jsbach one and is inert here. The two style checks it ran still work and should stay clean:

```bash
./scripts/ci_scripts/check_file_format.sh
./scripts/ci_scripts/check_license.sh
```

Nothing here builds or runs on its own; building and regression testing go through the ICON tree (see the parent `AGENTS.md`).

Regression test for this reduced version: build ICON (`config/dkrz/levante.intel --enable-openmp`), then run `run/checksuite.icon-land/test_jsb_standalone.config` with `EXP_OPTIONS+=serial,test/suite` (base run plus restart, OpenMP, nproma and MPI variants; 12 `ok` lines in the job log) and compare every output and restart file against the reference experiment with `cdo -s diffn` (must print nothing). Every change must keep the results bit-identical unless explicitly agreed otherwise.

## Repository layout

There is one repository: the reduced ICON tree, which contains this directory. The upstream workflow (jsbach repository with `land4icon-mpim`/`master`, paired merge requests with `icon-mpim`, submodule pointer updates) does not apply here. A change that needs a matching change in the host adapter `src/shared/mo_util_jsbach.f90` or in `ragnarok/` goes into the same commit.

## Known pitfalls

- **Reading another process's memory** is technically allowed and often wrong. Reading a variable from a process that runs further down the hierarchy may hand you the initialisation value if the variable is not in that process's `aggregate_` routine. Writing another process's variables can be silently discarded (overwritten by a later aggregation) or leak to siblings and parent through `INHERIT_` pointers. Do it only in exceptional, documented cases.
- `l_aggregate_all` does **not** mean "aggregate over all child tiles" in the naive sense — see Aggregation above.
- `isteptype=TSTEP_CONSTANT` variables appear to require an output file of their own; mixing them with time-dependent variables in one ICON output file produces corrupted fields.
- Stale preprocessed files: if a rebuild does not reflect your change, see the DSL consequences above.
- Checking out a jsbach branch that does not match ICON's `src/shared/mo_util_jsbach.f90` produces confusing compile errors — verify the adapter implementations match the ICON-Land version before debugging further.
- Wiki drift: the wiki was largely re-synced with the source tree in August 2026; known remaining drift: the usecases page misses `jsbach_forest_age_classes` and the `_tmx` variants. MPI-ESM/ECHAM as a second host is described *on purpose*, marked as historical — older ICON-Land versions can be run with it; do not "fix" those sections away.

## Commit messages

Do not add a `Co-Authored-By: Claude` (or other AI) trailer to commits. This repository feeds into ICON, whose contribution process gates merge requests on contributor-agreement status, and no CLA covering AI co-authorship is in place.

When a commit message references a merge request or an issue, write the **full URL** — `https://gitlab.dkrz.de/jsbach/jsbach/-/merge_requests/450` — not GitLab's `project!number` shorthand. Commit messages are read with `git log` and in mail, outside GitLab, where the shorthand resolves to nothing. In merge request descriptions the shorthand is the better choice and should be used instead.
