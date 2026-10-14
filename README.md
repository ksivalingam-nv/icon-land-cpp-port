# ICON-Land C++ port: reduced ICON scaffold and call tree

This public repository combines the current C++ port scaffold with source-derived HTML call-tree documentation. The application source is a snapshot of [the upstream reduced repository](https://gitlab.dkrz.de/m212005/icon-land-cpp-port) at commit `133173041863d2879900563b7a19476ec75d5b1c`. Original source files, executable permissions, symbolic links, license notices, and the seven dependency pins are preserved.

The C++ land-model port is at the preparation stage: the retained land physics and model/task framework are Fortran; `ragnarok/` supplies C++/Kokkos infrastructure and existing thermodynamics helpers.

- [Open the interactive call tree](https://ksivalingam-nv.github.io/icon-land-cpp-port/).
- [HTML file](icon-land-call-tree.html) and its [companion data file](icon-land-call-tree-data.js). Keep these together when opening locally. Native navigation remains available if JavaScript is blocked.
- [Call-tree generator](generate-icon-land-call-tree.py), [template](icon-land-call-tree-template.html), and [browser checks](check-icon-land-call-tree.cjs).

Regenerate the call tree from this checkout without compiling the application:

```bash
python3 generate-icon-land-call-tree.py --repository-url https://github.com/ksivalingam-nv/icon-land-cpp-port
```

The analysis records static calls, lexical function references, selected type-bound bindings, the 22 queued task entries, and the manually traced C++ bridge. It is not a runtime profile or a complete compiler-resolved call graph. Tracing limitations are included in the HTML.

The dependency submodules remain pinned to their upstream repositories. Their URLs are absolute so cloning from GitHub does not redirect the GitLab dependencies to GitHub. Access to those repositories is controlled by their respective hosts. The inherited GitLab CI configuration has not been exercised on GitHub. No application build or simulation was performed for this publication.

The following scaffold description and validation claims come from the upstream README.

This repository is a stripped-down copy of [ICON](http://icon-model.org) that contains only what is needed to run the land model ICON-Land (JSBACH) standalone, together with the C++/Kokkos framework `ragnarok`. Its purpose is the port of ICON-Land from Fortran to C++, starting with the technical infrastructure (model, tiles, processes and tasks, memory). A small code base makes that port simpler and faster to build and test.

It is a working scaffold, not a development line of ICON: it cannot be merged back into or synced with `icon-mpim` or the jsbach repository.

## What was done

Starting from `icon-mpim` master at `019924437b` and ICON-Land at jsbach `03f30edbb` (`icon-land-2026.04` plus 25 commits):

- **ICON sources** reduced to the dependency closure of the standalone land executable: about 226 of 949 files in `src/` and `support/`, roughly 164,000 of 860,000 lines. Atmosphere, ocean, waves, coupling and the NWP/AES physics are gone.
- **ICON-Land** moved from a submodule into this repository as the ordinary directory `externals/jsbach`, and reduced to the `jsbach_lite` usecase: surface energy balance, turbulence, soil and snow energy, hydrology, radiation, climatological phenology and the atmosphere/land interface. Carbon, assimilation, disturbance, land-cover change, forest age, QUINCY, PFTs and the HD river routing were removed (203 files and 121,000 lines down to 74 files and 40,000 lines).
- **ragnarok** kept without its microphysics package. Its ICON bridge is cut down to patch access, halo synchronisation, timers and logging, and ragnarok is initialised from the standalone land driver.
- **Everything else** pruned to match: externals (only cdi, mtime, fortran-support, iconmath, memman and kokkos remain), run setups, scripts, documentation, build wrappers (only DKRZ Levante) and CI (ragnarok unit jobs and the land standalone test).
- **MemMan** (`externals/memman`, release 2.2.0) is built unconditionally and allocates all var-list variables, including every ICON-Land variable added with `Add_var`; taken from the `icon-mpim` branch `feature-enable-memman-unconditionally`.
- **Build system** left as it is apart from the defaults: a plain `configure` now gives the land build with ragnarok. Options and preprocessor branches of removed components still exist but no longer work.

## Verification

Every reduction step reproduced the output and restart files of the standalone land test bit-identically with the Intel compiler on Levante, against a reference run made before the step. The final tree also builds and passes the test with NAG (without ragnarok, which NAG cannot compile), NVHPC on CPU and NVHPC on GPU. The pruned GitLab CI has not been run.

## Build and test

```bash
git submodule update --init
mkdir -p build/intel && cd build/intel
../../config/dkrz/levante.intel --enable-openmp
make -j24
make -j16 check-icon
```

The regression test is the standalone land experiment `run/checksuite.icon-land/test_jsb_standalone.config` (R2B4, 3 hours, with restart, OpenMP, nproma and MPI consistency runs). [AGENTS.md](./AGENTS.md) describes how to run it, how to compare results, and how the tree is organised; [externals/jsbach/AGENTS.md](./externals/jsbach/AGENTS.md) does the same for the reduced ICON-Land.

## License

ICON is available under a BSD 3-clause license; see [LICENSES/](./LICENSES) and [AUTHORS.TXT](./AUTHORS.TXT). ICON-Land carries its own copyright notice (MPI-M, MPI-BGC) and license files in `externals/jsbach`.
