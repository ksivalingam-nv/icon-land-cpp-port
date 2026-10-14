# ICON Ragnarok Contribution Guide

## Goals

C++/Kokkos development in ICON should aim to meet the following goals, listed in order of importance:

1. **Readability** - Favor simple, clear implementations and avoid unnecessary complexity or obscurity. Use explicit include paths whenever this improves clarity. A useful set of practices can be found in the [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines).
2. **Correctness** - Unit tests are mandatory, and including functional tests in a merge request is strongly encouraged. For unit-testing, ragnarok uses the [GoogleTest](https://google.github.io/googletest/) framework, which is either provided by the system (e.g., on Levante) or fetched automatically during the build.
3. **Performance portability** - New code must perform well across different target architectures.

ICON’s GitLab CI infrastructure supports both correctness and performance testing. If you are unsure how to validate your changes, feel free to reach out.

> **Note:** Ragnarok aims for a rigorous review process via GitLab merge requests. Ensure that your code follows the principles above, includes appropriate documentation, and is covered by automated tests before marking a merge request as 'Ready'.

## Building Ragnarok

Ragnarok uses the [CMake](https://cmake.org) build system and is automatically configured by ICON through [`configure.ac`](../configure.ac). To build ICON with Ragnarok, pass the `--enable-ragnarok` option to the appropriate configure wrapper. Ragnarok can also be built as a standalone project to allow isolated testing of physics components. Details on dependencies, project structure, compilation, and testing are provided in the project [README](README.md).

When ICON is configured with `--enable-ragnarok`, the [Kokkos](https://github.com/kokkos/kokkos) submodule is built automatically using the selected backend. For information about using or developing with Kokkos, refer to the official [Kokkos tutorials](https://github.com/kokkos/kokkos-tutorials).

Any configure wrapper using a CXX compiler with **C++20** support can be used to build ragnarok. Tested compiler stacks on Levante: Intel oneAPI 2025.3.2, GNU gcc 11.2, nvhpc 24.7.

There is an option to build only ragnarok, for the purpose of checking only the C++ code and running the unit-tests. For that, there are pre-defined build recipes for Levante and LUMI in [ragnarok/cmake/caches/](cmake/caches/) folder. Sample build script:

```sh
# inside icon-mpim root folder
cmake -C ragnarok/cmake/caches/levante_cpu_intel.cmake \
      -B build_only_ragnarok \
      -S ragnarok
cmake --build build_only_ragnarok --parallel
cd build_only_ragnarok && ctest
```

> **Note:** Ragnarok build recipes (also used by the gitlab CI pipelines) are strict in terms of warnings and treat them as compilation errors.

## Adding a New Physics Component to Ragnarok

To add a new physics component, create a directory named (`<component_name>`) in `ragnarok/` and place the component’s source and test files there.
Add the component to the build system by calling `add_subdirectory(<component_name>)` in the project's [CMakeLists.txt](CMakeLists.txt).
For reference, the `aes_thermodynamics` structure looks like this:
```
└── 📁ragnarok
    └── 📁aes_thermodynamics
        ├── CMakeLists.txt
        ├── thermo_constants.hpp
        ├── thermo_test.cpp
        ├── thermo.hpp
        ├── thermo.ipp
```

### Acceptance criteria for merging new code into ICON

- The code should be as self-explanatory as much as possible and approved via gitlab review process;
- The code should have documentation: public functions and interfaces must include inline documentation using Doxygen syntax `///`. Documentation for private functions is encouraged but optional. Place documentation close to the implementation where appropriate— for example, in `.cpp` rather than the `.hpp` file or in `.md` files in the source directory.
- The code must be tested:
    - automatic testing: add new unit- or ICON-ragnarok integration tests (via testbed or running full experiments in gitlab-ci) for the newly added code
    - signoff from scientific code owners
- The code it must run with acceptable performance on Levante CPU and Levante GPU. To replace existing code it must run with comparable performance to the Fortran code on Levante CPU, Levante GPU, and LUMI GPU in case of the AES physics, and on all systems included in the ICON test infrastructure in the case of shared model infrastructure.

## Extending the ICON testbed with integration tests

Integration tests that exercise ragnarok inside a fully constructed ICON model
can be added in [`src/testbed/mo_test_ragnarok.f90`](../src/testbed/mo_test_ragnarok.f90);
any C++ code they need goes into
[`ragnarok_testbed.cpp`](support/icon_bridge/ragnarok_testbed.cpp) with its
Fortran interface in
[`mo_ragnarok_testbed.f90`](support/icon_bridge/mo_ragnarok_testbed.f90).
Write a subroutine following the pattern below, call it from
`test_icon_bridge` or `run_ragnarok_tests`, and GitLab CI runs it via the `ragnarok:testbed` job family:

```fortran
  SUBROUTINE test_my_feature()
    CHARACTER(len=*), PARAMETER :: context = module_name//'::test_my_feature'
    REAL(wp) :: result, reference

    ! compare the ragnarok result against an ICON reference implementation
    IF (result /= reference) CALL finish(context, 'Test failed')
    CALL message(context, 'Test passed.')
  END SUBROUTINE test_my_feature
```

## Design concepts

### Multi-dimensional arrays in ragnarok

Ragnarok is using `Kokkos::View` to store multi-dimensional arrays. The current types which describe the agreed ICON conventions are defined in [ragnarok/common/types.hpp](common/types.hpp). **Please check these types before defining new ones.**

> **Note:** The default layout for the ragnarok code is row-major (aka `Kokkos::LayoutRight`). This is slightly counterintuitive since ICON Fortran uses column-major. This is just the interpretation of the data pointer (aka pointer arithmetics of the kokkos view indices), while the heap memory is serialized column-major wise (and will remain so while ICON has both C++ and Fortran code in parallel).  

To slice the data, use `Kokkos::subview`, which will create a strided view (see `StridedView3D` definitions in [ragnarok/common/types.hpp](common/types.hpp)).
> **Note:** Be aware that subviews can reference a non-contiguous array in memory which will lead to performance penalties. See example below:

```cpp
// assume the following 3D array
View3D<T> rho(rho_ptr, nblocks, nlev, nproma);

// retrieve the 2D array for a given block, jb = 5
// the result is still a continuous array in memory -> GOOD!
View2D<T> rho_block_5 = Kokkos::subview(rho, 5, Kokkos::ALL, Kokkos::ALL);

// assume the following 4D array: aes_phy_field%qtrc_phy
View4D<T> qtrc_phy(q_ptr, 6, nblocks, nlev, nproma);

//retrieve all the humidity fields for a given block, jb = 5
// the result has a gap in the slow dimension, OK for nlev, nproma loops, but slow for others
View3D<T> qx_block_5 = Kokkos::subview(qtrc_phy, Kokkos::ALL, 5, Kokkos::ALL, Kokkos::ALL);
```

#### Dealing with persistent data

1. Constant data and configurations

Since the `t_patch` is a user defined type, it cannot be sent via ISO-C bindings so then ragnarok reconstructs on the C++ side with only the fields that are needed.
Read more in the [docs](support/icon_bridge/docs/config_recon.md).

```cpp
#include "support/icon_bridge/icon_domain.hpp"

const icon::Patch& patch         = icon::get_patch(jg); // domain id
const int nproma                 = patch.cells.block_size;
const int nlev                   = patch.nlev;
const int nblocks                = patch.nblks_c;
```

Similar approach with namelists, ragnarok retrieves only the required fields.

```cpp
#include "support/icon_bridge/icon_config.hpp"

// pass the domain id
const TurbVdiffConfig config = icon::get_aes_vdf_config(jg);

printf("%d\n", config.solver_type);
```

> **Note:** If some fields are missing in ragnarok, please open a work item in gitlab and ask for it.  

2. Ragnarok gets its data pointers explicitly from ICON Fortran via Fortran ISO-C bindings interface

Contiguous arrays are passed as raw C pointers which are then wrapped into unmanaged `Kokkos::View`s. This avoids unnecessary memory copies (regardless if these are host or device pointers). Metadata describing the shape of the Fortran array may also have to be passed unless already known from the domain context, e.g., from the patch available in ragnarok.

```cpp

#include "common/types.hpp"

template <typename T>
void f(int nblocks, int nlev, int nproma, const T* rho_ptr) {
   View3D<T> rho(rho_ptr, nblocks, nlev, nproma);
}

```


> **Note:** Assumption: ICON always sends the pointer in the default execution space, so host pointer when ICON is compiled for CPU and device otherwise.
To ensure the interoperability with OpenACC, make sure to call `ACC HOST_DATA USE_DEVICE` as shown in the code snippet below. Moreover, ragnarok ensures the consistency of the GPU pointers between Fortran and C++ by using the same CUDA/HIP stream when constructing the default `Kokkos::ExecutionSpace` which is used in the Kokkos kernels, hence no CPU-GPU explicit synchronisation (e.g `ACC WAIT(1)`) is required.

```fortran
!$ACC HOST_DATA USE_DEVICE(rho)

! call Fortran ISO-C bindings routine implemented in ragnarok
CALL ragnarok_bindings(jg)

!$ACC END HOST_DATA
```

3. Ragnarok retrieves the data pointers from the MemoryManager - **NOT FULLY supported yet!**

Only the fields allocated with the memory manager can be retrieved via the memory manager.

#### Temporary arrays

To allocate multi-dimensional **local** variables, use the local memory pool. This is a pre-allocated region in heap memory in the default execution space (hence portable) and it acts like a stack: once the variable is out of scope, the "stack" pointer gets decremented without deallocating the buffer which can be reused by the next component.
A simple demonstrator is shown below but for full documentation, read the [docs](support/memory/docs/scratch_memory.md).

```cpp
#include "support/memory/scratch_view.hpp"
#include "common/types.hpp"

auto local_space      = ragnarok::get_local_space();
View2D<T> local_var = ragnarok::scratch_view(local_space, nlev, nproma);
```

### Portable loops

Ragnarok provides a few types of loops which infer the parallelisation strategy automatically and adapt to it, by transforming itself from a `for` loop to `Kokkos::parallel_for` loop.
For usage, see the following code snippet or the associated [docs](support/kernels/docs/icon_loops.md) for more complex use cases.

 ```cpp
#include "support/kernels/grid_range.hpp"
#include "support/kernels/icon_loops.hpp"
#include "support/icon_bridge/icon_domain.hpp"

const icon::Patch& patch           = icon::get_patch(jg);
const icon::GridRange grid_range   = patch.cells.prog_range;
const icon::LevelRange level_range{0, patch.nlev};

icon::block_loop(
   "block_loop", grid_range, block_lambda(const int jb) {

    const icon::IndexRange index_range = grid_range.index_range(jb);

    // extra code here ...

	icon::index_loop(
           "update_field", level_range, index_range, index_lambda(int jk, int jc) {
                 rho(jb, jk, jc) = T{0.};
         });
});
```

> **Note:** There is no longer needed to call `get_indices_*` because the `GridRange` already captures the information.

To ensure that the loops are **portable**, keep in mind a few GPU restrictions for the body of the loops:
- do not use `std::` functions  (including `std::cout`) or data structures (e.g. `std::vector`); C functions (e.g. `printf`) are safe;
- do not call recursive functions; the GPU stack is much smaller than the CPU one and will overflow rapidly;
- the only accessible data structures are the `Kokkos::View`, don't try to capture raw pointers or data structures. That will lead to either a Segmentation Fault or would imply an unnecessary memory copy which is expensive;
- do not create `Kokkos::View`s inside the loops. There is no way to dynamically allocate memory on the device.
- use only portable math operators (e.g. `Kokkos::fmax`, `Kokkos::pow` etc)
- if you ever need to use an `if cpu` statement, assume your code is wrong. There are some situations when the same code should execute on CPU and GPU during the same simulation run (e.g. calling ICON math for model initialisation on CPU and from inside GPU kernels), but these are limited.

### Timers in C++

Ragnarok timers appear in the ICON timer report. Usage is described in
[timer_usage.md](support/timer/docs/timer_usage.md), the model behind it in
[timer_design.md](support/timer/docs/timer_design.md). The main rules:
- define timers at namespace scope; a timer must never be constructed inside a parallel region;
- a timer states the `timers_level` at which it becomes active at construction. Timers around kernels use `timer::timers_level_sync_threshold`: from that level on every measurement synchronises the device first, which a proper kernel measurement on GPUs needs because the CPU has to wait for the completion of device work;
- place `start`/`stop` outside or inside a `block_loop` body, never inside an `index_loop` body, which may become a GPU kernel.

```cpp
#include "support/timer/timer.hpp"

namespace my_ns {

timer::Timer timer_bar("timing_function_bar");                                   // active from timers_level 1
timer::Timer timer_kernel("timing_kernel", timer::timers_level_sync_threshold);  // synchronised kernel timing

void foo() {
  // ...

  timer_bar.start();
  bar();
  timer_bar.stop();

  // ...
}

} // end namespace
```

## Coding Style

The ICON C++/Kokkos codebase follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html), with minor ICON-specific adaptations. We use [`pre-commit`](https://pre-commit.com) hooks to ensure that:

* code formatting follows the Google style via [`clang-format` tool](https://clang.llvm.org/docs/ClangFormat.html) and
* common technical issues are identified early through static analysis using [`clang-tidy` tool](https://clang.llvm.org/extra/clang-tidy/).

These tools automatically modify staged files. The same checks run in GitLab CI and are required for any merge request. Instructions for enabling `pre-commit` in ICON are available in the [ICON Contribution Guidelines](../CONTRIBUTING.md#code-contributions).
