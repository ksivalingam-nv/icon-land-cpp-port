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
#include "types.hpp"

#include <gtest/gtest.h>

#include "support/memory/scratch_view.hpp"

namespace {
using testing::Types;

template <typename T>
class TypesTest : public testing::Test {
 public:
  // workaround for nvcc + kokkos + gtest bug: cannot dispatch a KOKKOS lambda from a
  // non-public parent function when compiled with nvcc
  static void update_view_2d(ManagedView2D<T> view2d, int dim1, int dim2) {
    using RangePolicy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;

    Kokkos::parallel_for(
        "init2d", RangePolicy2D({0, 0}, {dim1, dim2}),
        KOKKOS_LAMBDA(const int i, const int j) { view2d(i, j) = static_cast<T>(i * dim2 + j); });
    Kokkos::fence();
  }

  // workaround for nvcc + kokkos + gtest bug: cannot dispatch a KOKKOS lambda from a
  // non-public parent function when compiled with nvcc
  static void update_view_3d(ManagedView3D<T> view3d, int dim1, int dim2, int dim3) {
    using RangePolicy3D = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

    Kokkos::parallel_for(
        "init3d", RangePolicy3D({0, 0, 0}, {dim1, dim2, dim3}), KOKKOS_LAMBDA(const int i, const int j, const int k) {
          view3d(i, j, k) = static_cast<T>(i * dim2 * dim3 + j * dim3 + k);
        });
    Kokkos::fence();
  }

  // same as above but for an unmanaged view
  static void update_view_2d(View2D<T, MemorySpace> view2d, int dim1, int dim2) {
    using RangePolicy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;

    Kokkos::parallel_for(
        "init_scratch_2d", RangePolicy2D({0, 0}, {dim1, dim2}),
        KOKKOS_LAMBDA(const int i, const int j) { view2d(i, j) = static_cast<T>(i * dim2 + j); });
    Kokkos::fence();
  }

  // Adding one to existing value for 2d strided view
  static void increment_view_2d(StridedRankedView<T, 2> strided_view2d, int dim1, int dim2) {
    using RangePolicy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;

    Kokkos::parallel_for(
        "increment_strided_view_2d", RangePolicy2D({0, 0}, {dim1, dim2}),
        KOKKOS_LAMBDA(const int i, const int j) { strided_view2d(i, j) += static_cast<T>(1); });
    Kokkos::fence();
  }

};  // end class

TYPED_TEST_SUITE_P(TypesTest);

TYPED_TEST_P(TypesTest, CheckManagedView2DLayout) {
  const int dim1 = 2;
  const int dim2 = 3;

  ManagedView2D<TypeParam> view2d("2d", dim1, dim2);
  Kokkos::View<TypeParam**, Kokkos::LayoutRight, Kokkos::HostSpace> hostView("host", dim1, dim2);
  TypesTest<TypeParam>::update_view_2d(view2d, dim1, dim2);
  Kokkos::deep_copy(hostView, view2d);

  // check the layout setups
  EXPECT_EQ(hostView(0, 2), TypeParam{2});
}

TYPED_TEST_P(TypesTest, CheckManagedView3DLayout) {
  const int dim1 = 2;
  const int dim2 = 3;
  const int dim3 = 2;

  ManagedView3D<TypeParam> view3d("3d", dim1, dim2, dim3);
  Kokkos::View<TypeParam***, Kokkos::LayoutRight, Kokkos::HostSpace> hostView("host", dim1, dim2, dim3);
  TypesTest<TypeParam>::update_view_3d(view3d, dim1, dim2, dim3);
  Kokkos::deep_copy(hostView, view3d);

  EXPECT_EQ(hostView(0, 1, 1), TypeParam{3});
}

TYPED_TEST_P(TypesTest, CheckUnmanagedView2DLayout) {
  const int dim1 = 2;
  const int dim2 = 3;

  // emulate ICON's action to allocate the pointer in device memory
  ManagedView2D<TypeParam> icon2d("2d", dim1, dim2);
  TypesTest<TypeParam>::update_view_2d(icon2d, dim1, dim2);

  // create the view based on the existing pointer in the execution space
  View2D<TypeParam, MemorySpace> view2d(icon2d.data(), dim1, dim2);

  // expect reuse of the pointer in the same memory space
  EXPECT_EQ(icon2d.data(), view2d.data());
}

TYPED_TEST_P(TypesTest, CheckUnmanagedView3DLayout) {
  const int dim1 = 2;
  const int dim2 = 3;
  const int dim3 = 2;

  // emulate ICON's action to allocate the pointer in device memory
  ManagedView3D<TypeParam> icon3d("3d", dim1, dim2, dim3);
  TypesTest<TypeParam>::update_view_3d(icon3d, dim1, dim2, dim3);

  // create the view based on the existing pointer in the execution space
  View3D<TypeParam, MemorySpace> view3d(icon3d.data(), dim1, dim2, dim3);

  // expect reuse of the pointer in the same memory space
  EXPECT_EQ(icon3d.data(), view3d.data());
}

TYPED_TEST_P(TypesTest, CheckHostView) {
  const int dim1          = 2;
  TypeParam array1D[dim1] = {TypeParam{0}, TypeParam{1}};
  auto view1d             = HostView1D<TypeParam>(array1D, dim1);

  // expect always the same pointer, regardless of the Kokkos backend
  EXPECT_TRUE(array1D == view1d.data());
}

TYPED_TEST_P(TypesTest, CheckScratchView2DLayout) {
  const int dim1 = 2;
  const int dim2 = 3;

  // Get local space:
  ragnarok::LocalScratchSpace<Kokkos::DefaultExecutionSpace> local_space =
      ragnarok::get_local_space<Kokkos::DefaultExecutionSpace>();

  // Create scratch view similar to how it's done in graupel.ipp
  View2D<TypeParam, MemorySpace> view2d = ragnarok::scratch_view(local_space, dim1, dim2);

  // Initialize the scratch view with some values
  TypesTest<TypeParam>::update_view_2d(view2d, dim1, dim2);

  // Verify the values were set correctly
  Kokkos::View<TypeParam**, Kokkos::LayoutRight, Kokkos::HostSpace> hostView("host", dim1, dim2);
  Kokkos::deep_copy(hostView, view2d);

  // Check some specific values
  EXPECT_EQ(hostView(0, 0), TypeParam{0});
  EXPECT_EQ(hostView(0, 2), TypeParam{2});
  EXPECT_EQ(hostView(1, 1), TypeParam{4});
}

TYPED_TEST_P(TypesTest, RankedViewBehaviour2D) {
  const int dim1 = 2;
  const int dim2 = 3;

  ManagedView2D<TypeParam> icon2d("2d", dim1, dim2);
  TypesTest<TypeParam>::update_view_2d(icon2d, dim1, dim2);

  RankedView<TypeParam, 2, MemorySpace> view2d(icon2d.data(), dim1, dim2);

  EXPECT_EQ(icon2d.data(), view2d.data());
}

TYPED_TEST_P(TypesTest, RankedViewBehaviour3D) {
  const int dim1 = 2;
  const int dim2 = 3;
  const int dim3 = 5;

  ManagedView3D<TypeParam> icon3d("3d", dim1, dim2, dim3);
  TypesTest<TypeParam>::update_view_3d(icon3d, dim1, dim2, dim3);

  RankedView<TypeParam, 3, MemorySpace> view3d(icon3d.data(), dim1, dim2, dim3);

  EXPECT_EQ(icon3d.data(), view3d.data());
}

TYPED_TEST_P(TypesTest, StridedRankedViewTest) {
  const int dim1 = 4;
  const int dim2 = 3;
  const int dim3 = 2;

  // All in execution Kokkos::DefaultExecutionSpace
  ManagedView3D<TypeParam> icon3d("3d", dim1, dim2, dim3);
  TypesTest<TypeParam>::update_view_3d(icon3d, dim1, dim2, dim3);

  StridedRankedView<TypeParam, 2> strided_slice = Kokkos::subview(icon3d, Kokkos::ALL, 1, Kokkos::ALL);

  EXPECT_EQ(strided_slice.extent(0), dim1);
  EXPECT_EQ(strided_slice.extent(1), dim3);

  TypesTest<TypeParam>::increment_view_2d(strided_slice, dim1, dim3);

  // Getting from Device to Host (icon3d)
  auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, icon3d);

  // Verify that only j = 1 was modified.
  for (int i = 0; i < dim1; ++i) {
    for (int j = 0; j < dim2; ++j) {
      for (int k = 0; k < dim3; ++k) {
        TypeParam expected = static_cast<TypeParam>(i * dim2 * dim3 + j * dim3 + k);
        if (j == 1) {
          expected += static_cast<TypeParam>(1);
        }
        EXPECT_EQ(host(i, j, k), expected) << "Mismatch at (" << i << "," << j << "," << k << ")";
      }
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(TypesTest, CheckManagedView2DLayout, CheckManagedView3DLayout, CheckUnmanagedView2DLayout,
                            CheckUnmanagedView3DLayout, CheckHostView, CheckScratchView2DLayout, RankedViewBehaviour2D,
                            RankedViewBehaviour3D, StridedRankedViewTest);

using MyTypes = ::testing::Types<float, double>;
INSTANTIATE_TYPED_TEST_SUITE_P(RagnarokTests, TypesTest, MyTypes);

}  // namespace
