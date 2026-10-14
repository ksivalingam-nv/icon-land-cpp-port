# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

import numpy as np
from chilitest.statistics import EnsembleStatistics, Statistics


class Test:
    def __init__(
        self,
        test_statistics: Statistics,
        ref_statistics: EnsembleStatistics,
        dimensions: list[str] | None = None,
    ):
        # Validate that test and reference statistics are compatible
        if test_statistics.variables != ref_statistics.variables:
            raise ValueError(
                f"Variable mismatch between test and reference statistics.\n"
                f"  Test variables:  {test_statistics.variables}\n"
                f"  Ref variables:   {ref_statistics.variables}"
            )
        if (
            test_statistics.spatial_average_type
            != ref_statistics.spatial_average_type
        ):
            raise ValueError(
                f"Spatial average type mismatch between test and reference statistics.\n"
                f"  Test spatial_average_type:  {test_statistics.spatial_average_type}\n"
                f"  Ref spatial_average_type:   {ref_statistics.spatial_average_type}"
            )
        if test_statistics.dimensions != ref_statistics.dimensions:
            raise ValueError(
                f"Dimensions mismatch between test and reference statistics.\n"
                f"  Test dimensions: {test_statistics.dimensions}\n"
                f"  Ref dimensions:  {ref_statistics.dimensions}"
            )

        self.exp_name = test_statistics.file_info.single_exp_name
        self.test_avgs = test_statistics.spatial_averages_by_exp
        self.ref_stats = ref_statistics.ensemble_statistics

        self.dimensions = (
            dimensions if dimensions is not None else test_statistics.dimensions
        )

    def _eps_scaled_ref(self, dim: str, var: str):
        """Return machine-epsilon absolute reference value scaled to the reference value.
        This is used as a fallback tolerance in case of zero positive and/or negative deviation.
        """
        return np.finfo(np.float64).eps * np.abs(
            self.ref_stats[dim][f"{var}_ref"]
        )

    def extreme_values(self, factor) -> None:
        factor = np.float64(factor)  # scaling factor for the tolerance range

        for dim in self.dimensions:
            for var, spatial_averages in self.test_avgs[dim][
                self.exp_name
            ].items():

                diff = spatial_averages - self.ref_stats[dim][f"{var}_ref"]

                max_dev = (
                    self.ref_stats[dim][f"{var}_max"]
                    - self.ref_stats[dim][f"{var}_ref"]
                )
                min_dev = (
                    self.ref_stats[dim][f"{var}_min"]
                    - self.ref_stats[dim][f"{var}_ref"]
                )
                tolerance_max = factor * np.maximum(
                    max_dev, self._eps_scaled_ref(dim, var)
                )
                tolerance_min = factor * np.minimum(
                    min_dev, -self._eps_scaled_ref(dim, var)
                )

                result = np.where(
                    diff >= 0,
                    diff <= tolerance_max,
                    np.abs(diff) <= np.abs(tolerance_min),
                )

                if np.all(result):
                    print(f"✓ Extreme-values test passed for {dim}/{var}")
                    continue

                failed = ~result
                idx = np.argwhere(failed)

                raise AssertionError(
                    f"Extreme-values test failed for {dim}/{var}\n"
                    f"Failed elements: {failed.sum()} / {failed.size}\n"
                    f"Indices:\n{idx}\n\n"
                    f"Observed values:\n{spatial_averages[failed]}\n"
                    f"Reference values:\n{self.ref_stats[dim][f'{var}_ref'][failed]}\n"
                    f"Differences:\n{diff[failed]}\n"
                    f"Allowed positive tolerance:\n{tolerance_max[failed]}\n"
                    f"Allowed negative tolerance:\n{tolerance_min[failed]}"
                )

    def z_score(self, tolerance) -> None:
        tolerance = np.float64(tolerance)

        for dim in self.dimensions:
            for var, spatial_averages in self.test_avgs[dim][
                self.exp_name
            ].items():

                x = spatial_averages
                mean = self.ref_stats[dim][f"{var}_mean"]
                std = np.maximum(
                    self.ref_stats[dim][f"{var}_std"],
                    self._eps_scaled_ref(dim, var),
                )

                z = np.abs(x - mean) / std
                passed = z < tolerance

                if np.all(passed):
                    print(f"✓ Z-score test passed for {dim}/{var}")
                    continue

                failed = ~passed
                idx = np.argwhere(failed)

                raise AssertionError(
                    f"Z-score test failed for {dim}/{var}\n"
                    f"Failed elements: {failed.sum()} / {failed.size}\n"
                    f"Indices:\n{idx}\n\n"
                    f"Observed values:\n{x[failed]}\n"
                    f"Mean:\n{mean[failed]}\n"
                    f"Std:\n{std[failed]}\n"
                    f"Z-scores:\n{z[failed]}\n"
                    f"Tolerance: {tolerance}"
                )
