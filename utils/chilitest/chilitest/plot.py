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

import matplotlib.pyplot as plt
import numpy as np
from chilitest.file_info import FileInfo
from chilitest.statistics import EnsembleStatistics, Statistics


def plot_spatial_average_differences(
    test_statistics: Statistics,
    ref_statistics: EnsembleStatistics,
    ref_file_info: FileInfo,
    dim: str,
    var_list: list[str] = None,
    figsize: tuple = (8, 4),
):
    test_stats = test_statistics.spatial_averages_by_exp[dim][
        ref_file_info.single_exp_name
    ]
    ref_stats = ref_statistics.spatial_averages_by_exp[dim]

    vars_to_plot = var_list if var_list is not None else list(test_stats.keys())

    for var in vars_to_plot:
        test_arr = test_stats[var]

        plt.figure(figsize=figsize)

        # Reference values
        ref_arr = ref_stats[ref_file_info.ensemble_ref_exp_name][var]

        differences = []

        # Plot ensemble members
        first_line = True
        for ref_name, ref in ref_stats.items():
            diff = ref[var] - ref_arr
            differences.append(diff)

            if ref_name != ref_file_info.ensemble_ref_exp_name:
                plt.plot(
                    diff,
                    color="0.7",
                    linewidth=0.8,
                    label="ensemble member" if first_line else None,
                )
                first_line = False

        # Plot test experiment difference
        test_diff = test_arr - ref_arr
        differences.append(test_diff)

        plt.plot(
            test_diff,
            color="black",
            linewidth=2.5,
            label="test experiment",
        )

        # Scale y-axis symmetrically around zero
        max_abs = max(np.max(np.abs(d)) for d in differences)
        plt.ylim(-max_abs, max_abs)

        plt.xlabel("Timestep index")
        plt.ylabel("$x - x_{\\mathrm{ref}}$")
        plt.title(f"{var} ({dim})")
        plt.legend()
        plt.tight_layout()
        plt.show()
