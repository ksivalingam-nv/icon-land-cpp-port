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

"""Central registry of the string values accepted in application configs.

Every config value that selects behaviour by name is declared here.
"""

# `spatial_average_type` config values: (dimension, average type) -> name of
# the `Statistics` method implementing it.
SPATIAL_AVERAGE_TYPES: dict[tuple[str, str], str] = {
    ("2d", "global_mean"): "global_mean_2d",
    ("3d", "surface_mean_per_p_level"): "surface_mean_per_p_level",
}

# `perturbation_method` config values accepted by `Ensemble`
PERTURBATION_METHODS: tuple[str, ...] = (
    "pinit_nml_params",
    "input_netcdf",
)
