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

import json
import os
from typing import Callable

import numpy as np
import xarray as xr
from chilitest.file_info import FileInfo
from chilitest.registry import SPATIAL_AVERAGE_TYPES


class Statistics:
    _io_precision_error = "Expected np.float64 precision for all variables"

    def __init__(
        self,
        file_info: FileInfo,
        variables: dict[str, list[str]],
        spatial_average_type: dict[str, str],
        dimensions: list[str] | None = None,
    ):
        self.file_info = file_info
        self.variables = variables
        self.spatial_average_type = spatial_average_type

        # The dimensions to process are already implied by the config's
        # `variables`/`spatial_average_type` keys -- derive them instead of
        # keeping a second source of truth. `dimensions` remains as an
        # override to restrict processing to a subset.
        if variables.keys() != spatial_average_type.keys():
            raise ValueError(
                "`variables` and `spatial_average_type` must declare the "
                f"same dimensions, got {sorted(variables)} vs "
                f"{sorted(spatial_average_type)}"
            )
        self.dimensions = (
            dimensions if dimensions is not None else list(variables)
        )

        self.spatial_average_functions: dict[str, Callable] = (
            self.get_average_functions()
        )

        self.dataset_by_exp: dict[str, dict[str, xr.Dataset]] = {}
        self.spatial_averages_by_exp: dict[
            str, dict[str, dict[str, np.ndarray]]
        ] = {}

        self._precision_conversion_error = RuntimeError(
            "Unexpected dtype change during xarray -> numpy conversion"
        )

    def get_average_functions(self) -> dict[str, Callable]:
        for dim in self.dimensions:
            avg = self.spatial_average_type[dim]

            if (dim, avg) not in SPATIAL_AVERAGE_TYPES:
                raise ValueError(
                    f"Unsupported combination: dim={dim}, average_type={avg}"
                    f" (supported: {sorted(SPATIAL_AVERAGE_TYPES)})"
                )

        return {
            dim: getattr(
                self,
                SPATIAL_AVERAGE_TYPES[(dim, self.spatial_average_type[dim])],
            )
            for dim in self.dimensions
        }

    def _ds_to_ndarr(self, ds: xr.Dataset, var: str) -> np.ndarray:
        arr = ds[var].to_numpy()
        if arr.dtype != ds[var].dtype:
            raise self._precision_conversion_error
        return arr

    def global_mean_2d(self, ds: xr.Dataset, var: str) -> np.ndarray:
        arr = self._ds_to_ndarr(ds, var)
        return arr.mean(axis=tuple(range(1, arr.ndim)))

    def surface_mean_per_p_level(self, ds: xr.Dataset, var: str) -> np.ndarray:
        arr = self._ds_to_ndarr(ds, var)
        return arr.mean(axis=tuple(range(2, arr.ndim)))

    def open_dataset(self):
        for dim in self.dimensions:
            ds = {
                name: xr.open_dataset(
                    self.file_info.output_file_by_exp(dim)[name]
                )
                for name in self.file_info.all_exp_names
            }

            assert all(
                ds_i[var].dtype == np.float64
                for ds_i in ds.values()
                for var in ds_i.data_vars
            ), self._io_precision_error

            self.dataset_by_exp[dim] = ds

    def compute_spatial_averages(self):
        for dim in self.dimensions:
            self.spatial_averages_by_exp[dim] = {
                name: {
                    var: self.spatial_average_functions[dim](ds, var)
                    for var in self.variables[dim]
                }
                for name, ds in self.dataset_by_exp[dim].items()
            }


class EnsembleStatistics(Statistics):
    def __init__(
        self,
        file_info: FileInfo,
        variables: dict[str, list[str]],
        spatial_average_type: dict[str, str],
        statistics: list[str],
        dimensions: list[str] | None = None,
    ):
        super().__init__(
            file_info=file_info,
            variables=variables,
            spatial_average_type=spatial_average_type,
            dimensions=dimensions,
        )

        self.statistics = statistics
        self.ensemble_statistics: dict[str, dict[str, np.ndarray]] = {}

    def open_dataset(self):
        if self.ensemble_statistics:
            raise RuntimeError(
                "Trying to open statistics from a second experiment (this EnsembleStatistics instance already contains experiment statistics)"
            )
        super().open_dataset()

    def generate_ensemble_statistics(self):
        if len(self.file_info.all_exp_names) <= 1:
            raise ValueError(
                "Please provide spatial averages calculated from ensemble experiments."
            )

        for dim in self.dimensions:
            self.ensemble_statistics[dim] = {}
            for var in self.variables[dim]:
                self.ensemble_statistics[dim][f"{var}_ref"] = (
                    self.spatial_averages_by_exp[dim][
                        self.file_info.ensemble_ref_exp_name
                    ][var]
                )

                x_name_list = [
                    self.spatial_averages_by_exp[dim][name][var]
                    for name in self.spatial_averages_by_exp[dim].keys()
                    if name != self.file_info.ensemble_ref_exp_name
                ]
                x_name = np.stack(x_name_list, axis=0)

                for stat in self.statistics:
                    np_func = getattr(np, stat)
                    self.ensemble_statistics[dim][f"{var}_{stat}"] = np_func(
                        x_name, axis=0
                    )

    def _write_meta_sidecar(self):
        meta = {
            "variables": self.variables,
            "spatial_average_type": self.spatial_average_type,
        }
        meta_path = self.file_info.ensemble_stats_meta_file()
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

    def save_ensemble_statistics(self):
        os.makedirs(self.file_info.stats_dir, exist_ok=True)
        for dim in self.dimensions:
            coordinates = {
                "times": self._ds_to_ndarr(
                    self.dataset_by_exp[dim][
                        self.file_info.ensemble_ref_exp_name
                    ],
                    "time",
                )
            }
            # 2d/global_mean averages need no coordinate beyond time
            if (
                dim == "3d"
                and self.spatial_average_type[dim] == "surface_mean_per_p_level"
            ):
                coordinates["plev"] = self._ds_to_ndarr(
                    self.dataset_by_exp[dim][
                        self.file_info.ensemble_ref_exp_name
                    ],
                    "plev",
                )

            np.savez(
                self.file_info.ensemble_stats_file(dim),
                **coordinates,
                **self.ensemble_statistics[dim],
            )

        self._write_meta_sidecar()

    @classmethod
    def load(
        cls, file_info: FileInfo, dimensions: list[str] | None = None
    ) -> "EnsembleStatistics":
        # Load sidecar metadata
        meta_path = file_info.ensemble_stats_meta_file()
        try:
            with open(meta_path, "r") as f:
                meta = json.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Statistics meta file not found at '{meta_path}'. "
                "This file is required since reference statistics are now "
                "self-describing. Regenerate reference statistics to create it."
            )

        variables = meta.get("variables", {})
        spatial_average_type = meta.get("spatial_average_type", {})

        instance = cls(
            file_info=file_info,
            variables=variables,
            spatial_average_type=spatial_average_type,
            statistics=[],
            dimensions=dimensions,
        )

        for dim in instance.dimensions:
            with np.load(file_info.ensemble_stats_file(dim)) as data:
                stats = {k: data[k] for k in data.files}

            assert all(
                arr.dtype == np.float64
                or np.issubdtype(arr.dtype, np.datetime64)
                for arr in stats.values()
            ), cls._io_precision_error

            instance.ensemble_statistics[dim] = stats

        return instance
