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
import logging
from datetime import datetime

import hiopy.configure as hc
import isodate
import numpy as np
import zarr


def init_data_request(
    dataset: zarr.Group,
    data_request: dict,
    simulation_params: dict,
    dict_of_names: dict,
):
    start = datetime.fromisoformat(simulation_params["start_date"])
    end = datetime.fromisoformat(simulation_params["end_date"])
    logging.debug(f"Simulation start date: {start}, end date: {end}")
    logging.debug(
        f"Total simulation duration in hours: {(end - start).total_seconds() / 3600}"
    )
    total_simulated_hours = (end - start).total_seconds() / 3600

    zarr_args = {}
    if dataset.metadata.zarr_format == 3:
        zarr_args["compressors"] = [
            zarr.codecs.BloscCodec(
                cname="zstd", clevel=6, shuffle=zarr.codecs.BloscShuffle.shuffle
            )
        ]
    else:
        import numcodecs

        zarr_args["compressors"] = numcodecs.Blosc(
            cname="zstd", clevel=6, shuffle=numcodecs.Blosc.SHUFFLE
        )
        zarr_args["chunk_key_encoding"] = {"name": "v2", "separator": "/"}

    for data_group in data_request:
        try:
            group_config = data_request[data_group]["settings"]
            group_vars = data_request[data_group]["variables"]
            group_grid_details = data_request[data_group]["grid_details"]

            for ts_needed in group_config["time_res_needed"]:
                parent_zoom_lev = np.max(group_grid_details["zoom_levels"])
                logging.debug(
                    f"Parent zoom level for {data_group} is {parent_zoom_lev}"
                )

                for zoom_lev in group_grid_details["zoom_levels"]:
                    parent_group = ""
                    hiopy_args = {
                        "attributes": {},
                        "yac_source_comp": group_config["yac_source_comp"],
                        "yac_source_grid": group_config["yac_source_grid"],
                    }
                    time_method = group_config["time_method"]
                    group_name = f"{ts_needed}_{time_method}_z{zoom_lev}"
                    zg = None
                    if group_name in dataset:
                        zg = zarr.open_group(
                            path=f"{group_name}", store=dataset.store
                        )
                        logging.debug(
                            f"opened existing group {group_name} for adding more variables."
                        )
                    else:
                        zg = dataset.create_group(
                            name=group_name,
                            attributes=dict(simulation_params["attribution"]),
                        )

                        hc.add_healpix_grid(zg, zoom_lev)
                        hc.add_time(
                            zg,
                            np.datetime64(
                                datetime.fromisoformat(
                                    simulation_params["start_date"]
                                ).replace(tzinfo=None)
                            ),
                            np.datetime64(
                                datetime.fromisoformat(
                                    simulation_params["end_date"]
                                ).replace(tzinfo=None)
                            ),
                            ts_needed,
                        )

                        if zoom_lev < parent_zoom_lev:
                            order = zoom_lev + 1
                            parent_group = f"{ts_needed}_{time_method}_z{order}"
                            zg.attrs["hiopy::parent"] = parent_group

                        logging.debug(
                            f"created group {group_name} ; parent group: {parent_group}"
                        )

                    ncells = 12 * 4**zoom_lev

                    cell_chunk = ncells

                    hours_per_chunk = simulation_params[
                        "simulated_hours_per_chunk"
                    ]
                    chunks_per_shard = simulation_params["chunks_per_shard"]

                    logging.debug(
                        f"{data_group}::{group_name} : Total simulated hours: {total_simulated_hours}, hours per chunk: {hours_per_chunk}, chunks per shard: {chunks_per_shard}"
                    )

                    if hours_per_chunk > total_simulated_hours:
                        logging.warning(
                            f"{data_group}::{group_name} : Upper limit of possible hours per chunk is {total_simulated_hours} instead of the configured {hours_per_chunk}. Selecting that for time_chunk."
                        )
                        hours_per_chunk = int(total_simulated_hours)

                    duration = isodate.parse_duration(ts_needed)
                    enforce_single_file_per_timestep = False
                    if not isinstance(duration, isodate.duration.Duration):
                        time_chunk = int(
                            simulation_params["simulated_hours_per_chunk"]
                            * 60
                            * 60
                            / duration.total_seconds()
                        )
                    else:
                        enforce_single_file_per_timestep = True
                        logging.warning(
                            f"{data_group}::{group_name} : Setting smaller chunk for duration {ts_needed} to avoid potential data loss from longer runs."
                        )

                    if zoom_lev == 8:
                        cell_chunk = int(ncells / 4)

                    if zoom_lev == 9:
                        cell_chunk = int(ncells / 16)

                    if zoom_lev >= 10:
                        cell_chunk = int(ncells / 64)

                    if zoom_lev >= 10 and time_chunk > 24:
                        enforce_single_file_per_timestep = True
                        logging.warning(
                            f"{data_group}::{group_name} : Selected chunking will lead to very large file sizes for zoomlev >= 10. Setting smaller chunk to avoid memory issues in runtime."
                        )

                    if enforce_single_file_per_timestep:
                        time_chunk = 1
                        chunks_per_shard = 1

                    chunk_shape = (time_chunk, cell_chunk)

                    is_3d = group_config.get("name_of_level", None)
                    if is_3d:
                        height_chunk = len(group_config["levels"])
                        if zoom_lev > 9:
                            logging.warning(
                                f"{data_group}::{group_name} : Selected chunking will lead to very large file sizes for 3D vars. Setting chunks_per_shard to 1 to avoid memory issues in runtime."
                            )
                            chunks_per_shard = 1
                        chunk_shape = (
                            chunk_shape[0],
                            height_chunk,
                            chunk_shape[1],
                        )

                    if is_3d and group_config["name_of_level"] not in zg:
                        collection_selection = group_config.get(
                            "indices_to_save", None
                        )
                        if collection_selection is not None:
                            collection_selection = [
                                int(idx) for idx in collection_selection
                            ]
                        hc.add_height(
                            zg,
                            group_config["name_of_level"],
                            group_config["positive_direction"],
                            group_config["levels"],
                            collection_selection,
                        )
                    if "oce" in hiopy_args["yac_source_comp"]:
                        if "ocean_frac_mask_sfc" not in zg:
                            hiopy_args["attributes"][
                                "hiopy::interpolation_stack"
                            ] = "[conservative: {partial_coverage: true}, fixed: {user_value: 0.0}]"

                            if zg.metadata.zarr_format == 3:
                                zarr_args["chunks_per_shard"] = 1
                            hc.add_variable(
                                zg,
                                name="ocean_frac_mask_sfc",
                                taxis=None,
                                chunk_shape=chunk_shape[-1:],
                                yac_name="valid_mask_sfc",
                                **{**hiopy_args, "zaxis": None},
                                **zarr_args,
                            )

                        if is_3d and group_config["frac_mask"] not in zg:
                            hiopy_args["attributes"][
                                "hiopy::interpolation_stack"
                            ] = "[conservative: {partial_coverage: true}, fixed: {user_value: 0.0}]"
                            frac_mask_yac_name = "valid_mask"
                            if "half" in group_config["frac_mask"]:
                                frac_mask_yac_name = "valid_mask_half"

                            if zg.metadata.zarr_format == 3:
                                zarr_args["chunks_per_shard"] = 1

                            hc.add_variable(
                                zg,
                                name=group_config["frac_mask"],
                                taxis=None,
                                chunk_shape=(
                                    len(group_config["levels"]),
                                    chunk_shape[-1],
                                ),
                                yac_name=frac_mask_yac_name,
                                zaxis=group_config.get("name_of_level", None),
                                **hiopy_args,
                                **zarr_args,
                            )

                    for variable in group_vars:
                        icon_var_name = variable
                        if (
                            "is_pressure_level" in group_config
                            and group_config["is_pressure_level"]
                        ):
                            icon_var_name = "pl::" + variable

                        hiopy_args["attributes"][
                            "hiopy::interpolation_stack"
                        ] = "[average, fixed: {user_value: NAN}]"
                        hiopy_args["attributes"]["hiopy::copy_metadata"] = True

                        if parent_group == "":
                            hiopy_args["attributes"][
                                "hiopy::interpolation_stack"
                            ] = "[conservative: {partial_coverage: true}, fixed: {user_value: NAN}]"

                        if zg.metadata.zarr_format == 3:
                            zarr_args["chunks_per_shard"] = chunks_per_shard

                        if (
                            variable in dict_of_names
                        ):  # use the user's overrides for names if set
                            icon_var_name = dict_of_names[variable]

                        if variable not in zg:
                            hc.add_variable(
                                zg,
                                name=variable,
                                time_method=time_method,
                                chunk_shape=chunk_shape,
                                frac_mask=group_config.get("frac_mask", None),
                                zaxis=group_config.get("name_of_level", None),
                                yac_name=icon_var_name,
                                **hiopy_args,
                                **zarr_args,
                            )
        except Exception as e:
            logging.error(
                f"Store initialisation failed for {data_request[data_group]} with error: {e}"
            )
            raise e

    zarr.consolidate_metadata(dataset.store)
