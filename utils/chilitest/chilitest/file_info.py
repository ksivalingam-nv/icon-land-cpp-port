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

import os
import subprocess
from pathlib import Path

from chilitest.ensemble import Ensemble


class FileInfo:
    _EXPERIMENTS_SUBDIR = "experiments"
    _SCRIPTS_SUBDIR = "scripts"
    _WORK_SUBDIR = "work"
    _LOG_SUBDIR = "log"
    _DATA_SUBDIR = "data"

    def __init__(
        self,
        mkexp_config_path: str,
        build_dir: str,
        stats_dir: str,
        ensemble: Ensemble | None = None,
        output_dirname: (
            str | None
        ) = None,  # Only necessary if user wants to give the `run_<>` folder to analyse manually
    ):
        self._ensemble_run = True if ensemble is not None else False

        self.mkexp_config_path = mkexp_config_path
        self.repo_root = self._get_repo_root_dir()
        self.mkexp_path = f"{self.repo_root}/utils/mkexp/mkexp"

        self.single_exp_name: str = Path(self.mkexp_config_path).stem
        self.ensemble_ref_exp_name: str = f"{self.single_exp_name}_ref"
        self.ensemble_exp_names: list = (
            [self.ensemble_ref_exp_name]
            + [
                f"{self.single_exp_name}_{i}"
                for i in range(1, ensemble.size + 1)
            ]
            if self._ensemble_run
            else None
        )
        self.all_exp_names: list = (
            self.ensemble_exp_names
            if self._ensemble_run
            else [self.single_exp_name]
        )

        self.stats_dir = stats_dir
        # Absolute: the path is baked into mkexp variables and runscripts
        # that execute with other working directories.
        self.build_dir = os.path.abspath(build_dir)
        self.exp_executable = f"{self.build_dir}/bin/icon"

        self.run_dir = f"{self.build_dir}/run"

        self.version_file = (
            f"{self.stats_dir}/{self.single_exp_name}_version.txt"
        )

        self.mkexp_script_dir_by_exp: dict = self._expand_to_dict(
            f"{self.build_dir}/{self._EXPERIMENTS_SUBDIR}/{{w_exp_name}}/{self._SCRIPTS_SUBDIR}"
        )
        self.mkexp_work_dir_by_exp: dict = self._expand_to_dict(
            f"{self.build_dir}/{self._EXPERIMENTS_SUBDIR}/{{w_exp_name}}/{self._WORK_SUBDIR}"
        )

        self.run_start_file_by_exp: dict = self._append(
            self.mkexp_script_dir_by_exp, "{w_exp_name}.run_start"
        )
        self.log_file_by_exp: dict = self._append(
            self.mkexp_script_dir_by_exp, "{w_exp_name}.run.log"
        )
        self.date_file_by_exp: dict = self._append(
            self.mkexp_script_dir_by_exp, "{w_exp_name}.date"
        )

        self._output_dirname = output_dirname

    def output_file_by_exp(self, dim: str) -> str:
        output_dirname = (
            self._output_dirname
            if self._output_dirname is not None
            else self._find_output_dirname()
        )
        return self._append(
            self.mkexp_work_dir_by_exp,
            f"{output_dirname}/{{w_exp_name}}_{dim}.nc",
        )

    def ensemble_stats_file(self, dim: str) -> str:
        return f"{self.stats_dir}/{self.single_exp_name}_stats_{dim}.npz"

    def ensemble_stats_meta_file(self) -> str:
        return f"{self.stats_dir}/{self.single_exp_name}_stats_meta.json"

    def mkexp_dir_args(self) -> dict[str, str]:
        exp_dir = f"{self.build_dir}/{self._EXPERIMENTS_SUBDIR}/$EXP_ID"
        return {
            "SCRIPT_DIR": f"{exp_dir}/{self._SCRIPTS_SUBDIR}",
            "LOG_DIR": f"{exp_dir}/{self._LOG_SUBDIR}",
            "WORK_DIR": f"{exp_dir}/{self._WORK_SUBDIR}",
            "DATA_DIR": f"{exp_dir}/{self._DATA_SUBDIR}",
        }

    @staticmethod
    def _get_repo_root_dir() -> str:
        # mkexp (utils/mkexp/mkexp) is not part of the chilitest package --
        # it is a sibling directory in the enclosing ICON checkout, only
        # reachable by locating that checkout via this file's own location.
        # That location is only meaningful for an editable install
        # ("pip install -e"): a regular install copies this file into
        # site-packages, severing the relationship to the checkout.
        current_file = Path(__file__).resolve()
        try:
            return subprocess.check_output(
                [
                    "git",
                    "-C",
                    str(current_file.parent),
                    "rev-parse",
                    "--show-toplevel",
                ],
                text=True,
                stderr=subprocess.PIPE,
            ).strip()
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Could not locate the enclosing ICON git checkout from "
                f"'{current_file}'. chilitest must be installed editable "
                f"(`pip install -e <path-to-chilitest-folder>`) so this "
                f"file stays inside the checkout it needs to find "
                f"utils/mkexp/mkexp in."
            ) from e

    def _expand_to_dict(self, template: str) -> dict[str, str]:
        return {
            name: template.format(w_exp_name=name)
            for name in self.all_exp_names
        }

    def _append(self, base: dict[str, str], filename: str) -> dict[str, str]:
        return {
            name: f"{base[name]}/{filename.format(w_exp_name=name)}"
            for name in self.all_exp_names
        }

    def _find_output_dirname(self) -> str:
        exp_dir_name = (
            f"{self.single_exp_name}_ref"
            if self._ensemble_run
            else self.single_exp_name
        )
        ref_work_dir_path = Path(self.mkexp_work_dir_by_exp[exp_dir_name])
        matches = [
            p
            for p in ref_work_dir_path.iterdir()
            if p.is_dir() and p.name.startswith("run_")
        ]

        if len(matches) == 0:
            raise FileNotFoundError("No folder matching 'run_*' found")
        if len(matches) > 1:
            raise RuntimeError(
                f"Expected exactly one 'run_*' folder in '{ref_work_dir_path}', found {len(matches)}: "
                f"{[p.name for p in matches]}"
            )

        return matches[0].name
