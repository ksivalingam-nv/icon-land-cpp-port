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

import getpass
import os
import subprocess
import time

from chilitest.ensemble import Ensemble
from chilitest.file_info import FileInfo
from numpy import nditer
from numpy.random import PCG64, Generator


class Experiment:
    def __init__(
        self,
        file_info: FileInfo,
        options: dict[str, str],
        ensemble: Ensemble = None,
        restart: bool = False,
    ):
        self.file_info = file_info
        self.ensemble = ensemble
        self.restart = restart

        # The source/build paths recorded in the generated run/SETUP.config
        # are stale when a build is relocated (e.g. restored from CI
        # artifacts), so pin them to the paths chilitest actually uses.
        self._mkexp_default_args = {
            "namelists.NAMELIST_atm.io_nml.lnetcdf_flt64_output": ".TRUE.",
            "abs_top_srcdir": self.file_info.repo_root,
            "abs_top_builddir": self.file_info.build_dir,
            **self.file_info.mkexp_dir_args(),
        }
        self._mkexp_default_exp_options = "rbf_coeffs/read"

        self.mkexp_args = self._get_mkexp_args(options)

    def generate_runfiles(self) -> None:
        os.environ["MKEXP_PATH"] = self.file_info.run_dir

        # EXP_OPTIONS given in <exp>.config must not be overriden
        mkexp_cmds = {
            name: [
                self.file_info.mkexp_path,
                self.file_info.mkexp_config_path,
                f"EXP_ID={name}",
                f"WRITE_COMPRESSION= ",
                *[
                    (
                        f"{option}+={value}"
                        if option == "EXP_OPTIONS"
                        else f"{option}={value}"
                    )
                    for option, value in self.mkexp_args.items()
                ],
            ]
            for name in self.file_info.all_exp_names
        }

        if self.ensemble is not None:
            if self.ensemble.perturbation_method == "pinit_nml_params":
                # By adding values to `pinit_seed` and `pinit_amplitude` (ICON builtin perturbation method),
                # a group of internal variables are perturbed with the *same* `perturbation_amplitude`.
                rng = Generator(PCG64(2394875))
                pinit_seeds = nditer(
                    rng.integers(low=1e10, high=1e11, size=self.ensemble.size)
                )

                for name, cmd in mkexp_cmds.items():
                    if not name.endswith("_ref"):
                        cmd.append(
                            f"namelists.NAMELIST_atm.initicon_nml.pinit_seed={next(pinit_seeds)}"
                        )
                        cmd.append(
                            f"namelists.NAMELIST_atm.initicon_nml.pinit_amplitude={self.ensemble.perturbation_amplitude}"
                        )

            elif self.ensemble.perturbation_method == "input_netcdf":
                # Perturb input netcdf files with the *same* `perturbation_amplitude`
                raise NotImplementedError
            else:
                # Unknown values are rejected by Ensemble against
                # registry.PERTURBATION_METHODS; reaching this branch means
                # the registry and the dispatch above drifted apart.
                raise RuntimeError(
                    "Perturbation method registered but not implemented: "
                    f"{self.ensemble.perturbation_method}"
                )

        for cmd in mkexp_cmds.values():
            subprocess.run(cmd, check=True)

        # New experiments will not run if `<exp_name>.date` files present in the run folders.
        # This happens when calling `generate_runfiles` more than once for the same experiment
        # (mkexp keeps old files by default)
        if not self.restart:
            for date_file in self.file_info.date_file_by_exp.values():
                if os.path.exists(date_file):
                    os.remove(date_file)

    def run(self) -> None:
        job_ids = []

        for (
            exp_name,
            run_start_file,
        ) in self.file_info.run_start_file_by_exp.items():
            cmd = ["sbatch", run_start_file]
            process = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )

            out, err = process.communicate()

            if process.returncode != 0:
                raise subprocess.CalledProcessError(
                    process.returncode, cmd, out, err
                )

            job_id = out.strip().split()[-1]
            job_ids.append(job_id)
            print(f"Submitted {exp_name} as job {job_id}")

        unfinished_jobs = set(job_ids)

        # getpass.getuser, unlike os.getlogin, also works without a
        # controlling terminal (batch jobs, CI runners).
        slurm_user = getpass.getuser()
        while unfinished_jobs:
            time.sleep(10)
            cmd = ["squeue", "-u", slurm_user, "-o", "%i"]

            result = subprocess.run(
                cmd, stdout=subprocess.PIPE, text=True, check=True
            )
            running_jobs = set(result.stdout.split())

            unfinished_jobs = set(job_ids) & running_jobs

        print("All jobs finished")

        self._check_experiment_success()

    def _get_mkexp_args(self, cli_options: str) -> dict[str, str]:
        # Merge input options with mkexp defaults, error if modified defaults
        # `EXP_ID` given in `generate_runfiles`, it is ensemble member dependent
        mkexp_args_overlap = cli_options.keys() & (
            self._mkexp_default_args.keys() | {"EXP_ID"}
        )
        if mkexp_args_overlap:
            raise ValueError(
                f"Attempting to modify the following default mkexp options: {sorted(mkexp_args_overlap)}"
            )

        # Add default `EXP_OPTIONS` to those provided by cli (if any)
        mkexp_args = self._mkexp_default_args.copy()
        cli_options_copy = cli_options.copy()
        if "EXP_OPTIONS" in cli_options:
            cli_options_copy[
                "EXP_OPTIONS"
            ] += f", {self._mkexp_default_exp_options}"
        else:
            cli_options_copy["EXP_OPTIONS"] = self._mkexp_default_exp_options

        mkexp_args.update(cli_options_copy)

        return mkexp_args

    def _check_experiment_success(self) -> None:
        def _experiment_succeeded(filename):
            marker = b"Experiment finished successfully"
            # Generous margin over any plausible run-log epilogue (SLURM
            # accounting footer, module-unload noise, ...); a line-count
            # cap missed the marker whenever more than ~100 lines followed
            # it.
            max_tail_bytes = 1024 * 1024

            with open(filename, "rb") as f:
                f.seek(0, 2)  # jump to end
                pos = f.tell()
                data = b""

                while pos > 0 and len(data) < max_tail_bytes:
                    read_size = min(8192, pos)  # read 8 KB chunks
                    pos -= read_size
                    f.seek(pos)
                    data = f.read(read_size) + data

            return marker in data

        failures = []
        for exp_name, log_file in self.file_info.log_file_by_exp.items():
            if not _experiment_succeeded(log_file):
                failures.append(f"Experiment {exp_name} failed ({log_file})")

        if failures:
            raise RuntimeError(
                "One or more experiments failed:\n" + "\n".join(failures)
            )
        print("All jobs succeeded")

    def save_version_info(self) -> None:
        version_info = subprocess.run(
            [self.file_info.exp_executable, "--version"],
            capture_output=True,
            text=True,
            check=True,
        )

        with open(self.file_info.version_file, "w") as f:
            f.write(version_info.stderr)
