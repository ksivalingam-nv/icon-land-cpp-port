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

from chilitest.application_config import ApplicationConfig, validate_config
from chilitest.experiment import Experiment
from chilitest.file_info import FileInfo
from chilitest.statistics import EnsembleStatistics, Statistics
from chilitest.test import Test

application_config = ApplicationConfig.load_with_cli(
    description="Run an ICON tolerance test: run the experiment with the "
    "build under test, postprocess the selected output variables, and test "
    "them against stored reference ensemble statistics."
)

validate_config(
    application_config,
    {
        "file_info": (FileInfo, {"ensemble"}),
        "experiment": (Experiment, {"file_info", "ensemble"}),
        "statistics": (Statistics, {"file_info"}),
        "test.extreme_values": (Test.extreme_values, set()),
    },
)

file_info = FileInfo(**application_config["file_info"])

experiment = Experiment(file_info=file_info, **application_config["experiment"])
experiment.generate_runfiles()
experiment.run()

test_stats = Statistics(file_info=file_info, **application_config["statistics"])
test_stats.open_dataset()
test_stats.compute_spatial_averages()

ref_stats = EnsembleStatistics.load(file_info=file_info)

test = Test(test_statistics=test_stats, ref_statistics=ref_stats)

test.extreme_values(**application_config["test"]["extreme_values"])
