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
from chilitest.ensemble import Ensemble
from chilitest.experiment import Experiment
from chilitest.file_info import FileInfo
from chilitest.statistics import EnsembleStatistics

application_config = ApplicationConfig.load_with_cli(
    description="Generate ICON tolerance test references: run a perturbed "
    "ensemble with the reference build, postprocess the selected output "
    "variables, and store the ensemble statistics for future tests."
)

validate_config(
    application_config,
    {
        "ensemble": (Ensemble, set()),
        "file_info": (FileInfo, {"ensemble"}),
        "experiment": (Experiment, {"file_info", "ensemble"}),
        "ensemble_statistics": (EnsembleStatistics, {"file_info"}),
    },
)

ensemble = Ensemble(**application_config["ensemble"])

file_info = FileInfo(ensemble=ensemble, **application_config["file_info"])

ensemble_experiment = Experiment(
    file_info=file_info,
    ensemble=ensemble,
    **application_config["experiment"],
)

ensemble_experiment.generate_runfiles()
ensemble_experiment.run()

stats = EnsembleStatistics(
    file_info=file_info, **application_config["ensemble_statistics"]
)

stats.open_dataset()
stats.compute_spatial_averages()
stats.generate_ensemble_statistics()

stats.save_ensemble_statistics()

ensemble_experiment.save_version_info()
