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
from chilitest.registry import PERTURBATION_METHODS


class Ensemble:
    def __init__(
        self,
        perturbation_method: str,
        size: int,  # This number does not include the reference member of the ensemble
        perturbation_amplitude: float,
    ) -> None:
        if perturbation_method not in PERTURBATION_METHODS:
            raise ValueError(
                f"Unknown perturbation method: {perturbation_method}"
                f" (available: {', '.join(PERTURBATION_METHODS)})"
            )
        self.perturbation_method = perturbation_method
        self.size = size
        self.perturbation_amplitude = np.float64(perturbation_amplitude)
