// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------
#ifndef RAGNAROK_COMMON_REAL_H_
#define RAGNAROK_COMMON_REAL_H_

namespace icon {

// working precision:
#ifdef __SINGLE_PRECISION
using Real = float;
#else
using Real = double;
#endif

}  // namespace icon

#endif  // RAGNAROK_COMMON_REAL_H_
