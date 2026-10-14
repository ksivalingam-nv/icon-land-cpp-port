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

#ifndef RAGNAROK_COMMON_ENTITY_H_
#define RAGNAROK_COMMON_ENTITY_H_

namespace icon {

// Grid entity a field or index range lives on. The integer values mirror the
// conventions used on the ICON Fortran side.
enum class Entity : int { cells = 1, edges = 2, vertices = 3 };

}  // namespace icon

#endif  // RAGNAROK_COMMON_ENTITY_H_
