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

///
/// @file
/// @brief initializes icon_bridge
///

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_BRIDGE_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_BRIDGE_H_

#include "icon_f2c.hpp"

namespace icon {

extern "C" void init_ragnarok_f2c(const f2c::FunTable* fun_table);

}  // namespace icon

#endif /* RAGNAROK_SUPPORT_ICON_BRIDGE_ICON_BRIDGE_H_ */
