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

#include "icon_bridge.hpp"

#include <cassert>

#include "finish.hpp"
#include "icon_domain.hpp"
#include "icon_f2c.hpp"
#include "logging.hpp"
#include "ragnarok.hpp"

namespace icon {

extern "C" void init_ragnarok_f2c(const f2c::FunTable* fun_table) {
  require(ragnarok::is_initialized());
  // basic bridge infrastructure:
  logging::init(fun_table->message);
  finish_settings::init(fun_table->finish);
  f2c::init(fun_table);
  // higher level states:
  domain::init();
}

}  // namespace icon
