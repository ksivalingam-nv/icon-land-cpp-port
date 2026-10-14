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

// Standalone implementation of icon::finish(). Used when RGK_ENABLE_STANDALONE=ON,
// i.e. when icon_bridge is not compiled and there is no Fortran runtime.
// In non-standalone builds, icon_bridge/finish.cpp provides this function instead.

#include "support/common/finish.hpp"

#include <cstdlib>
#include <iostream>

namespace icon {

[[noreturn]] void finish(std::string_view name, std::string_view text) {
  std::cerr << "FINISH: " << name;
  if (!text.empty()) std::cerr << ": " << text;
  std::cerr << '\n';
  std::abort();
}

[[noreturn]] void finish(std::string_view name) { finish(name, ""); }

}  // namespace icon
