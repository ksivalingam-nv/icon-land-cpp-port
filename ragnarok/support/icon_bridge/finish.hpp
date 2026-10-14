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

/// @file
/// @brief Fortran-bridge finish handler settings (finish, require, debug_require
///        are provided by support/common/finish.hpp).

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_FINISH_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_FINISH_H_

#include "support/common/finish.hpp"

namespace icon {

namespace finish_settings {

// Function pointer type for the finish handler used to bridge between C++ and
// Fortran. The handler receives a name and a text string, each as a raw char
// pointer with an explicit length (not null-terminated). It is expected to
// report the error and terminate the program.
using f2c_fin_fun_type = void (*)(const char* name, int name_len, const char* text, int text_len);

// Replaces the default finish handler with a custom one. This should be called
// once at program startup, typically by the Fortran side (mo_finish_utils.f90),
// to register its own error-reporting routine. Calling this function more than
// once is an error and will itself trigger a fatal abort.
void init(f2c_fin_fun_type f2c_finish);

bool is_initialized();

}  // namespace finish_settings

}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_FINISH_H_
