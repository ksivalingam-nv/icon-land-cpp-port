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

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_LOGGING_H_

/// @file
/// Logging interface for C++/Fortran integration: Fortran callback registration and message forwarding.
///

#define RAGNAROK_SUPPORT_ICON_BRIDGE_LOGGING_H_

#include <ostream>
#include <string_view>

namespace icon {

namespace logging {

using f2c_msg_fun_type = void (*)(const char* name, int name_len, const char* text, int text_len);
void init(f2c_msg_fun_type f2c_message);
bool is_initialized();
void set_stdio_proc(bool is_stdio);

}  // namespace logging

void message(std::string_view name, std::string_view text);

}  // namespace icon

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_LOGGING_H_
