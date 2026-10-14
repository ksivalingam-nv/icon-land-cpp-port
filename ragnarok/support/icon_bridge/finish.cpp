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

#include "finish.hpp"

#include <cstdlib>
#include <iostream>

extern "C" {
// This function is implemented in Fortran to write to a specified default Fortran unit.
void ragnarok_write_to_default_funit(const char* cname, int cname_len, const char* ctext, int ctext_len);
}

namespace {
// The default finish handler, which simply writes the error message and aborts the program.
// It will be replaced by the Fortran side at startup via finish_settings::init().
icon::finish_settings::f2c_fin_fun_type finish_impl = [](const char* name, int name_len, const char* text,
                                                         int text_len) {
  ragnarok_write_to_default_funit(name, name_len, text, text_len);
  std::abort();
};

bool init_state = false;

}  // namespace

namespace icon {

[[noreturn]] void finish(std::string_view name, std::string_view text) {
  finish_impl(name.data(), static_cast<int>(name.size()), text.data(), static_cast<int>(text.size()));
  std::abort();  // unreachable; tells the compiler this function never returns
}

[[noreturn]] void finish(std::string_view name) { finish(name, ""); }

namespace finish_settings {

void init(f2c_fin_fun_type f2c_fin_fun) {
  require(!init_state);
  finish_impl = f2c_fin_fun;
  init_state  = true;
}

bool is_initialized() { return init_state; }

}  // namespace finish_settings

}  // namespace icon
