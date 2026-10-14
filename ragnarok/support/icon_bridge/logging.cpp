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

#include "logging.hpp"

extern "C" {
void ragnarok_write_to_default_funit(const char* cname, int cname_len, const char* ctext, int ctext_len);
}

namespace {  // anon

bool is_stdio                                    = false;

// provide default implementation for use before init():
icon::logging::f2c_msg_fun_type f2c_message_impl = [](const char* name, int name_len, const char* text, int text_len) {
  if (is_stdio) {
    ragnarok_write_to_default_funit(name, name_len, text, text_len);
  }
};

bool init_state = false;

}  // namespace

namespace icon {

namespace logging {

void set_stdio_proc(bool stdio) { is_stdio = stdio; }

void init(f2c_msg_fun_type f2c_message_fun) {
  if (init_state) return;
  f2c_message_impl = f2c_message_fun;
  init_state       = true;
}

bool is_initialized() { return init_state; }

}  // namespace logging

// this is the message interface we expose to the icon namespace:
void message(std::string_view name, std::string_view text) {
  f2c_message_impl(name.data(), static_cast<int>(name.size()), text.data(), static_cast<int>(text.size()));
}

}  // namespace icon
