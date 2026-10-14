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
/// @brief Program termination utilities and assert-like functionality.

#ifndef RAGNAROK_SUPPORT_COMMON_FINISH_H_
#define RAGNAROK_SUPPORT_COMMON_FINISH_H_

#include <string>
#include <string_view>

#include "source_location.hpp"

namespace icon {

// Terminates the program with an error message consisting of a short label
// (name) and a longer description (text). The message is passed to the active
// finish handler, which by default writes the message to stderr and aborts.
// When compiled with icon_bridge, the handler can be replaced at startup via
// finish_settings::init() (see icon_bridge/finish.hpp).
[[noreturn]] void finish(std::string_view name, std::string_view text);

// Same as above, but without an error description. Useful when the label
// alone is sufficient to identify the problem.
[[noreturn]] void finish(std::string_view name);

namespace finish_utils {

// Returns a string describing the source location of the caller, in the form
// "file: <file>, line: <line>, fun: <function>". The location defaults to the
// call site, so you can pass the result directly to finish() to include the
// exact position in the source code in the error message.
inline std::string full_str(const rgk::source_location location = rgk::source_location::current()) {
  return "file: " + std::string(location.file_name()) + ", line: " + std::to_string(location.line()) +
         ", fun: " + location.function_name();
}

}  // namespace finish_utils

// assert-like functionality:
// the caller is not supposed to provide the source_location argument,
// the default location argument gets evaluated at the call site

inline void require(bool cond, const rgk::source_location location = rgk::source_location::current()) {
  if (!cond) finish(finish_utils::full_str(location), "require condition is false");
}
inline void require(bool cond, std::string_view reason,
                    const rgk::source_location location = rgk::source_location::current()) {
  if (!cond) finish(finish_utils::full_str(location), reason);
}

// alternative of `require` that is only active in debug mode:
#ifndef NDEBUG
// debug mode: enable checks
inline void debug_require(bool cond, const rgk::source_location location = rgk::source_location::current()) {
  if (!cond) finish(finish_utils::full_str(location), "debug_require condition is false");
}
inline void debug_require(bool cond, std::string_view reason,
                          const rgk::source_location location = rgk::source_location::current()) {
  if (!cond) finish(finish_utils::full_str(location), reason);
}
#else
// production mode: provide no-op dummy implementations of checks
inline void debug_require([[maybe_unused]] bool) {}
inline void debug_require([[maybe_unused]] bool, [[maybe_unused]] std::string_view) {}
#endif

}  // namespace icon

#endif  // RAGNAROK_SUPPORT_COMMON_FINISH_H_
