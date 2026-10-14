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

#ifndef RAGNAROK_SUPPORT_COMMON_SOURCE_LOCATION_H_
#define RAGNAROK_SUPPORT_COMMON_SOURCE_LOCATION_H_

#include <cstdint>

// detect if the standard library natively provides <source_location>
#if __has_include(<version>)
#include <version>
#endif

#if defined(__cpp_lib_source_location) && __has_include(<source_location>)
// --- environment has native C++20 source location ---
#include <source_location>
namespace rgk {
// we use rgk as namespace to access the appropriate implementation for our environment
using source_location = std::source_location;
}  // namespace rgk

#else
// --- environment does not have native source location ---
// fall back to our implementation using compiler built-ins
#include <string_view>

namespace rgk {
struct source_location {
 public:
  // __builtin_COLUMN() is supported by Clang, ICX (oneAPI), and GCC >= 11.
  // classic icc (__INTEL_COMPILER) supports __builtin_FILE/LINE/FUNCTION but
  // not __builtin_COLUMN(), so it correctly falls through to col = 0.
  // nvc++ defines __GNUC__ for compatibility but does NOT support
  // __builtin_COLUMN(), so we explicitly exclude it via __NVCOMPILER.
  [[nodiscard]] static constexpr source_location current(const char* file         = __builtin_FILE(),
                                                         const char* func         = __builtin_FUNCTION(),
                                                         std::uint_least32_t line = __builtin_LINE(),
#if defined(__clang__) || defined(__INTEL_LLVM_COMPILER) || \
    (defined(__GNUC__) && __GNUC__ >= 11 && !defined(__NVCOMPILER))
                                                         std::uint_least32_t col = __builtin_COLUMN()
#else
                                                         std::uint_least32_t col = 0
#endif
                                                             ) noexcept {
    source_location loc;
    loc.m_file = file;
    loc.m_func = func;
    loc.m_line = line;
    loc.m_col  = col;
    return loc;
  }

  constexpr source_location() noexcept = default;

  [[nodiscard]] constexpr std::uint_least32_t line() const noexcept { return m_line; }
  [[nodiscard]] constexpr std::uint_least32_t column() const noexcept { return m_col; }
  [[nodiscard]] constexpr const char* file_name() const noexcept { return m_file; }
  [[nodiscard]] constexpr const char* function_name() const noexcept { return m_func; }

 private:
  const char* m_file         = "";
  const char* m_func         = "";
  std::uint_least32_t m_line = 0;
  std::uint_least32_t m_col  = 0;
};
}  // namespace rgk
#endif

#endif  // RAGNAROK_SUPPORT_COMMON_SOURCE_LOCATION_H_
