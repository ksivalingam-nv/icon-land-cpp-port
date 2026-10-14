<!--
ICON

---------------------------------------------------------------
Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
Contact information: icon-model.org

See AUTHORS.TXT for a list of authors
See LICENSES/ for license information
SPDX-License-Identifier: BSD-3-Clause
---------------------------------------------------------------
-->

# Friend-Gated Details Pattern

**Applies to:** C++17 and later

---

## Intent

Compiler-enforced encapsulation of internal state and helpers using a
non-instantiable class with private `static inline` members and inline
friend function definitions. Replaces `namespace detail`, which relies
on naming convention alone.

---

## Structure

```cpp
namespace mylib {

// Forward-declare the public API
void init(SomeData d);
void foo(int x);
bool bar();

class Details {
    Details() = delete;                     // non-instantiable

    static inline SomeData state;           // private state
    static void helper() { /* ... */ }      // private helper

    friend void init(SomeData d) { state = d; }
    friend void foo(int x)  { state.value = x; helper(); }
    friend bool bar()       { return state.flag; }
};

} // namespace mylib
```

**Key properties:**
- Private members are compiler-enforced inaccessible (unlike `namespace detail`)
- Only the listed friend functions get access — explicit and auditable
- `Details` itself has no public interface; callers never name it
- Private static methods (e.g. `helper()`) are inaccessible from outside

---

## Language Mechanics

(References in brackets are ISO C++ standard section identifiers.)

**Friend definitions** — defined inside a class body, they are injected
into the enclosing namespace `[class.friend]`, `[namespace.memdef]/3`.
However, they are only findable via argument-dependent lookup (ADL),
which requires at least one argument whose type belongs to the same
namespace. Functions with no such arguments (e.g. `bool bar()`) would
be invisible to callers despite being defined. The forward declarations
above the class make the names visible to normal lookup; the friend
definitions then match those existing declarations rather than
introducing new hidden ones.

**Static inline members** — C++17 allows defining `static inline` data
members in-class, eliminating the need for a separate `.cpp` file
`[class.static.data]/3`.

---

## Comparison with Alternatives

| | `namespace detail` | This pattern | C++20 modules |
|---|---|---|---|
| Enforcement | Convention | Compiler | Compiler |
| Granularity | Namespace | Per-function | Per-symbol |
| Single header | Yes | Yes | No |
| Min. standard | C++98 | C++17 | C++20 |

---

## When to Use / Not Use

**Use** for small, stable subsystems where inlining matters and
`namespace detail` enforcement is desired — especially in multi-compiler
HPC environments without reliable C++20 module support.

**Avoid** for large subsystems where the friend list would grow unwieldy
(prefer splitting into multiple smaller Details classes), or where
`static inline` state creates thread-safety concerns.

---

## References

- ISO C++17 standard: `[class.friend]`, `[class.static.data]/3`, `[namespace.memdef]/3`
- Stroustrup, *The C++ Programming Language* 4ed, §20.5
