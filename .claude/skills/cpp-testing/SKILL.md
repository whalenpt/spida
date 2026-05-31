---
name: cpp-testing
description: >
  Testing workflow for modern C++ (C++17/20/23) using GoogleTest/GoogleMock with
  CMake/CTest. Use when writing/updating/fixing C++ tests, configuring GoogleTest/CTest,
  diagnosing failing or flaky tests, or adding coverage/sanitizers. For language style and
  idioms, use the companion cpp-coding-standards skill. Triggers: "write a test", "TEST_F",
  "gmock", "ctest", "flaky test", "coverage", "ASan/UBSan/TSan".
---

# C++ Testing (Agent Skill)

Agent-focused testing workflow for modern C++ (C++17/20/23) using GoogleTest/GoogleMock with CMake/CTest.

> **Project conventions win.** If the project defines its own build, test, coverage, or
> test-registration commands (e.g. `CLAUDE.md`, `CMakePresets.json`, a `Makefile`, an
> `add_<project>_test()` helper), those are authoritative — use them verbatim. Treat the
> generic invocations here as a fallback only for projects with no such convention. In
> particular: if the project is preset-driven, build and test through `cmake --preset` /
> `ctest --preset`, never raw `cmake -B`/`--build`/`ctest --test-dir`.

## When to Use

- Writing new C++ tests or fixing existing tests
- Designing unit/integration test coverage for C++ components
- Adding test coverage, CI gating, or regression protection
- Configuring CMake/CTest workflows for consistent execution
- Investigating test failures or flaky behavior
- Enabling sanitizers for memory/race diagnostics

### When NOT to Use

- Implementing new product features without test changes
- Large-scale refactors unrelated to test coverage or failures
- Performance tuning without test regressions to validate
- Non-C++ projects or non-test tasks
- Language style, idioms, or API design — use the companion `cpp-coding-standards` skill

## Core Concepts

- **TDD loop**: red → green → refactor (tests first, minimal fix, then cleanups).
- **Isolation**: prefer dependency injection and fakes over global state.
- **Test layout**: `tests/unit`, `tests/integration`, `tests/testdata`.
- **Mocks vs fakes**: mock for interactions, fake for stateful behavior.
- **CTest discovery**: use `gtest_discover_tests()` for stable test discovery.
- **CTest labels**: tag tests `unit`/`integration` so CI can run fast subsets first.
- **Assertions on results**: `EXPECT_THAT` + matchers for expressive checks; `EXPECT_THROW` for exceptions; inspect `std::expected` directly for error-as-value APIs.
- **CI signal**: run subset first, then full suite with `--output-on-failure`.

## TDD Workflow

Follow the RED → GREEN → REFACTOR loop:

1. **RED**: write a failing test that captures the new behavior
2. **GREEN**: implement the smallest change to pass
3. **REFACTOR**: clean up while tests stay green

```cpp
// tests/add_test.cpp
#include <gtest/gtest.h>

int Add(int a, int b); // Provided by production code.

TEST(AddTest, AddsTwoNumbers) { // RED
  EXPECT_EQ(Add(2, 3), 5);
}

// src/add.cpp
int Add(int a, int b) { // GREEN
  return a + b;
}

// REFACTOR: simplify/rename once tests pass
```

## Code Examples

### Basic Unit Test (gtest)

```cpp
// tests/calculator_test.cpp
#include <gtest/gtest.h>

int Add(int a, int b); // Provided by production code.

TEST(CalculatorTest, AddsTwoNumbers) {
    EXPECT_EQ(Add(2, 3), 5);
}
```

### Fixture (gtest)

```cpp
// tests/user_store_test.cpp
// Pseudocode stub: replace UserStore/User with project types.
#include <gtest/gtest.h>
#include <memory>
#include <optional>
#include <string>

struct User { std::string name; };
class UserStore {
public:
    explicit UserStore(std::string /*path*/) {}
    void Seed(std::initializer_list<User> /*users*/) {}
    std::optional<User> Find(const std::string &/*name*/) { return User{"alice"}; }
};

class UserStoreTest : public ::testing::Test {
protected:
    void SetUp() override {
        store = std::make_unique<UserStore>(":memory:");
        store->Seed({{"alice"}, {"bob"}});
    }

    std::unique_ptr<UserStore> store;
};

TEST_F(UserStoreTest, FindsExistingUser) {
    auto user = store->Find("alice");
    ASSERT_TRUE(user.has_value());
    EXPECT_EQ(user->name, "alice");
}
```

### Mock (gmock)

```cpp
// tests/notifier_test.cpp
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <string>

class Notifier {
public:
    virtual ~Notifier() = default;
    virtual void Send(const std::string &message) = 0;
};

class MockNotifier : public Notifier {
public:
    MOCK_METHOD(void, Send, (const std::string &message), (override));
};

class Service {
public:
    explicit Service(Notifier &notifier) : notifier_(notifier) {}
    void Publish(const std::string &message) { notifier_.Send(message); }

private:
    Notifier &notifier_;
};

TEST(ServiceTest, SendsNotifications) {
    MockNotifier notifier;
    Service service(notifier);

    EXPECT_CALL(notifier, Send("hello")).Times(1);
    service.Publish("hello");
}
```

### Matchers (`EXPECT_THAT`)

Prefer matchers over chains of `EXPECT_EQ` — failures print the expected *shape*, not just
a boolean. Matchers ship with gmock and work in plain gtest tests.

```cpp
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::ElementsAre;
using ::testing::HasSubstr;
using ::testing::Optional;
using ::testing::UnorderedElementsAre;

TEST(MatchersTest, ReadsExpressively) {
    const std::vector<int> v{1, 2, 3};
    EXPECT_THAT(v, ElementsAre(1, 2, 3));            // ordered
    EXPECT_THAT(v, UnorderedElementsAre(3, 1, 2));   // set-like
    EXPECT_THAT(std::string{"alice@example.com"}, HasSubstr("@"));
    EXPECT_THAT(std::optional<int>{5}, Optional(5));
}
```

### Parameterized Tests (`TEST_P`)

Use a value-parameterized suite to cover a table of cases without copy-pasting tests.

```cpp
#include <algorithm>
#include <gtest/gtest.h>

struct ClampCase { int input; int lo; int hi; int want; };

class ClampTest : public ::testing::TestWithParam<ClampCase> {};

TEST_P(ClampTest, ClampsIntoRange) {
    const auto &c = GetParam();
    EXPECT_EQ(std::clamp(c.input, c.lo, c.hi), c.want);
}

INSTANTIATE_TEST_SUITE_P(Bounds, ClampTest, ::testing::Values(
    ClampCase{-1, 0, 10, 0},
    ClampCase{ 5, 0, 10, 5},
    ClampCase{99, 0, 10, 10}
));
```

### Testing Exceptions and `std::expected`

Match the assertion to the error strategy of the code under test (see the Error Handling
section of `cpp-coding-standards`).

```cpp
#include <expected>
#include <gtest/gtest.h>

// Exception-based API: assert the type that propagates
TEST(ConnectTest, ThrowsOnUnreachable) {
    EXPECT_THROW(connect(Endpoint{"203.0.113.1"}), NetworkError);
}

// std::expected API (C++23): inspect value/error directly, no try/catch
TEST(ParseTest, RejectsEmpty) {
    const auto r = to_int("");
    ASSERT_FALSE(r.has_value());
    EXPECT_EQ(r.error(), ParseError::empty);
}

TEST(ParseTest, ParsesDecimal) {
    const auto r = to_int("42");
    ASSERT_TRUE(r.has_value());   // ASSERT: stop if false, *r below would be UB
    EXPECT_EQ(*r, 42);
}
```

### Death Tests

Verify that a precondition/contract violation aborts. Suite name must end in `DeathTest`.

```cpp
#include <gtest/gtest.h>

TEST(VectorAccessDeathTest, AbortsOnOutOfRange) {
    std::vector<int> v;
    EXPECT_DEATH(checked_at(v, 0), "index out of range");  // matches stderr regex
}
```

### Skipping at Runtime

```cpp
TEST(GpuTest, RunsKernel) {
    if (!gpu_available()) GTEST_SKIP() << "no GPU on this runner";
    // ...
}
```

### CMake/CTest Quickstart

Provision GoogleTest from a package manager (Conan 2 / vcpkg) when one is in use, with an
in-tree submodule fallback. This is the shape most projects ship; reach for `FetchContent`
only when there is no package manager and no vendored submodule.

```cmake
# CMakeLists.txt (excerpt)
cmake_minimum_required(VERSION 3.20)
project(example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)        # default to C++23; pin lower per project policy
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

# Preferred: package-manager-provided (Conan 2 / vcpkg), else vendored submodule.
find_package(GTest CONFIG QUIET)
if(GTest_FOUND)
    message(STATUS "Found GTest package")
else()
    message(STATUS "Building googletest from submodule")
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    add_subdirectory(external/googletest)   # requires clone --recurse-submodules
endif()

# Fallback ONLY when neither a package nor a submodule is available:
#   include(FetchContent)
#   set(GTEST_VERSION v1.17.0)   # pin per project policy
#   FetchContent_Declare(googletest
#     URL https://github.com/google/googletest/archive/refs/tags/${GTEST_VERSION}.zip)
#   FetchContent_MakeAvailable(googletest)

add_executable(example_tests
  tests/calculator_test.cpp
  src/calculator.cpp
)
# gtest_main provides main(); do not write a custom one.
target_link_libraries(example_tests PRIVATE GTest::gtest GTest::gmock GTest::gtest_main)

include(GoogleTest)
# Tag for fast CI subsets: ctest -L unit
gtest_discover_tests(example_tests PROPERTIES LABELS "unit")
```

> If the project provides a registration helper (e.g. `add_spida_test(name name.cpp)` that
> links `GTest::gtest_main` and the project library automatically), use it instead of a raw
> `add_executable` + `target_link_libraries` + `gtest_discover_tests` block, and follow the
> project's `TEST(FEATURE_TEST, CASE_NAME)` naming convention.

```bash
# Generic (no preset convention):
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure

# Preset-driven project (preferred when CMakePresets.json exists):
cmake --preset <configure-preset>
cmake --build --preset <build-preset> --parallel
ctest --preset <test-preset> --output-on-failure
```

## Running Tests

```bash
ctest --test-dir build --output-on-failure
ctest --test-dir build -R ClampTest                      # by test name regex
ctest --test-dir build -L unit                           # by label
ctest --test-dir build -R "UserStoreTest.*" --output-on-failure
# Preset equivalent: ctest --preset <test-preset> -L unit --output-on-failure
```

```bash
./build/example_tests --gtest_filter=ClampTest.*
./build/example_tests --gtest_filter=UserStoreTest.FindsExistingUser
./build/example_tests --gtest_repeat=100 --gtest_shuffle    # surface order-dependent flakes
```

## Debugging Failures

1. Re-run the single failing test with a gtest filter.
2. Add scoped logging around the failing assertion.
3. Run it under GDB: `gdb --args ./build/example_tests --gtest_filter=Suite.Case`.
4. Re-run with sanitizers enabled (ASan/UBSan).
5. Expand to the full suite once the root cause is fixed.

## Coverage

**Prefer target-level instrumentation** so coverage flags stay scoped to the code under
test and don't leak into dependencies. Some projects deliberately set coverage flags
**globally** (e.g. `add_compile_options(--coverage -O0 -g -fprofile-update=atomic)`) to
guarantee consistent instrumentation across all translation units, including when tests
run in parallel — when the project file does this, it is authoritative; don't fight it.

Note `-fprofile-update=atomic`: it prevents profile-counter races when tests execute in
parallel, at a small runtime cost. Include it (target- or global-level) whenever coverage
runs concurrently.

```cmake
option(ENABLE_COVERAGE "Enable coverage flags" OFF)

if(ENABLE_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # -fprofile-update=atomic: safe counters under parallel test execution
    target_compile_options(example_tests PRIVATE --coverage -O0 -g -fprofile-update=atomic)
    target_link_options(example_tests PRIVATE --coverage)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(example_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    target_link_options(example_tests PRIVATE -fprofile-instr-generate)
  endif()
endif()
```

GCC + gcov + lcov:

```bash
cmake -S . -B build-cov -DENABLE_COVERAGE=ON
cmake --build build-cov -j
ctest --test-dir build-cov

# Capture, then filter noise. Two independent fixes for version skew (use either/both):
#   (a) point lcov at the gcov matching the compiler:    --gcov-tool gcov-13
#   (b) tolerate format mismatches from a newer toolchain: --ignore-errors mismatch,gcov
lcov --capture --directory build-cov --gcov-tool gcov-13 \
     --output-file coverage.info --ignore-errors mismatch,gcov
lcov --remove coverage.info '/usr/*' '*/external/*' '*/test/*' '*/build/*' \
     --output-file coverage.filtered.info --ignore-errors unused
lcov --list coverage.filtered.info
```

Clang + llvm-cov:

```bash
cmake -S . -B build-llvm -DENABLE_COVERAGE=ON -DCMAKE_CXX_COMPILER=clang++
cmake --build build-llvm -j
LLVM_PROFILE_FILE="build-llvm/default.profraw" ctest --test-dir build-llvm
llvm-profdata merge -sparse build-llvm/default.profraw -o build-llvm/default.profdata
llvm-cov report build-llvm/example_tests -instr-profile=build-llvm/default.profdata
```

**Gating:** pick a project coverage threshold and fail CI below it rather than tracking a
trend by eye. `lcov --list`/`--summary` and `llvm-cov report` print rates a CI step can
parse. Keep coverage on a single, consistent build configuration — mixing
debug/release or compilers skews numbers.

## Sanitizers

```cmake
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(ENABLE_TSAN "Enable ThreadSanitizer" OFF)

if(ENABLE_ASAN)
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address)
endif()
if(ENABLE_UBSAN)
  add_compile_options(-fsanitize=undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=undefined)
endif()
if(ENABLE_TSAN)
  add_compile_options(-fsanitize=thread)
  add_link_options(-fsanitize=thread)
endif()
```

**Combination rules:** ASan and UBSan run together; TSan must run **alone** (incompatible
with ASan). Use separate build dirs per configuration. In a Conan/preset project, model
each sanitizer build as its own profile/preset rather than ad-hoc flags on an existing
build dir. Make failures hard-fail in CI:

```bash
cmake -S . -B build-asan -DENABLE_ASAN=ON -DENABLE_UBSAN=ON
cmake --build build-asan -j
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ctest --test-dir build-asan --output-on-failure

cmake -S . -B build-tsan -DENABLE_TSAN=ON
cmake --build build-tsan -j
TSAN_OPTIONS=halt_on_error=1 ctest --test-dir build-tsan --output-on-failure
```

Prefer sanitizers over Valgrind when both are available (far faster, better diagnostics);
keep Valgrind as the fallback for targets where sanitizers aren't supported.

## Flaky Tests Guardrails

- Never use `sleep` for synchronization; use condition variables or latches.
- Make temp directories unique per test and always clean them.
- Avoid real time, network, or filesystem dependencies in unit tests.
- Use deterministic seeds for randomized inputs.
- Run `--gtest_shuffle --gtest_repeat=N` locally to catch order-dependence early.

## Best Practices

### DO

- Keep tests deterministic and isolated
- Prefer dependency injection over globals
- Use `ASSERT_*` for preconditions, `EXPECT_*` for multiple checks
- Use `EXPECT_THAT` + matchers for collections, strings, and optionals
- Separate unit vs integration tests in CTest labels or directories
- Run sanitizers in CI for memory and race detection

### DON'T

- Don't depend on real time or network in unit tests
- Don't use sleeps as synchronization when a condition variable can be used
- Don't over-mock simple value objects
- Don't use brittle string matching for non-critical logs
- Don't dereference an `std::expected`/`std::optional` after an `EXPECT_*` — use `ASSERT_*` so a failure stops before the UB

### Common Pitfalls

- **Using fixed temp paths** → Generate unique temp directories per test and clean them.
- **Relying on wall clock time** → Inject a clock or use fake time sources.
- **Flaky concurrency tests** → Use condition variables/latches and bounded waits.
- **Hidden global state** → Reset global state in fixtures or remove globals.
- **Over-mocking** → Prefer fakes for stateful behavior and only mock interactions.
- **Missing sanitizer runs** → Add ASan/UBSan and TSan builds in CI (separately).
- **Coverage on debug-only builds** → Ensure coverage targets use consistent flags.
- **lcov/gcov version skew** → Pass `--gcov-tool` for the matching gcov and/or `--ignore-errors mismatch,gcov`.
- **Coverage counter races under `-j`** → Build coverage with `-fprofile-update=atomic`.

## Optional Appendix: Fuzzing / Property Testing

Only use if the project already supports LLVM/libFuzzer or a property-testing library.

- **libFuzzer**: best for pure functions with minimal I/O.
- **RapidCheck**: property-based tests to validate invariants.

Minimal libFuzzer harness (pseudocode: replace ParseConfig):

```cpp
#include <cstddef>
#include <cstdint>
#include <string>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    std::string input(reinterpret_cast<const char *>(data), size);
    // ParseConfig(input); // project function
    return 0;
}
```

## Quick Reference Checklist

Before marking C++ test work complete:

- [ ] Used the project's own build/test/registration conventions where they exist (presets, `add_<project>_test()`, naming)
- [ ] New behavior driven by a failing test first (RED → GREEN → REFACTOR)
- [ ] Tests are deterministic — no real time, network, or fixed temp paths
- [ ] `ASSERT_*` guards anything dereferenced afterward; `EXPECT_*` for independent checks
- [ ] Collections/strings/optionals checked with `EXPECT_THAT` + matchers
- [ ] Table-driven cases use `TEST_P` rather than copy-pasted tests
- [ ] Error paths tested per strategy: `EXPECT_THROW` for exceptions, value/error checks for `std::expected`
- [ ] Interactions mocked, stateful behavior faked (not over-mocked)
- [ ] Registered via `gtest_discover_tests()` (or project helper) with `unit`/`integration` labels
- [ ] Full suite green via `ctest --output-on-failure` (preset or `--test-dir`)
- [ ] Passes under ASan+UBSan; TSan run separately for concurrent code
- [ ] Coverage built with `-fprofile-update=atomic` when run in parallel
- [ ] Coverage measured on one consistent config and gated against a threshold
