# SPIDA — Claude Code project guide

Spectral Integration and Differentiation Algorithms (C++17 library).

## Build

Dependencies are managed with Conan 2. The CI image (`ghcr.io/whalenpt/spida:v1.4`) has
Conan pre-configured and packages pre-cached, so `--build=missing` is instant there.

```bash
# Install deps + configure + build (Release)
conan install . -of build/Release -s build_type=Release
cmake --preset conan-release
cmake --build --preset conan-release --parallel

# Same for Coverage (RelWithDebInfo + gcov instrumentation)
conan install . -of build/RelWithDebInfo -s build_type=RelWithDebInfo
cmake --preset conan-relwithdebinfo -DSPIDA_TEST=ON -DSPIDA_COVERAGE=ON
cmake --build --preset conan-relwithdebinfo --parallel
```

## Run tests

```bash
ctest --preset conan-release --output-on-failure
ctest --preset conan-relwithdebinfo --output-on-failure
```

## Generate coverage report

```bash
lcov --capture --directory build/RelWithDebInfo --output-file coverage.info \
     --ignore-errors mismatch,inconsistent,inconsistent,gcov
lcov --remove coverage.info '/usr/*' '*/external/*' '/opt/conan/*' '*/test/*' '*/build/*' \
     --output-file coverage.filtered.info --ignore-errors unused,inconsistent
lcov --list coverage.filtered.info
```

## Code layout

```
src/          Production headers and sources (transforms, grids, solvers, shapes)
test/         GoogleTest suites — one file per feature area
external/     Bundled third-party code (kissfft, nayukidct, pwutils, json11)
```

## Adding tests

Use the `add_spida_test()` helper in `test/CMakeLists.txt`:

```cmake
add_spida_test(mynewtest mynewtest.cpp)
```

This links `GTest::gtest_main` and `SPIDA::spida` automatically.
There is **no custom `main()`** — `GTest::gtest_main` provides it.

Test naming convention: `TEST(FEATURE_TEST, CASE_NAME)` — see existing files.

## Agents

| Agent | Purpose |
|-------|---------|
| `spida-test-writer` | Writes/improves GoogleTest coverage (used by nightly CI) |
| `gtest-writer` | Generic GoogleTest scaffolding |
| `cmake-auditor` | CMakeLists.txt review |
| `coverage-reporter` | Parses lcov and reports untested branches |
| `cpp-reviewer` | C++23 / Core Guidelines review |
| `conan-packager` | conanfile.py authoring |

## CI workflows

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| `cmake.yml` | push/PR to main, develop | Build + test on Linux, Windows, macOS |
| `coverage.yml` | push/PR to main | lcov → Codecov upload |
| `nightly-test-agent.yml` | 02:00 UTC daily | Claude Code writes tests, opens PR |

The nightly agent requires `ANTHROPIC_API_KEY` in repository secrets.
