# SPIDA — Claude Code project guide

Spectral Integration and Differentiation Algorithms (C++23 library).

## Git Conventions
- Default branch is `develop`, not `main`. Always use `develop` in CI configs, workflow files, and git commands unless explicitly told otherwise.

## Code Layout

```
src/          Production headers and sources (transforms, grids, solvers, shapes)
test/         GoogleTest suites — one file per feature area
external/     Bundled third-party code (kissfft, nayukidct, pwutils) — git submodules
demos/        Demo programs
```

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

Available CMake presets: `conan-release`, `conan-relwithdebinfo`. Always verify preset names against `CMakePresets.json` before using them in CI or build commands.

The build is **static-only**: `BUILD_SHARED_LIBS` is forced `OFF`. Do not introduce shared-library patterns.

## Run Tests

```bash
ctest --preset conan-release --output-on-failure
ctest --preset conan-relwithdebinfo --output-on-failure
```

## Generate Coverage Report

```bash
lcov --capture --directory build/RelWithDebInfo --output-file coverage.info \
     --ignore-errors mismatch,inconsistent,inconsistent,gcov
lcov --remove coverage.info '/usr/*' '*/external/*' '/opt/conan/*' '*/test/*' '*/build/*' \
     --output-file coverage.filtered.info --ignore-errors unused,inconsistent
lcov --list coverage.filtered.info
```

## CMake & Build Conventions
- This is a C++23 project; verify C++23 enforcement when modifying CMakeLists.txt.
- When fixing CMake install/export errors, prefer adjusting target visibility (PUBLIC/PRIVATE/INTERFACE) carefully — INTERFACE libraries and PUBLIC linkage frequently break the install export set. Verify the export set still resolves before declaring done.
- After CMake changes, clear the build cache (`rm -rf build/`) before reconfiguring to avoid stale cache issues.
- Do NOT add target aliases that may collide with submodule-provided targets; check submodules first.
- **Before any structural CMake change** (target type, linkage visibility, name, merge/split): map every target that depends on it transitively, every test/demo that includes its headers, and every `install(TARGETS ...)`/`export(TARGETS ...)` rule that references it. Propose the change only after this discovery pass, noting explicitly which dependents need updates.

## Submodules
- `kissfft` and other libraries under `external/` are git submodules. Clone requires `--recurse-submodules`.
- Conan also provides `kissfft` independently — this dual-provider setup is the main source of target alias collision risk. Always check which provider a consumer is actually using before touching target names.

## Conan / Dependencies
- Required packages in `conanfile.py` — do not remove any of these: `spdlog`, `gtest`, `kissfft`, `boost`.
- Do not use overly broad git config like `safe.directory '*'`; scope it to the specific project path.

## Adding Tests

Use the `add_spida_test()` helper in `test/CMakeLists.txt`:

```cmake
add_spida_test(mynewtest mynewtest.cpp)
```

This links `GTest::gtest_main` and `SPIDA::spida` automatically. There is **no custom `main()`** — `GTest::gtest_main` provides it.

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

## CI Workflows

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| `cmake.yml` | push/PR to `main`, `develop` | Build + test on Linux, Windows, macOS |
| `coverage.yml` | push/PR to `main` | lcov → Codecov upload |
| `nightly-test-agent.yml` | 02:00 UTC daily | Claude Code writes tests, opens PR |

The nightly agent requires `ANTHROPIC_API_KEY` in repository secrets.

## Verification Before Declaring Done
- After any CMake/build change, run a clean configure + build + test cycle locally before declaring the task complete.
- When fixing CI failures, double-check branch names, package names, and preset names against the actual repo state.
