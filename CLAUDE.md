## Build

Dependencies are managed with Conan 2. The build works in any environment with Conan 2
and a C++23 toolchain — no container or prebuilt image is required.

### MANDATORY build sequence — follow exactly, do not deviate

**Release:**
```bash
conan install . -of build/Release --build=missing -s build_type=Release
cmake --preset conan-release -DSPIDA_TEST=ON -DSPIDA_DEMOS=ON
cmake --build --preset conan-release --parallel
```

**Coverage (RelWithDebInfo + gcov instrumentation):**
```bash
conan install . -of build/RelWithDebInfo --build=missing -s build_type=RelWithDebInfo
cmake --preset conan-relwithdebinfo -DSPIDA_TEST=ON -DSPIDA_COVERAGE=ON
cmake --build --preset conan-relwithdebinfo --parallel
```

### Build rules (ENFORCED — no exceptions without explicit user instruction)
- ALWAYS use `--build=missing` on `conan install`. Do not assume any package is pre-cached.
- ALWAYS configure and build via `cmake --preset`. NEVER run raw `cmake -B <dir>` or `cmake --build <dir>`.
- NEVER build, configure, or test an individual target unless the user explicitly asks for it. Build the whole preset.
- NEVER guess or infer a CMake target name. Valid target names come only from `cmake --build --preset <preset> --target help` or from explicit declarations in `CMakeLists.txt`. If you don't have a confirmed name, do not invent one.
- NEVER mix directory casing: it is `build/Release` and `build/RelWithDebInfo` exactly. `build/release` is wrong on Linux.
- The only valid presets are `conan-release` and `conan-relwithdebinfo`. Verify against `CMakePresets.json` before using any preset name.
- If ANY step fails, STOP and report that exact failure. Do NOT try alternative targets, directories, build types, or ad-hoc `cmake`/`make` invocations to work around it.
- The build is **static-only**: `BUILD_SHARED_LIBS` is forced `OFF`. Do not introduce shared-library patterns.
- After any CMake change, run `rm -rf build/` before reconfiguring to avoid stale cache issues.

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
- `cmake_minimum_required` appears **only in the top-level `CMakeLists.txt`** (currently `VERSION 3.21`). Subdirectory `CMakeLists.txt` files must NOT declare their own `cmake_minimum_required` — the top-level declaration governs the entire build.
- This is a C++23 project; verify C++23 enforcement when modifying CMakeLists.txt.
- When fixing CMake install/export errors, prefer adjusting target visibility (PUBLIC/PRIVATE/INTERFACE) carefully — INTERFACE libraries and PUBLIC linkage frequently break the install export set. Verify the export set still resolves before declaring done.
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
