# Build + Test + Coverage (Canonical Coverage Pipeline)

Produces `coverage.filtered.info` — the artifact consumed by the review-test-parallel
command. Run this BEFORE review-test-parallel.

Execute the pipeline exactly as defined in CLAUDE.md. Do not deviate. Do not build
individual targets. Do not guess names or directories.

---

1. Clean environment check (only if needed)

```bash
rm -rf build/
```

Only do this if:
- build directory is inconsistent
- CMake errors reference missing presets or cache mismatch
- a structural CMake change was just made (CLAUDE.md requires a clean reconfigure)
- user explicitly requests a clean rebuild

Note: this removes both build/Release and build/RelWithDebInfo.

---

2. Coverage build (RelWithDebInfo + gcov instrumentation)

```bash
conan install . -of build/RelWithDebInfo --build=missing -s build_type=RelWithDebInfo
cmake --preset conan-relwithdebinfo -DSPIDA_TEST=ON -DSPIDA_COVERAGE=ON
cmake --build --preset conan-relwithdebinfo --parallel
```

---

3. Run tests (populates .gcda counters)

```bash
ctest --preset conan-relwithdebinfo --output-on-failure
```

---

4. Capture + filter coverage

```bash
lcov --capture --directory build/RelWithDebInfo --output-file coverage.info \
     --ignore-errors mismatch,inconsistent,inconsistent,gcov
lcov --remove coverage.info '/usr/*' '*/external/*' '/opt/conan/*' '*/test/*' '*/build/*' \
     --output-file coverage.filtered.info --ignore-errors unused,inconsistent
lcov --list coverage.filtered.info
```

---

EXECUTION RULES (STRICT)
- NEVER use cmake -B, cmake --build <dir>, or any non-preset invocation
- NEVER build individual targets
- NEVER guess target names (propagatortest, spidatest are invalid unless explicitly confirmed
  via `cmake --build --preset conan-relwithdebinfo --target help`)
- NEVER change build directory casing (build/RelWithDebInfo only)
- ALWAYS use Conan 2 (conan install ...)
- ALWAYS use --build=missing
- ALWAYS use the CMake preset conan-relwithdebinfo
- Use the lcov --ignore-errors flags above verbatim — they absorb gcov/lcov version skew
- If ANY step fails: STOP immediately and report the error verbatim
- Do NOT attempt alternative builds, workarounds, or partial compilation

---

EXPECTED OUTPUT VALIDATION
A successful run must include:
- Conan dependency resolution complete
- CMake configure succeeds using the conan-relwithdebinfo preset
- Full target build completes
- ctest runs without build failures
- coverage.filtered.info exists and reports nonzero line coverage in `lcov --list`

---

PURPOSE
Single source of the coverage artifact (coverage.filtered.info) for:
- the review-test-parallel command (which consumes it read-only, no rebuild)
- CI coverage upload parity
- consistent Conan + CMake preset usage with deterministic lcov filtering
