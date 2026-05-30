# Build and Test (Canonical Workflow)

Execute the full deterministic build pipeline exactly as defined in CLAUDE.md.

Do not deviate. Do not build individual targets. Do not guess names or directories.

---

1. Clean environment check (only if needed)

If there is any suspicion of stale state or previous failed builds:

rm -rf build/

Only do this if:
- build directory is inconsistent
- CMake errors reference missing presets or cache mismatch
- user explicitly requests a clean rebuild

---

2. Release Build (PRIMARY)

conan install . -of build/Release --build=missing -s build_type=Release
cmake --preset conan-release -DSPIDA_TEST=ON -DSPIDA_DEMOS=ON
cmake --build --preset conan-release --parallel

---

3. Run Tests (Release)

ctest --preset conan-release --output-on-failure

---

4. RelWithDebInfo + Coverage Build

conan install . -of build/RelWithDebInfo --build=missing -s build_type=RelWithDebInfo
cmake --preset conan-relwithdebinfo -DSPIDA_TEST=ON -DSPIDA_COVERAGE=ON
cmake --build --preset conan-relwithdebinfo --parallel

---

5. Run Tests (Coverage Build)

ctest --preset conan-relwithdebinfo --output-on-failure

---

EXECUTION RULES (STRICT)

- NEVER use cmake -B, cmake --build <dir>, or any non-preset invocation
- NEVER build individual targets
- NEVER guess target names (propagatortest, spidatest are invalid unless explicitly confirmed)
- NEVER change build directory casing (build/Release and build/RelWithDebInfo only)
- ALWAYS use Conan 2 (conan install ...)
- ALWAYS use --build=missing
- ALWAYS use CMake presets (conan-release, conan-relwithdebinfo)
- If ANY step fails: STOP immediately and report the error verbatim
- Do NOT attempt alternative builds, workarounds, or partial compilation

---

EXPECTED OUTPUT VALIDATION

A successful run must include:
- Conan dependency resolution complete
- CMake configure succeeds using preset
- Full target build completes
- ctest runs without build failures

---

PURPOSE

Single canonical build path for:
- reproducibility
- CI parity
- elimination of ad-hoc target guessing
- consistent Conan + CMake preset usage
