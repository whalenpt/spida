Goal: review test coverage and test quality across the repo using parallel subagents.

## Step 0 — Coverage artifact (consume only, never build)

This command does not build. If `coverage.filtered.info` exists, use it.
If it does not exist, STOP and tell the user to run the `build-test-coverage` command
first, then re-run this one.

Parallel agents must never build, run ctest, or touch build/RelWithDebInfo — they would
collide on the build directory and the static-only preset.

## Step 1 — Fan out (one subagent per module, in parallel)

Use the spida-test-writer subagent in read-only mode (analysis only — no file changes,
no building, no ctest).

Each agent reviews one src/ folder together with its matching tests under test/:
1. src/grid
2. src/helper
3. src/propagator
4. src/pwutils
5. src/rkstiff
6. src/shape
7. src/transform

Each agent should:
- only read its assigned src/ folder, the corresponding test file(s) in test/, and the
  relevant entries in coverage.filtered.info
- not modify any files
- map the module's public surface (functions, classes, templates, error/std::expected
  paths, exceptions) to existing test coverage, using the coverage data for line/branch
  gaps and static reading for what coverage can't show
- identify coverage gaps: untested functions/branches, boundary and empty inputs,
  error/exception paths, numeric edge cases (NaN/inf, zero-size, large N)
- review existing tests for improvements: loose assertions (prefer EXPECT_THAT +
  matchers), missing fixtures, copy-pasted cases that should be TEST_P, missing
  EXPECT_THROW / std::expected error checks, nondeterminism or flakiness, over-mocking,
  reliance on real time/filesystem
- propose concrete new tests as stubs, named TEST(FEATURE_TEST, CASE_NAME) and wired via
  add_spida_test() with GTest::gtest_main (no custom main)

## Step 2 — Synthesize (after all agents finish)

- deduplicate findings across modules
- rank by risk: untested numeric/critical paths and error handling first
- group output by module:
  (a) coverage gaps
  (b) existing-test improvements
  (c) proposed new test stubs
- flag shared/helper code (e.g. src/helper, src/pwutils) exercised indirectly so its
  tests aren't proposed redundantly across modules
- end with a short prioritized list of the highest-value tests to add next
