Use the cpp-reviewer subagent.

Split the repo into independent modules:

1. src/grid
2. src/helper
3. src/propagator
4. src/pwutils
5. src/rkstiff
6. src/shape
7. src/transform
8. test/
9. demos/

Run one cpp-reviewer per module in parallel.

Each agent should:
- only read its assigned folder
- not modify files
- focus on correctness + C++ safety + build issues

After all finish:
- deduplicate findings
- rank severity
- suggest fixes grouped by module
