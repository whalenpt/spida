# CI Fix Workflow
1. Read the failing job log via gh run view or the provided URL
2. Identify failure category: missing preset, package, branch name, or build error
3. Default branch is `develop` — verify before editing workflow files
4. Apply minimal fix, then run `cmake --preset <name> && cmake --build --preset <name>` locally if possible
5. Summarize the root cause in 1-2 sentences
