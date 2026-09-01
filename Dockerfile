# syntax=docker/dockerfile:1

# -----------------------------
# Bumping this image
#
# FROM digest, SNAPSHOT, and CLAUDE_CODE_VERSION are mutually constrained.
# Move them together, in this order:
#
#   1. docker pull ubuntu:26.04
#   2. docker image inspect --format '{{.Created}}' ubuntu:26.04
#      -> set FROM to the pulled digest; set SNAPSHOT to the next midnight
#         after Created. SNAPSHOT must never predate the base image.
#   3. docker run --rm <new digest> sh -c 'apt-get update && \
#        apt-get install -y --no-install-recommends ca-certificates curl gnupg && \
#        install -d -m 0755 /etc/apt/keyrings && \
#        curl -fsSL https://downloads.claude.ai/keys/claude-code.asc \
#          -o /etc/apt/keyrings/claude-code.asc && \
#        echo "deb [signed-by=/etc/apt/keyrings/claude-code.asc] \
#          https://downloads.claude.ai/claude-code/apt/stable stable main" \
#          > /etc/apt/sources.list.d/claude-code.list && \
#        apt-get update && apt-cache madison claude-code'
#      -> new CLAUDE_CODE_VERSION
#   4. Rebuild. A mismatched SNAPSHOT trips the assertion layer; a wall of
#      unmet dependencies means SNAPSHOT predates the base image.
# -----------------------------

# Pinned by digest: a floating tag moves forward while SNAPSHOT below does not,
# and the base image's packages must not be newer than the snapshot or apt is
# asked to downgrade and refuses.
FROM ubuntu:26.04@sha256:2260313b31c8c011cd2eebe728008efac1b3982be73eb71348ea2648d2c0e09b

# ARG, not ENV: noninteractive would otherwise leak into every shell in the
# running container and break anything that legitimately wants a prompt.
ARG DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8

# -----------------------------
# Bootstrap: ca-certificates before switching to the HTTPS snapshot
#
# apt reaches archive.ubuntu.com over plain HTTP, but once the Snapshot field
# takes effect it redirects to https://snapshot.ubuntu.com/. With no CA store
# that fetch can't succeed, and apt quietly stays on the live archive.
# This one package therefore comes from the live archive; it is re-pinned
# from the snapshot below.
# -----------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------
# Pin archive to a snapshot for reproducible builds
# -----------------------------
ARG SNAPSHOT=20260818T000000Z
RUN rm -f /etc/apt/sources.list.d/ubuntu.sources && \
    cat >/etc/apt/sources.list.d/ubuntu.sources <<EOF
Types: deb
URIs: http://archive.ubuntu.com/ubuntu
Suites: resolute resolute-updates resolute-security
Components: main universe
Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg
Snapshot: ${SNAPSHOT}
EOF

# -----------------------------
# Assert the snapshot took
#
# apt silently ignores source options it doesn't understand, so a bad field
# or a failed HTTPS fetch degrades to the live archive with no error.
# List filenames encode the host apt actually fetched from.
# -----------------------------
RUN set -eux; \
    apt-get update > /tmp/apt-update.log 2>&1 || { cat /tmp/apt-update.log; exit 1; }; \
    cat /tmp/apt-update.log; \
    ls /var/lib/apt/lists/; \
    if ! ls /var/lib/apt/lists/ | grep -q "^snapshot\.ubuntu\.com_ubuntu_${SNAPSHOT}_"; then \
        echo "FATAL: snapshot ${SNAPSHOT} not in effect; apt is using the live archive"; \
        exit 1; \
    fi; \
    rm -f /tmp/apt-update.log

# -----------------------------
# System dependencies
#
# No per-package version pins: SNAPSHOT is the single source of truth, and
# pins that must independently agree with it only add ways to break.
# ca-certificates is listed plainly rather than with --reinstall, which fails
# outright when the snapshot's version differs from the installed one.
# -----------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        gdb \
        gh \
        git \
        gnupg \
        lcov \
        less \
        ninja-build \
        openssh-client \
        pkg-config \
        python3 \
        python3-pip \
        ripgrep \
        tree \
        valgrind && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------
# Node.js
#
# Not from the Ubuntu archive: its `nodejs` package omits npm entirely, and
# its separate `npm` package drags in node-gyp -> libnode-dev -> libssl-dev,
# which is unsatisfiable under the SNAPSHOT pin. NodeSource's nodejs bundles
# npm and corepack and has no such dependency graph.
#
# Added via manual key+source-list steps rather than their curl-|-bash setup
# script — same reasoning and pattern as the Claude Code repo below.
#
# The Pin-Priority is load-bearing: both repos sit at 500 by default and the
# archive's nodejs would otherwise outrank NodeSource's. 1001 rather than
# 990 so the pin also permits a downgrade if a rebuild ever sees the archive
# version installed first.
# -----------------------------
ARG NODE_MAJOR=22
RUN set -eux; \
    install -d -m 0755 /etc/apt/keyrings; \
    curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key \
        | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg; \
    echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_${NODE_MAJOR}.x nodistro main" \
        > /etc/apt/sources.list.d/nodesource.list; \
    printf 'Package: nodejs\nPin: origin deb.nodesource.com\nPin-Priority: 1001\n' \
        > /etc/apt/preferences.d/nodesource; \
    apt-get update; \
    apt-get install -y --no-install-recommends nodejs; \
    rm -rf /var/lib/apt/lists/*

# -----------------------------
# Claude Code
#
# Installed from Anthropic's signed APT repository, pinned: that repo has no
# snapshot service, and the stable channel moves roughly weekly.
# -----------------------------
ARG CLAUDE_CODE_VERSION=2.1.211-1
RUN set -eux; \
    install -d -m 0755 /etc/apt/keyrings; \
    curl -fsSL https://downloads.claude.ai/keys/claude-code.asc \
        -o /etc/apt/keyrings/claude-code.asc; \
    gpg --show-keys --with-colons /etc/apt/keyrings/claude-code.asc \
        | grep -q '31DDDE24DDFAB679F42D7BD2BAA929FF1A7ECACE'; \
    echo "deb [signed-by=/etc/apt/keyrings/claude-code.asc] https://downloads.claude.ai/claude-code/apt/stable stable main" \
        > /etc/apt/sources.list.d/claude-code.list; \
    apt-get update; \
    apt-get install -y --no-install-recommends "claude-code=${CLAUDE_CODE_VERSION}"; \
    rm -rf /var/lib/apt/lists/*

# -----------------------------
# Disable Claude's auto-updater
#
# The package manager controls the installed version.
# -----------------------------
ENV DISABLE_AUTOUPDATER=1

# -----------------------------
# Conan
# -----------------------------
RUN python3 -m pip install --break-system-packages --no-cache-dir conan==2.28.1

# -----------------------------
# Match the container user to the host UID/GID
#
# Files written into the bind-mounted workspace land as this UID.
# Build with --build-arg UID=$(id -u) --build-arg GID=$(id -g) if yours differ.
# Placed here, after the apt layers, so changing it doesn't invalidate them —
# it only needs to precede the first chown ubuntu:ubuntu below.
# -----------------------------
ARG UID=1000
ARG GID=1000
RUN if [ "$UID" != "1000" ] || [ "$GID" != "1000" ]; then \
        groupmod -g "$GID" ubuntu && usermod -u "$UID" -g "$GID" ubuntu; \
    fi

# -----------------------------
# Conan home at a neutral, user-independent path
# -----------------------------
ENV CONAN_HOME=/opt/conan
RUN mkdir -p /opt/conan && chown ubuntu:ubuntu /opt/conan

# -----------------------------
# Claude configuration
#
# Claude normally keeps ~/.claude and ~/.claude.json.
# Point CLAUDE_CONFIG_DIR at ~/.claude so the config can be persisted together
# in a mounted volume. The directory must exist and be owned by ubuntu, or a
# named volume mounted there is created root-owned and Claude can't write.
# -----------------------------
ENV CLAUDE_CONFIG_DIR=/home/ubuntu/.claude
RUN install -d -o ubuntu -g ubuntu -m 0700 /home/ubuntu/.claude

# -----------------------------
# Workspace
# -----------------------------
RUN mkdir -p /workspace && chown ubuntu:ubuntu /workspace
WORKDIR /workspace

# -----------------------------
# Git safety
# -----------------------------
RUN git config --system --add safe.directory '*'

# -----------------------------
# Conan setup
# -----------------------------
USER ubuntu
RUN conan --version && \
    (conan remote add conancenter https://center2.conan.io || true) && \
    conan profile detect --force

# -----------------------------
# Pre-cache dependencies
#
# conan cache clean strips sources and build trees, which are usually
# gigabytes on a C++ dependency graph.
# -----------------------------
WORKDIR /tmp/spida-deps
COPY --chown=ubuntu:ubuntu conanfile.py .
COPY --chown=ubuntu:ubuntu conan.lock .
RUN conan install . \
        --build=missing \
        -s build_type=Release \
        --lockfile=conan.lock && \
    conan install . \
        --build=missing \
        -s build_type=RelWithDebInfo \
        --lockfile=conan.lock && \
    conan cache clean "*" --source --build --temp && \
    cd / && \
    rm -rf /tmp/spida-deps

# -----------------------------
# Smoke test
#
# Fail at build time rather than mid-session if a launcher is broken.
# The Node major is asserted because a silently-losing apt pin is exactly
# how the archive's npm-less nodejs would get in.
# -----------------------------
RUN set -eux; \
    claude --version; \
    conan --version; \
    cmake --version; \
    gdb --version; \
    valgrind --version; \
    node --version; \
    npm --version; \
    test "$(node -p 'process.versions.node.split(".")[0]')" = "${NODE_MAJOR}"

# -----------------------------
# Final workspace
# -----------------------------
WORKDIR /workspace
CMD ["sleep", "infinity"]
