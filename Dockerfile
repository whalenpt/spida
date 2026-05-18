FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------
# Bootstrap: get ca-certificates before switching to HTTPS snapshot
# -----------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates curl && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------
# Pin archive to a snapshot for reproducible builds
# -----------------------------
ARG SNAPSHOT=20260501T000000Z
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
# System dependencies (pinned via snapshot)
# -----------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential=12.12ubuntu2 \
        clang-tools \
        cmake=4.2.3-2ubuntu2 \
        curl \
        gcc-15 \
        g++-15 \
        gdb \
        git \
        lcov \
        less \
        ninja-build \
        python3=3.14.3-0ubuntu2 \
        python3-pip \
        rsync \
        tree \
        unzip \
        vim \
        zip && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-15 100 \
        --slave /usr/bin/g++ g++ /usr/bin/g++-15 \
        --slave /usr/bin/gcov gcov /usr/bin/gcov-15 && \
    rm -rf /var/lib/apt/lists/*

# Sanity-check the toolchain version
RUN gcc --version | grep -q ' 15\.' && g++ --version | grep -q ' 15\.'

RUN python3 -m pip install --break-system-packages conan==2.28.1

# -----------------------------
# Node.js + Claude Code CLI
# -----------------------------
RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash - && \
    apt-get install -y --no-install-recommends nodejs && \
    npm install -g @anthropic-ai/claude-code && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------
# Conan home at a neutral, user-independent path
# -----------------------------
ENV CONAN_HOME=/opt/conan
RUN mkdir -p /opt/conan && chown ubuntu:ubuntu /opt/conan

# -----------------------------
# Workspace
# -----------------------------
RUN mkdir -p /workspace && chown ubuntu:ubuntu /workspace
WORKDIR /workspace
USER ubuntu

# -----------------------------
# Conan setup (CONAN_HOME=/opt/conan from ENV above)
# -----------------------------
RUN conan --version && \
    (conan remote add conancenter https://center2.conan.io || true) && \
    conan profile detect --force

# -----------------------------
# Pre-cache dependencies
# -----------------------------
WORKDIR /tmp/spida-deps
COPY --chown=ubuntu:ubuntu conanfile.py .
RUN conan install . --build=missing -s build_type=Release && \
    conan install . --build=missing -s build_type=RelWithDebInfo && \
    rm -rf /tmp/spida-deps

# -----------------------------
# Git safety (system-wide, applies to all users including root in CI)
# -----------------------------
USER root
RUN git config --system --add safe.directory '*'
USER ubuntu
WORKDIR /workspace
CMD ["/bin/bash"]
