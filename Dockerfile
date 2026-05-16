FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------
# Bootstrap: get ca-certificates before switching to HTTPS snapshot
# -----------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates && \
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
        cmake=4.2.3-2ubuntu2 \
        git \
        g++ \
        lcov \
        ninja-build \
        python3=3.14.3-0ubuntu2 \
        python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --break-system-packages conan==2.28.1

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
