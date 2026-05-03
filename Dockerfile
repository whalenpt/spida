FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------
# System dependencies
# -----------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential=12.12ubuntu2 \
        cmake=4.2.3-2ubuntu2 \
        git=1:2.53.0-1ubuntu1 \
        python3=3.14.3-0ubuntu2 \
        python3-pip=25.1.1+dfsg-1ubuntu2 \
        ca-certificates=20260223 && \
    rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --break-system-packages conan==2.28.1

RUN mkdir -p /workspace
WORKDIR /workspace
USER ubuntu

# -----------------------------
# Conan setup (runs as dev — profile and cache land in /home/dev/.conan2)
# -----------------------------
RUN conan --version && \
    conan remote add conancenter https://center2.conan.io || true && \
    conan profile detect --force

# -----------------------------
# Pre-cache dependencies
# -----------------------------
WORKDIR /tmp/spida-deps
COPY --chown=ubuntu:ubuntu conanfile.py .
RUN conan install . --build=missing -s build_type=Release && \
    rm -rf /tmp/spida-deps

# -----------------------------
# Git safety
# -----------------------------
RUN git config --global --add safe.directory /workspace

WORKDIR /workspace
CMD ["/bin/bash"]
