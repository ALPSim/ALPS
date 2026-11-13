# CI build environments for the platform matrix in ci.yml (job "platforms"), one image per
# distribution: ghcr.io/alpsim/alps-ci:<name>, built by ci-images.yml. Each image holds the
# compilers of its matrix entries, ALPS's native dependencies and a Python environment with
# the build tools, so jobs start building at once. Boost is built from source (pinned in
# .github/dependencies.json) where the distribution's is older than ALPS's minimum.
ARG BASE
FROM ${BASE}
ARG PACKAGES=""
ARG BOOST=""
LABEL org.opencontainers.image.source=https://github.com/ALPSim/ALPS
LABEL org.opencontainers.image.description="ALPS CI build environment"
COPY .github/ci-images/install.sh /tmp/ci/.github/ci-images/install.sh
COPY .github/scripts/prepare_dependencies.py .github/scripts/build-boost.sh /tmp/ci/.github/scripts/
COPY .github/dependencies.json /tmp/ci/.github/dependencies.json
RUN PACKAGES="${PACKAGES}" BOOST="${BOOST}" bash /tmp/ci/.github/ci-images/install.sh && rm -rf /tmp/ci
ENV PATH=/opt/ci-venv/bin:${PATH} \
    CMAKE_PREFIX_PATH=/opt/boost \
    CMAKE_GENERATOR=Ninja
