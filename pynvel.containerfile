ARG PYTHON_VERSION=3.11
ARG OS_TAG=slim-bookworm

FROM python:$PYTHON_VERSION-$OS_TAG

RUN apt update && apt install -y \
    git \
    gcc \
    g++ \
    gfortran \
    wget \
  && rm -rf /var/lib/apt/lists/* \
  && apt clean

RUN python -m pip install --upgrade pip
RUN python -m pip install --no-cache-dir ninja meson patchelf meson-python setuptools setuptools_scm click numpy pandas cython toml build wheel pytest

ARG PYNVEL_VERSION=main

# Clone, build, install PyNVEL
WORKDIR /
RUN git clone https://github.com/forest-modeling/pynvel.git --recurse-submodules

WORKDIR /pynvel

RUN git checkout "${PYNVEL_VERSION}"
RUN python -m pip install --no-build-isolation .

ENTRYPOINT ["pynvel"]