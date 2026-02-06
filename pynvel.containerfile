# ---- Stage 1: Build ----
ARG PYTHON_VERSION=3.11
ARG OS_TAG=slim-bookworm

FROM python:$PYTHON_VERSION-$OS_TAG AS builder

WORKDIR /pynvel

RUN apt update && apt install -y --no-install-recommends \
    #build-essential \
    git \
    gcc \
    g++ \
    gfortran \
    wget \
  && rm -rf /var/lib/apt/lists/* \
  && apt clean

# Set environment variables for the virtual environment location
ENV VIRTUAL_ENV=/pynvel/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN python -m venv /pynvel/venv
RUN python -m pip install --no-cache-dir --upgrade pip && \
  python -m pip install --no-cache-dir ninja meson patchelf meson-python setuptools setuptools_scm click numpy pandas cython toml build wheel pytest

ARG PYNVEL_VERSION=main

# Clone, build, install PyNVEL
RUN git clone https://github.com/forest-modeling/pynvel.git  /pynvel_git --recurse-submodules

WORKDIR /pynvel_git

RUN git checkout "${PYNVEL_VERSION}"
RUN python -m pip install --no-build-isolation .

# # Cleanup, remove build tools and source
# RUN apt remove -y git gcc g++ wget \
#   && rm -rf /var/lib/apt/lists/* \
#   && apt clean

RUN python -m pip uninstall -y ninja meson patchelf meson-python setuptools setuptools_scm cython build wheel pytest

RUN rm -r /pynvel_git


# ---- Stage 2: Production ----
FROM python:$PYTHON_VERSION-$OS_TAG

RUN apt update && apt install -y --no-install-recommends libgfortran5 \
  && rm -rf /var/lib/apt/lists/* \
  && apt clean

# Copy the cleaned virtual environment from the build stage
COPY --from=builder /pynvel/venv /pynvel/venv

# Set environment variables for the virtual environment location
ENV VIRTUAL_ENV=/pynvel/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN python -m pip install --no-cache-dir click numpy pandas toml

ENTRYPOINT ["pynvel"]

