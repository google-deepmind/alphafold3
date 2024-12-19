# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

FROM nvidia/cuda:12.6.0-base-ubuntu22.04 AS hmmer_build

ARG DEBIAN_FRONTEND=noninteractive
ARG HMMER_VERSION=3.4

# gcc, g++, make are required for compiling hmmer
RUN apt-get update && apt-get install --no-install-recommends -y \
        wget \
        gcc \
        g++ \
        make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Build HMMER
RUN mkdir /hmmer_build /hmmer \
    && wget http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz --directory-prefix /hmmer_build \
    && cd /hmmer_build && tar zxf hmmer-*.tar.gz && rm hmmer-*.tar.gz \
    && cd ./hmmer-* \
    && ./configure --prefix /hmmer \
    && make -j8 \
    && make install \
    && cd ./easel && make install \
    && rm -R /hmmer_build


FROM nvidia/cuda:12.6.0-base-ubuntu22.04 AS alphafold_build

ARG DEBIAN_FRONTEND=noninteractive
ARG ALPHAFOLD_VERSION=3.0.0

# Get latest package listing, install software-properties-common, git, wget,
# compilers and libraries.
# git is required for pyproject.toml toolchain's use of CMakeLists.txt.
# gcc, g++, make are required for compiling AlphaFold 3 libaries.
# zlib is a required dependency of AlphaFold 3.
# Use deadsnakes/ppa apt repository to install specific Python versions
RUN apt-get update && apt-get install --no-install-recommends -y \
        software-properties-common \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update && apt-get install --no-install-recommends -y \
        git \
        wget \
        gcc \
        g++ \
        make \
        zlib1g-dev \
        zstd \
        python3.11 \
        python3-pip \
        python3.11-venv \
        python3.11-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create Python 3.11 virtual environment for the AlphaFold build
RUN python3.11 -m venv /alphafold3_venv

# Update pip to the latest version. Not necessary in Docker, but good to do when
# this is used as a recipe for local installation since we rely on new pip
# features for secure installs.
RUN /alphafold3_venv/bin/pip3 install --no-cache-dir --upgrade pip

# Add the AlphaFold 3 source code from GitHub to the image and
# set the working directory to there
ADD https://github.com/google-deepmind/alphafold3.git#v${ALPHAFOLD_VERSION} /app/alphafold

WORKDIR /app/alphafold

# Install the Python dependencies AlphaFold 3 needs and build AlphaFold itself.
RUN /alphafold3_venv/bin/pip3 install --no-cache-dir -r dev-requirements.txt \
    && /alphafold3_venv/bin/pip3 install --no-cache-dir --no-deps --verbose .

    # Build chemical components database (this binary was installed by pip).
RUN /alphafold3_venv/bin/build_data

# Remove build-related content from the AlphaFold directory, ready for copy in the
# "final" stage
RUN [ "/bin/bash", \
    "-c", \
    "rm -rf /app/alphafold/{CMakeLists.txt,dev-requirements.txt,docker,pyproject.toml,requirements.txt,src}" ]


FROM nvidia/cuda:12.6.0-base-ubuntu22.04 AS final

ARG DEBIAN_FRONTEND=noninteractive

COPY --from=hmmer_build /hmmer /hmmer
COPY --from=alphafold_build /alphafold3_venv /alphafold3_venv
COPY --from=alphafold_build /app/alphafold /app/alphafold

# Package installs without dev or build related requirements
RUN apt-get update && apt-get install --no-install-recommends -y \
        software-properties-common \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update && apt-get install --no-install-recommends -y \
        git \
        wget \
        zlib1g \
        zstd \
        python3.11 \
        python3-pip \
        python3.11-venv \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# To work around a known XLA issue causing the compilation time to greatly
# increase, the following environment variable setting XLA flags must be enabled
# when running AlphaFold 3. Note that if using CUDA capability 7 GPUs, it is
# necessary to set the following XLA_FLAGS value instead:
# ENV XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
# (no need to disable gemm in that case as it is not supported for such GPU).
ENV XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
# Memory settings used for folding up to 5,120 tokens on A100 80 GB.
ENV XLA_PYTHON_CLIENT_PREALLOCATE=true
ENV XLA_CLIENT_MEM_FRACTION=0.95

ENV PATH="/hmmer/bin:/alphafold3_venv/bin:${PATH}"

WORKDIR /app/alphafold

CMD ["python3", "run_alphafold.py", "--help"]
