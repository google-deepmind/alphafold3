#!/bin/sh
#

set -e

# conda create -n alphafold3 python=3.12
# conda activate alphafold3

pip install uv

# Install the exact dependency tree using uv and cache the build artifacts.
# --frozen: do not update the lockfile during build.
# --all-groups: install development/test dependencies defined in pyproject.toml.
# --no-editable: install as a static package.
# If using this as a recipe for local installation, we recommend removing the
# --frozen and --no-editable flags.
UV_LINK_MODE=copy uv sync --frozen --all-groups --no-editable

# Build chemical components database (this binary was installed by uv sync).
uv run build_data
