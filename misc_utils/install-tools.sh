#!/usr/bin/env bash
# Local setup for the sketching / ANI benchmark tools. No root, no conda.
# Place this script in misc_utils/. Everything is installed into THIRD_PARTY
# (default: third_party/ in the parent directory of misc_utils/):
#   - bioconda tools: installed with micromamba (one small binary, fast solver) into third_party/env
#   - HyperGen:       cloned from GitHub and built with cargo into third_party/bin
# Tools that are already on your PATH are detected and skipped.
#
# Usage (from any directory):  bash misc_utils/install-tools.sh [THIRD_PARTY]
# Linux x86_64 only (Dashing 2 is Linux only).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
THIRD_PARTY_DIR="${1:-$PROJECT_ROOT_DIR/third_party}"
mkdir -p "$THIRD_PARTY_DIR" && THIRD_PARTY_DIR="$(cd "$THIRD_PARTY_DIR" && pwd)"  # absolute: it is written into activate.sh
MICROMAMBA_DIR="$THIRD_PARTY_DIR/micromamba"
MICROMAMBA_BINARY="$MICROMAMBA_DIR/bin/micromamba"
MICROMAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-64/latest"
CONDA_ENV_DIR="$THIRD_PARTY_DIR/env"
SKMER_ENV_DIR="$THIRD_PARTY_DIR/env-skmer"  # skmer 3.x pins seqtk 1.3 (zlib < 1.3), which python 3.13 in env cannot have
LOCAL_BIN_DIR="$THIRD_PARTY_DIR/bin"
ACTIVATE_FILE="$THIRD_PARTY_DIR/activate.sh"

HYPERGEN_REPO_URL="https://github.com/wh-xu/Hyper-Gen.git"
HYPERGEN_SOURCE_DIR="$THIRD_PARTY_DIR/src/Hyper-Gen"
MIN_SOURMASH_VERSION="4.4.0"
MIN_BADREAD_VERSION="0.4.1"  # first release with the pacbio2021 error model
MIN_DASHING2_VERSION="2.1.20" # current bioconda release
MIN_SKMER_VERSION="3.3.0"    # installed into its own env; in the main env the solver can only pick skmer 1.1.0 (2018)

mkdir -p "$THIRD_PARTY_DIR" "$LOCAL_BIN_DIR"

# Put the local dirs first so tools installed by earlier runs are found.
# bin/ before env/bin: the skmer launcher in bin/ must win over an old skmer in env/bin.
export PATH="$LOCAL_BIN_DIR:$CONDA_ENV_DIR/bin:$PATH"

# "command name:bioconda package" pairs.
BIOCONDA_TOOLS=(
  "mash:mash"
  "bindash:bindash"
  "jellyfish:kmer-jellyfish"
  "seqtk:seqtk"
)

packages_to_install=()
hypergen_needs_build=false
skmer_needs_env=false


# is_new_enough <command> <min version> <version flag>: first x.y[.z] in "<command> <flag>" output must be >= min
is_new_enough() {
  command -v "$1" >/dev/null 2>&1 || return 1
  local installed_version
  installed_version="$("$1" "$3" 2>&1 | grep -oE '[0-9]+(\.[0-9]+)+' | head -n1)"
  [ -n "$installed_version" ] || return 1
  [ "$(printf '%s\n%s\n' "$2" "$installed_version" | sort -V | head -n1)" = "$2" ]
}

ensure_micromamba() {
  if [ ! -x "$MICROMAMBA_BINARY" ]; then
    mkdir -p "$MICROMAMBA_DIR"
    curl -L "$MICROMAMBA_URL" | tar -xj -C "$MICROMAMBA_DIR" bin/micromamba
  fi
}


echo "== Checking tools =="

for tool_entry in "${BIOCONDA_TOOLS[@]}"; do
  command_name="${tool_entry%%:*}"
  package_name="${tool_entry#*:}"
  if command -v "$command_name" >/dev/null 2>&1; then
    echo "[ok]      $command_name"
  else
    echo "[missing] $command_name"
    packages_to_install+=("$package_name")
  fi
done

# "command name:minimum version:flag that prints the version" (command name = bioconda package); older versions get upgraded.
for tool_entry in "sourmash:$MIN_SOURMASH_VERSION:--version" "badread:$MIN_BADREAD_VERSION:--version" "dashing2:$MIN_DASHING2_VERSION:--help"; do
  IFS=: read -r command_name min_version version_flag <<< "$tool_entry"
  if is_new_enough "$command_name" "$min_version" "$version_flag"; then
    echo "[ok]      $command_name"
  else
    echo "[missing] $command_name >= $min_version"
    packages_to_install+=("$command_name>=$min_version")
  fi
done

if is_new_enough skmer "$MIN_SKMER_VERSION" -h; then
  echo "[ok]      skmer"
else
  echo "[missing] skmer >= $MIN_SKMER_VERSION (own env)"
  skmer_needs_env=true
fi

if python3 -c "import numpy" >/dev/null 2>&1; then
  echo "[ok]      python3 + numpy"
else
  echo "[missing] python3 + numpy"
  packages_to_install+=("numpy")
fi

if command -v hyper-gen >/dev/null 2>&1; then
  echo "[ok]      hyper-gen"
else
  echo "[missing] hyper-gen"
  hypergen_needs_build=true
  if ! command -v cargo >/dev/null 2>&1; then
    echo "[missing] cargo (needed to build hyper-gen)"
    packages_to_install+=("rust")
  fi
fi


if [ "${#packages_to_install[@]}" -gt 0 ]; then
  echo
  echo "== Installing with micromamba: ${packages_to_install[*]} =="
  ensure_micromamba

  if [ -d "$CONDA_ENV_DIR/conda-meta" ]; then
    micromamba_subcommand="install"
  else
    micromamba_subcommand="create"
  fi

  "$MICROMAMBA_BINARY" "$micromamba_subcommand" -y -p "$CONDA_ENV_DIR" \
    -c conda-forge -c bioconda "${packages_to_install[@]}"
fi


if [ "$skmer_needs_env" = true ]; then
  echo
  echo "== Installing skmer >= $MIN_SKMER_VERSION into $SKMER_ENV_DIR =="
  ensure_micromamba

  if [ -d "$SKMER_ENV_DIR/conda-meta" ]; then
    micromamba_subcommand="install"
  else
    micromamba_subcommand="create"
  fi

  "$MICROMAMBA_BINARY" "$micromamba_subcommand" -y -p "$SKMER_ENV_DIR" \
    -c conda-forge -c bioconda "skmer>=$MIN_SKMER_VERSION"

  # Launcher in bin/: runs skmer with its own pinned mash, jellyfish and seqtk first on PATH
  printf '#!/bin/bash\nPATH="%s/bin:$PATH" exec "%s/bin/skmer" "$@"\n' "$SKMER_ENV_DIR" "$SKMER_ENV_DIR" > "$LOCAL_BIN_DIR/skmer"
  chmod +x "$LOCAL_BIN_DIR/skmer"
fi


if [ "$hypergen_needs_build" = true ]; then
  echo
  echo "== Building HyperGen =="

  if [ ! -d "$HYPERGEN_SOURCE_DIR" ]; then
    git clone "$HYPERGEN_REPO_URL" "$HYPERGEN_SOURCE_DIR"
  fi

  # --root puts the binary in third_party/bin
  cargo install --path "$HYPERGEN_SOURCE_DIR" --root "$THIRD_PARTY_DIR"
fi


# Script that adds the local tools to PATH (also use it inside sbatch scripts).
cat > "$ACTIVATE_FILE" <<EOF
export PATH="$LOCAL_BIN_DIR:$CONDA_ENV_DIR/bin:\$PATH"
EOF


echo
echo "== Final check =="

all_tools_found=true
for command_name in mash sourmash dashing2 bindash hyper-gen skmer jellyfish seqtk badread python3; do
  if command -v "$command_name" >/dev/null 2>&1; then
    echo "[ok]      $command_name"
  else
    echo "[MISSING] $command_name"
    all_tools_found=false
  fi
done

if ! python3 -c "import numpy" >/dev/null 2>&1; then
  echo "[MISSING] numpy"
  all_tools_found=false
fi

if [ "$all_tools_found" = false ]; then
  exit 1
fi

echo
echo "All tools found. Run: source $ACTIVATE_FILE"