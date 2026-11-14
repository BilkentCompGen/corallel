#!/usr/bin/env bash
# Removes every tool output and log of the simulation experiments of the pipeline, keeping the simulated data
# and its metadata in WORK_DIR/genomes:
#   fa-sim/        ref.fa, div_*.fa, sv_*.fa, filenames.txt, shortnames.txt
#   fq-sim/        ref.fa, mut_*.fa
#   fq-sim/fq_*/   ref.fq, mut.fq, filenames.txt, shortnames.txt, min_count.txt, sim-data-out.txt (Badread log)
# It also removes fa-sim/, fq-sim/ and sim-results.tsv from WORK_DIR/results/latest. Without -f it only lists
# what it would delete.
#
# Usage: bash cleanup-sim.sh [-f]   (runs from any directory)
#   WORK_DIR comes from corallel-config.sh next to this script, as in corallel-pipeline.sh (default: <repo>).

SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
WORK_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="${SCRIPT_DIR}/corallel-config.sh"

# only WORK_DIR is taken from the config
[ -f "$CONFIG_FILE" ] && WORK_DIR=$(source "$CONFIG_FILE" > /dev/null 2>&1; echo "$WORK_DIR")
WORK_DIR=$(realpath -m "${WORK_DIR}")

force=false
[ "$1" = "-f" ] && force=true

cd "${WORK_DIR}/genomes" || exit 1

mapfile -t targets < <(
    if [ -d fa-sim ]; then
        find fa-sim -mindepth 1 -maxdepth 1 ! -name ref.fa ! -name 'div_*.fa' ! -name 'sv_*.fa' ! -name filenames.txt ! -name shortnames.txt
    fi
    if [ -d fq-sim ]; then
        find fq-sim -mindepth 2 -maxdepth 2 -path 'fq-sim/fq_*/*' ! -name ref.fq ! -name mut.fq ! -name filenames.txt ! -name shortnames.txt ! -name min_count.txt ! -name sim-data-out.txt
    fi
    for f in fa-sim fq-sim sim-results.tsv; do
        [ -e "../results/latest/$f" ] && echo "../results/latest/$f"
    done
)

if [ ${#targets[@]} -eq 0 ]; then
    echo "Nothing to clean in $WORK_DIR"
    exit 0
fi

if [ "$force" = true ]; then
    rm -rf -- "${targets[@]}"
    echo "Deleted ${#targets[@]} files and folders in $WORK_DIR"
else
    printf '%s\n' "${targets[@]}"
    echo "Dry run: the ${#targets[@]} files and folders above would be deleted from $WORK_DIR. Add -f to delete them."
fi
