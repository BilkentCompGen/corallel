#!/usr/bin/env bash
# Removes every tool output and log of the simulation experiments from a pipeline WORK_DIR,
# keeping the simulated data and its metadata:
#   fa-sim/        ref.fa, div_*.fa, sv_*.fa, filenames.txt, shortnames.txt
#   fq-sim/        mut_*.fa
#   fq-sim/fq_*/   ref.fq, mut.fq, filenames.txt, shortnames.txt, min_count.txt, sim-data-out.txt (Badread log)
# It also removes WORK_DIR/sim-results.tsv. Without -f it only lists what it would delete.
#
# Usage: bash cleanup-sim.sh [WORK_DIR] [-f]   (WORK_DIR defaults to the current directory)

work_dir=.
force=false
for arg in "$@"; do
    case $arg in
        -f) force=true ;;
        *)  work_dir=$arg ;;
    esac
done
cd "$work_dir" || exit 1

mapfile -t targets < <(
    if [ -d fa-sim ]; then
        find fa-sim -mindepth 1 -maxdepth 1 \
            ! -name ref.fa ! -name 'div_*.fa' ! -name 'sv_*.fa' ! -name filenames.txt ! -name shortnames.txt
    fi
    if [ -d fq-sim ]; then
        find fq-sim -mindepth 2 -maxdepth 2 -path 'fq-sim/fq_*/*' \
            ! -name ref.fq ! -name mut.fq ! -name filenames.txt ! -name shortnames.txt \
            ! -name min_count.txt ! -name sim-data-out.txt
    fi
    if [ -e sim-results.tsv ]; then
        echo sim-results.tsv
    fi
)

if [ ${#targets[@]} -eq 0 ]; then
    echo "Nothing to clean in $PWD"
    exit 0
fi

if [ "$force" = true ]; then
    rm -rf -- "${targets[@]}"
    echo "Deleted ${#targets[@]} files and folders in $PWD"
else
    printf '%s\n' "${targets[@]}"
    echo "Dry run: the ${#targets[@]} files and folders above would be deleted from $PWD. Add -f to delete them."
fi
