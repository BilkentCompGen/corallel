# Usage: bash corallel-pipeline.sh   (no arguments; runs from any directory)
#   Settings come from corallel-config.sh next to this script, including WORK_DIR, where inputs are read and
#   results written (default: <repo>/genomes; set a full path in the config). Tools come from <repo>/third_party.
#
# Required files:
#   fa-primates
#       ├── hg002v1.1.fasta (+ .fai)
#       ├── mGorGor1.dip.cur.20231122.fasta (+ .fai)
#       ├── mPanPan1.dip.cur.20231122.fasta (+ .fai)
#       ├── mPanTro3.dip.cur.20231122.fasta (+ .fai)
#       ├── mPonAbe1.dip.cur.20231205.fasta (+ .fai)
#       ├── mPonPyg2.dip.cur.20231122.fasta (+ .fai)
#       ├── mSymSyn1.analysis-dip.20240514.fasta (+ .fai)
#       ├── mSymSyn1.analysis-dip.20240514.fasta.fai
#       ├── filenames.txt
#       └── shortnames.txt
#
#   The content of the input can be filled via following command:
#       filenames.txt:
#           find . -maxdepth 1 -name "*.fasta" | sort -V | sed 's|^\./||' > filenames.txt
#   The content of the shortname is as follows:
#       shortnames.txt:
#           Human
#           Gorilla
#           Bonobo
#           Chimpanzee
#           Sum Orang
#           Bor Orang
#           Siamang
#

RUN_CORALLEL_TUMOR="true"
RUN_MASH_TUMOR="true"
RUN_SOURMASH_TUMOR="true"
RUN_DASHING2_TUMOR="true"
RUN_ANI="true"
RUN_CORALLEL_PRIM="true"
RUN_MASH_PRIM="true"
RUN_SOURMASH_PRIM="true"
RUN_DASHING2_PRIM="true"
RUN_CORALLEL_FQ="true"
RUN_CORALLEL_FQ_CC="true"
RUN_SIM_DATA="true"
RUN_CORALLEL_SIM="true"
RUN_MASH_SIM="true"
RUN_SOURMASH_SIM="true"
RUN_DASHING2_SIM="true"
RUN_BINDASH_SIM="true"
RUN_HYPERGEN_SIM="true"
RUN_CORALLEL_SIM_FQ="true"
RUN_MASH_SIM_FQ="true"
RUN_DASHING2_SIM_FQ="true"
RUN_BINDASH_SIM_FQ="true"
RUN_SKMER_SIM_FQ="true"

SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
PROJECT_ROOT_DIR="$(dirname "$SCRIPT_DIR")"

CORALLEL="${PROJECT_ROOT_DIR}/corallel"
PHYLOWIZARD="${PROJECT_ROOT_DIR}/phylowizard.py"
PLOT_ANI="${SCRIPT_DIR}/plot-mut-vs-dist.py"
MASH="mash"
CONVERT="${SCRIPT_DIR}/convert.py"
SOURMASH="sourmash"
DASHING2="dashing2"
CMP_TREE="${SCRIPT_DIR}/cmp_trees.py"
PLOT_DIST="${SCRIPT_DIR}/plot_tree_disruption.py"
BINDASH="bindash"
HYPERGEN="hyper-gen"
SKMER="skmer"
BADREAD="badread"
SIM_PY="${SCRIPT_DIR}/sim.py"

SIM_REF="chm13v2.0.chr20.fa"
SIM_ERROR_MODEL="pacbio2021"  # badread error/qscore model: pacbio2021 (HiFi-like) or random (uniform sub/ins/del)
SIM_IDENTITY_SD=0.05          # spread of per-read identity, in percentage points (0 = every read exactly 1 - error)
SIM_SEEDS="1 2 3"
SIM_DIVERGENCES="0.0001 0.0005 0.001 0.005 0.01 0.02"
SIM_SV_DIVERGENCE=0.001
SIM_INDEL_RATE=0.0001
SIM_SV_COUNTS="0 100 250 500 1000"
SIM_READ_DIVERGENCE=0.001
SIM_COVERAGES="5 10 15 20 30"
SIM_ERROR_RATES="0.001 0.005 0.01 0.02 0.05"
SIM_FIXED_COVERAGE=20
SIM_FIXED_ERROR=0.001
SIM_KMERS="21 31"
SIM_LEVELS_FA="4 5 6 7"
SIM_LEVELS_FQ="4 5 6"
SIM_SKETCH_SIZE=1000000
SIM_THREADS=${SLURM_CPUS_PER_TASK:-64}

PRIMATES_FQ_DIR=.

WORK_DIR="${PROJECT_ROOT_DIR}/genomes"
THIRD_PARTY="${PROJECT_ROOT_DIR}/third_party"

MUTATION_COUNT_PHY="mutation.count.phy"

CONFIG_FILE="$(dirname "${BASH_SOURCE[0]}")/corallel-config.sh"

# activate the tools before reading the config, so the config can use them (e.g. $(which sourmash))
source "${THIRD_PARTY}/activate.sh"

if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    echo "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE"
fi

mkdir -p "$WORK_DIR" && cd "$WORK_DIR" || exit 1
WORK_DIR=$PWD
echo "Working in $WORK_DIR"

SIM_RESULTS="$PWD/sim-results.tsv"
[ -f "${SIM_RESULTS}" ] || printf "sample\ttool\tparam\tdistance\n" > ${SIM_RESULTS}
[ -f "${SIM_REF}" ] && SIM_SCALED=$(( $(grep -v '>' ${SIM_REF} | tr -d '\n' | wc -c) / SIM_SKETCH_SIZE ))

# Corallel Tumor
if [ "$RUN_CORALLEL_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f corallel-tumor-fa-out.txt

    for l in 4 5 6 7 8 9; do

        /bin/time -v ${CORALLEL} fa \
            -i filenames.txt \
            -s shortnames.txt \
            -t 4 \
            -l "$l" \
            -p tumor \
            -v >> corallel-tumor-fa-out.txt 2>&1

        python3 ${PHYLOWIZARD} tumor.set.evol.lvl${l}.phy --normalize >> corallel-tumor-fa-out.txt 2>&1

    done

    cd $WORK_DIR
fi

# Mash Tumor
if [ "$RUN_MASH_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f mash-tumor-fa-out.txt

    for s in 1000 5000 50000 500000 5000000; do

        /bin/time -v ${MASH} sketch \
            -o tumor.mash.${s}.msh *.fa \
            -s ${s} \
            -p 4 >> mash-tumor-fa-out.txt 2>&1

        /bin/time -v ${MASH} dist \
            tumor.mash.${s}.msh tumor.mash.${s}.msh > tumor.mash.${s}.dist

        python3 ${CONVERT} --mash tumor.mash.${s}.dist tumor.mash.${s}.phy

        if [ -f "tumor.mash.${s}.phy" ]; then
            python3 ${PHYLOWIZARD} tumor.mash.${s}.phy --normalize >> mash-tumor-fa-out.txt 2>&1
        else
            echo "File tumor.mash.${s}.phy not found!"
        fi
    done

    for s in 1000 5000 50000 500000 5000000; do
        if [ -f "tumor.mash.$s.nj.newick" ]; then

            newick=$(cat "tumor.mash.$s.nj.newick" | sed -E "s/.fa//g")
            node_names=$(echo "$newick" | grep -oE 'node[0-9]+' | sort -uV)
            declare -A mash_node_to_tumor
            tumor_id=1

            for node in $node_names; do
                mash_node_to_tumor["$node"]="tumor$tumor_id"
                tumor_id=$((tumor_id + 1))
            done

            for node in "${!mash_node_to_tumor[@]}"; do
                tumor=${mash_node_to_tumor[$node]}
                newick=$(echo "$newick" | sed -E "s/\b$node\b/$tumor/g")
            done

            echo "$newick" > "tumor.mash.$s.nj.renamed.newick"
        fi
    done

    cd $WORK_DIR
fi

# Sourmash tumor
if [ "$RUN_SOURMASH_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f sourmash-tumor-fa-out.txt

    /bin/time -v ${SOURMASH} compute -k 21 *.fa > sourmash-tumor-fa-out.txt 2>&1
    /bin/time -v ${SOURMASH} compare -p 8 *.sig -o sourmash.tumor >> sourmash-tumor-fa-out.txt 2>&1

    cd $WORK_DIR
fi

# Dashing2 Tumor
if [ "$RUN_DASHING2_TUMOR" = "true" ]; then

    cd fa-tumor-human

    rm -f dashing2-tumor-fa-out.txt

    for s in 1000 5000 50000 500000 5000000; do

        /bin/time -v ${DASHING2} sketch \
            -F filenames.txt \
            -p 4 \
            -S${s} \
            --cmpout "tumor.dashing2.s${s}.phy" >> dashing2-tumor-fa-out.txt 2>&1

        python3 ${CONVERT} --dashing2 tumor.dashing2.s${s}.phy tumor.dashing2.s${s}.norm.phy

        if [ -f "tumor.dashing2.s${s}.norm.phy" ]; then
            python3 ${PHYLOWIZARD} tumor.dashing2.s${s}.norm.phy --normalize >> dashing2-tumor-fa-out.txt 2>&1
        else
            echo "tumor.dashing2.s${s}.norm.phy not found!"
        fi

    done

    cd $WORK_DIR
fi

# Corallel-ANI tumor
if [ "$RUN_ANI" = "true" ]; then

    mkdir -p mut-vs-ani-plots
    cd mut-vs-ani-plots

    if [ -f "${MUTATION_COUNT_PHY}" ]; then
        for l in 4 5 6 7 8 9; do
            python3 ${PLOT_ANI} ${MUTATION_COUNT_PHY} ../fa-tumor-human/tumor.set.jaccard.lvl${l}.phy --label "Corallel (${l})" --ylabel "Jaccard Distance" -o count-corallel-j-${l}.pdf
            python3 ${PLOT_ANI} ${MUTATION_COUNT_PHY} ../fa-tumor-human/tumor.set.evol.lvl${l}.phy --label "Corallel (${l})" -o count-corallel-p-${l}.pdf
        done
    fi

    cd $WORK_DIR
fi

# Corallel primates
if [ "$RUN_CORALLEL_PRIM" = "true" ]; then

    cd fa-primates

    rm -f corallel-primates-fa-out.txt

    # Corallel
    for l in 4 5 6 7; do
        /bin/time -v ${CORALLEL} fa \
            -i filenames.txt \
            -s shortnames.txt \
            -t 7 \
            -l "$l" \
            -p primates \
            -v >> corallel-primates-fa-out.txt 2>&1;

        python3 ${PHYLOWIZARD} primates.set.evol.lvl${l}.phy;
        python3 ${CMP_TREE} ${PRIMATES_GROUND} primates.set.evol.lvl${l}.nj.newick > primates.set.evol.lvl${l}.cmp.txt
    done

    cd $WORK_DIR
fi

# Mash primates
if [ "$RUN_MASH_PRIM" = "true" ]; then

    cd fa-primates

    rm -f mash-primates-fa-out.txt

    /bin/time -v ${MASH} sketch \
        -o primates.mash.1000.msh *.fasta \
        -s 1000 \
        -p 7 >> mash-primates-fa-out.txt 2>&1

    /bin/time -v mash dist primates.mash.1000.msh primates.mash.1000.msh > primates.mash.1000.dist
    python3 ${CONVERT} --mash primates.mash.1000.dist primates.mash.1000.phy
    python3 ${PHYLOWIZARD} primates.mash.1000.phy

    cd $WORK_DIR
fi

# Mash primates
if [ "$RUN_SOURMASH_PRIM" = "true" ]; then

    cd fa-primates

    rm -f sourmash-primates-fa-out.txt

    /bin/time -v ${SOURMASH} compute -k 21 *.fasta > sourmash-primates-fa-out.txt 2>&1
    /bin/time -v ${SOURMASH} compare -p 8 *.sig -o sourmash.primates >> sourmash-primates-fa-out.txt 2>&1

    cd $WORK_DIR
fi

# Dashing2 primates
if [ "$RUN_DASHING2_PRIM" = "true" ]; then

    cd fa-primates

    rm -f dashing2-primates-fa-out.txt

    /bin/time -v ${DASHING2} sketch \
        -F filenames.txt \
        -p 6 \
        --cmpout "primates.dashing2.default.phy" >> dashing2-primates-fa-out.txt 2>&1

    python3 ${CONVERT} --dashing2 primates.dashing2.default.phy primates.dashing2.default.norm.phy
    python3 ${PHYLOWIZARD} primates.dashing2.default.norm.phy >> dashing2-primates-fa-out.txt 2>&1

    cd $WORK_DIR
fi

# Corallel primates fq
if [ "$RUN_CORALLEL_FQ" = "true" ]; then

    mkdir -p fq-primates

    cd $PRIMATES_FQ_DIR

    rm -f corallel-primates-fq-out.txt

    for l in 4 5 6; do

        /bin/time -v ${CORALLEL} fq \
            -i filenames.txt \
            -s shortnames.txt \
            -l "$l" \
            -t 32 \
            -r 11 \
            -p primates.fq \
            -v >> corallel-primates-fq-out.txt 2>&1

        python3 ${PHYLOWIZARD} primates.fq.set.evol.lvl${l}.phy >> corallel-primates-fq-out.txt 2>&1

    done

    mv *.phy $WORK_DIR/fq-primates
    mv *.newick $WORK_DIR/fq-primates
    mv corallel-primates-fq-out.txt $WORK_DIR/fq-primates


    cd $WORK_DIR
fi

# Corallel primates fq
if [ "$RUN_CORALLEL_FQ_CC" = "true" ]; then

    mkdir -p fq-primates-cc

    cd $PRIMATES_FQ_DIR

    for l in 4 5 6; do
        /bin/time -v ${CORALLEL} fq \
            -i filenames.txt \
            -s shortnames.txt \
            -o binnames${l}.txt \
            -l "$l" \
            -t 32 \
            -r 6 \
            -p primates.fq \
            -v >> corallel-primates-fq-out.txt 2>&1
    done

    for l in 4 5 6; do
        /bin/time -v ${CORALLEL} ld \
            -i binnames${l}.txt \
            -s shortnames.txt \
            -t 6 \
            -r 6 \
            -l "$l" \
            -p primates.fq -v
        
        rm -f distances.lvl${l}.txt

        for cc in {0..200..4}; do
            python3 ${PHYLOWIZARD} primates.fq.cc${cc}.set.evol.lvl${l}.phy >> corallel-primates-fq-out.txt 2>&1
            python3 ${CMP_TREE} ${PRIMATES_GROUND} primates.fq.cc${cc}.set.evol.lvl${l}.nj.newick >> distances.lvl${l}.txt
            echo "" >> distances.lvl${l}.txt
        done

        python3 ${PLOT_DIST} distances.lvl${l}.txt -o tree_disruption.lvl${l}.pdf
    done

    mv *.phy $WORK_DIR/fq-primates-cc
    mv *.newick $WORK_DIR/fq-primates-cc
    mv corallel-primates-fq-out.txt $WORK_DIR/fq-primates-cc
    mv distances.lvl*.txt $WORK_DIR/fq-primates-cc
    mv tree_disruption.lvl*.pdf $WORK_DIR/fq-primates-cc
    

    cd $WORK_DIR
fi

# Simulated assemblies and reads
if [ "$RUN_SIM_DATA" = "true" ]; then

    [ -f "${SIM_REF}" ] || { echo "ERROR: ${SIM_REF} not found in $PWD (set SIM_REF in ${CONFIG_FILE})" >&2; exit 1; }

    mkdir -p fa-sim fq-sim
    ln -sf "$(realpath ${SIM_REF})" fa-sim/ref.fa

    for seed in ${SIM_SEEDS}; do

        for divergence in ${SIM_DIVERGENCES}; do
            python3 ${SIM_PY} fa-sim/ref.fa fa-sim/div_${divergence}_0_${seed}.fa ${divergence} 0 0 ${seed}
        done

        for sv_count in ${SIM_SV_COUNTS}; do
            python3 ${SIM_PY} fa-sim/ref.fa fa-sim/sv_${SIM_SV_DIVERGENCE}_${sv_count}_${seed}.fa ${SIM_SV_DIVERGENCE} ${SIM_INDEL_RATE} ${sv_count} ${seed}
        done
    done

    (cd fa-sim && echo ref.fa && ls div_*.fa sv_*.fa | sort -V) > fa-sim/filenames.txt
    sed 's/\.fa$//' fa-sim/filenames.txt > fa-sim/shortnames.txt

    for seed in ${SIM_SEEDS}; do

        mutant=$PWD/fq-sim/mut_${SIM_READ_DIVERGENCE}_${seed}.fa
        python3 ${SIM_PY} fa-sim/ref.fa ${mutant} ${SIM_READ_DIVERGENCE} 0 0 ${seed}

        (for coverage in ${SIM_COVERAGES}; do echo "${coverage} ${SIM_FIXED_ERROR}"; done
         for error in ${SIM_ERROR_RATES}; do echo "${SIM_FIXED_COVERAGE} ${error}"; done) | sort -u | while read coverage error; do

            mkdir -p fq-sim/fq_${coverage}_${error}_${seed}
            cd fq-sim/fq_${coverage}_${error}_${seed}

            for genome in ref:$(realpath ../../fa-sim/ref.fa):${seed} mut:${mutant}:$((seed + 100)); do
                IFS=: read name fasta read_seed <<< "${genome}"
                # badread is single-threaded: run SIM_THREADS copies, each making 1/SIM_THREADS of the
                # coverage with its own seed, then join them (the reads therefore depend on SIM_THREADS).
                # Adapters, chimeras, junk/random reads and glitches are off, so sequencing error is the only variable.
                pids=()
                for part in $(seq ${SIM_THREADS}); do
                    ${BADREAD} simulate \
                        --reference ${fasta} \
                        --quantity $(awk "BEGIN {print ${coverage} / ${SIM_THREADS}}")x \
                        --length 15000,7000 \
                        --identity $(awk "BEGIN {print 100 * (1 - ${error})}"),100,${SIM_IDENTITY_SD} \
                        --error_model ${SIM_ERROR_MODEL} \
                        --qscore_model ${SIM_ERROR_MODEL} \
                        --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
                        --start_adapter_seq "" --end_adapter_seq "" \
                        --seed $(( read_seed * 1000 + part )) > ${name}.part${part}.fq 2> ${name}.part${part}.log &
                    pids+=($!)
                done
                for i in "${!pids[@]}"; do
                    wait ${pids[i]} || { echo "ERROR: badread failed for ${name} in $PWD:" >&2; tail -n 3 ${name}.part$((i + 1)).log >&2; }
                done
                cat $(seq -f "${name}.part%g.fq" ${SIM_THREADS}) > ${name}.fq
                cat $(seq -f "${name}.part%g.log" ${SIM_THREADS}) >> sim-data-out.txt
                rm -f ${name}.part*
            done

            printf "ref.fq\nmut.fq\n" > filenames.txt
            printf "ref\nmut\n" > shortnames.txt
            echo $(( coverage / 2 > 2 ? coverage / 2 : 2 )) > min_count.txt

            cd $WORK_DIR
        done
    done
fi

# Corallel simulated assemblies
if [ "$RUN_CORALLEL_SIM" = "true" ]; then

    cd fa-sim

    rm -f corallel-sim-fa-out.txt

    for l in ${SIM_LEVELS_FA}; do

        /bin/time -v ${CORALLEL} fa \
            -i filenames.txt \
            -s shortnames.txt \
            -t ${SIM_THREADS} \
            -l "$l" \
            -p sim.fa \
            --vec \
            -v >> corallel-sim-fa-out.txt 2>&1

        awk -v param=lvl$l 'NR == FNR {name[FNR] = $1; next} FNR == 2 {for (i = 3; i <= NF; i++) print name[i - 1] "\tcorallel\t" param "\t" $i}' shortnames.txt sim.fa.vec.evol.lvl${l}.phy >> ${SIM_RESULTS}
    done

    cd $WORK_DIR
fi

# Mash simulated assemblies
if [ "$RUN_MASH_SIM" = "true" ]; then

    cd fa-sim

    rm -f mash-sim-fa-out.txt

    for k in ${SIM_KMERS}; do

        /bin/time -v ${MASH} sketch \
            -o sim.mash.k${k}.msh $(cat filenames.txt) \
            -k ${k} \
            -s ${SIM_SKETCH_SIZE} \
            -p ${SIM_THREADS} >> mash-sim-fa-out.txt 2>&1

        /bin/time -v ${MASH} dist \
            -p ${SIM_THREADS} \
            sim.mash.k${k}.msh sim.mash.k${k}.msh > sim.mash.k${k}.dist 2>> mash-sim-fa-out.txt

        awk -v param=k$k '$1 == "ref.fa" && $2 != "ref.fa" {sub(/\.fa$/, "", $2); print $2 "\tmash\t" param "\t" $3}' sim.mash.k${k}.dist >> ${SIM_RESULTS}
    done

    cd $WORK_DIR
fi

# Sourmash simulated assemblies
if [ "$RUN_SOURMASH_SIM" = "true" ]; then

    cd fa-sim

    rm -f sourmash-sim-fa-out.txt

    for k in ${SIM_KMERS}; do

        for fasta in $(cat filenames.txt); do
            /bin/time -v ${SOURMASH} sketch dna \
                -p k=${k},scaled=${SIM_SCALED} \
                --name ${fasta%.fa} \
                -o sim.sourmash.k${k}.${fasta%.fa}.sig ${fasta} >> sourmash-sim-fa-out.txt 2>&1
        done

        /bin/time -v ${SOURMASH} compare \
            --ani \
            -p ${SIM_THREADS} \
            --csv sim.sourmash.k${k}.csv \
            $(sed "s/\.fa$/.sig/; s/^/sim.sourmash.k${k}./" filenames.txt) >> sourmash-sim-fa-out.txt 2>&1

        # sourmash writes CRLF line ends: strip the \r, or it sticks to the last sample name
        awk -F, -v param=k$k '{sub(/\r$/, "")} NR == 1 {for (i = 1; i <= NF; i++) {name[i] = $i; if ($i == "ref") ref_row = i + 1}} NR == ref_row {for (i = 1; i <= NF; i++) if (name[i] != "ref") print name[i] "\tsourmash\t" param "\t" 1 - $i}' sim.sourmash.k${k}.csv >> ${SIM_RESULTS}
    done

    cd $WORK_DIR
fi

# Dashing2 simulated assemblies
if [ "$RUN_DASHING2_SIM" = "true" ]; then

    cd fa-sim

    rm -f dashing2-sim-fa-out.txt

    for k in ${SIM_KMERS}; do

        /bin/time -v ${DASHING2} sketch \
            -F filenames.txt \
            -p ${SIM_THREADS} \
            -k ${k} \
            -S ${SIM_SKETCH_SIZE} \
            --mash-distance \
            --cmpout sim.dashing2.k${k}.dist >> dashing2-sim-fa-out.txt 2>&1

        awk -v param=k$k '/^#Sources/ {for (i = 2; i <= NF; i++) {name[i] = $i; sub(/\.fa$/, "", name[i])}} $1 == "ref.fa" {for (i = 2; i <= NF; i++) if ($i != "-" && name[i] != "ref") print name[i] "\tdashing2\t" param "\t" $i}' sim.dashing2.k${k}.dist >> ${SIM_RESULTS}
    done

    cd $WORK_DIR
fi

# BinDash simulated assemblies
if [ "$RUN_BINDASH_SIM" = "true" ]; then

    cd fa-sim

    rm -f bindash-sim-fa-out.txt

    for k in ${SIM_KMERS}; do

        /bin/time -v ${BINDASH} sketch \
            --kmerlen=${k} \
            --sketchsize64=$(( SIM_SKETCH_SIZE / 64 )) \
            --nthreads=${SIM_THREADS} \
            --outfname=sim.bindash.k${k}.sketch $(cat filenames.txt) >> bindash-sim-fa-out.txt 2>&1

        /bin/time -v ${BINDASH} dist \
            --nthreads=${SIM_THREADS} \
            --mthres=10 \
            --outfname=sim.bindash.k${k}.dist \
            sim.bindash.k${k}.sketch >> bindash-sim-fa-out.txt 2>&1

        awk -v param=k$k '$1 == "ref.fa" || $2 == "ref.fa" {other = ($1 == "ref.fa") ? $2 : $1; if (other != "ref.fa") {sub(/\.fa$/, "", other); print other "\tbindash\t" param "\t" $3}}' sim.bindash.k${k}.dist >> ${SIM_RESULTS}
    done

    cd $WORK_DIR
fi

# HyperGen simulated assemblies
if [ "$RUN_HYPERGEN_SIM" = "true" ]; then

    cd fa-sim

    rm -f hypergen-sim-fa-out.txt
    mkdir -p hypergen-ref && ln -sf ../ref.fa hypergen-ref/ref.fa
    # HyperGen at its defaults (1-in-1500 sampling, 4096-d vectors): ~1M hashes per 4096-d vector gave ANI = 100 for everything

    for k in ${SIM_KMERS}; do

        /bin/time -v ${HYPERGEN} sketch \
            -p hypergen-ref \
            -o sim.hypergen.ref.k${k}.sketch \
            -t ${SIM_THREADS} \
            -k ${k} >> hypergen-sim-fa-out.txt 2>&1

        /bin/time -v ${HYPERGEN} sketch \
            -p . \
            -o sim.hypergen.k${k}.sketch \
            -t ${SIM_THREADS} \
            -k ${k} >> hypergen-sim-fa-out.txt 2>&1

        /bin/time -v ${HYPERGEN} dist \
            -r sim.hypergen.ref.k${k}.sketch \
            -q sim.hypergen.k${k}.sketch \
            -o sim.hypergen.k${k}.ani \
            -t ${SIM_THREADS} \
            -a 0 >> hypergen-sim-fa-out.txt 2>&1

        awk -v param=k$k '{sub(/.*\//, "", $2); sub(/\.fa$/, "", $2)} $2 != "ref" {print $2 "\thypergen\t" param "\t" 1 - $3 / 100}' sim.hypergen.k${k}.ani >> ${SIM_RESULTS}
    done

    cd $WORK_DIR
fi

# Corallel simulated reads
if [ "$RUN_CORALLEL_SIM_FQ" = "true" ]; then

    for sample_dir in fq-sim/*/; do

        cd ${sample_dir}

        sample=$(basename ${sample_dir})
        rm -f corallel-sim-fq-out.txt

        for l in ${SIM_LEVELS_FQ}; do

            /bin/time -v ${CORALLEL} fq \
                -i filenames.txt \
                -s shortnames.txt \
                -t ${SIM_THREADS} \
                -l "$l" \
                --min-cc $(cat min_count.txt) \
                -p sim.${sample} \
                -v >> corallel-sim-fq-out.txt 2>&1

            awk -v sample=${sample} -v param=lvl$l 'NR == 2 {print sample "\tcorallel\t" param "\t" $3}' sim.${sample}.set.evol.lvl${l}.phy >> ${SIM_RESULTS}
        done

        cd $WORK_DIR
    done
fi

# Mash simulated reads
if [ "$RUN_MASH_SIM_FQ" = "true" ]; then

    for sample_dir in fq-sim/*/; do

        cd ${sample_dir}

        sample=$(basename ${sample_dir})

        /bin/time -v ${MASH} dist \
            -p ${SIM_THREADS} \
            -k 21 \
            -s ${SIM_SKETCH_SIZE} \
            -r \
            -m $(cat min_count.txt) \
            ref.fq mut.fq > sim.mash.dist 2> mash-sim-fq-out.txt

        awk -v sample=${sample} '{print sample "\tmash\tk21\t" $3}' sim.mash.dist >> ${SIM_RESULTS}

        cd $WORK_DIR
    done
fi

# Dashing2 simulated reads
if [ "$RUN_DASHING2_SIM_FQ" = "true" ]; then

    for sample_dir in fq-sim/*/; do

        cd ${sample_dir}

        sample=$(basename ${sample_dir})

        /bin/time -v ${DASHING2} sketch \
            -p ${SIM_THREADS} \
            -k 21 \
            -S ${SIM_SKETCH_SIZE} \
            -m $(cat min_count.txt) \
            --mash-distance \
            --cmpout sim.dashing2.dist \
            ref.fq mut.fq > dashing2-sim-fq-out.txt 2>&1

        awk -v sample=${sample} '/^#Sources/ {for (i = 2; i <= NF; i++) if ($i == "mut.fq") mut_col = i} $1 == "ref.fq" {print sample "\tdashing2\tk21\t" $mut_col}' sim.dashing2.dist >> ${SIM_RESULTS}

        cd $WORK_DIR
    done
fi

# BinDash simulated reads
if [ "$RUN_BINDASH_SIM_FQ" = "true" ]; then

    for sample_dir in fq-sim/*/; do

        cd ${sample_dir}

        sample=$(basename ${sample_dir})

        /bin/time -v ${BINDASH} sketch \
            --kmerlen=21 \
            --sketchsize64=$(( SIM_SKETCH_SIZE / 64 )) \
            --nthreads=${SIM_THREADS} \
            --outfname=sim.bindash.sketch ref.fq mut.fq > bindash-sim-fq-out.txt 2>&1

        /bin/time -v ${BINDASH} dist \
            --nthreads=${SIM_THREADS} \
            --mthres=10 \
            --outfname=sim.bindash.dist \
            sim.bindash.sketch >> bindash-sim-fq-out.txt 2>&1

        awk -v sample=${sample} '($1 == "ref.fq" && $2 == "mut.fq") || ($1 == "mut.fq" && $2 == "ref.fq") {print sample "\tbindash\tk21\t" $3}' sim.bindash.dist >> ${SIM_RESULTS}

        cd $WORK_DIR
    done
fi

# Skmer simulated reads
if [ "$RUN_SKMER_SIM_FQ" = "true" ]; then

    for sample_dir in fq-sim/*/; do

        cd ${sample_dir}

        sample=$(basename ${sample_dir})
        rm -rf skmer-lib skmer-input
        # .fastq names: skmer 1.x only reads *.fastq
        mkdir -p skmer-input && ln -sf ../ref.fq skmer-input/ref.fastq && ln -sf ../mut.fq skmer-input/mut.fastq

        /bin/time -v ${SKMER} reference skmer-input \
            -p ${SIM_THREADS} \
            -s ${SIM_SKETCH_SIZE} \
            -l skmer-lib \
            -o sim.skmer > skmer-sim-fq-out.txt 2>&1

        awk -F'\t' -v sample=${sample} 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "mut") mut_col = i} $1 == "ref" {print sample "\tskmer\tk31\t" $mut_col}' sim.skmer.txt >> ${SIM_RESULTS}

        cd $WORK_DIR
    done
fi

cd $WORK_DIR
