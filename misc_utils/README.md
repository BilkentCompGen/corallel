# misc_utils

`corallel-pipeline.sh` reads settings from `corallel-config.sh` next to it. Without a config,
`WORK_DIR` is the repository root.

## Layout

```
WORK_DIR/
├── genomes/                    # inputs
│   ├── fa-tumor-human/         # tumor assemblies       (you provide)
│   │   ├── *.fa (+ .fai)
│   │   ├── filenames.txt
│   │   └── shortnames.txt
│   ├── fa-primates/            # primate assemblies     (you provide)
│   │   ├── *.fasta (+ .fai)
│   │   ├── filenames.txt
│   │   └── shortnames.txt
│   ├── fq-primates/            # primate reads          (you provide, or set PRIMATES_FQ_DIR)
│   │   ├── reads
│   │   ├── filenames.txt
│   │   └── shortnames.txt
│   ├── fa-sim/                 # simulated assemblies   (created by the pipeline)
│   └── fq-sim/                 # simulated reads        (created by the pipeline)
└── results/
    └── latest/                 # all outputs, one folder per input folder above
        ├── fa-tumor-human/
        ├── fa-primates/
        ├── fq-primates/
        ├── fq-primates-cc/
        ├── mut-vs-ani-plots/
        ├── fa-sim/
        ├── fq-sim/
        └── sim-results.tsv
```

## Other files

Set in the config. Give full paths; a plain file name is looked up in `genomes/`.

| Variable             | File                                       |
|----------------------|--------------------------------------------|
| `SIM_REF_FA`         | reference for simulated assemblies         |
| `SIM_REF_FQ`         | reference for simulated reads              |
| `MUTATION_COUNT_PHY` | tumor mutation counts                      |
| `PRIMATES_GROUND`    | primate ground-truth tree (newick)         |
