#!/usr/bin/env python3
"""Print the simulation results of corallel-pipeline.sh (sim-results.tsv) as LaTeX tables.

Each cell is the mean over seeds 1-3 of estimate / truth:
  assembly divergence   truth = the simulated divergence
  assembly SVs          relative to the same tool at 0 SVs
  reads (coverage, error)  truth = 0.001 (SIM_READ_DIVERGENCE)
A cell with no usable value is "--"; a mean over fewer than 3 seeds is marked "*".
Missing rows and unusable values (nan, inf, empty) are listed on stderr.

Usage: python3 sim-tables.py sim-results.tsv > tables.tex
"""
import math
import sys

SEEDS = ['1', '2', '3']
FA_TOOLS = [('corallel', f'lvl{l}') for l in (4, 5, 6, 7)] + \
           [(t, k) for t in ('mash', 'dashing2', 'bindash', 'sourmash', 'hypergen') for k in ('k21', 'k31')]
FQ_TOOLS = [('corallel', f'lvl{l}') for l in (4, 5, 6)] + \
           [('mash', 'k21'), ('dashing2', 'k21'), ('bindash', 'k21'), ('skmer', 'k31')]
DIVERGENCES = ['0.0001', '0.0005', '0.001', '0.005', '0.01', '0.02', '0.05']
SV_COUNTS = ['100', '1000', '10000']   # each relative to 0 SVs
COVERAGES = ['5', '10', '15', '20', '30']
ERRORS = ['0.001', '0.005', '0.01', '0.02', '0.05']
READ_TRUTH = 0.001

results = {}
with open(sys.argv[1], newline='\n') as f:   # split on \n only: a stray \r must not break a line
    next(f)
    for line in f:
        fields = [x.strip() for x in line.split('\t')]
        if len(fields) == 4:
            results[tuple(fields[:3])] = fields[3]   # a later row (rerun) replaces an earlier one

problems = {}


def mean(samples, tool, param):
    values = []
    for sample in samples:
        raw = results.get((sample, tool, param))
        try:
            value = float(raw)
        except (TypeError, ValueError):
            value = math.nan
        if math.isfinite(value):
            values.append(value)
        else:
            problems[(sample, tool, param)] = 'no row' if raw is None else repr(raw)
    return (sum(values) / len(values), len(values)) if values else (None, 0)


def table(title, header, columns, tools, samples_of, truth_of):
    print(f'% {title}')
    print('\\begin{tabular}{l' + 'r' * len(columns) + '}')
    print(' & '.join(['Tool'] + header) + ' \\\\ \\hline')
    for tool, param in tools:
        cells = []
        for column in columns:
            value, n = mean(samples_of(column), tool, param)
            truth = truth_of(column, tool, param)
            if value is None or not truth:
                cells.append('--')
            else:
                cells.append(f'{value / truth:.2f}' + ('' if n == len(SEEDS) else '*'))
        print(' & '.join([f'{tool} {param}'] + cells) + ' \\\\')
    print('\\end{tabular}\n')


percent = lambda x: f'{float(x) * 100:g}\\%'

table('Assemblies, divergence: estimate / true divergence',
      [percent(d) for d in DIVERGENCES], DIVERGENCES, FA_TOOLS,
      lambda d: [f'div_{d}_0_{s}' for s in SEEDS],
      lambda d, tool, param: float(d))
table('Assemblies, SVs: estimate / same tool at 0 SVs',
      SV_COUNTS, SV_COUNTS, FA_TOOLS,
      lambda n: [f'sv_0.001_{n}_{s}' for s in SEEDS],
      lambda n, tool, param: mean([f'sv_0.001_0_{s}' for s in SEEDS], tool, param)[0])
table('Reads, coverage (0.1% error): estimate / 0.001',
      [f'{c}$\\times$' for c in COVERAGES], COVERAGES, FQ_TOOLS,
      lambda c: [f'fq_{c}_0.001_{s}' for s in SEEDS],
      lambda c, tool, param: READ_TRUTH)
table('Reads, error (20x): estimate / 0.001',
      [percent(e) for e in ERRORS], ERRORS, FQ_TOOLS,
      lambda e: [f'fq_20_{e}_{s}' for s in SEEDS],
      lambda e, tool, param: READ_TRUTH)

if problems:
    print(f'{len(problems)} missing or unusable results:', file=sys.stderr)
    for (sample, tool, param), reason in sorted(problems.items()):
        print(f'  {sample}\t{tool}\t{param}\t{reason}', file=sys.stderr)
