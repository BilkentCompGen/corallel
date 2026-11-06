# usage: python3 sim.py ref.fa out.fa divergence indel_rate sv_count seed
import sys
import numpy as np

ref_path, out_path = sys.argv[1], sys.argv[2]
divergence, indel_rate = float(sys.argv[3]), float(sys.argv[4])
sv_count, seed = int(sys.argv[5]), int(sys.argv[6])

DELETION, INSERTION, INVERSION, DUPLICATION = 0, 1, 2, 3
N = 4
BASES = np.frombuffer(b'ACGTN', np.uint8)
COMPLEMENT = np.array([3, 2, 1, 0, N], np.uint8)
ENCODE = np.full(256, N, np.uint8)
ENCODE[np.frombuffer(b'ACGTacgt', np.uint8)] = [0, 1, 2, 3, 0, 1, 2, 3]

rng = np.random.default_rng(seed)

chromosomes = []
for line in open(ref_path):
    if line.startswith('>'):
        chromosomes.append((line[1:].split()[0], []))
    else:
        chromosomes[-1][1].append(line.strip())

genome_length = sum(len(line) for _, lines in chromosomes for line in lines)

with open(out_path, 'w') as out:
    for name, lines in chromosomes:
        sequence = ENCODE[np.frombuffer(''.join(lines).encode(), np.uint8)]
        length = sequence.size

        snp_positions = np.flatnonzero((rng.random(length) < divergence) & (sequence != N))
        sequence[snp_positions] = (sequence[snp_positions] + rng.integers(1, 4, snp_positions.size)) % 4

        indel_count = rng.poisson(indel_rate * length)
        events = [(position, int(rng.integers(1, 11)), int(rng.integers(2)))
                  for position in rng.integers(0, length, indel_count)]

        chromosome_sv_count = rng.poisson(sv_count * length / genome_length)
        events += [(position, int(10 ** rng.uniform(np.log10(50), 4)), int(rng.integers(4)))
                   for position in rng.integers(0, length, chromosome_sv_count)]

        pieces, cursor = [], 0
        for position, size, kind in sorted(events):
            if position < cursor:
                continue
            pieces.append(sequence[cursor:position])
            if kind == DELETION:
                cursor = position + size
            elif kind == INSERTION:
                pieces.append(rng.integers(0, 4, size).astype(np.uint8))
                cursor = position
            elif kind == INVERSION:
                pieces.append(COMPLEMENT[sequence[position:position + size][::-1]])
                cursor = position + size
            else:
                pieces += [sequence[position:position + size], sequence[position:position + size]]
                cursor = position + size
        pieces.append(sequence[cursor:])

        mutated = BASES[np.concatenate(pieces)].tobytes().decode()
        out.write(f'>{name}\n')
        out.writelines(mutated[i:i + 80] + '\n' for i in range(0, len(mutated), 80))
