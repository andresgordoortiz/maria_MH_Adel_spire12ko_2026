"""Compare regional presence with an explicit search over every native window."""
import importlib.util
from fractions import Fraction
import math
from pathlib import Path
import re

spec = importlib.util.spec_from_file_location('regional_motifs', Path('analysis/regional_motifs.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
for sequence in ('ACGTAACCGGTT', 'A'*70+'CGT'+'A'*120, 'CGT'+'A'*100+'CGT'):
    for pattern in (re.compile('CGT'), re.compile('AA'), re.compile('C[CG]')):
        for mode in ('start', 'end', 'scaled'):
            extent, width = 7, 5
            expected = False
            for j in range(extent):
                start = j if mode == 'start' else len(sequence)-extent-width+j if mode == 'end' else Fraction(j*len(sequence), extent)
                tract = sequence[max(0, math.ceil(start)):max(0, min(len(sequence), math.floor(start+width)))]
                expected |= pattern.search(tract) is not None
            assert module.present(sequence, [pattern], mode, extent, width) == expected
assert module.present('GGGACGTTT', [re.compile('TTT'), re.compile('ACG')], 'start', 5, 5)
assert module.present('CCCAAAAAAA', [re.compile('AAA')], 'scaled', 2, 3)
print('PASS: fixed-region presence matches explicit windows, both edges, scaled gaps and motif unions.')
