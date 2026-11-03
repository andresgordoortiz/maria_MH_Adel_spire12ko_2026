"""Recount the two reported associations by explicit native windows."""
import csv
import hashlib
from pathlib import Path
import re

root = Path('reports/regional_motif_presence')
assert (root/'COMPLETE').is_file()
for row in csv.DictReader((root/'input_checksums.csv').open()):
    with Path(row['file']).open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == row['sha256'], row['file']
rows = list(csv.DictReader((root/'regional_presence.csv').open()))
assert len(rows) == 2*91*29*3
for kind, rbp, region, mode, extent, positive, background in [
    ('EX', 'PCBP2', 'R1', 'end', 50, 48, 7654),
    ('Alt5', 'PPRC1', 'R5', 'start', 250, 16, 1530),
]:
    for group, expected in [('up', positive), ('bg', background)]:
        row, = [r for r in rows if (r['dataset'], r['event_type'], r['RBP'], r['region'], r['group'])
                == ('spire', kind, rbp, region, group)]
        patterns = [re.compile(p) for p in row['variants'].split(';')]
        n, found = 0, 0
        with Path(row['fasta']).open() as handle:
            for header in handle:
                assert header.startswith('>')
                sequence = next(handle).strip()
                n += 1
                present = False
                for j in range(extent):
                    start = j if mode == 'start' else len(sequence)-extent-50+j
                    tract = sequence[max(0, start):max(0, min(len(sequence), start+50))]
                    for pattern in patterns:
                        if pattern.search(tract) is not None:
                            present = True
                            break
                    if present:
                        break
                found += present
        assert found == int(row['positive']) == expected and n == int(row['n']), (kind, group, found, n)
print('PASS: all presence-input hashes and independent foreground/background window recounts for PCBP2 and PPRC1.')
