#!/usr/bin/env python3
"""Validate completed Matt feature tables against their exact input records."""
import collections
import csv
import pathlib
import sys
from prepare_references import file_hash

if len(sys.argv) != 3:
    raise SystemExit('Usage: python analysis/validate_matt.py ANALYSIS_DIR MATT_DIR')
source, target = map(pathlib.Path, sys.argv[1:])
dataset = (source / 'DATASET').read_text().strip()
assert dataset in ('spire', 'fmn2')
csv.field_size_limit(sys.maxsize)
with (target / 'input_checksums.csv').open() as f:
    for check in csv.DictReader(f):
        assert file_hash(source / 'external_inputs' / check['file']) == check['sha256']
for kind, suffix in [('EX', 'efeatures'), ('IR', 'ifeatures')]:
    with (source / 'external_inputs' / f'{dataset}_{kind}_matt.tab').open() as f:
        expected = {r['EVENT']: r for r in csv.DictReader(f, delimiter='\t')}
    seen = set()
    counts = collections.Counter()
    with (target / kind / f'{dataset}_{kind}_matt_with_{suffix}.tab').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            event = row['EVENT']
            assert event in expected and event not in seen, event
            assert all(row[k] == v for k, v in expected[event].items()), event
            assert float(row['EXON_LENGTH' if kind == 'EX' else 'INTRON_LENGTH']) > 0
            counts[row['DATASET']] += 1
            seen.add(event)
    assert seen and counts['ndiff'] > 0
    assert (target / kind / 'summary.pdf').stat().st_size > 0
    print(kind, dict(counts))
