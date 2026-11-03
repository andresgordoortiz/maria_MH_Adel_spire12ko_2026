"""Verify the actual deliverable: hashes, table counts, archives and local links."""
import csv
import hashlib
from html.parser import HTMLParser
from pathlib import Path
import sys
import zipfile

root = Path(sys.argv[1] if len(sys.argv) > 1 else 'results/publication')
assert (root / 'PACKAGE_COMPLETE').is_file()
for row in csv.DictReader((root / 'output_checksums.csv').open()):
    path = root / row['file']
    assert path.resolve().is_relative_to(root.resolve())
    with path.open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == row['sha256'], path
    assert path.stat().st_size == int(row['bytes'])
for name in ('Supplementary_tables.zip', 'Figures.zip', 'Figures_SPIRE.zip', 'Figures_FMN2.zip', 'Figures_Convergence.zip'):
    with zipfile.ZipFile(root / name) as archive:
        assert archive.testzip() is None
        assert all(not p.startswith('/') and '..' not in Path(p).parts for p in archive.namelist())
rows = list(csv.DictReader((root / 'tables/S1_SPIRE_samples.csv').open()))
assert [r['sample'] for r in rows] == ['SPIRE_Cont_1', 'SPIRE_Cont_2', 'SPIRE_Cont_7',
                                      'SPIRE_DKO_4', 'SPIRE_DKO_5', 'SPIRE_DKO_6']
assert {r['experimental_unit'] for r in rows} == {'independent biological replicate'}
rows = list(csv.DictReader((root / 'source_data/Convergence_shared_sample_PSI.csv').open()))
assert len(rows) == 110
assert {r['sample'] for r in rows if r['dataset'].startswith('SPIRE')} == {
    'Cont_1', 'Cont_2', 'Cont_7', 'DKO_4', 'DKO_5', 'DKO_6'}
rows = list(csv.DictReader((root / 'tables/figure_manifest.csv').open()))
assert len(rows) == 13
assert {r['figure'] for r in rows} == {'SPIRE_Figure_1', 'SPIRE_Figure_2', *(f'SPIRE_S{i}' for i in range(1, 6)),
                                     'FMN2_Figure_1', *(f'FMN2_S{i}' for i in range(1, 5)), 'Convergence_Figure_1'}
rows = list(csv.DictReader((root / 'source_data/SPIRE_division_event_annotations.csv').open()))
assert len({r['EVENT'] for r in rows}) == 119 and len({r['GENE'] for r in rows}) == 75
assert len({r['GOALL'] for r in rows}) == 4
assert all(r['GENE'] not in ('Spire1', 'Spire2') for r in rows)
for study, events, samples in [('SPIRE', 1547, 6), ('FMN2', 275, 4)]:
    rows = list(csv.DictReader((root / f'source_data/{study}_clustered_PSI.csv').open()))
    assert len(rows) == events * samples
    assert len({r['EVENT'] for r in rows}) == events
    assert {int(r['row']) for r in rows} == set(range(1, events + 1))
    for kind, labels, regions, bins in [('EX', 30, 8, 68), ('IR', 15, 5, 88)]:
        rows = list(csv.DictReader((root / f'source_data/{study}_RBP_{kind}_heatmap.csv').open()))
        assert len(rows) == labels * regions * 2
        assert len({r['RBP'] for r in rows}) == labels
        assert all(r['position'] == {'up': '1', 'down': '2'}[r['direction']] for r in rows)
        assert {int(r['row']) for r in rows} == set(range(1, labels + 1))
        rows = list(csv.DictReader((root / f'source_data/{study}_GC_{kind}_profiles.csv').open()))
        assert len(rows) == bins * 3
        rows = list(csv.DictReader((root / f'source_data/{study}_GC_{kind}_difference.csv').open()))
        assert len(rows) == bins * 2 and all(float(r['se']) >= 0 for r in rows)
rows = list(csv.DictReader((root / 'tables/S7_intron_GC_audit_summary.csv').open()))
assert len(rows) == 6 and all(float(r['max_absolute_disagreement']) < 1e-12 for r in rows)
for path in (root / 'text').glob('*.md'):
    assert 'animal' not in path.read_text().lower(), path


class Links(HTMLParser):
    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src') and not value.startswith(('http:', 'https:', '#', 'data:')):
                assert (root / value).is_file(), value


Links().feed((root / 'index.html').read_text())
print('PASS: output hashes, ZIP integrity, sample labels, all division events, row clustering and sequence GC reconciliation and gallery links.')
