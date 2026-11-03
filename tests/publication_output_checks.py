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
for name in ('Supplementary_tables.zip', 'Figures.zip'):
    with zipfile.ZipFile(root / name) as archive:
        assert archive.testzip() is None
        assert all(not p.startswith('/') and '..' not in Path(p).parts for p in archive.namelist())
rows = list(csv.DictReader((root / 'tables/S1_SPIRE_samples.csv').open()))
assert [r['sample'] for r in rows] == ['SPIRE_Cont_1', 'SPIRE_Cont_2', 'SPIRE_Cont_7',
                                      'SPIRE_DKO_4', 'SPIRE_DKO_5', 'SPIRE_DKO_6']
assert {r['experimental_unit'] for r in rows} == {'independent animal'}
rows = list(csv.DictReader((root / 'source_data/FigureS3D_shared_sample_PSI.csv').open()))
assert len(rows) == 110
assert {r['sample'] for r in rows if r['dataset'].startswith('SPIRE')} == {
    'Cont_1', 'Cont_2', 'Cont_7', 'DKO_4', 'DKO_5', 'DKO_6'}
rows = list(csv.DictReader((root / 'source_data/Figure1E_event_selection.csv').open()))
assert [r['GENE'] for r in rows] == ['Cit', 'Gen1', 'Nek4']
assert all(float(r['loo_min_signed_pp']) >= 10 for r in rows)
assert all(r['GENE'] not in ('Spire1', 'Spire2') for r in rows)


rows = list(csv.DictReader((root / 'tables/figure_manifest.csv').open()))
assert {r['figure'] for r in rows} == {'Figure_1', 'Figure_2', *(f'Figure_S{i}' for i in range(1, 8))}
rows = list(csv.DictReader((root / 'source_data/Figure1D_division_PSI.csv').open()))
assert len(rows) == 90 and len({r['GENE'] for r in rows}) == 15
assert all(0 <= float(r['PSI']) <= 100 for r in rows)
rows = list(csv.DictReader((root / 'source_data/Figure2A_feature_effects.csv').open()))
assert len(rows) == 24 and len({r['feature'] for r in rows}) == 6
assert {r['group1'] for r in rows} == {'up', 'down'}
for kind, labels, regions, bins in [('EX', 30, 8, 68), ('IR', 15, 5, 88)]:
    rows = list(csv.DictReader((root / f'source_data/RBP_{kind}_heatmap.csv').open()))
    assert len(rows) == labels * regions * 4
    assert len({r['RBP'] for r in rows}) == labels
    positions = {('spire', 'up'): '1', ('spire', 'down'): '2', ('fmn2', 'up'): '3', ('fmn2', 'down'): '4'}
    assert all(r['position'] == positions[(r['dataset'], r['direction'])] for r in rows)
    assert {r['dataset'] for r in rows} == {'spire', 'fmn2'}
    rows = list(csv.DictReader((root / f'source_data/GC_{kind}_profiles.csv').open()))
    assert len(rows) == bins * 6
    assert {r['category'] for r in rows} == {'up', 'down', 'ndiff'}


class Links(HTMLParser):
    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src') and not value.startswith(('http:', 'https:', '#', 'data:')):
                assert (root / value).is_file(), value


Links().feed((root / 'index.html').read_text())
print('PASS: output hashes, ZIP integrity, sample labels, candidate ranking, feature/RBP/GC coverage and gallery links.')
