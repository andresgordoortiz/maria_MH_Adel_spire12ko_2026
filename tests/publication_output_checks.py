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
rows = list(csv.DictReader((root / 'source_data/Figure2D_shared_sample_PSI.csv').open()))
assert len(rows) == 110
assert {r['sample'] for r in rows if r['dataset'].startswith('SPIRE')} == {
    'Cont_1', 'Cont_2', 'Cont_7', 'DKO_4', 'DKO_5', 'DKO_6'}
rows = list(csv.DictReader((root / 'source_data/Figure1E_event_selection.csv').open()))
assert [r['GENE'] for r in rows] == ['Cit', 'Gen1', 'Nek4']
assert all(float(r['loo_min_signed_pp']) >= 10 for r in rows)
assert all(r['GENE'] not in ('Spire1', 'Spire2') for r in rows)


class Links(HTMLParser):
    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src') and not value.startswith(('http:', 'https:', '#', 'data:')):
                assert (root / value).is_file(), value


Links().feed((root / 'index.html').read_text())
print('PASS: output hashes, ZIP integrity, sample labels, example ranking and gallery links.')
