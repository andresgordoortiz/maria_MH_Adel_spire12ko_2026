"""Coordinate, classification, manifest and checksum checks for the handoff."""
import contextlib
import csv
import io
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from prepare_rnamotifs import convert_row, file_hash, prepare

row = dict(chr='chr1', strand='+', firstExonStart=101, firstExonEnd=150,
           exonStart=201, exonEnd=250, secondExonStart=301, secondExonEnd=350)
assert convert_row(row, 'EX', 1, 'up') == [1, 1, 'chr1', '+', 150, 200, 250, 300, 2]
assert convert_row({**row, 'strand': '-'}, 'EX', 2, 'down') == [
    2, 2, 'chr1', '-', 150, 200, 250, 300, -2]
ir = {**row, 'exonStart': 151, 'exonEnd': 300}
for strand in ('+', '-'):
    assert convert_row({**ir, 'strand': strand}, 'IR', 3, 'ndiff') == [
        3, 3, 'chr1', strand, 100, 150, 300, 350, 0, 'RI']
for invalid, kind in (({**row, 'firstExonEnd': 201}, 'EX'),
                      ({**row, 'strand': '?'}, 'EX'),
                      ({**ir, 'exonStart': 152}, 'IR')):
    try:
        convert_row(invalid, kind, 1, 'up')
        raise AssertionError('Invalid coordinates or strand accepted')
    except ValueError:
        pass

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    inputs, tables = root / 'external_inputs', root / 'tables'
    inputs.mkdir(); tables.mkdir()
    hashes = []
    for kind in ('EX', 'IR'):
        for group in ('up', 'down', 'ndiff'):
            path = inputs / f'spire_{kind}_{group}.txt'
            with path.open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=list(row), delimiter='\t')
                writer.writeheader(); writer.writerow(ir if kind == 'IR' else row)
            hashes.append([path.name, file_hash(path)])
            with (tables / f'rmaps_{kind}_{group}_event_manifest.csv').open('w') as f:
                f.write(f'EVENT,GENE,DATASET\n{kind}_{group},Gene,{group}\n')
    with (inputs / 'input_checksums.csv').open('w', newline='') as f:
        writer = csv.writer(f); writer.writerow(['file', 'sha256']); writer.writerows(hashes)
    with contextlib.redirect_stdout(io.StringIO()):
        prepare(root, 'spire', root / 'converted')
    lines = (root / 'converted/spire_EX.txt').read_text().splitlines()
    assert [line.split(';')[-1] for line in lines] == ['2', '-2', '0']
    assert (root / 'converted/COMPLETE').exists()
    try:
        prepare(root, 'spire', root / 'converted')
        raise AssertionError('Existing output overwritten')
    except FileExistsError:
        pass
    path.write_text(path.read_text() + '\n')
    try:
        prepare(root, 'spire', root / 'bad')
        raise AssertionError('Modified input accepted')
    except ValueError:
        assert not (root / 'bad').exists()

print('RNAmotifs handoff checks passed')
