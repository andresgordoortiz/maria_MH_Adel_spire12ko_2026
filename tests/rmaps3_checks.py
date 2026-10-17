"""Coordinate conversion, checksum and output-boundary checks; standard library only."""
import csv
import importlib.util
from pathlib import Path
import tempfile

spec = importlib.util.spec_from_file_location('rmaps3', 'analysis/run_rmaps3.py')
rmaps3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rmaps3)


def row(strand='+', start=200, end=249):
    return dict(chr='chr1', strand=strand, firstExonStart=100, firstExonEnd=199,
                exonStart=start, exonEnd=end, secondExonStart=250, secondExonEnd=349)


for strand in ('+', '-'):
    assert rmaps3.convert_row(row(strand), 'EX') == ['chr1', strand, 199, 249, 99, 199, 249, 349]
    assert rmaps3.convert_row(row(strand), 'IR') == ['chr1', strand, 99, 349, 99, 199, 249, 349]
for kind, strand in [('Alt5', '+'), ('Alt3', '-')]:
    assert rmaps3.convert_row(row(strand), kind) == ['chr1', strand, 99, 249, 99, 199, 249, 349]
for kind, strand in [('Alt3', '+'), ('Alt5', '-')]:
    assert rmaps3.convert_row(row(strand), kind) == ['chr1', strand, 199, 349, 249, 349, 99, 199]
assert rmaps3.convert_row(row(start=200, end=200), 'EX')[2:4] == [199, 200]
for bad, kind in [(row('?'), 'EX'), (row(start=199), 'EX'),
                  (row(start=201), 'IR'), (row(start=201), 'Alt5'),
                  (row(end=248), 'Alt3')]:
    try:
        rmaps3.convert_row(bad, kind)
        raise AssertionError('Invalid geometry accepted')
    except ValueError:
        pass

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    analysis = root / 'analysis'
    (analysis / 'external_inputs').mkdir(parents=True)
    (analysis / 'tables').mkdir()
    checks = []
    for kind in rmaps3.EVENTS:
        for index, group in enumerate(('up', 'down', 'ndiff')):
            value = row()
            for key in value:
                if key not in ('chr', 'strand'):
                    value[key] += index * 1000
            path = analysis / f'external_inputs/test_{kind}_{group}.txt'
            rmaps3.write_rows(path, list(value), [list(value.values())])
            checks.append([path.name, rmaps3.file_hash(path)])
            rmaps3.write_rows(analysis / f'tables/rmaps_{kind}_{group}_event_manifest.csv',
                             ['EVENT', 'GENE', 'DATASET'], [[f'{kind}{index}', 'Gene', group]])
    rmaps3.write_rows(analysis / 'external_inputs/input_checksums.csv', ['file', 'sha256'], checks)
    fasta = root / 'mm10.fa'
    Path(str(fasta) + '.fai').write_text('chr1\t100000\t0\t60\t61\n')
    sources, counts = rmaps3.prepare(analysis, 'test', root / 'output', fasta)
    assert len(sources) == 24 and len(counts) == 12
    assert len(list((root / 'output/inputs').glob('*.txt'))) == 12
    try:
        rmaps3.prepare(analysis, 'test', root / 'output', fasta)
        raise AssertionError('Existing output overwritten')
    except FileExistsError:
        pass
    with (analysis / 'external_inputs/test_EX_up.txt').open('a') as handle:
        handle.write('\n')
    try:
        rmaps3.prepare(analysis, 'test', root / 'bad-checksum', fasta)
        raise AssertionError('Changed input accepted')
    except ValueError as error:
        assert 'checksum' in str(error)
    assert not (root / 'bad-checksum').exists()
print('rMAPS3 conversion (+/- for all types), one-base intervals, invalid geometry, checksums and overwrite checks passed.')
