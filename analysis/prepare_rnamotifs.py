#!/usr/bin/env python3
"""Convert this analysis's mm10 rMAPS exports to native RNAmotifs SE/RI input."""
import argparse
import csv
import hashlib
from pathlib import Path


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def convert_row(row, kind, index, group):
    left_start, left_end, start, end, right_start, right_end = (
        int(row[key]) for key in ('firstExonStart', 'firstExonEnd', 'exonStart',
                                 'exonEnd', 'secondExonStart', 'secondExonEnd'))
    if row['strand'] not in ('+', '-') or not row['chr'].startswith('chr'):
        raise ValueError('Invalid chromosome or strand')
    if not 0 < left_start <= left_end < start <= end < right_start <= right_end:
        raise ValueError('Invalid or overlapping event coordinates')
    if kind == 'IR' and (start != left_end + 1 or end != right_start - 1):
        raise ValueError('Retained intron must join the flanking exon boundaries')
    # VASTDB is 1-based closed; RNAmotifs uses genomic-order, 0-based boundaries.
    coordinates = (left_end, start - 1, end, right_start - 1) if kind == 'EX' else (
        left_start - 1, left_end, right_start - 1, right_end)
    return [index, index, row['chr'], row['strand'], *coordinates,
            {'up': 2, 'down': -2, 'ndiff': 0}[group], *(['RI'] if kind == 'IR' else [])]


def prepare(analysis_dir, dataset, output_dir):
    if output_dir.exists():
        raise FileExistsError(f'Output must be new: {output_dir}')
    inputs = analysis_dir / 'external_inputs'
    with (inputs / 'input_checksums.csv').open() as f:
        expected = {r['file']: r['sha256'] for r in csv.DictReader(f)}
    converted, manifest, checks, counts = {}, [], [], []
    for kind in ('EX', 'IR'):
        rows, seen = [], set()
        for group in ('up', 'down', 'ndiff'):
            path = inputs / f'{dataset}_{kind}_{group}.txt'
            digest = file_hash(path)
            if digest != expected[path.name]:
                raise ValueError(f'Input checksum mismatch: {path}')
            events_path = analysis_dir / 'tables' / f'rmaps_{kind}_{group}_event_manifest.csv'
            with path.open() as f, events_path.open() as g:
                coordinates = list(csv.DictReader(f, delimiter='\t'))
                events = list(csv.DictReader(g))
            if len(coordinates) != len(events):
                raise ValueError(f'Coordinate/manifest count mismatch: {path}')
            for row, event in zip(coordinates, events):
                if not event['EVENT'] or event['EVENT'] in seen or event['DATASET'] != group:
                    raise ValueError(f'Duplicate event or mismatched group: {events_path}')
                seen.add(event['EVENT'])
                index = len(rows) + 1
                rows.append(convert_row(row, kind, index, group))
                manifest.append([kind, index, event['EVENT'], event['GENE'], group])
            for source in (path, events_path):
                checks.append([str(source.resolve()), file_hash(source)])
            counts.append([kind, group, len(coordinates)])
        converted[kind] = rows
    output_dir.mkdir(parents=True)
    for kind, rows in converted.items():
        with (output_dir / f'{dataset}_{kind}.txt').open('w', newline='') as f:
            csv.writer(f, delimiter=';', lineterminator='\n').writerows(rows)
    for name, header, rows in (
        ('event_manifest.csv', ['type', 'row_id', 'EVENT', 'GENE', 'direction'], manifest),
        ('source_checksums.csv', ['file', 'sha256'], checks),
        ('counts.csv', ['type', 'direction', 'events'], counts),
    ):
        with (output_dir / name).open('w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
    (output_dir / 'COMPLETE').write_text('Converted EX/IR only; all ndiff events retained.\n')
    for kind, group, n in counts:
        print(f'{dataset} {kind} {group}: {n}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis_dir', type=Path)
    parser.add_argument('dataset', choices=('spire', 'fmn2'))
    parser.add_argument('output_dir', type=Path)
    args = parser.parse_args()
    prepare(args.analysis_dir, args.dataset, args.output_dir)
