#!/usr/bin/env python3
"""Run local rMAPS3 in its Pixi environment on the checked final splicing exports."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys


EVENTS = {'EX': 'se', 'IR': 'ri', 'Alt5': 'a5ss', 'Alt3': 'a3ss'}


def file_hash(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def convert_row(row, kind):
    left, le, start, end, right, re = (int(row[k]) for k in (
        'firstExonStart', 'firstExonEnd', 'exonStart', 'exonEnd',
        'secondExonStart', 'secondExonEnd'))
    strand = row['strand']
    if strand not in ('+', '-') or not 0 < left <= le < start <= end < right <= re:
        raise ValueError(f'Invalid strand or exon geometry: {row}')
    if kind == 'EX':
        coords = start - 1, end, left - 1, le, right - 1, re
    elif kind == 'IR':
        if start != le + 1 or end != right - 1:
            raise ValueError('Retained intron must join both flanking exons')
        coords = left - 1, re, left - 1, le, right - 1, re
    elif kind in ('Alt5', 'Alt3'):
        extends_left = (kind == 'Alt5') == (strand == '+')
        if extends_left:
            if start != le + 1:
                raise ValueError('Alternative exon segment must adjoin the left exon')
            coords = left - 1, end, left - 1, le, right - 1, re
        else:
            if end != right - 1:
                raise ValueError('Alternative exon segment must adjoin the right exon')
            coords = start - 1, re, right - 1, re, left - 1, le
    else:
        raise ValueError(f'Unsupported event type: {kind}')
    return [row['chr'], strand, *coords]


def write_rows(path, header, rows):
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t' if path.suffix == '.txt' else ',')
        writer.writerow(header)
        writer.writerows(rows)


def prepare(analysis_dir, dataset, output, fasta):
    if output.exists():
        raise FileExistsError(f'Output must be new: {output}')
    with (analysis_dir / 'external_inputs/input_checksums.csv').open() as handle:
        expected = {r['file']: r['sha256'] for r in csv.DictReader(handle)}
    with Path(str(fasta) + '.fai').open() as handle:
        lengths = {r[0]: int(r[1]) for r in csv.reader(handle, delimiter='\t')}
    converted, manifest, sources, counts = {}, [], {}, []
    for kind in EVENTS:
        seen, geometries = set(), set()
        for group in ('up', 'down', 'ndiff'):
            source = analysis_dir / f'external_inputs/{dataset}_{kind}_{group}.txt'
            mapping = analysis_dir / f'tables/rmaps_{kind}_{group}_event_manifest.csv'
            digest = file_hash(source)
            if digest != expected[source.name]:
                raise ValueError(f'Input checksum mismatch: {source}')
            with source.open() as a, mapping.open() as b:
                rows = list(csv.DictReader(a, delimiter='\t'))
                events = list(csv.DictReader(b))
            if len(rows) != len(events) or not rows:
                raise ValueError(f'Empty or mismatched coordinate/event files: {source}')
            values = []
            for row, event in zip(rows, events):
                if not event['EVENT'] or event['EVENT'] in seen or event['DATASET'] != group:
                    raise ValueError(f'Duplicate event or incorrect direction: {event}')
                seen.add(event['EVENT'])
                value = convert_row(row, kind)
                if tuple(value) in geometries:
                    raise ValueError(f'Duplicate event geometry: {event}')
                geometries.add(tuple(value))
                if value[0] not in lengths or max(value[2:]) + 300 > lengths[value[0]]:
                    raise ValueError(f'Coordinates outside mm10 reference: {event}')
                values.append(value)
                manifest.append([kind, group, len(values), event['EVENT'], event['GENE']])
            converted[f'{dataset}_{kind}_{group}.txt'] = values
            counts.append([kind, group, len(values)])
            sources[str(source)] = digest
            sources[str(mapping)] = file_hash(mapping)
    (output / 'inputs').mkdir(parents=True)
    for filename, rows in converted.items():
        kind = filename.split('_')[1]
        middle = ['exonStart', 'exonEnd', 'firstExonStart', 'firstExonEnd', 'secondExonStart', 'secondExonEnd']
        if kind == 'IR':
            middle = ['riExonStart', 'riExonEnd', 'upstreamExonStart', 'upstreamExonEnd', 'downstreamExonStart', 'downstreamExonEnd']
        elif kind in ('Alt5', 'Alt3'):
            middle = ['longExonStart', 'longExonEnd', 'shortExonStart', 'shortExonEnd', 'flankingExonStart', 'flankingExonEnd']
        write_rows(output / 'inputs' / filename, ['chr', 'strand', *middle], rows)
    write_rows(output / 'event_manifest.csv', ['type', 'direction', 'input_row', 'EVENT', 'GENE'], manifest)
    write_rows(output / 'event_counts.csv', ['type', 'direction', 'events'], counts)
    shutil.copy2(analysis_dir / 'external_inputs/input_checksums.csv', output / 'input_checksums.csv')
    return sources, counts


def validate_result(output, kind, motifs):
    regions = 5 if kind == 'IR' else 8
    for direction in ('up', 'dn'):
        with (output / f'pVal.{direction}.vs.bg.RNAmap.txt').open() as handle:
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader)
            rows = list(reader)
        if len(header) != regions + 1 or len(rows) != motifs:
            raise ValueError(f'Unexpected rMAPS3 result dimensions: {output}')
        for row in rows:
            if len(row) != len(header) or not all(math.isfinite(float(p)) and 0 <= float(p) <= 1 for p in row[1:]):
                raise ValueError(f'Invalid rMAPS3 p-values: {output}')
    if len(list((output / 'maps').glob('*.pdf'))) != motifs:
        raise ValueError(f'Missing motif-map PDFs: {output}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis_dir', type=Path)
    parser.add_argument('dataset', choices=('spire', 'fmn2'))
    parser.add_argument('output_dir', type=Path)
    parser.add_argument('rmaps_repo', type=Path)
    args = parser.parse_args()
    analysis_dir, output, repo = (p.resolve() for p in (args.analysis_dir, args.output_dir, args.rmaps_repo))
    fasta = repo / 'genomedata/mm10/mm10.fa'
    motifs = repo / 'data/knownMotifs.human.mouse.txt'
    for path in (fasta, Path(str(fasta) + '.fai'), motifs, repo / 'cli.py', repo / 'pixi.lock'):
        if not path.is_file():
            raise FileNotFoundError(path)
    sources, counts = prepare(analysis_dir, args.dataset, output, fasta)
    for path in (fasta, Path(str(fasta) + '.fai'), motifs, repo / 'pixi.lock', repo / 'pixi.toml', Path(__file__).resolve()):
        sources[str(path)] = file_hash(path)
    for path in [repo/'cli.py', *sorted((repo/'rmaps_core').glob('*.py')), *sorted((repo/'legacy').glob('motifMap*_MP.py'))]:
        sources[str(path)] = file_hash(path)
    with motifs.open() as handle:
        motif_count = sum(1 for line in handle if line.strip()) - 1
    env = dict(os.environ)
    env.setdefault('PYTHON_CPU_COUNT', '4')
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    record = {'dataset': args.dataset, 'python': sys.version, 'executable': sys.executable,
        'rmaps_commit': subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip(),
        'rmaps_version': subprocess.check_output([sys.executable, str(repo / 'cli.py'), '--version'], env=env, text=True).strip(),
        'counter': 'corrected_half_open_windows' if (repo/'rmaps_core/corrected_motif_counter.py').exists() else 'unmodified_beta',
        'stat_method': 'mannwhitney_greater',
        'source_sha256': sources, 'converted_sha256': {p.name: file_hash(p) for p in (output / 'inputs').glob('*.txt')},
        'PYTHON_CPU_COUNT': env['PYTHON_CPU_COUNT'], 'motifs': motif_count, 'commands': {}}
    for kind, event in EVENTS.items():
        cmd = [sys.executable, str(repo / 'cli.py'), 'motif-map', event,
            '--known-motifs', str(motifs), '--motifs', 'NA', '--fasta-root', str(repo / 'genomedata'),
            '--genome', 'mm10', '--output', str(output / kind),
            '--up', str(output / f'inputs/{args.dataset}_{kind}_up.txt'),
            '--down', str(output / f'inputs/{args.dataset}_{kind}_down.txt'),
            '--background', str(output / f'inputs/{args.dataset}_{kind}_ndiff.txt'),
            '--intron', '250', '--exon', '50', '--window', '50', '--step', '1',
            '--sig-fdr', '0.05', '--sig-delta-psi', '0.1', '--stat-method', 'mannwhitney_greater', '--keep-temp']
        record['commands'][kind] = cmd
        (output / 'provenance.json').write_text(json.dumps(record, indent=2) + '\n')
        print(f'Starting {args.dataset} {kind}: {motif_count} motifs; all ndiff background events retained', flush=True)
        with (output / f'{kind}.cli.log').open('w') as log:
            subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        validate_result(output / kind, kind, motif_count)
        windows = output / kind / 'window_pvalues'
        windows.mkdir()
        for path in (output / kind / 'temp').glob('*.pVal.*.txt'):
            shutil.move(str(path), windows / path.name)
        # Per-event count matrices are reproducible and very large; preserve window p-values and FASTAs.
        shutil.rmtree(output / kind / 'temp')
        (output / kind / 'COMPLETE').write_text('CLI exited successfully; motif tables and PDFs validated.\n')
        print(f'Completed {args.dataset} {kind}', flush=True)
    (output / 'COMPLETE').write_text('All four rMAPS3 motif analyses complete; original event groups and sample labels retained.\n')


if __name__ == '__main__':
    main()
