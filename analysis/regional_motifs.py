"""Count RBP motif presence over the existing, fixed rMAPS search regions."""
import concurrent.futures
import csv
import hashlib
import json
from pathlib import Path
import re
import sys


LAYOUTS = {
    'EX': [('UpstreamExon', 'end', 50), ('UpstreamExonIntron', 'start', 250),
           ('UpstreamIntron', 'end', 250), ('TargetExon', 'start', 50),
           ('TargetExon', 'end', 50), ('DownstreamIntron', 'start', 250),
           ('DownstreamExonIntron', 'end', 250), ('DownstreamExon', 'start', 50)],
    'IR': [('UpstreamExon', 'start', 50), ('UpstreamExon', 'end', 50),
           ('Intron', 'scaled', 250), ('downstreamExon', 'start', 50), ('downstreamExon', 'end', 50)],
    'Alt5': [('LongExon', 'start', 50), ('LongExon', 'end', 50),
             ('SplicingRegion', 'start', 50), ('SplicingRegion', 'end', 50),
             ('FirstIntron', 'start', 250), ('FlankinfExonIntron', 'end', 250),
             ('FlankingExon', 'start', 50), ('FlankingExon', 'end', 50)],
    'Alt3': [('FlankingExon', 'start', 50), ('FlankingExon', 'end', 50),
             ('FirstIntron', 'start', 250), ('secondIntron', 'end', 250),
             ('SplicingRegion', 'start', 50), ('SplicingRegion', 'end', 50),
             ('downstreamExon', 'start', 50), ('downstreamExon', 'end', 50)],
}


def present(sequence, patterns, mode, extent, window=50):
    """Union of motif variants and searched windows; each sequence contributes once."""
    if mode != 'scaled':
        n = len(sequence)
        tract = sequence[:extent + window - 1] if mode == 'start' else sequence[max(0, n-extent-window):n-1]
        return any(pattern.search(tract) for pattern in patterns)
    n = len(sequence)
    for pattern in patterns:
        cursor = 0
        while (hit := pattern.search(sequence, cursor)) is not None:
            first = max(0, ((hit.end()-window)*extent + n-1)//n)
            last = min(extent-1, hit.start()*extent//n)
            if first <= last:
                return True
            cursor = hit.start()+1
    return False


def count_file(task):
    dataset, kind, group, path, slots, motifs = task
    patterns = {name: [re.compile(m) for m in variants] for name, variants in motifs.items()}
    # Non-scaled tracts need one search for the union of catalogue variants.
    unions = {name: [re.compile('|'.join('(?:'+m+')' for m in variants))] for name, variants in motifs.items()}
    counts = {(index, name): 0 for index, _, _ in slots for name in motifs}
    n = 0
    with Path(path).open() as handle:
        for header in handle:
            sequence = next(handle, '').strip().upper()
            if not header.startswith('>') or not sequence or re.search('[^ACGTN]', sequence):
                raise ValueError(f'Invalid FASTA record: {path}, {n+1}')
            n += 1
            for index, mode, extent in slots:
                for name in motifs:
                    counts[index, name] += present(sequence, patterns[name] if mode == 'scaled' else unions[name], mode, extent)
    return [dict(dataset=dataset, event_type=kind, group=group, region=f'R{index+1}',
                 RBP=name, n=n, positive=k, variants=';'.join(motifs[name]), fasta=path)
            for (index, name), k in counts.items()]


def run(output):
    if output.exists():
        raise ValueError(f'Choose a new output directory: {output}')
    output.mkdir(parents=True)
    tasks, inputs = [], set()
    for dataset in ('spire', 'fmn2'):
        root = Path(f'reports/{dataset}_final/rmaps3_outputs')
        if not (root/'COMPLETE').is_file():
            raise ValueError(f'Incomplete rMAPS inputs: {root}')
        inputs.add(root/'provenance.json')
        for kind, layout in LAYOUTS.items():
            path = root/kind/'pVal.up.vs.bg.RNAmap.txt'
            inputs.add(path)
            with path.open() as handle:
                names = [row['RBP'] for row in csv.DictReader(handle, delimiter='\t')]
            assert len(names) == len(set(names)) == 114
            motifs = {}
            for name in names:
                rbp, pattern = name.split('.', 1)
                assert re.fullmatch(r'[ACGT\[\]]+', pattern)
                assert len(re.sub(r'\[[ACGT]+\]', 'A', pattern)) <= 50
                motifs.setdefault(rbp, []).append(pattern)
            assert len(motifs) == 91
            for group in ('up', 'dn', 'bg'):
                for region in dict.fromkeys(row[0] for row in layout):
                    path = root/kind/'fasta'/f'{group}.{region}.fasta'
                    assert path.is_file(), path
                    inputs.add(path)
                    slots = [(i, mode, extent) for i, (name, mode, extent) in enumerate(layout) if name == region]
                    tasks.append((dataset, kind, group, str(path), slots, motifs))
    rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        for result in pool.map(count_file, tasks):
            rows.extend(result)
            print(f'Counted {result[0]["dataset"]} {result[0]["event_type"]} {result[0]["group"]} {result[0]["fasta"]}', flush=True)
    with (output/'regional_presence.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    inputs.add(Path(__file__).relative_to(Path.cwd()) if Path(__file__).is_absolute() else Path(__file__))
    with (output/'input_checksums.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['file', 'sha256'])
        for path in sorted(inputs):
            with path.open('rb') as source:
                writer.writerow([str(path), hashlib.file_digest(source, 'sha256').hexdigest()])
    (output/'COMPLETE').write_text('Fixed-region RBP presence counted; no window P selection.\n')
    print(f'Complete: {len(rows)} RBP/region/group counts in {output}', flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python3 analysis/regional_motifs.py NEW_OUTPUT_DIRECTORY')
    run(Path(sys.argv[1]))
