"""Audit completed rMAPS3 outputs; run in the existing rMAPS3 Pixi environment."""
import csv, json, math, sys
from fractions import Fraction
import re

from scipy.stats import hypergeom, mannwhitneyu

assert math.isclose(float(hypergeom.sf(0,10,2,2)),17/45)
from pathlib import Path

dataset = sys.argv[1]
assert dataset in ('spire', 'fmn2')
root = Path(f'results/{dataset}_final/rmaps3_outputs')
assert (root / 'COMPLETE').is_file()
se = ['UpstreamExon_3prime', 'UpstreamExonIntron', 'UpstreamIntron',
      'TargetExon_5prime', 'TargetExon-3prime', 'DownstreamIntron',
      'DownstreamExonIntron', 'DownstreamExon_5prime']
report = {}
for kind in ('EX', 'IR', 'Alt5', 'Alt3'):
    regions = se if kind == 'EX' else [f'R{i}' for i in range(1, 6 if kind == 'IR' else 9)]
    checked = 0
    for direction in ('up', 'dn'):
        with (root / kind / f'pVal.{direction}.vs.bg.RNAmap.txt').open() as handle:
            summary = list(csv.DictReader(handle, delimiter='\t'))
        files = list((root / kind / 'window_pvalues').glob(f'*.pVal.{direction}.vs.bg.txt'))
        assert len(files) == len(summary) == 114
        byname = {r['RBP']: r for r in summary}
        assert len(byname) == 114
        for path in files:
            key = path.name.removesuffix(f'.pVal.{direction}.vs.bg.txt')
            mins = dict.fromkeys(regions, 1.)
            with path.open() as handle:
                windows = list(csv.DictReader(handle, delimiter='\t'))
            assert len(windows) == {'EX': 1200, 'IR': 450, 'Alt5': 800, 'Alt3': 800}[kind]
            positions = set()
            for row in windows:
                p = float(list(row.values())[-1])
                assert math.isfinite(p) and 0 <= p <= 1
                position = row['Region'], row['position']
                assert position not in positions
                positions.add(position)
                mins[row['Region']] = min(mins[row['Region']], p)
            assert all(abs(float(a)-b) < 1e-14
                       for a, b in zip(list(byname[key].values())[1:], mins.values())), (kind, key)
            checked += len(windows)
    sequences = []
    for group, label in [('up', 'up'), ('dn', 'down'), ('bg', 'ndiff')]:
        with (root / 'inputs' / f'{dataset}_{kind}_{label}.txt').open() as handle:
            reader = csv.reader(handle, delimiter='\t')
            next(reader)
            headers = ['>' + ':'.join(row) for row in reader]
        fastas = sorted((root / kind / 'fasta').glob(f'{group}.*.fasta'))
        assert len(fastas) == {'EX': 7, 'IR': 3, 'Alt5': 5, 'Alt3': 5}[kind]
        for path in fastas:
            all_n, n_bases, bases = 0, 0, 0
            with path.open() as handle:
                for expected in headers:
                    assert handle.readline().strip() == expected, path
                    sequence = handle.readline().strip()
                    assert sequence, path
                    n = sequence.count('N')
                    all_n += n == len(sequence)
                    n_bases += n
                    bases += len(sequence)
                assert handle.readline() == '', path
            sequences.append(dict(file=path.name, records=len(headers), bases=bases,
                                  all_N_records=all_n, N_base_fraction=n_bases/bases))
    pdfs = list((root / kind / 'maps').glob('*.pdf'))
    pngs = list((root / kind / 'maps').glob('*.png'))
    assert len(pdfs) == len(pngs) == 114
    for path in pdfs:
        with path.open('rb') as handle:
            assert handle.read(4) == b'%PDF', path
    for path in pngs:
        with path.open('rb') as handle:
            assert handle.read(8) == b'\x89PNG\r\n\x1a\n', path
    report[kind] = dict(window_pvalues_checked=checked, summary_tables_match_recomputed_minima=True,
                        FASTA_headers_match_inputs=True, FASTA_regions=sequences,
                        motif_PDF_and_PNG_pairs=114)
(root / 'output_audit.json').write_text(json.dumps(report, indent=2) + '\n')
print(dataset, 'PASS:', sum(r['window_pvalues_checked'] for r in report.values()),
      'window p-values; all regional minima, FASTA record orders/counts and 456 PDF/PNG pairs checked.')
print('All-N sequence records:', sum(f['all_N_records'] for r in report.values() for f in r['FASTA_regions']))

layouts = {
 'EX': [('UpstreamExon','end',50),('UpstreamExonIntron','start',250),
        ('UpstreamIntron','end',250),('TargetExon','start',50),('TargetExon','end',50),
        ('DownstreamIntron','start',250),('DownstreamExonIntron','end',250),('DownstreamExon','start',50)],
 'IR': [('UpstreamExon','start',50),('UpstreamExon','end',50),('Intron','scaled',250),
        ('downstreamExon','start',50),('downstreamExon','end',50)],
 'Alt5': [('LongExon','start',50),('LongExon','end',50),('SplicingRegion','start',50),
          ('SplicingRegion','end',50),('FirstIntron','start',250),('FlankinfExonIntron','end',250),
          ('FlankingExon','start',50),('FlankingExon','end',50)],
 'Alt3': [('FlankingExon','start',50),('FlankingExon','end',50),('FirstIntron','start',250),
          ('secondIntron','end',250),('SplicingRegion','start',50),('SplicingRegion','end',50),
          ('downstreamExon','start',50),('downstreamExon','end',50)],
}
report = {}
for kind, layout in layouts.items():
    candidates = []
    for direction in ('up','dn'):
        with (root/kind/f'pVal.{direction}.vs.bg.RNAmap.txt').open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                candidates.extend((float(value),row['RBP'],region,direction,index)
                                  for index,(region,value) in enumerate(list(row.items())[1:]))
    pvalue,motif_id,region,direction,index = min(candidates)
    if kind == 'EX':
        region = ['UpstreamExon_3prime','UpstreamExonIntron','UpstreamIntron',
                  'TargetExon_5prime','TargetExon-3prime','DownstreamIntron',
                  'DownstreamExonIntron','DownstreamExon_5prime'][index]
    with (root/kind/'window_pvalues'/f'{motif_id}.pVal.{direction}.vs.bg.txt').open() as handle:
        rows = [r for r in csv.DictReader(handle,delimiter='\t') if r['Region']==region]
    best = min(rows,key=lambda r:float(list(r.values())[-1]))
    assert float(list(best.values())[-1]) == pvalue
    position = int(best['position'])
    sequence_region,mode,extent = layout[index]
    motif = re.compile(motif_id.split('.',1)[1])
    counts = {}
    for group in (direction,'bg'):
        values = []
        with (root/kind/'fasta'/f'{group}.{sequence_region}.fasta').open() as handle:
            for header in handle:
                assert header.startswith('>')
                sequence = next(handle).strip()
                start = (position if mode=='start' else len(sequence)-extent-50+position
                         if mode=='end' else Fraction(position*len(sequence),250))
                values.append(sum(start <= hit.start() and hit.end() <= start+50
                                  for hit in motif.finditer(sequence)))
        counts[group] = values
    calculated = float(mannwhitneyu(counts[direction],counts['bg'],alternative='greater').pvalue)
    assert math.isclose(calculated,pvalue,rel_tol=1e-10,abs_tol=0), (kind,pvalue,calculated)
    report[kind] = dict(motif=motif_id,region=region,position=position,direction=direction,
        reported_p=pvalue,independently_recomputed_p=calculated,
        foreground_events=len(counts[direction]),background_events=len(counts['bg']),
        foreground_total_hits=sum(counts[direction]),background_total_hits=sum(counts['bg']),
        foreground_events_with_hits=sum(x>0 for x in counts[direction]),
        background_events_with_hits=sum(x>0 for x in counts['bg']))
    if all(x<=1 for values in counts.values() for x in values):
        # Exact conditional permutation for binary counts; diagnostic only.
        report[kind]['exact_binary_allocation_p'] = float(hypergeom.sf(
            sum(counts[direction])-1,len(counts[direction])+len(counts['bg']),
            sum(counts[direction])+sum(counts['bg']),len(counts[direction])))
(root/'real_window_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(dataset,'PASS: strongest real window independently counted and tested in all four event types.')
print(json.dumps(report,indent=2))
