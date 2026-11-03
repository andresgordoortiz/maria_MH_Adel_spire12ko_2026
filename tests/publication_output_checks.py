"""Verify the actual deliverable: hashes, table counts, archives and local links."""
import csv
import hashlib
import math
from html.parser import HTMLParser
from pathlib import Path
import sys
import statistics
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
assert len(rows) == 15
assert {r['figure'] for r in rows} == {'SPIRE_Figure_1', 'SPIRE_Figure_2', *(f'SPIRE_S{i}' for i in range(1, 7)),
                                     'FMN2_Figure_1', *(f'FMN2_S{i}' for i in range(1, 6)), 'Convergence_Figure_1'}
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
        rows = list(csv.DictReader((root / f'source_data/{study}_RBP_{kind}_regional_scores.csv').open()))
        assert len({r['region_label'] for r in rows}) == regions
        assert all('exon' in r['region_label'].lower() or 'intron' in r['region_label'].lower() for r in rows)
        rows = list(csv.DictReader((root / f'source_data/{study}_GC_{kind}_profiles.csv').open()))
        assert len(rows) == bins * 3
        rows = list(csv.DictReader((root / f'source_data/{study}_GC_{kind}_difference.csv').open()))
        assert len(rows) == bins * 2 and all(float(r['se']) >= 0 for r in rows)
    original_go = {}
    expected_top = set()
    for ontology in ('BP', 'CC', 'MF'):
        for background in ('tested_N10', 'paper_OrgDb'):
            path = root / f'tables/S6_{study.lower()}_GO_{background}_all_{ontology}_all_terms.csv'
            rows = list(csv.DictReader(path.open()))
            expected_top.update(r['ID'] for r in sorted(rows, key=lambda r: (float(r['pvalue']), r['ID']))[:3])
            original_go.update(((r['ID'], background), r) for r in rows)
    for selection in ('top_terms', 'manuscript_themes'):
        rows = list(csv.DictReader((root / f'source_data/{study}_GO_{selection}.csv').open()))
        ids = {r['ID'] for r in rows}
        assert len(rows) == 2 * len(ids)
        expected_ids = expected_top if selection == 'top_terms' else {
            'GO:0051301', 'GO:0051321', 'GO:0007059', 'GO:0007051', 'GO:0030036',
            'GO:0005635', 'GO:0006338', 'GO:0016607', 'GO:0008380'}
        assert ids == expected_ids
        for r in rows:
            original = original_go.get((r['ID'], r['background']))
            if original:
                assert r['tested'] == 'TRUE'
                for field in ('pvalue', 'p.adjust', 'FoldEnrichment', 'Count'):
                    assert math.isclose(float(r[field]), float(original[field]), rel_tol=1e-12)
            else:
                assert r['tested'] == 'FALSE' and r['p.adjust'] == 'NA'
                assert int(r['background_term_genes']) < 10 or int(r['background_term_genes']) > 500
measurements = list(csv.DictReader((root / 'tables/S7_six_feature_measurements.csv').open()))
for study in ('SPIRE', 'FMN2'):
    features = {r['feature'] for r in measurements}
    for feature in features:
        rows = list(csv.DictReader((root / f'source_data/{study}_{feature}_display_summary.csv').open()))
        assert len(rows) == 3
        for r in rows:
            values = [float(m['value']) for m in measurements
                      if (m['dataset'], m['feature'], m['DATASET']) == (study.lower(), feature, r['DATASET'])]
            assert int(r['n']) == len(values)
            assert math.isclose(float(r['median']), statistics.median(values), rel_tol=1e-12)
            assert math.isclose(float(r['mean']), statistics.fmean(values), rel_tol=1e-12)
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
print('PASS: output integrity, sample/event labels, feature summaries, GO selections/backgrounds, RBP regions, clustering, GC and gallery links.')
