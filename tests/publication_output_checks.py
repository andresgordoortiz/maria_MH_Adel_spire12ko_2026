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
assert len(rows) == 16
assert {r['figure'] for r in rows} == {'SPIRE_Figure_1', 'SPIRE_Figure_2', *(f'SPIRE_S{i}' for i in range(1, 8)),
                                     'FMN2_Figure_1', *(f'FMN2_S{i}' for i in range(1, 6)), 'Convergence_Figure_1'}
rows = list(csv.DictReader((root / 'source_data/SPIRE_division_event_annotations.csv').open()))
assert len({r['EVENT'] for r in rows}) == 127 and len({r['GENE'] for r in rows}) == 79
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
    revised = [r for r in csv.DictReader((root / 'tables/S6_GO_revised_all_tests.csv').open())
               if r['dataset'] == study.lower() and r['selection'] == 'selected']
    expected_top = set()
    for ontology in ('BP', 'CC', 'MF'):
        for background in ('tested_N10', 'paper_OrgDb'):
            candidates = [r for r in revised if (r['ontology'], r['background']) == (ontology, background)]
            expected_top.update(r['ID'] for r in sorted(candidates, key=lambda r: (float(r['pvalue']), r['ID']))[:3])
    lookup = {(r['ID'], r['background']): r for r in revised}
    for selection in ('top_terms', 'manuscript_themes'):
        rows = list(csv.DictReader((root / f'source_data/{study}_GO_{selection}.csv').open()))
        ids = {r['ID'] for r in rows}
        assert len(rows) == 2 * len(ids)
        expected_ids = expected_top if selection == 'top_terms' else {
            'GO:0051301', 'GO:0051321', 'GO:0007059', 'GO:0007051', 'GO:0030036',
            'GO:0005635', 'GO:0006338', 'GO:0016607', 'GO:0008380'}
        assert ids == expected_ids
        for r in rows:
            original = lookup[r['ID'], r['background']]
            for field in ('pvalue', 'q_BH', 'FoldEnrichment', 'Count'):
                assert math.isclose(float(r[field]), float(original[field]), rel_tol=1e-12)
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
# Independently reproduce every primary BH family from its nominal P values.
def check_bh(rows, pfield, qfield):
    rows = [r for r in rows if r[pfield] not in ('NA', '') and r[qfield] not in ('NA', '')]
    order = sorted(range(len(rows)), key=lambda i: float(rows[i][pfield]))
    running = 1.0
    for rank in range(len(rows), 0, -1):
        r = rows[order[rank-1]]
        running = min(running, float(r[pfield])*len(rows)/rank)
        assert math.isclose(float(r[qfield]), running, rel_tol=1e-10, abs_tol=1e-14), r


go = list(csv.DictReader((root / 'tables/S6_GO_revised_all_tests.csv').open()))
for study in ('spire', 'fmn2'):
    families = []
    for background in ('tested_N10', 'paper_OrgDb'):
        for selection in ('selected', 'consistent', 'separated'):
            family = [r for r in go if (r['dataset'], r['background'], r['selection']) == (study, background, selection)]
            check_bh(family, 'pvalue', 'q_BH')
            families.append({r['ID'] for r in family})
    assert all(f == families[0] for f in families)
    assert all(float(r['q_BH']) > .05 for r in go if (r['dataset'], r['background'], r['selection']) == (study, 'tested_N10', 'selected'))
    motifs = [r for r in csv.DictReader((root / 'tables/S8_RBP_regional_presence_exact_BH.csv').open()) if r['dataset'] == study]
    assert len(motifs) == 5278
    check_bh(motifs, 'pvalue', 'q_BH')
    hits = {(r['RBP'], r['event_type'], r['region'], r['direction']) for r in motifs if float(r['q_BH']) <= .05}
    assert hits == ({('PCBP2', 'EX', 'R1', 'up'), ('PPRC1', 'Alt5', 'R5', 'up')} if study == 'spire' else set())
    matt = [r for r in csv.DictReader((root / 'tables/S7_Matt_all_feature_tests.csv').open()) if r['dataset'] == study]
    check_bh(matt, 'p.value', 'p.BH')
for name, pfield, qfield in [('S2_SPIRE_all_gene_expression', 'pvalue', 'padj'),
                            ('S6_GO_expression_all_tests', 'pvalue', 'q_BH'),
                            ('S6_GO_event_count_matched_themes', 'pvalue', 'q_BH'),
                            ('S5_overlap_inference_BH', 'p.value', 'q_BH')]:
    check_bh(list(csv.DictReader((root / f'tables/{name}.csv').open())), pfield, qfield)
mapping = list(csv.DictReader((root / 'tables/S6_GO_resolved_symbol_mapping.csv').open()))
assert len({r['SYMBOL'] for r in mapping}) == len(mapping)
assert {r['method'] for r in mapping} == {'direct symbol', 'unique alias'}
aliases = {}
for r in csv.DictReader((root / 'tables/S6_GO_alias_candidates.csv').open()):
    aliases.setdefault(r['SYMBOL'], set()).add(r['ENTREZID'])
assert all(aliases[r['SYMBOL']] == {r['ENTREZID']} for r in mapping if r['method'] == 'unique alias')
union, = csv.DictReader((root / 'tables/S6_spire_division_union_counts.csv').open())
assert int(union['observed']) == 79 and int(union['background_annotated']) == 780
assert math.isclose(float(union['expected_uniform']), 76.3624496899815, rel_tol=1e-12)
main = list(csv.DictReader((root / 'source_data/SPIRE_main_GO.csv').open()))
assert len(main) == 7 and all(r['background'] == 'tested_N10' and r['selection'] == 'selected' for r in main)
assert all(float(r['q_BH']) > .05 for r in main)
opportunities = [r for r in csv.DictReader((root / 'tables/S6_GO_gene_event_opportunities.csv').open()) if r['dataset'] == 'spire']
for r in csv.DictReader((root / 'source_data/SPIRE_GO_event_opportunity_rates.csv').open()):
    genes = [g for g in opportunities if min(math.floor(math.log2(int(g['tested_events']))), 7) == int(r['stratum'])]
    assert int(r['genes']) == len(genes)
    assert int(r['selected_genes']) == sum(g['selected'] == 'TRUE' for g in genes)
    assert math.isclose(float(r['fraction']), int(r['selected_genes'])/len(genes), rel_tol=1e-12)
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
print('PASS: output integrity, sample/event labels, feature summaries, all BH families, GO mappings/backgrounds, exact motifs, clustering, GC and gallery links.')
