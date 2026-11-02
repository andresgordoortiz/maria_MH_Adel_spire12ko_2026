"""Render the current figure galleries from the checked publication manifest."""
import csv
from html import escape
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GROUPS = {
    '01_SPIRE': 'SPIRE corregido', '02_FMN2': 'FMN2',
    '03_COMPARACION': 'SPIRE frente a FMN2', '04_RBP': 'Motivos RBP',
    '05_CONTROLES': 'Muestras y control de calidad',
}
TITLES = {
    'Figure2_spire_fastq_splicing_panels': 'SPIRE · Figura 2 completa: splicing y Matt',
    'Figure2_fmn2_splicing_panels': 'FMN2 · Figura 2 completa: splicing y Matt',
    'Fig3A_protein_impact': 'Impacto funcional previsto de los eventos',
    'DGE_volcano': 'SPIRE · Expresión diferencial',
    'Target_gene_expression': 'SPIRE · Expresión de Spire1, Spire2 y Fmn2',
    'Comparison_overlap_counts': 'Genes y eventos compartidos',
    'Comparison_delta_PSI': 'Dirección y magnitud de los cambios de splicing',
    'Comparison_shared_event_sample_PSI': 'PSI de cada muestra en los eventos compartidos',
    'Comparison_Matt_feature_signatures': 'Características de secuencia compartidas · Matt',
    'Comparison_event_rates': 'Comparación dentro del universo de eventos común',
    'Fig2H_EX_RBP_heatmap': '2H · Motivos RBP en exones',
    'FigS6A_EX_RBP_tracks': 'S6A · Perfiles regionales de motivos en exones',
    'FigS6B_EX_RBP_boxplots': 'S6B · Motivos en exones por dirección del cambio',
    'FigS6C_IR_RBP_tracks': 'S6C · Perfiles regionales de motivos en intrones',
    'FigS6D_IR_RBP_heatmap': 'S6D · Motivos RBP en intrones',
    'FigS6E_IR_RBP_boxplots': 'S6E · Motivos en intrones por dirección del cambio',
    'Expression_PCA': 'SPIRE · PCA de expresión',
    'Splicing_PCA': 'SPIRE · PCA de splicing',
    'Target_exon_support': 'SPIRE · Soporte de exones en los loci diana',
    'Sample_correspondence': 'Correspondencia entre etiquetas antiguas y FASTQ',
    'Old_new_splicing': 'Comparación diagnóstica: análisis anterior y corregido',
}
STYLE = '''body{font:17px/1.5 system-ui;color:#222;max-width:1200px;margin:32px auto;padding:0 20px}
a{color:#174e93}nav{display:flex;flex-wrap:wrap;gap:18px;margin:24px 0}.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(290px,1fr));gap:20px}
article{border:1px solid #ddd;border-radius:6px;padding:14px}h3{font-size:17px;margin:0 0 12px}img{width:100%;height:260px;object-fit:contain}
section{margin:40px 0}.note{padding:15px;background:#f2f5f8}details{border-top:1px solid #ddd;padding:16px 0}summary{cursor:pointer;font-weight:600;margin-bottom:16px}footer{border-top:1px solid #ddd;margin-top:40px;padding-top:20px}'''


def url(path, page):
    return escape(os.path.relpath(ROOT / path, page.parent), quote=True)


def card(row, page):
    title = escape(TITLES.get(row['figure'], row['figure'].replace('_', ' ')))
    png, pdf = url(row['png'], page), url(row['pdf'], page)
    return f'<article><h3>{title}</h3><a href="{png}"><img src="{png}" alt="{title}" loading="lazy"></a><p><a href="{pdf}">PDF</a> · <a href="{png}">PNG a tamaño completo</a></p></article>'


def page_start(title):
    return f'<!doctype html><html lang="es"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{escape(title)}</title><style>{STYLE}</style><h1>{escape(title)}</h1>'


def main():
    rows = list(csv.DictReader((ROOT / 'results/figures_manifest.csv').open()))
    assert rows and len({r['png'] for r in rows}) == len(rows)
    for row in rows:
        assert row['group'] in GROUPS
        for kind in ('png', 'pdf'):
            path = (ROOT / row[kind]).resolve()
            assert path.is_relative_to((ROOT / 'results').resolve()) and path.is_file(), path
    principal = [r for r in rows if r['section'] == 'principales']
    page = ROOT / 'results/index.html'
    parts = [page_start('Figuras actuales · SPIRE y FMN2'),
        '<p class="note"><strong>Versión correcta: SPIRE Cont_1/2/7 frente a DKO_4/5/6.</strong> Todas las figuras de esta carpeta pertenecen al análisis actual. Las versiones anteriores están archivadas fuera de results.</p>',
        '<p>Empieza por las figuras principales de abajo. Cada imagen abre el PNG completo; el enlace PDF sirve para guardar o compartir. Los paneles individuales y las complementarias están en la galería de cada grupo.</p>',
        '<nav>' + ''.join(f'<a href="#{g}">{escape(t)}</a>' for g, t in GROUPS.items()) + '</nav>']
    if (ROOT / 'results/publication/PACKAGE_COMPLETE').is_file():
        parts.insert(1, '<p class="note"><a href="publication/index.html"><strong>Figuras finales para el manuscrito</strong></a>: dos principales, cuatro suplementarias, métodos, resultados y tablas S1–S8.</p>')
    for group, title in GROUPS.items():
        selected = [r for r in principal if r['group'] == group]
        if group == '05_CONTROLES':
            selected = [r for r in rows if r['figure'] in ('Expression_PCA','Splicing_PCA','Target_exon_support')]
        parts.append(f'<section id="{group}"><h2>{escape(title)}</h2><p><a href="{group}/index.html">Ver todas las figuras de este grupo</a></p><div class="gallery">')
        parts.extend(card(r, page) for r in selected)
        parts.append('</div></section>')
    parts.append('<footer><p>GO y motivos RBP: resultados exploratorios; consultar los límites de interpretación.</p><p><a href="../docs/CONCLUSIONS.md">Conclusiones finales</a> · <a href="../reports/index.html">Informes completos y tablas</a> · <a href="figures_manifest.csv">Listado de figuras actuales</a></p></footer></html>')
    page.write_text('\n'.join(parts))
    for group, title in GROUPS.items():
        page = ROOT / 'results' / group / 'index.html'
        selected = [r for r in rows if r['group'] == group]
        parts = [page_start(title), '<p><a href="../index.html">Volver a las figuras principales</a></p><p>Todas estas figuras son actuales. Los paneles individuales desglosan las figuras completas; las complementarias amplían el análisis.</p>']
        sections = sorted({r['section'] for r in selected}, key=lambda x: (x != 'principales', x))
        for section in sections:
            figures = [r for r in selected if r['section'] == section]
            title2 = section.replace('_', ' ').replace('/', ' · ').capitalize()
            block = '<div class="gallery">' + ''.join(card(r, page) for r in figures) + '</div>'
            if section == 'principales':
                parts.append('<h2>Principales</h2>' + block)
            else:
                parts.append(f'<details><summary>{escape(title2)} ({len(figures)} figuras actuales)</summary>{block}</details>')
        parts.append('</html>')
        page.write_text('\n'.join(parts))
    print(f'Published {len(principal)} principal figures and {len(rows)} total current PNG/PDF pairs.')


if __name__ == '__main__':
    main()
