"""Package an already completed publication assembly; standard library only."""
import csv
import hashlib
from html import escape
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def package(root):
    if not (root / 'COMPLETE').is_file():
        raise ValueError('Publication R assembly must complete before packaging')
    figures = list(csv.DictReader((root / 'tables/figure_manifest.csv').open()))
    if {r['figure'] for r in figures} != {'Figure_1', 'Figure_2', 'Figure_S1', 'Figure_S2', 'Figure_S3', 'Figure_S4'}:
        raise ValueError('Expected two main and four supplementary figures')
    for row in figures:
        for extension in ('pdf', 'png'):
            path = root / 'figures' / (row['figure'] + '.' + extension)
            if not path.is_file() or path.stat().st_size < 1000:
                raise ValueError(f'Missing or empty figure: {path}')
    inventory = []
    for path in sorted([*root.glob('tables/*.csv'), *root.glob('source_data/*.csv')]):
        if path.name == 'table_inventory.csv':
            continue
        with path.open(newline='') as handle:
            rows = csv.reader(handle)
            columns = next(rows)
            count = sum(1 for _ in rows)
        inventory.append(dict(file=str(path.relative_to(root)), rows=count,
                              columns=';'.join(columns), sha256=sha256(path)))
    with (root / 'tables/table_inventory.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['file', 'rows', 'columns', 'sha256'])
        writer.writeheader()
        writer.writerows(inventory)
    manuscript = [root / 'text' / (name + '.md') for name in ('RESULTS', 'METHODS', 'LEGENDS', 'TABLES')]
    if not all(path.is_file() for path in manuscript):
        raise ValueError('Incomplete manuscript text')
    pandoc = shutil.which('pandoc')
    if pandoc:
        for extension in ('docx', 'html'):
            subprocess.run([pandoc, '--standalone', '--metadata', 'title=SPIRE RNA analysis',
                            *map(str, manuscript), '-o', str(root / ('Manuscript.' + extension))], check=True)
    for name, paths in (
        ('Supplementary_tables.zip', [*root.glob('tables/*.csv'), *root.glob('source_data/*.csv'),
                                      *root.glob('text/*.md'), root / 'sessionInfo.txt']),
        ('Figures.zip', list(root.glob('figures/*'))),
    ):
        temporary = root / (name + '.tmp')
        with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for path in sorted(paths):
                archive.write(path, arcname=str(path.relative_to(root)))
        with zipfile.ZipFile(temporary) as archive:
            if archive.testzip() is not None:
                raise ValueError(f'Corrupt archive: {temporary}')
        temporary.replace(root / name)
    page = ['<!doctype html><html lang="en"><meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            '<title>SPIRE · Publication figures and analysis</title>',
            '<style>body{font:17px/1.55 system-ui;color:#222;max-width:1100px;margin:36px auto;padding:0 22px}'
            'a{color:#235887}img{width:100%;height:auto}article{border-top:1px solid #ddd;margin:40px 0;padding-top:22px}'
            '.note{background:#f3f2f5;padding:16px 22px}h1{font-size:30px}h2{font-size:23px}nav{display:flex;gap:18px;flex-wrap:wrap}</style>',
            '<h1>SPIRE loss: expression and splicing</h1>',
            '<p class="note">Three independent animals per condition. Selective RNA changes and division-associated '
            'candidates; limited shared FMN2 targeting. The RNA results do not establish a global molecular phenocopy.</p>',
            '<nav><a href="Figures.zip">Download all figures</a><a href="Supplementary_tables.zip">Supplementary tables and source data</a>'
            '<a href="text/RESULTS.md">Results</a><a href="text/METHODS.md">Methods</a><a href="text/LEGENDS.md">Legends</a></nav>']
    if pandoc:
        page.append('<p><a href="Manuscript.docx">Editable manuscript (Word)</a> · '
                    '<a href="Manuscript.html">Read the manuscript</a></p>')
    page.append('<p>Final-size layouts: 180 mm wide; embedded-font PDFs and 450 dpi PNGs. '
                'Figure numbering is local to the RNA analysis and can be changed when integrated with microscopy.</p>')
    for row in figures:
        name, title = escape(row['figure']), escape(row['title'])
        page.append(f'<article><h2>{name.replace("_", " ")} · {title}</h2>'
                    f'<p><a href="figures/{name}.pdf">PDF</a> · <a href="figures/{name}.png">PNG</a></p>'
                    f'<a href="figures/{name}.png"><img src="figures/{name}.png" alt="{title}" loading="lazy"></a></article>')
    page.extend(['<h2>Tables and provenance</h2><p><a href="text/TABLES.md">Table guide and field definitions</a> · '
                 '<a href="tables/table_inventory.csv">File inventory and SHA-256</a> · '
                 '<a href="tables/input_checksums.csv">Input checksums</a> · '
                 '<a href="sessionInfo.txt">R environment</a></p><ul>'])
    for row in inventory:
        if Path(row['file']).name.startswith(tuple(f'S{i}_' for i in range(1, 9))):
            page.append(f'<li><a href="{escape(row["file"])}">{escape(Path(row["file"]).stem)}</a> '
                        f'({row["rows"]:,} rows)</li>')
    page.append('</ul><p>Experimental methods and the SPIRE data accession still require author completion. '
                'Per-event betAS FPR and exploratory motif p-values must not be described as multiple-testing FDR.</p></html>')
    (root / 'index.html').write_text('\n'.join(page))
    outputs = sorted(p for p in root.rglob('*') if p.is_file() and p.name not in ('output_checksums.csv', 'PACKAGE_COMPLETE'))
    with (root / 'output_checksums.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['file', 'bytes', 'sha256'])
        writer.writerows((str(p.relative_to(root)), p.stat().st_size, sha256(p)) for p in outputs)
    (root / 'PACKAGE_COMPLETE').write_text('Archive integrity and complete output inventory verified.\n')
    print(f'Packaged {len(figures)} figures and {len(inventory)} data files in {root}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python3 analysis/package_publication.py PUBLICATION_DIRECTORY')
    package(Path(sys.argv[1]))
