#!/usr/bin/env python3
"""Download published mm10 annotations and stream the union of tested events."""
import csv
import fcntl
import gzip
import hashlib
import json
import pathlib
import re
import sys
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
REF = ROOT / 'data/reference'
URLS = {
    'EVENT_INFO-mm10.tab.gz': 'https://vastdb.crg.eu/downloads/mm10/EVENT_INFO-mm10.tab.gz',
    'PROT_IMPACT-mm10-v3.tab.gz': 'https://vastdb.crg.eu/downloads/mm10/PROT_IMPACT-mm10-v3.tab.gz',
    'gencode.vM10.annotation.gtf.gz': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz',
    'GRCm38.primary_assembly.genome.fa.gz': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/GRCm38.primary_assembly.genome.fa.gz',
}
EXPECTED_SHA256 = {
    'EVENT_INFO-mm10.tab.gz': 'b5e43aeaadfd7a870c9b1dda46dd7e591501bd20f72560f727fcd122b31d4cf4',
    'PROT_IMPACT-mm10-v3.tab.gz': '4f8166c3e6a3c602a3ec947fb3b4415f8e47159e9b2c8d350d7ecf08fb7c668c',
    'gencode.vM10.annotation.gtf.gz': 'a47ccb2cd2a6648c16a016d36f9091581c57c5c19a889e8bfd47d87617bf634a',
    'GRCm38.primary_assembly.genome.fa.gz': '09de3c56747c69948d5ca39125caa5bebeac41e0de2e4a67bac947bfc80864f0',
}

def file_hash(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()

def gc(seq):
    seq = seq.upper()
    if seq in ('NA', '-', ''):
        return ''
    return 100 * (seq.count('G') + seq.count('C')) / len(seq) if seq else ''

def parse_coord(value):
    m = re.fullmatch(r'([^:]+):(\d+)-(\d+)', value)
    if not m:
        return '', '', ''
    return m[1], int(m[2]), int(m[3])

def main():
    if len(sys.argv) < 2:
        raise SystemExit('Usage: python analysis/prepare_references.py RESULTS1/tables/all_events.csv [RESULTS2/tables/all_events.csv ...]')
    REF.mkdir(parents=True, exist_ok=True)
    lock = (REF / '.prepare.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    sources = []
    for name, url in URLS.items():
        path = REF / name
        if not path.exists():
            print('Downloading', url, flush=True)
            temporary = path.with_suffix(path.suffix + '.part')
            with urllib.request.urlopen(url, timeout=120) as response, temporary.open('wb') as f:
                while chunk := response.read(1024 * 1024):
                    f.write(chunk)
            temporary.replace(path)
        h = file_hash(path)
        if h != EXPECTED_SHA256[name]:
            raise RuntimeError('Reference checksum changed: ' + name)
        sources.append(dict(file=str(path.relative_to(ROOT)), url=url, sha256=h))
        if name.endswith(('.gtf.gz', '.fa.gz')):
            target = path.with_suffix('')
            if not target.exists():
                print('Decompressing', name, flush=True)
                with gzip.open(path, 'rb') as f, target.with_suffix(target.suffix + '.part').open('wb') as dest:
                    while chunk := f.read(1024 * 1024):
                        dest.write(chunk)
                target.with_suffix(target.suffix + '.part').replace(target)
    (REF / 'sources.json').write_text(json.dumps(sources, indent=2) + '\n')
    csv.field_size_limit(sys.maxsize)
    events = set()
    for table in sys.argv[1:]:
        with open(table, newline='') as f:
            events.update(r['EVENT'] for r in csv.DictReader(f))
    genes = {}
    with gzip.open(REF / 'gencode.vM10.annotation.gtf.gz', 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] == 'gene':
                attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
                genes.setdefault(attrs['gene_name'], attrs['gene_id'])
    columns = ['EVENT', 'GENE', 'COORD_o', 'REF_CO', 'CO_C1', 'CO_A', 'CO_C2', 'LE_o', 'LE_n',
               'SCAFFOLD', 'START', 'END', 'STRAND', 'GENEID_ENSEMBL', 'UPSTRM_EX_BORDER',
               'DOSTRM_EX_BORDER', 'TARGET_LENGTH', 'TARGET_GC', 'C1_GC', 'C2_GC']
    found = set()
    temporary = REF / 'event_annotations.tsv.part'
    with gzip.open(REF / 'EVENT_INFO-mm10.tab.gz', 'rt') as f, temporary.open('w', newline='') as out:
        writer = csv.DictWriter(out, columns, delimiter='\t')
        writer.writeheader()
        for r in csv.DictReader(f, delimiter='\t'):
            if r['EVENT'] not in events:
                continue
            assert r['EVENT'] not in found, 'Duplicate reference event: ' + r['EVENT']
            found.add(r['EVENT'])
            scaffold, start, end = parse_coord(r['CO_A'])
            c1, c2 = parse_coord(r['CO_C1']), parse_coord(r['CO_C2'])
            strand = r['REF_CO'][-1:]
            # Matt borders are genomic left/right even on the minus strand.
            left, right = (c1[2], c2[1]) if strand == '+' else (c2[2], c1[1])
            writer.writerow({**{k: r[k] for k in columns if k in r},
                'SCAFFOLD': scaffold, 'START': start, 'END': end, 'STRAND': strand,
                'GENEID_ENSEMBL': genes.get(r['GENE'], ''),
                'UPSTRM_EX_BORDER': left, 'DOSTRM_EX_BORDER': right,
                'TARGET_LENGTH': len(r['Seq_A']) if r['Seq_A'] not in ('NA', '-', '') else '', 'TARGET_GC': gc(r['Seq_A']),
                'C1_GC': gc(r['Seq_C1']), 'C2_GC': gc(r['Seq_C2'])})
    with (REF / 'missing_reference_events.txt').open('w') as f:
        f.write('\n'.join(sorted(events - found)))
    temporary.replace(REF / 'event_annotations.tsv')
    print('Annotated', len(found), 'tested events;', len(events - found), 'missing IDs recorded separately', flush=True)

if __name__ == '__main__':
    assert parse_coord('chrX:10-20') == ('chrX', 10, 20)
    assert parse_coord('NA') == ('', '', '')
    assert gc('ACGT') == 50
    main()
