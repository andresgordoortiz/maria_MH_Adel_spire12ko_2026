"""Read-only checks of PSI/count pairing and exon-position support bias."""
from pathlib import Path
from collections import defaultdict
import csv,gzip,math,re,statistics,json,hashlib
out=Path('results/spire_identity_check/tables')
def write(name,rows):
 assert rows
 with (out/f'{name}.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def number(x):
 try:return float(x)
 except ValueError:return float('nan')
def coverage(q):
 p=q.rsplit('@',1)[-1].split(',')
 if len(p)!=2:return [float('nan')]*2
 return [number(x) for x in p]
assert coverage('N,N@3.25,6.75')==[3.25,6.75]
assert math.isnan(number('NA'))
genes={};ambiguous=set()
with gzip.open('data/reference/gencode.vM10.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  p=line.rstrip('\n').split('\t')
  if p[2]!='gene':continue
  m=re.search('gene_name "([^"]+)"',p[8])
  if not m:continue
  name=m[1];value=(p[0],p[6])
  if name in genes and genes[name]!=value:ambiguous.add(name)
  genes[name]=value
for name in ambiguous:genes.pop(name,None)
exons=defaultdict(dict);targets=[];bad_examples=[]
raw=Path('final_input_tables/spiredko_INCLUSION_LEVELS_FULL-mm10.tab')
with raw.open() as f:
 reader=csv.DictReader(f,delimiter='\t');samples=[s for s in reader.fieldnames[6:] if not s.endswith('-Q')]
 stat=[dict(sample=s,finite_PSI_and_counts_N10=0,PSI_count_mismatch_gt_006pp=0,max_PSI_count_difference_pp=0.0,total_corrected_support_overlapping_events=0.0,finite_count_cells=0) for s in samples]
 assert len(samples)==6
 for r in reader:
  totals=[]
  for j,s in enumerate(samples):
   inc,exc=coverage(r[s+'-Q']);total=inc+exc;psi=number(r[s]);totals.append(total)
   if math.isfinite(total) and inc>=0 and exc>=0:
    stat[j]['total_corrected_support_overlapping_events']+=total;stat[j]['finite_count_cells']+=1
   if math.isfinite(psi) and math.isfinite(total) and total>=10 and inc>=0 and exc>=0:
    diff=abs(psi-100*inc/total)
    stat[j]['finite_PSI_and_counts_N10']+=1
    stat[j]['max_PSI_count_difference_pp']=max(stat[j]['max_PSI_count_difference_pp'],diff)
    if diff>.06:
     stat[j]['PSI_count_mismatch_gt_006pp']+=1
     if len(bad_examples)<30:bad_examples.append(dict(sample=s,EVENT=r['EVENT'],PSI=psi,counts_PSI=100*inc/total,difference_pp=diff))
   if r['GENE'] in ['Spire1','Spire2','Fmn2'] and r['COMPLEX'] in ['S','C1','C2','C3','MIC']:
    targets.append(dict(GENE=r['GENE'],EVENT=r['EVENT'],COORD=r['COORD'],COMPLEX=r['COMPLEX'],sample=s,PSI=r[s],corrected_total=total))
  if r['COMPLEX']!='S' or r['GENE'] not in genes:continue
  m=re.fullmatch(r'([^:]+):(\d+)-(\d+)',r['COORD'])
  if not m or m[1]!=genes[r['GENE']][0]:continue
  # Identical exon coordinates are represented once; annotation modules are excluded.
  exons[r['GENE']].setdefault(r['COORD'],(int(m[2]),totals))
write('PSI_count_pairing_QC',stat)
if bad_examples:write('PSI_count_pairing_discrepancy_examples',bad_examples)
write('target_exon_support',targets)
position=[]
for gene,coordinates in exons.items():
 if gene in ['Spire1','Spire2'] or len(coordinates)<6:continue
 rows=sorted(coordinates.values(),reverse=genes[gene][1]=='-')
 anchor=[statistics.median(row[1][j] for row in rows if math.isfinite(row[1][j])) for j in [0,1,3,4]] if all(any(math.isfinite(row[1][j]) for row in rows) for j in [0,1,3,4]) else []
 if len(anchor)!=4 or min(anchor)<20:continue
 for j,s in enumerate(samples):
  five=[row[1][j] for row in rows[:3] if math.isfinite(row[1][j])]
  three=[row[1][j] for row in rows[-3:] if math.isfinite(row[1][j])]
  if len(five)!=3 or len(three)!=3:continue
  a,b=statistics.median(five),statistics.median(three)
  position.append(dict(GENE=gene,sample=s,exon_events=len(rows),five_prime_median_support=a,three_prime_median_support=b,log2_five_to_three_support=math.log2((a+.5)/(b+.5))))
write('exon_position_support_by_gene',position)
summary=[]
for s in samples:
 vals=[r['log2_five_to_three_support'] for r in position if r['sample']==s]
 summary.append(dict(sample=s,genes=len(vals),median_log2_five_to_three_support=statistics.median(vals)))
write('exon_position_support_summary',summary)
print(json.dumps(dict(pairing=stat,position_bias=summary),indent=2))
