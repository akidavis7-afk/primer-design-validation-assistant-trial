import argparse, pandas as pd
from pathlib import Path
from src.core import read_fasta_text, analyze_pair
p=argparse.ArgumentParser(); p.add_argument('--target',required=True); p.add_argument('--primers',required=True); p.add_argument('--output',default='outputs'); a=p.parse_args()
target=read_fasta_text(Path(a.target).read_text())
df=pd.read_csv(a.primers)
rows=[]
for _,r in df.iterrows(): rows.append({'pair_id':r['pair_id'],**analyze_pair(target,r['forward'],r['reverse'])})
out=Path(a.output); out.mkdir(exist_ok=True); pd.DataFrame(rows).to_csv(out/'primer_report.csv',index=False)
print(out/'primer_report.csv')
