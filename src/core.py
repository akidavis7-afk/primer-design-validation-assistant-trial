from __future__ import annotations
import re
from dataclasses import dataclass, asdict

COMP=str.maketrans('ACGTacgt','TGCAtgca')
def rc(s:str)->str: return s.translate(COMP)[::-1].upper()
def gc(s): return round(100*(s.upper().count('G')+s.upper().count('C'))/max(len(s),1),2)
def tm(s):
    s=s.upper(); return 2*(s.count('A')+s.count('T'))+4*(s.count('G')+s.count('C'))
def longest_complement(a,b):
    b=rc(b); best=0
    for shift in range(-len(b),len(a)+1):
        run=0
        for i,ch in enumerate(a.upper()):
            j=i-shift
            if 0<=j<len(b) and ch==b[j]: run+=1; best=max(best,run)
            else: run=0
    return best

def read_fasta_text(text): return ''.join(x.strip() for x in text.splitlines() if not x.startswith('>')).upper()

def analyze_pair(target, forward, reverse):
    target=target.upper(); forward=forward.upper(); reverse=reverse.upper()
    fpos=target.find(forward); rbind=rc(reverse); rpos=target.find(rbind)
    amp=None
    if fpos>=0 and rpos>=0 and rpos>fpos: amp=rpos+len(rbind)-fpos
    warnings=[]
    for label,s in [('forward',forward),('reverse',reverse)]:
        if not 18<=len(s)<=35: warnings.append(f'{label}: unusual length')
        if not 35<=gc(s)<=65: warnings.append(f'{label}: GC% outside 35-65')
        if re.search(r'(A{5,}|C{5,}|G{5,}|T{5,})',s): warnings.append(f'{label}: homopolymer >=5')
        if longest_complement(s,s)>=5: warnings.append(f'{label}: possible hairpin/self-dimer')
    if longest_complement(forward,reverse)>=5: warnings.append('pair: possible heterodimer')
    if amp is None: warnings.append('pair: binding sites not found in valid orientation')
    return {'forward_length':len(forward),'forward_gc':gc(forward),'forward_tm':tm(forward),
            'reverse_length':len(reverse),'reverse_gc':gc(reverse),'reverse_tm':tm(reverse),
            'forward_position':fpos+1 if fpos>=0 else None,'reverse_position':rpos+1 if rpos>=0 else None,
            'amplicon_size':amp,'warning_count':len(warnings),'warnings':' | '.join(warnings)}
