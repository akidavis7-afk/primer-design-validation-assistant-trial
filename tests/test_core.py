from src.core import gc,tm,rc,analyze_pair

def test_metrics():
    assert rc('ATGC')=='GCAT'
    assert gc('GGCC')==100
    assert tm('ATGC')==12

def test_pair():
    target='AAAATGCGTACGTACGTACGTTTTTACGATCGATCGATCGATCG'
    r=analyze_pair(target,'ATGCGTACGTACGTACGT','CGATCGATCGATCGATCG')
    assert r['forward_position']==4
