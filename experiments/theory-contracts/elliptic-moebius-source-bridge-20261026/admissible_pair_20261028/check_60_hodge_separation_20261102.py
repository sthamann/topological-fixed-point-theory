import importlib.util
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[4]
def load(name,path):
 s=importlib.util.spec_from_file_location(name,root/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
o=load('original','experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py')
h=load('hodge','experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/hodge_complement_20261028/hodge_complement_channel.py')
r=o.source_rays();U,six0,_=o.reflection_actions(r);six,_=o.outer_mark_action(six0);mapping=o.find_bell_partition_intertwiner(U,six);Cs=h.hodge_amplitudes()
M=np.stack([x.reshape(-1) for x in U],axis=1);Gc=M.conj().T@M
if np.any(Gc.imag) or not np.array_equal(Gc.real,np.rint(Gc.real)): raise ValueError('Gram must be real integral')
G=Gc.real.astype(np.int64);I=np.eye(60,dtype=np.int64)
if any(np.trace(a.T@b) != 6*(i==j) for i,a in enumerate(Cs) for j,b in enumerate(Cs)): raise ValueError('Hodge gram')
if any(not np.array_equal(x.conj().T@x,np.eye(10)) for x in U): raise ValueError('not unitary')
if np.any(G@(G-8*I)@(G-16*I)@(G-96*I)) or np.trace(G)!=600 or np.trace(G@G)!=14400 or np.trace(G@G@G)!=944640: raise ValueError('spectral certificate')
print('Gram polynomial and trace moments passed; normalized top5 sum = 4/15')
for mark in range(6):
 chart,_=h.duad_chart(mark,mapping)
 if not np.array_equal(chart.T@chart,np.eye(10)): raise ValueError('chart is not orthogonal')
 val=0
 for u in U:
  v=chart.T@u@chart
  val+=sum(int(round(abs(np.trace(c.T@v))**2)) for c in Cs)
 if val != 144: raise ValueError('chart overlap differs')
 print('mark',mark,'integer squared-overlap',val,'p',str(__import__('fractions').Fraction(val,3600)))
print('Ky Fan bound: P5 probability <= 4/15 under arbitrary unitary input/output chart; Hodge P5 probability 1')
