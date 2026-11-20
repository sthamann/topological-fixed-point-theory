#!/usr/bin/env python3
"""Check syndrome incidence symmetry's induced physical symplectic lift."""
import sys,json
from pathlib import Path
from qec_origin import word_pauli,syndrome,span,symplectic
from syndrome_spread import invertible_binary_maps,derived_subgroup

rows=[word_pauli(word) for word in ('XZZXI','IXZZX','XIXZZ','ZXIXZ')]
S=span(rows)
lines=[];error_of={}
for i in range(5):
 es=(1<<i,1<<(i+5),(1<<i)|(1<<(i+5)))
 lines.append(frozenset(syndrome(e,rows,5) for e in es))
 for e in es:error_of[syndrome(e,rows,5)]=e
G={g for g in invertible_binary_maps() if {frozenset(g[x] for x in line) for line in lines}==set(lines)}
core=derived_subgroup(derived_subgroup(G))
logical_X=31;logical_Z=31<<5
images={};code_preserved=0
for g in G:
 columns=[error_of[g[syndrome(1<<i,rows,5)]] for i in range(10)]
 def M(v):
  out=0
  for i in range(10):
   if (v>>i)&1:out^=columns[i]
  return out
 assert all(symplectic(columns[i],columns[j],5)==symplectic(1<<i,1<<j,5) for i in range(10) for j in range(10))
 assert {M(v) for v in S}==S;code_preserved+=1
 # coefficients [logical X, logical Z] modulo stabilizer, columns MX,MZ
 images[g]=tuple((symplectic(M(v),logical_Z,5),symplectic(M(v),logical_X,5)) for v in (logical_X,logical_Z))
assert len({images[g] for g in G})==6
assert {images[g] for g in core}=={((1,0),(0,1))}
print(json.dumps({'spread_automorphisms':len(G),'physical_symplectic_lifts':code_preserved,
 'full_logical_symplectic_image_order':len(set(images.values())),
 'perfect_core_order':len(core),'perfect_core_logical_symplectic_image':list({images[g] for g in core}),
 'interpretation':'A5 core acts trivially on logical Pauli labels. Syndrome symmetry is not a faithful icosahedral action on the logical qubit. Clifford phase lifts and Pauli corrections are separate and cannot create a faithful projective A5 inside the logical Clifford group.'},indent=2))
