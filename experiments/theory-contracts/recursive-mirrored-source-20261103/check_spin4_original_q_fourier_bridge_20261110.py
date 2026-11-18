"""Embed the actual variable-H Q Fourier matrix and test transverse mixing."""
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np

base = Path(__file__).resolve().parent
repo = base.parents[2]
with contextlib.redirect_stdout(io.StringIO()):
    s=runpy.run_path(str(base/'check_aps_higgs_normal_jet_20261110.py'))
runtime=repo/'experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime'
sys.path.insert(0,str(runtime))
from assemble import Assembler

assembler=Assembler(runtime.parent/'domain/native-diagnostic-v1.json')
item=assembler.assemble({'schema':'tfpt.field-configuration.v1',
                        'domain_sha256':assembler.domain['_sha256'],'genes':[2,1,1]})
g = item['gamma96']
inner_sign = g.diagonal().real
I=np.eye(len(g))
gam=s['gamma']; g5=s['gamma5']; gn=gam[3]; Jm=1j*gn@g5
Ax, Ay=(1j*gn@ga for ga in gam[:2])
modes=item['modes']; alpha=item['alpha']; phi=item['phi_fourier']
zero=np.zeros_like(g)
phi_full=np.block([[phi.get(a-b,zero) for b in modes] for a in modes])
ps=np.diag([a+alpha for a in modes])
diag=np.kron(ps,g)
h_plus=item['h_N']
assert all(e<1e-14 for e in item['kinematic_residuals'].values())
g_full=np.kron(np.eye(len(modes)),g)
h_minus=-g_full@h_plus@g_full
B=np.kron(Ax,np.kron(ps,I))+np.kron(Jm,phi_full)
errors=[]
for c in (1,-1):
    for b,h in ((1,h_plus),(-1,h_minus)):
        P=(np.eye(4)+c*g5)@(np.eye(4)+b*Ax)/4
        u=P[:,int(np.argmax(np.diag(P).real))]
        u/=np.linalg.norm(u)
        v=-1j*Jm@u
        W=np.vstack([np.diag(np.tile(np.where(inner_sign>0,u[k],v[k]),len(modes)))
                     for k in range(4)])
        errors.append(float(np.max(abs(B@W-W@h))))
assert max(errors)<2e-13

# Apply the same actual first harmonic to two coupled, non-collinear momenta.
px0, px1 = alpha, alpha+1
transverse = .2
A0=(px0*Ax+transverse*Ay)/np.hypot(px0,transverse)
A1=(px1*Ax+transverse*Ay)/np.hypot(px1,transverse)
Q0,Q1=np.kron(A0,g),np.kron(A1,g)
M=np.kron(Jm,phi[1])
idfull=np.eye(len(Q0))
leak=(idfull-Q1)@M@(idfull+Q0)/4
leak_norm=float(np.linalg.norm(leak,ord=2))
cosine=(px0*px1+transverse**2)/(np.hypot(px0,transverse)*np.hypot(px1,transverse))
predicted=np.linalg.norm(phi[1],ord=2)*np.sqrt((1-cosine)/2)
assert abs(leak_norm-predicted)<1e-13 and leak_norm>1e-3
print(json.dumps({'actual_domain':item['domain_id'],'modes':modes,
    'original_candidate_id':item['candidate_id'],'bridge_max_error':max(errors),
    'transverse_tangential_component':transverse,'sector_mixing_operator_norm':leak_norm,
    'exact_formula_numerically_evaluated':float(predicted)},indent=2))
