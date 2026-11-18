"""Three decisive exact role tests in the already given native root/source chart.

No new source, metric, coupling or state is chosen. Reuse the existing
original eightfold-source root builder; never run its writing main().
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import importlib.util
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[4]
path=ROOT/'experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/dynamic_reconstruction_20261027/check_a7_e7_clock.py'
spec=importlib.util.spec_from_file_location('existing_a7_e7_source',path)
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def run():
    roots=frozenset(tuple(int(2*x) for x in r) for r in old.e8_roots())
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    axes=((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))
    corners=[frozenset(r for r in roots if
                      r[5]*a[0]==r[6]*a[1]==r[7]*a[2]) for a in axes]
    need(all(len(c)==72 for c in corners),'four E6 corners')
    normals=[]
    for r in sorted(roots):
        if any(r[:5]) or sum(x!=0 for x in r)!=2:
            continue
        if next(x for x in r if x)!=2:
            continue
        normals.append(r)
    need(len(normals)==6,'six unoriented D3 root normals')
    edges=[]
    for i,j in combinations(range(4),2):
        need(len(corners[i]&corners[j])==40,'common D5')
        matches=[n for n in normals if all(dot(n,r)==0 for r in corners[i]|corners[j])]
        need(len(matches)==1,'unique E7 edge')
        target={r for r in roots if dot(matches[0],r)==0}
        closed=set(corners[i]|corners[j]);front=set(closed)
        while front:
            add={tuple(x+y for x,y in zip(a,b)) for a in front for b in closed}
            front=(add&roots)-closed;closed.update(front)
        need(closed==target and len(closed)==126,'two E6 generate precisely E7')
        edges.append({'vertices':[i,j],'normal':list(matches[0]),
                      'intersection_roots':40,'union_roots':104,'closure_roots':126})
    need(len({tuple(x['normal']) for x in edges})==6,'whole tetrahedron edge dictionary')

    # A root-current cut is an inclusion of a subalgebra, not a quotient map.
    # Removed root alpha has a nonzero retained Cartan projection, so
    # P[e_alpha,e_-alpha] != [Pe_alpha,Pe_-alpha]=0, with either bracket sign.
    alpha6=[F(1),0,0,0,0,F(1),0,0]
    p6=alpha6[:5]+[F(1,3)]*3
    alpha7=[F(1),F(1),0,0,0,0,0,0]
    p7=[F(x)-F(1,4) for x in alpha7]
    need(sum(x*x for x in p6)==F(4,3),'nonzero E6 projected coroot')
    need(sum(x*x for x in p7)==F(3,2),'nonzero E7 projected coroot')

    rho2=(7,5,3,1,-1,-3,-5,-7)
    neutral={r for r in roots if sum(r)==0}
    need(len(neutral)==126,'actual Q neutral E7')
    commuting={r for r in neutral if dot(rho2,r)==0}
    need(len(commuting)==8 and sp.Matrix(list(commuting)).rank()==4,'four A1 zero-mode roots')
    need(all(dot(a,b) in (0,8,-8) for a in commuting for b in commuting),'orthogonal A1 factors')
    counter=(2,-2,0,0,0,0,0,0)
    need(counter in neutral and -F(dot(rho2,counter),4)==-1,'neutral is not source-time invariant')
    for a in (F(1,4),F(3,4)):
        all_commuting={r for r in roots if
                       -F(dot(rho2,r),4)+(2*a-1)*F(sum(r),4)==0}
        need(all_commuting==commuting,'same horizontal source-time centralizer in both quarters')
    result={'scope':'given native E8 lattice and original eightfold source H8 only',
            'tetrahedron':{'E6_vertices':4,'E7_edges':edges},
            'root_cut_projection':{'E6_nonzero_projected_coroot_norm2':str(sum(x*x for x in p6)),
                                   'E7_nonzero_projected_coroot_norm2':str(sum(x*x for x in p7)),
                                   'is_Lie_homomorphism':False},
            'actual_H8':{'Q_neutral_E7_roots':126,
                         'ordinary_horizontal_commuting_roots':8,
                         'horizontal_root_type':'4A1',
                         'neutral_counterexample':'[H8,J_(e0-e1,0)]=-J_(e0-e1,0)',
                         'all_affine_commuting_modes_classified':False},
            'physical_E8_E7_E6_gauge_or_breaking_rule_derived':False}
    Path(__file__).with_name('role-check-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    run()
