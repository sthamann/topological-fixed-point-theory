"""Quantum readout of the seven marks, and blind reconstruction from native TFPT."""
import argparse
from pathlib import Path
import sys
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"fortsetzung"))
from originalvergleich import load_bridge
from apartment import twice_projector
from itertools import combinations, permutations

def native(root):
    bridge=load_bridge(root)
    src,_=bridge.inherited(root)
    q=bridge.finite_data(src)["q"]
    gens=bridge.generators()
    channel=s.zeros(16)
    gamma=[]
    for v in src.W16:
        word=bridge.monomial(v,gens)
        col=word.vec()/2
        if any(v) and not q[v]:
            gamma.append(word)
        eigen=s.Rational(2,3) if q[v] else s.Rational(1,3)
        channel+=(eigen if any(v) else 1)*col*col.H
    return channel,gamma

def matrix_units(d):
    return [s.eye(d)[:,i]*s.eye(d)[j,:] for j in range(d) for i in range(d)]

def reconstruct(root):
    phi,gamma=native(root)
    D=s.eye(4).vec()*s.eye(4).vec().T/4
    reflection=3*s.eye(16)+4*D-6*phi
    units=matrix_units(4)
    conj=lambda A:(reflection*A.vec()).reshape(4,4).T
    # Linear intertwiner equations S(E)B = B E^T for every matrix unit.
    constraints=s.Matrix.vstack(*(s.kronecker_product(s.eye(4),conj(E))
                    -s.kronecker_product(E,s.eye(4)) for E in units))
    kernel=constraints.nullspace()
    assert len(kernel)==1
    B=kernel[0].reshape(4,4).T
    assert B.T==-B and B*B.H==s.eye(4) and B*B.conjugate()==-s.eye(4)
    assert all(conj(u*v)==conj(v)*conj(u) for u in units for v in units)
    assert all(conj(u)==B*u.T*B.H for u in units)
    print("Native reconstructed form, modulo scalar:",B.tolist())
    print("Intertwiner dimension: 1; reversed products: 256/256; antiunitary square -I")
    forms=[B]+[s.I*g*B for g in gamma]
    assert len(forms)==6 and all(j.T==-j for j in forms)
    pf=lambda m:m[0,1]*m[2,3]-m[0,2]*m[1,3]+m[0,3]*m[1,2]
    for i,U in enumerate(forms):
        assert pf(U)==pf(B)
        for j,V in enumerate(forms):
            assert U*V.H+V*U.H==2*int(i==j)*s.eye(4)
            if i!=j:
                assert s.expand(pf(U+V)-2*pf(B))==0
    # A rotated unit normal yields the corresponding spinor transfer;
    # all coefficients live in Q(i,sqrt(35)) for this double-chart probe.
    new=(-B+s.sqrt(35)*forms[1])/6
    new_reflection=s.Matrix.hstack(*[(new*E.T*new.H).vec() for E in units])
    original_rankfive=(s.eye(16)+reflection)/2-D
    new_rankfive=(s.eye(16)+new_reflection)/2-D
    assert s.simplify(s.trace(original_rankfive*new_rankfive))==s.Rational(10,9)
    print("Six native skew forms: full Clifford identities; Pfaffian quadratic form is sum of squares")
    print("Two-spinor-chart overlap from native matrices:",s.simplify(s.trace(original_rankfive*new_rankfive)))
    return phi

def averages():
    # Coordinates on V=1^perp in R^7 use b_i=e_i-e_7, i<7.
    frame=s.eye(7)[:,:6]-s.eye(7)[:,6]*s.ones(1,6)
    metric=frame.T*frame
    verts=[s.eye(7)[:,i]-s.ones(7,1)/7 for i in range(7)]
    pairs=list(combinations(range(6),2))
    def induced_on_planes(P):
        # d exterior^2(P) has eigenvalue one on n wedge n^perp when
        # P projects on a line and zero on exterior^2(n^perp).
        def entry(i,j,k,l):
            return (int(j==l)*P[i,k]-int(i==l)*P[j,k]
                   +int(i==k)*P[j,l]-int(j==k)*P[i,l])
        return s.Matrix([[entry(i,j,k,l) for k,l in pairs] for i,j in pairs])
    projs=[]
    for i,vec in enumerate(verts):
        line=vec[:6,:]*vec[:6,:].T*metric/(vec.dot(vec))
        tang=s.eye(6)-line
        fam=s.Matrix(twice_projector(i+1))[1:,1:]*frame/2
        fam=fam[:6,:]
        assert fam.rank()==3 and tang*fam==fam and (tang-fam).rank()==2
        proj=induced_on_planes(line)
        assert proj*proj==proj and proj.rank()==5
        projs.append(proj)
    assert sum(projs,s.zeros(15))==s.Rational(7,3)*s.eye(15)
    assert {s.trace(p*q) for i,p in enumerate(projs) for q in projs[i+1:]}=={s.Rational(10,9)}
    assert s.Matrix.hstack(*(p.vec() for p in projs)).rank()==7
    chans=[(2*s.eye(15)-p)/3 for p in projs]
    avg=sum(chans,s.zeros(15))/7
    assert avg==s.Rational(5,9)*s.eye(15)
    responses={}
    for name,p in [('fresh',s.Rational(1,7)),('shared',s.Rational(1,4)),('frozen',s.S.One)]:
        second=sum((V*U*(p/7 if i==j else (1-p)/42)
                    for i,U in enumerate(chans) for j,V in enumerate(chans)),s.zeros(15))
        assert second==second[0,0]*s.eye(15)
        responses[name]=second[0,0]
    assert list(responses.values())==[s.Rational(25,81),s.Rational(101,324),s.Rational(1,3)]
    # Propagate the target flag as classical memory without retaining labels
    # at the output. Channel multiplication uses chronological order.
    from flaggen import transfers
    kernels=transfers()
    blocks=[s.eye(15) for _ in range(21)]
    for tick in range(1,5):
        new=[s.zeros(15) for _ in blocks]
        for u,adj in kernels.items():
            for j in range(21):
                prev=[i for i in range(21) if adj[i,j]]
                if prev:
                    new[j]+=chans[u-1]*sum((blocks[i] for i in prev),s.zeros(15))/8
        blocks=new
        mixed=sum(blocks,s.zeros(15))/21
        spectrum=mixed.eigenvals()
        print("Unrecorded shared-flag tick",tick, "exact spectrum",spectrum)
        if tick==3:
            low=s.Rational(24427,139968)
            high=s.Rational(195479,1119744)
            assert spectrum=={low:7,high:8}
            sector=(high*s.eye(15)-mixed)/(high-low)
            assert sector*sector==sector
            compare_planes(frame,pairs,sector,high-low)
    print("Seven rank-five channel projectors: span dimension 7; distinct trace overlaps 10/9")
    print("Traceless responses: one tick",avg[0,0],"; two ticks",responses)

def compare_planes(frame,pairs,sector,gap):
    from apartment import LINES
    reps={}
    for p in permutations(range(7)):
        table=frozenset(tuple(sorted(p[k-1] for k in line)) for line in LINES)
        reps.setdefault(table,p)
    def wedge(p):
        perm=s.eye(7)[:,list(p)]
        R=(perm*frame)[:6,:]
        return s.Matrix([[R[i,k]*R[j,l]-R[i,l]*R[j,k] for k,l in pairs] for i,j in pairs])
    planes=[]
    for p in reps.values():
        pinv=tuple(p.index(k) for k in range(7))
        planes.append(wedge(p)*sector*wedge(pinv))
    assert len({tuple(mat) for mat in planes})==len(reps)==30
    overlaps={s.trace(sector*mat) for mat in planes}
    closest=max(x for x in overlaps if x!=7)
    print("Different full three-tick maps on 30 labeled incidence tables:",len(planes))
    print("Rank-seven overlap traces:",sorted(overlaps),"; gap:",gap)
    print("Closest maps, squared HS superoperator distance:",gap**2*(14-2*closest))

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--repo",type=Path,default=Path(__file__).resolve().parents[4])
    reconstruct(parser.parse_args().repo)
    averages()
