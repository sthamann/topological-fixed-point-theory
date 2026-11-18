"""Original affine reflected Higgs response, with and without a cut."""
from itertools import product
from pathlib import Path
import json
import sys
from functools import lru_cache

import numpy as np
import sympy as s

repo = Path(sys.argv[1])
contract = repo/'experiments/theory-contracts/recursive-mirrored-source-20261103'
sys.path.insert(0, str(contract))
from check_higgs_gram_complement_20261111 import (adjoint_sign, add, charged_slots,
                                                  dot, epsilon, family_vertex,
                                                  neg, vertex_response)
from e8_boundary_response import current_response


def tail(*coordinates):
    return (0,)*5+coordinates


def hermitian(root):
    return ((root,1),(neg(root),adjoint_sign(root)))


A = hermitian(tail(1,-1,0)) + hermitian(tail(-1,-1,0))
B = hermitian(tail(0,1,-1))
I = (((0,)*8,1),)
LEFT = s.Matrix([-2,1,0])/s.sqrt(5)
RIGHT = s.eye(3)[:,0]


def finite_word(p,b,h,a,q,t):
    # Keep only actual insertions, in the same radial and reflection order.
    marked=[(neg(b),1/t),(h,s.S.One),(a,t),(q,s.S.Zero)]
    marked=[(root,z) for root,z in marked if any(root)]
    roots=(neg(p),)+tuple(root for root,z in marked)
    z=(s.S.One,)+tuple(z/(1+z) for root,z in marked)
    pref=adjoint_sign(p)*s.prod((1+z)**-2 for root,z in marked)
    if any(b):
        pref*=adjoint_sign(b)/t**2
    return roots,z,pref


def direct(p,b,h,a,q,t):
    roots,z,pref=finite_word(p,b,h,a,q,t)
    return pref*vertex_response(roots,z)


def staged(p,b,h,a,q,t):
    # First cut fixes total lattice momentum. The factors below are the
    # exact all-mode Heisenberg contractions at the two intermediate cuts.
    if add(p,b)!=add(add(q,a),h):
        return s.S.Zero
    return (epsilon(b,p)*epsilon(a,q)*epsilon(h,add(q,a))*t**(dot(b,p)+dot(a,q))
            *(1-t)**(dot(h,a)-dot(b,h))*(1-t*t)**(-dot(b,a)))


def coefficient(slot,h,bra,ket,t):
    left,right=slot
    out=s.zeros(4)
    for i,j in product(range(4),repeat=2):
        p,p_sign=family_vertex(left,i)
        q,q_sign=family_vertex(right,j)
        p=neg(p)
        for (b,bc),(a,ac) in product(bra,ket):
            d=direct(p,b,h,a,q,t)
            u=staged(p,b,h,a,q,t)
            assert d==u,(p,b,h,a,q,d,u)
            out[i,j]+=p_sign*q_sign*adjoint_sign(p)*bc*ac*d
    return out


def zero_mode(slot,h,bra,ket):
    left,right=slot
    out=s.zeros(4)
    for i,j in product(range(4),repeat=2):
        p,p_sign=family_vertex(left,i)
        q,q_sign=family_vertex(right,j)
        p=neg(p)
        for (b,bc),(a,ac) in product(bra,ket):
            x,y=add(p,b),add(q,a)
            if x==add(y,h) and dot(x,x)==dot(y,y)==2:
                out[i,j]+=p_sign*q_sign*adjoint_sign(p)*bc*ac*epsilon(b,p)*epsilon(a,q)*epsilon(h,y)
    return out


@lru_cache(maxsize=None)
def same_leg(root,term,conjugate,t):
    states=[family_vertex(root,i) for i in range(3)]
    out=s.zeros(3)
    for i,j in product(range(3),repeat=2):
        p,ps=states[i];q,qs=states[j]
        if conjugate:
            ps*=adjoint_sign(p);qs*=adjoint_sign(q);p,q=neg(p),neg(q)
        for (b,bc),(a,ac) in product(term,term):
            d=direct(p,b,(0,)*8,a,q,t)
            assert d==staged(p,b,(0,)*8,a,q,t)
            out[i,j]+=ps*qs*s.conjugate(bc)*ac*d
    return out


def normalize(g):
    m=np.array(g,dtype=complex)
    assert np.allclose(m,m.conj().T)
    eig,v=np.linalg.eigh(m)
    assert eig.min()>0
    return (v*eig**(-.5))@v.conj().T


def branch(t):
    slots=charged_slots()
    roots=[(0,0,0,0,-1)+tuple(sign if k==j else 0 for k in range(3))
           for j in range(3) for sign in (-1,1)]
    rows={}
    matrices={}
    exact={}
    lengths={'e':(8,5,3),'u':(7,3,0),'d':(7,5,2)}
    seed=float(1/(6*s.pi)+3/(256*s.pi**4))
    base=(seed*(1-seed))**.5
    for tag,rows_by_color in slots.items():
        values={}
        matrices[tag]={}
        exact[tag]={}
        for name,bra,ket in [('empty',I,I),('AA',A,A),('AB',A,B),('BA',B,A),('BB',B,B)]:
            null_values=[]
            exact[tag][name]=[]
            for h in roots:
                phase=adjoint_sign(h) if tag=='u' else 1
                h=neg(h) if tag=='u' else h
                outputs=[phase*coefficient(slot,h,bra,ket,t) for slot in rows_by_color]
                assert all(m in (outputs[0],-outputs[0]) for m in outputs)
                exact[tag][name].append(outputs[0][:3,:3])
                null_values.append(s.simplify((LEFT.T*outputs[0][:3,:3]*RIGHT)[0]))
                if len(null_values)==1:
                    test=outputs[0][:3,:3]
                    left,right=rows_by_color[0]
                    l,r=(normalize(same_leg(left,bra,True,t)),normalize(same_leg(right,ket,False,t)))
                    word=np.diag(base**(np.array(lengths[tag])/2))
                    sv=np.sort(np.linalg.svd(word@l@np.array(test,dtype=float)@r@word,compute_uv=False))
                    zm=phase*zero_mode(rows_by_color[0],h,bra,ket)
                    matrices[tag][name]={'singulars':sv.tolist(),'exact_rank_before_norm':test.rank(),
                                         'grade_one_only_null_overlap':str(s.simplify((LEFT.T*zm[:3,:3]*RIGHT)[0]))}
            values[name]=[str(x) for x in null_values]
        rows[tag]=values
    ratios={}
    for name in ('AB','BA'):
        e,u,d=(np.array(matrices[tag][name]['singulars']) for tag in ('e','u','d'))
        if min(e[0],u[0],d[0])>1e-13:
            ratios[name]={'quarks':(u/d).tolist(),'leptons':(e[:-1]/e[1:]).tolist(),'e_over_u':float(e[0]/u[0])}
    z,x1,x3,x5=s.symbols('z h1 h3 h5')
    # Compare the two orders coherently, first at the actual seed root,
    # then for the three available original cubic Higgs directions.
    coherent={}
    for tag in ('e','d','u'):
        mixed=[a+z*b for a,b in zip(exact[tag]['AB'],exact[tag]['BA'])]
        triplet=x1*mixed[0]+x3*mixed[2]+x5*mixed[4]
        coherent[tag]={'original_h1_det':str(s.factor(mixed[0].det())),
                       'cubic_triplet_det':str(s.factor(triplet.det())),
                       'det_at_fixed_equal_triplet_and_orders':str(triplet.det().subs({z:1,x1:1,x3:1,x5:1}))}
        assert triplet.det().subs(z,1).expand()==0
    fixed_profiles={}
    for phase_left,phase_right in product((s.I,-s.I),repeat=2):
        norm_l=A+tuple((root,phase_left*c) for root,c in B)
        norm_r=A+tuple((root,phase_right*c) for root,c in B)
        singulars={};ranks={}
        for tag in ('e','d','u'):
            p=exact[tag]
            def trip(name):
                return p[name][0]+p[name][2]+p[name][4]
            total=(trip('AA')+phase_right*trip('AB')+s.conjugate(phase_left)*trip('BA')
                   +s.conjugate(phase_left)*phase_right*trip('BB'))
            left,right=slots[tag][0]
            l,r=(normalize(same_leg(left,norm_l,True,t)),normalize(same_leg(right,norm_r,False,t)))
            word=np.diag(base**(np.array(lengths[tag])/2))
            ranks[tag]=total.rank()
            singulars[tag]=np.sort(np.linalg.svd(word@l@np.array(total,dtype=complex)@r@word,compute_uv=False))
        val={'ranks':ranks}
        if min(ranks.values())==3:
            e,u,d=(singulars[tag] for tag in ('e','u','d'))
            val.update(quarks=(u/d).tolist(),leptons=(e[:-1]/e[1:]).tolist(),e_over_u=float(e[0]/u[0]))
            assert not np.allclose(u/d,[55/117,34/(47*seed),3/(26*seed**2)],rtol=.01)
        fixed_profiles[f'{phase_left},{phase_right}']=val
    return {'null_overlaps_all_six_higgs_roots':rows,'original_h1_gram_normalized':matrices,
            'h1_full_rank_ratios':ratios,'coherent_order_tests':coherent,'fixed_quarter_phase_controls':fixed_profiles}


def ward_probe():
    # Same nonzero AB term from the reflected original H1 matrix element.
    left,right=charged_slots()['e'][0]
    p=neg(family_vertex(left,1)[0]); q=family_vertex(right,0)[0]
    h=(0,0,0,0,-1,-1,0,0)
    t=s.Rational(4,5)
    p_sign=family_vertex(left,1)[1]*adjoint_sign(p)
    terms=[]
    for (b,bc),(a,ac) in product(A,B):
        roots,z,pref=finite_word(p,b,h,a,q,t)
        if vertex_response(roots,z)!=0:
            terms.append((roots,z,pref*bc*ac*p_sign))
    assert terms
    y=(s.Rational(-1,3),)*3+(s.Rational(1,2),)*2+(0,)*3
    w,u=s.Rational(7,6),s.Rational(6,5)
    baseline=sum(pref*vertex_response(roots,z) for roots,z,pref in terms)
    assert baseline != 0
    pole=lambda x:sum(dot(y,r)/(x-z0) for r,z0 in zip(terms[0][0],terms[0][1]))
    one=sum(pref*current_response(roots,z,((y,u),)) for roots,z,pref in terms)/baseline
    two=sum(pref*current_response(roots,z,((y,u),(y,w))) for roots,z,pref in terms)/baseline
    assert one==pole(u) and s.factor((two-pole(u)*pole(w))*(u-w)**2)==s.Rational(5,6)
    # Extract the finite Sugawara term by subtracting the exact double pole
    # and taking the coincident regular part in each Cartan direction.
    stress=s.S.Zero
    for roots,z,pref in terms:
        vacuum=vertex_response(roots,z)
        for j in range(8):
            v=tuple(int(k==j) for k in range(8))
            mixed=current_response(roots,z,((v,u),(v,w)))/vacuum-1/(u-w)**2
            pu=sum(r[j]/(u-x) for r,x in zip(roots,z))
            pw=sum(r[j]/(w-x) for r,x in zip(roots,z))
            assert mixed==pu*pw
            stress+=pref*vacuum*pu**2/(2*baseline)
    assert stress!=0
    return {'connected_hypercharge_level':'5/6','affine_stress_over_same_vertex':str(stress),
            'seam_word_absolute_vertex_in_chosen_frame':str(baseline)}


def fuchs_alternative():
    a,b,z=s.symbols('a b z',real=True)
    T=s.Matrix([[0,1,0],[0,2,0],[0,2,1]])
    F=s.Matrix([[0,0,a],[0,0,b],[0,0,0]])
    N=-s.I*s.Matrix([[0,0,b],[0,0,a],[0,0,0]])
    L=s.Matrix([[1,0],[2,0],[0,s.sqrt(5)]])/s.sqrt(5)
    R=s.eye(3)[:,1:]
    out={}
    for name,mark in [('F',F),('N',N)]:
        B=-(mark+mark.T)
        D=T+z*B
        assert (LEFT.T*B*RIGHT)[0]==0
        complement=L.T*D*R
        direct=s.factor(D.det())
        cut=s.factor((LEFT.T*D*RIGHT-LEFT.T*D*R*complement.inv()*L.T*D*RIGHT)[0])
        frame=LEFT.row_join(L)
        assert s.simplify(direct-frame.det()*complement.det()*cut)==0
        assert s.expand(direct).coeff(z,2)!=0
        out[name]={'direct_det':str(direct),'staged_null_schur':str(cut)}
    return out


if __name__=='__main__':
    output={str(t):branch(t) for t in (s.Rational(4,5),s.Rational(3,4))}
    output['same_word_ward_probe']=ward_probe()
    output['fuchs_alternative']=fuchs_alternative()
    print(json.dumps(output,indent=2))
