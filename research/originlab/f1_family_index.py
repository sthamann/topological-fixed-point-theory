"""Four-mark family Toeplitz index and conjugate mirror obstruction.

The symbol is a DECLARED MODEL, never the original weighted P1/R113 source.
"""
from itertools import permutations
import sympy as sp

def projectors():
    q = sp.ones(4)/4
    return q, sp.eye(4)-q

def shift(n):
    if type(n) is not int or n < 3:
        raise ValueError("at least three Hardy levels required")
    s = sp.zeros(n)
    for j in range(n-1):
        s[j+1,j] = 1
    return s

def shift_power(n, power):
    s = shift(n)
    return s**power if power >= 0 else s.T**(-power)

def finite_toeplitz(n, common, relative):
    q,p = projectors()
    return sp.kronecker_product(q,shift_power(n,common))+sp.kronecker_product(p,shift_power(n,relative))

def index_data(common, relative):
    if type(common) is not int or type(relative) is not int:
        raise ValueError("integer windings required")
    return dict(common=common,relative=relative,det_winding=common+3*relative,
                reference_index=-common,relative_index=-3*relative,
                total_index=-common-3*relative,
                kernel=max(-common,0)+3*max(-relative,0),
                cokernel=max(common,0)+3*max(relative,0))

def paired_completion(n):
    """Exact finite unitary with two edge couplings, not physical gauge mass."""
    q,p=projectors()
    s=shift(n)
    a=sp.zeros(n);b=sp.zeros(n)
    a[0,0]=1;b[n-1,n-1]=1
    shared=sp.kronecker_product(q,sp.eye(n))
    return sp.BlockMatrix([
        [shared+sp.kronecker_product(p,s),sp.kronecker_product(p,a)],
        [sp.kronecker_product(p,b),shared+sp.kronecker_product(p,s.T)]
    ]).as_explicit()

def radial_time_covariance(n=6):
    """L0 number operator commutators, not a physical source identification."""
    q,p=projectors()
    h=sp.kronecker_product(sp.eye(4),sp.diag(*range(n)))
    replies=[]
    for k0,k3 in ((0,0),(1,0),(0,1),(-2,1)):
        d=finite_toeplitz(n,k0,k3)
        expected=(k0*sp.kronecker_product(q,shift_power(n,k0))
                  +k3*sp.kronecker_product(p,shift_power(n,k3)))
        assert h*d-d*h==expected
        replies.append(dict(k0=k0,k3=k3,commutes_with_L0=expected==sp.zeros(4*n)))
    assert [v["commutes_with_L0"] for v in replies]==[True,False,False,False]
    return replies

def finite_mark_blindness():
    """Two S4-equivariant unit-determinant-winding loops agreeing at all four marks."""
    q,p=projectors()
    z,f=sp.symbols("z f",nonzero=True)
    reference_only=z*q+p
    relative_with_counterflow=z*f**(-3)*q+f*p
    assert sp.simplify(reference_only.det()-z)==0
    assert sp.simplify(relative_with_counterflow.det()-z)==0
    assert relative_with_counterflow.subs(f,1)==reference_only
    assert relative_with_counterflow.subs(f,-1)==-reference_only
    # f(e^(i theta))=1 outside (pi/8,3pi/8), and
    # exp(i*(8 theta-pi)) inside this arc: continuous, winding one.
    marks=[sp.Integer(0),sp.pi/2,sp.pi,3*sp.pi/2]
    assert all(not (sp.pi/8 < angle < 3*sp.pi/8) for angle in marks)
    assert sp.simplify(sp.exp(sp.I*(8*sp.pi/4-sp.pi)))==-1
    left=index_data(1,0)
    right=index_data(-2,1)
    assert left["det_winding"]==right["det_winding"]==1
    assert left["relative_index"]==0 and right["relative_index"]==-3
    return dict(status="EXACT_FINITE_MARK_NONIDENTIFIABILITY",
                same_four_mark_matrices=True,same_full_determinant_function="z",
                separated_relative_indices=[0,-3],
                explicit_unmeasured_angle="pi/4",
                required_original_observation="A source-certified continuous four-mark family loop and its sector-resolved winding, not four samples or determinant alone")


def controls():
    q,p=projectors()
    mark_nonidentifiability=finite_mark_blindness()
    radial_covariance=radial_time_covariance()
    assert q*q==q and p*p==p and q*p==sp.zeros(4)
    assert q.rank()==1 and p.rank()==3
    for perm in permutations(range(4)):
        m=sp.zeros(4)
        for j,i in enumerate(perm):m[i,j]=1
        assert m*q==q*m and m*p==p*m
    z=sp.symbols("z",nonzero=True)
    cases={}
    for k0,k3 in [(0,0),(1,0),(0,1),(-2,1),(2,-1)]:
        d=index_data(k0,k3)
        assert sp.simplify((z**k0*q+z**k3*p).det()-z**(k0+3*k3))==0
        assert d["kernel"]-d["cokernel"]==d["total_index"]
        assert d["total_index"]+index_data(-k0,-k3)["total_index"]==0
        cases[f"{k0},{k3}"]=d
    assert cases["0,1"]["det_winding"]==3
    assert cases["-2,1"]["det_winding"]==1
    assert cases["-2,1"]["relative_index"]==-3
    assert cases["-2,1"]["reference_index"]==2
    finite=[]
    for n in (3,4,6):
        t=finite_toeplitz(n,-2,1)
        assert 4*n-t.rank()==5 and 4*n-t.T.rank()==5
        v=paired_completion(n)
        assert v.T*v==sp.eye(8*n) and v*v.T==sp.eye(8*n)
        finite.append(dict(N=n,zero_modes_each_end=5,unitary=True))
    choices=[dict(k0=1-3*k3,k3=k3,variation=abs(1-3*k3)+3*abs(k3))
             for k3 in range(-3,4)]
    assert min(choices,key=lambda d:d["variation"])["k3"]==0
    return dict(status="EXACT_FINITE_FAMILY_INDEX_OBSTRUCTION",
        s4_projector_ranks=[1,3],cases=cases,finite_controls=finite,
        unit_determinant_winding_with_three_relative=dict(k0=-2,k3=1),
        primitive_winding_alone_selects_no_family=True,
        conjugate_pair_total_index=0,radial_time_covariance=radial_covariance,
        finite_mark_nonidentifiability=mark_nonidentifiability,
        physical_source_identity_proved=False,
        original_r113_values=None,physical_mirror_decoupling_proved=False,
        next_original_observation="Same-source weighted P1 four-mark loop, determinant winding, conjugate boundary domains, positive time, whole selector and absolute R113 reference")

if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2))
