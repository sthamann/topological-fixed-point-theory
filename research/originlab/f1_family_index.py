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
    # Along |z|=|f|=1, these inverses are adjoints and both
    # reflected positive Grams are I despite different windings.
    inverse_a=z**(-1)*q+p
    inverse_b=z**(-1)*f**3*q+f**(-1)*p
    assert (inverse_a*reference_only-sp.eye(4)).applyfunc(sp.simplify)==sp.zeros(4)
    assert (inverse_b*relative_with_counterflow-sp.eye(4)).applyfunc(sp.simplify)==sp.zeros(4)
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
                same_positive_reflected_gram="I4",
                separated_relative_indices=[0,-3],
                explicit_unmeasured_angle="pi/4",
                required_original_observation="A source-certified continuous four-mark family loop and its sector-resolved winding, not four samples or determinant alone")


def analytic_one_sheet_obstruction():
    """Conditional theorem: common analytic disk polarization plus unit determinant winding.

    A holomorphic scalar block with no boundary zeros has winding equal
    to its number of disk zeros. This premise is NOT implied by the
    previously established positive P1 Euclidean time semigroup.
    """
    q,p=projectors()
    z=sp.symbols("z",nonzero=True)
    # Sectorwise analytic functions have nonnegative zero counts.
    solutions=[(n0,n3) for n0 in range(5) for n3 in range(5)
               if n0+3*n3==1]
    assert solutions==[(1,0)]
    # One explicit analytic inner solution: degree one only on Q.
    a=sp.Rational(1,3)
    b=(z-a)/(1-a*z)
    inner=b*q+p
    assert sp.factor(inner.det()-b)==0
    assert sp.factor(b*b.subs(z,1/z)-1)==0
    # The proposed three-relative-family winding has a second-order
    # pole on the common channel and cannot be analytic in the disk.
    mixed=z**(-2)*q+z*p
    assert sp.factor(mixed.det()-z)==0
    assert (z*z*mixed).applyfunc(lambda x: sp.limit(x,z,0))==q
    # Full continuous SU(4) mixing does not preserve the rank-three P.
    h=sp.diag(1,-1,0,0)
    assert p*h-h*p!=sp.zeros(4)
    return dict(status="CONDITIONAL_EXACT_ONE_SHEET_NO_GO",
                assumptions=["fixed S4 action throughout the loop",
                             "both scalar family blocks holomorphic in one disk",
                             "boundary invertibility",
                             "primitive P1 winding equals determinant winding one"],
                allowed_zero_counts=dict(common=1,relative=0),
                relative_family_index_from_one_sheet=0,
                mixed_family_symbol_requires_reference_pole_order_at_least=2,
                full_SU4_mixing_breaks_relative_projector=True,
                original_physical_analyticity_proved=False)


def based_su4_unwinding():
    """Based contraction keeping all four marked values and det fixed.

    Mixes Q with relative channels at intermediate stages, so it does
    NOT preserve the fixed S4 commutant. This is the physical caveat.
    """
    f,c,s=sp.symbols("f c s",nonzero=True)
    h=sp.Matrix([[c/f,s],[-s,c*f]])
    at_one_inverse=sp.Matrix([[c,-s],[s,c]])
    b=(h*at_one_inverse).applyfunc(sp.expand)
    norm=c*c+s*s
    assert sp.factor(b.det()-norm*norm)==0
    bstar=b.T.applyfunc(lambda x: x.subs(f,1/f))
    assert (bstar*b-norm*norm*sp.eye(2)).applyfunc(sp.simplify)==sp.zeros(2)
    assert b.subs({c:1,s:0})==sp.diag(1/f,f)
    assert b.subs({c:0,s:1})==sp.eye(2)
    assert b.subs(f,1)==norm*sp.eye(2)

    def pair(matrix,j):
        out=sp.eye(4)
        for x,i in enumerate((0,j)):
            for y,k in enumerate((0,j)):
                out[i,k]=matrix[x,y]
        return out

    product=sp.eye(4)
    for j in (1,2,3):
        product=product*pair(b,j)
    assert product.subs({c:1,s:0})==sp.diag(f**(-3),f,f,f)
    assert product.subs({c:0,s:1})==sp.eye(4)
    assert product.subs(f,1)[0,0]==norm**3
    # The existing witness f equals 1 at 1,i,-1,-i, and winds once.
    return dict(status="EXACT_SU4_BASED_UNWINDING",
                same_four_mark_values_throughout=True,
                same_full_determinant_throughout=True,
                same_unitary_gram_throughout=True,
                destroys_fixed_S4_block_equivariance=True,
                result="relative winding 3 not protected if common/relative mixing is allowed",
                physical_mixing_admissibility_proved=False)


def certified_winding_sampling():
    """Finite samples can certify a winding only with a genuine source bound.

    If |u'|<=L, min_j |u(theta_j)|>=m and max angular gap h
    satisfies L*h<m, each true arc and its polygon chord lie in
    the convex nonzero disk D(u(theta_j),L*h). Thus they are homotopic.
    No L or m is assumed known for the original weighted R113 source.
    """
    assert sp.pi/4 < 1     # u(z)=z, eight equal samples, L=m=1
    assert 3*sp.pi/10 < 1  # u(z)=z^3, twenty samples, L=3, m=1
    assert not (sp.pi/2 < 1)  # four mark-only samples do not meet this bound
    return dict(status="CONDITIONAL_CERTIFIED_WINDING_SAMPLING",
                sufficient_condition="L * maximum_angular_gap < minimum_sample_modulus_lower_bound",
                test_examples=[dict(function="z",samples=8,certified=True),
                               dict(function="z^3",samples=20,certified=True),
                               dict(function="z",samples=4,certified=False)],
                source_derivative_bound_supplied=False,
                original_winding_certified=False)


def controls():
    q,p=projectors()
    mark_nonidentifiability=finite_mark_blindness()
    radial_covariance=radial_time_covariance()
    analytic_obstruction=analytic_one_sheet_obstruction()
    sector_unwinding=based_su4_unwinding()
    sampling=certified_winding_sampling()
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
        analytic_one_sheet_obstruction=analytic_obstruction,
        sector_unwinding_without_S4=sector_unwinding,
        certified_winding_sampling=sampling,
        physical_source_identity_proved=False,
        original_r113_values=None,physical_mirror_decoupling_proved=False,
        next_original_observation="Same-source weighted P1 four-mark loop, determinant winding, conjugate boundary domains, positive time, whole selector and absolute R113 reference")

if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2))
