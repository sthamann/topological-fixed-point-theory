"""Original signed family selector, Pin normal modes and signed compiler closure."""
from contextlib import redirect_stdout
from itertools import product
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    gate = runpy.run_path(str(base / "check_tfpt_gate_a_family_closure_20261110.py"))
    cocycle_path = repo / "experiments/theory-contracts/canonical-seam-fermionization-20261025/check_same_source_affine_defect.py"
    original_cocycle = runpy.run_path(str(cocycle_path))["eps_ext"] if cocycle_path.is_file() else None

sp, q, Q, vertex = (gate[k] for k in ("sp","q","Q","vertex"))
signed = gate["signed"]
I3 = sp.eye(3)
ladder = sp.sqrt(2)*sp.Matrix([[0,0,0],[1,0,0],[0,1,0]])
spin = sp.diag(-1,1,-1)
direct = ladder**2/2
outer = sp.diag(1,0,1)
assert spin*ladder*spin==-ladder and spin*direct*spin==direct
assert direct[2,0]==1 and direct.rank()==1
old = {}
norms = json.loads((base/"ramond-full-return-results-20261110.json").read_text())
packet_norms = gate["old"]["m"].code()[1]


def charge_label(momentum):
    half = int(not momentum[0].is_Integer)
    integers = [v-sp.Rational(half,2) for v in momentum]
    assert all(v.is_Integer for v in integers)
    return tuple(map(int,integers)),half


def phase(a,b):
    # Fixed ordered CAR section E_(n+c*s)=E_n X^c in the original
    # canonical-seam-fermionization contract. Both legs are E8 charges.
    n,c = charge_label(a)
    m,d = charge_label(b)
    exponent = sum(n[i]*m[j] for i in range(8) for j in range(i))
    value = sp.I**(c*sum(m))*(-1)**(exponent%2)
    if c and d:
        carry = sum(i*(n[i]+m[i]) for i in range(8))
        value *= -(-1)**(carry%2)
    if original_cocycle is not None:
        assert complex(value)==original_cocycle((n,c),(m,d))
    return value
for alpha, rows in signed["controls"].items():
    quarter = sp.Rational(alpha)
    rows_by_sign = {}
    for sign in (-1,1):
        lower,higher = [r for r in rows if r["sign"]==sign]
        p = [q+sp.Matrix(list(map(sp.sympify,r["charge_old"]))) for r in (lower,higher)]
        assert [int(sp.exp(2*sp.pi*sp.I*x[0])) for x in p]==[-1,-1]
        assert p[1]==p[0]-2*q
        midpoint = (p[0]+p[1])/2
        top, bottom = sorted(p,key=lambda x:-q.dot(x))
        assert int(sp.exp(2*sp.pi*sp.I*midpoint[0]))==1
        first = vertex({():Q.one},-q,top)
        second = vertex(first,-q,midpoint)
        assert first=={((1,4),):-Q.one}
        assert second in ({():2*Q.one},{():-2*Q.one})
        # Both images have the same original H and gauge transport. The
        # phases in this directed pairing use the original minus-2q mark.
        i = 2 if sign==-1 else 1
        g = sp.Matrix(norms[alpha]["full_return_Gram"])
        lo = sp.sympify(lower["exact_ground_amplitudes_up_to_old_cocycle"][i])
        hi = sp.sympify(higher["exact_ground_amplitudes_up_to_old_cocycle"][0])
        energy = 32+2*sign*quarter
        lower_charge,higher_charge = (v-q for v in p)
        direct_phase = phase(-2*q,top)
        direct_word_phase = sp.conjugate(phase(higher_charge,q))*direct_phase*phase(lower_charge,q)
        assert direct_word_phase==1
        ladder_phase = phase(-q,top)*phase(-q,midpoint)*Q.to_sympy(second[()])/2
        assert ladder_phase==direct_phase
        ratio = sp.simplify(ladder_phase/direct_phase)
        normalization = sp.sqrt(g[0,0]*g[i,i]*packet_norms[0]*packet_norms[i])
        amplitude = direct_word_phase*sp.conjugate(hi)*lo*sp.exp(-energy/8)/normalization
        rows_by_sign[str(sign)] = {"packet_pair":[0,i],
                                  "direct_old_mark":[float(sp.re(amplitude)),float(sp.im(amplitude))],
                                  "two_steps_with_retained_middle":[float(sp.re(ratio*amplitude)),float(sp.im(ratio*amplitude))],
                                  "two_step_phase_relative_to_direct":str(ratio),
                                  "two_steps_with_old_R_only_filter":[0.,0.]}
    old[alpha] = rows_by_sign

# Within the three-state space any orthogonal projector retaining the
# two ends is diagonal in this basis. The actual Q_Theta could act on
# a larger space and mix the middle with orthogonal companion modes.
selector_cases = []
for bits in product((0,1),repeat=3):
    p = sp.diag(*bits)
    if p*outer!=outer:
        continue
    assert p*direct*p==direct
    step_defect = (p*ladder**2*p - (p*ladder*p)**2)/2
    assert step_defect==(1-bits[1])*direct
    first_loss = p*ladder.H*(I3-p)*ladder*p
    assert first_loss==sp.diag(2*(1-bits[1]),0,0)
    selector_cases.append({"retained_diagonal":list(bits),
                           "direct_mark":1,
                           "two_step_mark":int(((p*ladder*p)**2)[2,0]/2),
                           "lost_first_step_norm_squared":int(first_loss[0,0])})
assert len(selector_cases)==2
s = sp.Symbol("s",real=True,nonnegative=True)
mid_mix = sp.sqrt(s*(1-s))
P_big = sp.Matrix([[1,0,0,0],[0,s,0,mid_mix],[0,0,1,0],[0,mid_mix,0,1-s]])
L_big = sp.diag(ladder,0)
assert (P_big**2-P_big).applyfunc(sp.simplify)==sp.zeros(4)
assert sp.simplify((L_big*P_big*L_big)[2,0]/2)==s
assert sp.simplify(((2*P_big-sp.eye(4))[1,1]+1)/2)==s
# This is the original family rotation, not a postulated physical clock.
Jx = (ladder+ladder.H)/2
assert Jx**3==Jx
angle = sp.Symbol("theta",real=True)
rotate = lambda theta: I3-sp.I*sp.sin(theta)*Jx+(sp.cos(theta)-1)*Jx**2
quarter_rotation = outer*rotate(sp.pi/2)*outer
twice_projected = quarter_rotation**2
whole_half = outer*rotate(sp.pi)*outer
assert twice_projected==quarter_rotation
assert (whole_half-twice_projected).rank()==1
assert rotate(2*sp.pi)==I3 and quarter_rotation.rank()==1

# Original four labels (c,b): native Higgs sends (-c,-b), linear Green
# reflection sends (-c,b). The only reflection-even/Higgs-odd sign is b.
labels = tuple(product((1,-1),(1,-1)))  # first entry is c/i
profiles = []
for values in product((1,-1),repeat=4):
    chi = dict(zip(labels,values,strict=True))
    if all(chi[-c,b]==chi[c,b] and chi[-c,-b]==-chi[c,b] for c,b in labels):
        profiles.append(list(values))
assert profiles==[[b for c,b in labels],[-b for c,b in labels]]
sigma1 = sp.Matrix([[0,1],[1,0]])
sigma2 = sp.Matrix([[0,-sp.I],[sp.I,0]])
sigma3 = sp.diag(1,-1)
assert sigma1*sigma2*sigma1==-sigma2
assert sigma1*sigma3*sigma1==-sigma3
assert sigma3*sigma2 != sigma2*sigma3

# In the original paired Pin chart the b matrices are tau_i and
# K=I tensor H_b - tau2 tensor sigma1 rho Phi. The b-swap i*tau2
# commutes with both displayed terms for the supplied native blocks.
pin_swap = sp.I*sigma2
assert pin_swap**2==-sp.eye(2) and pin_swap.conjugate()==pin_swap
assert pin_swap*sigma3==-sigma3*pin_swap and pin_swap*sigma2==sigma2*pin_swap
internal = gate["native"]
phi_h = internal["inner"](np.array([0.,1.]))
phi_m = internal["majorana"]
pair_generator_on_normal_zero = -np.kron(np.array(sigma2,dtype=complex),phi_m+.55*phi_h)
physical_pin_swap = np.kron(np.array(pin_swap,dtype=complex),np.eye(96))
b_sign = np.kron(np.diag([1,-1]),np.eye(96))
assert np.linalg.norm(physical_pin_swap@pair_generator_on_normal_zero
                      -pair_generator_on_normal_zero@physical_pin_swap)<1e-14
assert np.linalg.norm(b_sign@pair_generator_on_normal_zero-pair_generator_on_normal_zero@b_sign)>1

# Can one repair the independent b-sign with a constant native internal
# grading C that is J-real, preserves Gamma and the SM representation?
# Both Higgs and Majorana terms flip b. In particular C must anticommute
# with Majorana on the isolated gauge-neutral nu_R and its J partner.
rep = internal["native"]
v = [rep.IDX[("nu","R",0,k)] for k in range(3)]
mr = phi_m[np.ix_(v,[k+rep.NP for k in v])]
assert np.all(mr.imag==0) and np.array_equal(mr,mr.T)
m_exact = sp.Matrix([[sp.Rational(float(x.real)) for x in row] for row in mr])
square = m_exact**2
char = square.charpoly().as_poly()
assert m_exact.det()!=0 and char.gcd(char.diff()).degree()==0
# A Hermitian C block obeying C M+M conjugate(C)=0 commutes with M^2.
# The simple spectrum of real M^2 makes C real diagonal in its real
# orthonormal eigenbasis. Invertibility of M then forces C=0 there.


# Original NS normal domain H_b=i sigma3 d_w, w in [-1,1].
# The signed normal coordinate is the deck-odd section. It preserves
# the two boundary lines but does not commute with the original H_b.
n = sp.Symbol("n",integer=True,nonzero=True)
w = sp.Symbol("w",real=True)
coeff = sp.simplify(sp.integrate(w*sp.cos(n*sp.pi*(w+1)/2),(w,-1,1))/2)
assert sp.simplify(coeff-2*((-1)**n-1)/(sp.pi*n)**2)==0
norm_square = sp.simplify(2*16/sp.pi**4*(1-sp.Rational(1,16))*sp.zeta(4))
energy_square = sp.simplify(2*4/sp.pi**2*(1-sp.Rational(1,4))*sp.zeta(2))
assert (norm_square,energy_square)==(sp.Rational(1,3),1)
x, weights = leggauss(96)
for k in range(-9,10):
    numerical = np.dot(weights,x*np.cos(k*np.pi*(x+1)/2))/2
    expected = float(coeff.subs(n,k)) if k else 0.
    np.testing.assert_allclose(numerical,expected,atol=3e-15)


def normal_pair(u,cutoff):
    ns = np.arange(-cutoff,cutoff+1)
    ns = ns[ns%2!=0]
    e = ns*np.pi/2
    power = 16/(ns*np.pi)**4
    return float(power@np.exp(u*e-np.logaddexp(0.,2*e)))


normal = {}
for u in (.3,.5,.7):
    normal[str(u)] = normal_pair(u,128)
    np.testing.assert_allclose(normal[str(u)],normal_pair(u,64),atol=1e-15)
assert len(set(normal.values()))==3 and min(normal.values())>0

# Double the original normal interval. The mirror R_w: w->2-w and the
# old normal deck J: w->-w differ by four in their two orders. Thus
# their lifts commute for periodic transverse seam, anticommute for
# antiperiodic. On R_w=b subspace the signs at +1,-1 are b,eta*b.
y = sp.Symbol("y")
mirror = lambda w: 2-w
deck = lambda w: -w
assert mirror(deck(y))-deck(mirror(y))==4
alternatives = {
    "original_mobius_normal_line": {"normal_momenta":"n*pi/2","zero_multiplicity":1,
                                    "next_abs_normal_energy":"pi/2"},
    "klein_double_periodic_transverse": {"normal_momenta":"n*pi/2","zero_multiplicity":2,
                                          "mirror_and_deck_commute":True,
                                          "fixed_b_recovers_original_boundary_signs":True},
    "klein_double_antiperiodic_transverse": {"normal_momenta":"(n+1/2)*pi/2","zero_multiplicity":0,
                                              "min_abs_normal_energy":"pi/4",
                                              "mirror_and_deck_commute":False,
                                              "fixed_b_boundary_signs":"(b,-b)"},
}
# Independent numerical boundary checks for the reflected periodic
# extension of the original b eigenfunction, on the doubled interval.
for b,k in product((-1,1),range(-4,5)):
    def psi(w):
        return np.array([np.exp(-1j*k*np.pi*(w+1)/2),
                         b*np.exp(1j*k*np.pi*(w+1)/2)])/2
    s1 = np.array(sigma1,dtype=float)
    for w in (-1.,-.47,.25,1.):
        np.testing.assert_allclose(psi(w),b*s1@psi(2-w),atol=2e-15)
        np.testing.assert_allclose(psi(w),psi(w+4),atol=2e-15)

# Original Fuchs and quarter gates in the actual common signed frame.
sys.path.insert(0,str(repo/"verification"))
from v117_monodromy_weyl_a3 import M0_EXACT as M, U_EXACT as U, enumerate_group

# Common signed frame of the earlier Möbius/Fuchs test.
M6,U6 = sp.diag(M,-M),sp.diag(U,U)
swap = sp.BlockMatrix([[sp.zeros(3),I3],[I3,sp.zeros(3)]]).as_explicit()
sign = sp.diag(I3,-I3)
assert M**3==I3 and U**4==I3
assert all(sp.simplify(g*sign-sign*g)==sp.zeros(6) for g in (M6,U6))
assert swap*sign*swap==-sign and swap*U6*swap==U6 and swap*M6*swap==-M6
# U has three simple eigenspaces; its spectral projectors and the
# nonzero M links generate all nine matrix units on either sign block.
eigenvalues = (1,sp.I,-sp.I)
projectors = [sp.simplify(sp.prod((U-v*I3)/(u-v) for v in eigenvalues if v!=u)) for u in eigenvalues]
assert projectors==[sp.diag(*[int(i==j) for i in range(3)]) for j in range(3)]
assert all((projectors[i]*M*projectors[j])[i,j]!=0 for i in range(3) for j in range(3) if i!=j)
assert sign==M6**3
# The swap makes the off-diagonal matrix units from either diagonal
# M3 block. The group is still finite: the swap merely interchanges the
# two signs of (g, +/-g), and the native 3x3 subgroup is enumerated below.
size = len(enumerate_group([U,M]))
assert size>0
print(json.dumps({
    "old_signed_selector_with_both_outer_lines_retained":selector_cases,
    "old_same_time_reference_phase_sensitive_probe":old,
    "if_the_actual_ladder_and_endpoints_are_identified":{
        "two_step_over_whole_word_for_a_fixed_inner_selector":"s=<mid|P|mid>",
        "range":"0<=s<=1",
        "nested_retention":"1>=s_Theta>=s_str>=s_adm>=0",
        "Q_Theta_midpoint_expectation":"2*s_Theta-1",
        "s_invariant_under_a_static_selector_on_the_single_midpoint_for_all_packet_pairs":True,
    },
    "same_family_rotation_iterated_with_and_without_internal_projection":{
        "single_quarter_compressed_rank":quarter_rotation.rank(),
        "two_quarters_projected_vs_full_defect_rank":(whole_half-twice_projected).rank(),
        "full_rotation_recurs_after_four_quarters":True,
        "iterated_compressed_quarter_is_stationary_after_one_step":True,
    },
    "native_pin_label_reflection_even_higgs_odd_signs":profiles,
    "stationary_native_Paired_Pin_companion":{
        "swap_in_paired_b_frame":"i*tau2",
        "squares_to_minus_one":True,
        "commutes_with_fixed_original_paired_Green_generator_and_Higgs_mark":True,
        "b_sign_alone_does_not_commute_with_that_generator":True,
        "no_invertible_compensating_constant_gauge_grading_J_real_internal_sign":"nu_R Majorana block invertible with simple squared spectrum",
    },
    "actual_QTheta_Pstr_Padm":None,
    "original_deck_odd_normal_coordinate":{
        "normal_branch": "H_b=i sigma3 partial_w on b-fixed original interval",
        "coeff_to_n_neq0":"2*((-1)^n-1)/(pi*n)^2",
        "null_to_null":0,
        "norm_squared":str(norm_square),
        "commutator_with_H_b_norm_squared_on_null":str(energy_square),
        "conditional_NS_reflected_normal_gram_by_time_sum":normal,
        "flat_local_scalar_compensator_commuting_with_full_D":"none on connected cover",
        "interpretation":"Odd normal coordinate as geometric test section; no native Ramond coupling assumed.",
    },
    "geometry_alternatives":alternatives,
    "existing_fuchs_compiler_gate_algebra":{
        "block_diagonal_word_algebra_dimension":18,
        "same_with_actual_sign_exchange_available_dimension":36,
        "M_U_group_order":size,
        "signed_M6_U6_group_order":2*size,
        "signed_group_with_exchange_order":4*size,
        "physical_implementation_of_added_exchange_gate":None,
    },
},indent=2))
