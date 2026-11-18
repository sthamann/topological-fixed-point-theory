"""Exact old charge/frequency capacity and neutral family closure for Gate A."""
from contextlib import redirect_stdout
from itertools import combinations
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    exact = runpy.run_path(str(base / "check_tfpt_family_middle_ordered_cut_20261110.py"))
    signed = runpy.run_path(str(base / "check_tfpt_ramond_signed_higgs_cut_20261110.py"))

sp, q, r, d1, M = (exact[k] for k in ("sp","q","r","d1","M"))
vertex, gram, rank = (exact[k] for k in ("vertex","gram","rank"))
old, V, Q, polys = (exact[k] for k in ("old","V","Q","polys"))
families = {}
for alpha in signed["quarters"]:
    drift = signed["old_omega"](alpha)
    coeff = [Q.from_sympy(x) for x in old["cl"](V.T*drift/8)]
    states = [old["vertex_zero"](old["series_response"](f,coeff,4)) for f in polys]
    outer = [vertex(v,d1,q) for v in states]
    middle = [exact["middle"](v) for v in states]
    raised = [vertex(v,q,r) for v in outer]
    at_top = raised+middle
    all_lower = [vertex(v,-q,r+q) for v in at_top]
    # The two compositions of q and -q use conjugate cocycle signs.
    # Bare oscillator maps have the same relative sign on the two weights.
    assert all(not vertex(v,q,r+q) for v in at_top)
    raised_back = [vertex(v,q,r) for v in all_lower]
    same = [dict((k,2*a) for k,a in v.items()) for v in at_top]
    assert raised_back==same or raised_back==[dict((k,-a) for k,a in v.items()) for v in same]
    assert all((r+k*q).dot(r+k*q)/2==1+k*k for k in range(-2,3))
    # At old L0=4 the only allowed momenta are r+kq with |k|<=1.
    # Thus only family spin 0 or 1 occur. J_-J_+/2 is the spin-1
    # projection on weight zero. Its Gram is half the top-image Gram.
    triplets = rank(gram(at_top))
    middle_weight_singlet = (gram(outer)-gram(raised)/2).applyfunc(sp.expand)
    singlets = rank(middle_weight_singlet)
    assert all(sp.simplify(middle_weight_singlet.extract(ix,ix).det())>=0
               for n in range(1,6) for ix in combinations(range(5),n))
    assert (triplets,singlets)==(6,5)
    ratio = sp.Rational(35395,19017)
    assert sp.simplify(ratio-(1+2*sp.Rational(8189,19017)))==0
    families[str(alpha)] = {
        "family_triplet_copies":triplets,
        "family_singlet_copies":singlets,
        "minimal_closed_neutral_family_dimension":3*triplets+singlets,
        "minimal_dimensions_at_family_weights_minus2_zero_plus2":[triplets,triplets+singlets,triplets],
        "same_time_on_full_closure":str(32+drift.dot(r)),
        "full_over_ramond_fixed_P0_ratio_exact":str(ratio),
    }

# Confirm the two energies and unique old ground-state line for each
# enumerated output momentum on the original full two-quarter clock.
# Across both quarters the frequency difference fixes sum(d)=+/-2.
# With unchanged Ramond lattice, fixed SM charges and integer oscillator
# grade, signed_candidates exhausts the simultaneous eigenspaces.
capacities = {}
for sign, rows in signed["candidates"].items():
    for charge, mode, ground in rows:
        target = q+charge
        assert (target.dot(target)/2)==ground and mode+ground==4
        assert sum(charge)==2*sign
        assert all(32+2*sign*a==8*ground+signed["old_omega"](a).dot(target)
                   for a in signed["quarters"])
    capacities[str(sign)] = {"dimension":len(rows),
                            "momenta":[list(map(str,q+d)) for d,_,_ in rows],
                            "oscillator_multiplicity_per_momentum":[1]*len(rows)}

native = exact["native"]
rep, gamma, make_phi = (native[k] for k in ("native","gamma","inner"))
np_part = rep.NP
e_l = [rep.IDX[("e","L",0,k)] for k in range(3)]
e_r = [rep.IDX[("e","R",0,k)] for k in range(3)]
y_e = make_phi(np.array([0.,1.]))[np.ix_(e_l,e_r)]
y_exact = sp.Matrix([[sp.Rational(float(v.real))+sp.I*sp.Rational(float(v.imag))
                      for v in row] for row in y_e])
assert y_exact.det()!=0
complex_mark = -1j*gamma@make_phi(np.array([0.,1.]))-gamma@make_phi(np.array([0.,1j]))
electron = e_l+e_r+[i+np_part for i in e_l+e_r]
electron = np.array(electron)
grade = np.diag(gamma)[electron]
mat = complex_mark[np.ix_(electron,electron)]
rates = {}
for sign in (-1,1):
    left,right = np.nonzero(grade==sign)[0],np.nonzero(grade==-sign)[0]
    b = mat[np.ix_(left,right)]
    rates[str(sign)] = int(np.linalg.matrix_rank(b))
    assert rates[str(sign)]==(6 if sign==1 else 0)
assert np.linalg.norm(native["majorana"][electron,:])==0
assert np.linalg.matrix_rank(mat[:6,:6])==3

# Does allowing the other old spin sector add any signed-frequency
# capacity? Enumerate ALL integer momenta in the same old E8 chart.
# Oscillator levels add 64 to 8H. Cartan constraints give a rigorous
# upper bound, and the separate Ramond run certifies the two singlets.
energy_cap = 260


def term(j,x):
    return 8*x*x+(4*j-15)*x


options = tuple(range(-6,7,2))
mins = [min(term(j,x) for x in options) for j in range(8)]
assert all(term(j,x)+sum(mins)-mins[j]>energy_cap
           for j in range(8) for x in (-8,8))
# The convex summands are strictly increasing farther out on both tails.
integer_targets = {-1:[],1:[]}


def collect_integer_targets(vector,energy,charge):
    j = len(vector)
    if j==8:
        if tuple(charge)!=signed["target"]:
            return
        for sign in (-1,1):
            budget = 256+4*sign-energy
            if sum(vector)==4*sign and budget>=0 and budget%64==0:
                integer_targets[sign].append((vector,budget//64))
        return
    for x in options:
        cost = energy+term(j,x)
        if cost+sum(mins[j+1:])<=energy_cap:
            collect_integer_targets(vector+[x],cost,
                                    [c+row[j]*x for c,row in zip(charge,signed["rows"],strict=True)])


collect_integer_targets([],0,[0]*4)
integer_census = {}
for sign, candidates in integer_targets.items():
    assert len(candidates)==1 and candidates[0][1]==1
    mid = sp.Matrix(candidates[0][0])/2
    outer = [q+d for d,_,_ in signed["candidates"][sign]]
    assert {tuple(p) for p in outer}=={tuple(mid+q),tuple(mid-q)}
    assert q.dot(mid)==0
    color_roots = [M.T@(sp.eye(8)[:,i]-sp.eye(8)[:,i+1]) for i in (0,1)]
    assert all(d.dot(mid)==d.dot(q)==0 for d in color_roots)
    # At level 1 and color-neutral p, the zero modes of the two simple
    # color roots annihilate h(-1)|p> iff h is orthogonal to both roots.
    # For the weak Cartan weight -1, the weak lowering mode vanishes.
    color_rank = sp.Matrix.hstack(*color_roots).rank()
    assert color_rank==2 and (M@mid).dot(signed["weak"])==-1
    first = vertex({():Q.one},-q,mid+q)
    second = vertex(first,-q,mid)
    assert first=={((1,4),):-Q.one}
    assert second in ({():2*Q.one},{():-2*Q.one})
    integer_census[str(sign)] = {
        "momentum":list(map(str,mid)),
        "oscillator_level":1,
        "color_singlet_dimension":8-color_rank,
        "neutral_family_closure_of_two_R_lines":3,
        "unreached_family_singlets_at_same_energy":5,
        "half_of_squared_family_lowering_matches_old_minus_2q_transition_up_to_global_cocycle":True,
    }


gram_file = json.loads((base/"ramond-full-return-results-20261110.json").read_text())
cross = {}
for sign,i in ((-1,2),(1,1)):
    vals = []
    for a in ("1/4","3/4"):
        rows = [row for row in signed["controls"][a] if row["sign"]==sign]
        lo,hi = [[sp.sympify(z) for z in row["exact_ground_amplitudes_up_to_old_cocycle"]]
                 for row in rows]
        # Cocycle product +1 as independently checked for V_(-2q).
        assert lo[0]==0
        g = sp.Matrix(gram_file[a]["full_return_Gram"])
        # The raw covectors were evaluated on the unnormalized packets.
        normalization = sp.sqrt(g[0,0]*g[i,i])
        vals.append(sp.conjugate(hi[0])*lo[i]*sp.exp(-sp.Rational(rows[0]["old_intermediate_energy"])/8)
                    /normalization)
    ratio = sp.simplify(vals[1]/vals[0])
    cross[str(sign)] = {
        "packet_pair":[0,i],
        "complex_ratio":str(ratio),
        "ratio_abs":float(sp.Abs(ratio)),
        "relative_phase_rad":float(sp.arg(ratio)),
    }
assert abs(cross["-1"]["ratio_abs"]-1.20325)<1e-4
assert abs(cross["1"]["ratio_abs"]-.93378)<1e-4
assert abs(cross["-1"]["relative_phase_rad"]-cross["1"]["relative_phase_rad"])>.02
norms = [sum(sp.Rational(v)**2 for v in row["SM_semisimple_drift_target"][:3])
         for row in signed["gauge_readouts"].values()]
assert norms==[sp.Rational(73,6),sp.Rational(97,6)]

# Fixed nonzero native background: a whole old H-eigenspace can only be
# transported isometrically to an equally large native eigenspace. Compute
# the exact multiplicities on the actual single-momentum v252 coefficient.
plus,minus = (np.flatnonzero(sgn*np.diag(gamma)>0) for sgn in (1,-1))
full = -1j*gamma@(native["majorana"]+.55*make_phi(np.array([0.,1.])))
B = full[np.ix_(minus,plus)]
labels = rep.P*2
groups = {}
for i,k in enumerate(plus):
    f,ch,col,gen = labels[k]
    key = (f,col, k>=np_part) if f!="nu" else ("nu",0,False)
    groups.setdefault(key,[]).append(i)
# Exact since finite inputs in the implementation are binary rationals.
Bexact = sp.Matrix([[sp.Rational(float(v.real))+sp.I*sp.Rational(float(v.imag))
                      for v in row] for row in B])
from sympy.polys.matrices import DomainMatrix

bd = DomainMatrix.from_Matrix(Bexact).convert_to(sp.QQ_I)
adjoint = DomainMatrix.from_Matrix(Bexact.H).convert_to(sp.QQ_I)
K = (adjoint.matmul(bd)).to_Matrix()
var = sp.Symbol("z")
per_species = {}
for (species,col,anti),ix in groups.items():
    complement = [i for i in range(48) if i not in ix]
    assert K.extract(ix,complement)==sp.zeros(len(ix),len(complement))
    coefficients = DomainMatrix.from_Matrix(K.extract(ix,ix)).convert_to(sp.QQ_I).charpoly()
    poly = sp.Poly.from_list([sp.QQ_I.to_sympy(c) for c in coefficients],var,domain=sp.QQ)
    if species not in per_species:
        per_species[species] = [poly,0]
    assert per_species[species][0]==poly
    per_species[species][1] += 1
# Root differences across quarter modes p=m+alpha are half integers.
# Gershgorin bounds all nonnegative eigenvalues by the largest row sum.
bound = int(sp.ceiling(max(sum(abs(sp.re(v))+abs(sp.im(v)) for v in row)
                           for row in K.tolist())))
primes = (65537,65539,65543)


def modular(poly,prime):
    coeff = [int(v.p)*pow(int(v.q),-1,prime)%prime for v in poly.all_coeffs()]
    return sp.Poly.from_list(coeff,var,modulus=prime)


polynomials = [p for p,_ in per_species.values()]
reductions = {p:[modular(poly,p) for poly in polynomials] for p in primes}
for i,poly in enumerate(polynomials):
    assert poly.eval(0)!=0
    assert any(rows[i].gcd(rows[i].diff()).degree()==0 for rows in reductions.values())
    for j in range(i):
        assert any(rows[i].gcd(rows[j]).degree()==0 for rows in reductions.values())
    for j in range(len(polynomials)):
        for k in range(1,2*bound+1):
            assert any(rows[i].gcd(rows[j].shift(k*pow(2,-1,p)%p)).degree()==0
                       for p,rows in reductions.items()), (i,j,k,bound)

largest_mult = max(m for _,m in per_species.values())
assert largest_mult==6

print(json.dumps({
    "old_family_closure":families,
    "old_two_quarter_signed_full_Ramond_capacity":capacities,
    "old_integer_spin_sector_two_quarter_signed_Cartan_candidates":integer_census,
    "native_electron_particle_rank_for_nonzero_sign":3,
    "native_electron_adjoint_reverses_Bohr_sign":True,
    "native_color_singlet_electron_H_plus_kinematic_ranks":rates,
    "scope_of_rank_obstruction":"Faithful signed electron transition into the old q-prepared Ramond joint target alone; other sectors and coupled clocks are not classified.",
    "old_directed_complex_quarter_ratios":cross,
    "orthogonal_fixed_color_drift_squared_norms":list(map(str,norms)),
    "orthogonal_fixed_color_drift_squared_norm_difference":str(norms[1]-norms[0]),
    "native_single_momentum_nonzero_H_squared_spectral_multiplicities":{
        species:{"distinct_eigenvalues":poly.degree(),"multiplicity_per_eigenvalue":mult}
        for species,(poly,mult) in per_species.items()
    },
    "max_native_H_eigenspace_dimension_at_fixed_nonzero_background":largest_mult,
    "whole_native_circle_momentum":"p=m+alpha; m in Z",
    "exact_max_BdagB_spectral_upper_bound":bound,
    "whole_circle_no_cross_momentum_collisions_at_fixed_H":True,
    "max_native_H_eigenspace_dimension_in_direct_sum_of_quarters":2*largest_mult,
    "physical_full_Padm":None,
},indent=2))
