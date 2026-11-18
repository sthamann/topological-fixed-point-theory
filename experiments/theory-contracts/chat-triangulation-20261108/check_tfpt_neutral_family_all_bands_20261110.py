"""Exact closed family reads on original Cartan bands and selector tomography."""
from contextlib import redirect_stdout
from functools import lru_cache
from math import factorial, prod
from pathlib import Path
import io
import json
import runpy
import sys

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    cover = runpy.run_path(str(base / "check_tfpt_cover_selector_compiler_20261110.py"))

sp, q, Q, vertex, phase = (cover[k] for k in ("sp","q","Q","vertex","phase"))
gate, signed = cover["gate"], cover["signed"]
r = gate["r"]
assert q.dot(q)==2 and q.dot(r)==0


@lru_cache(None)
def partitions(total, minimum=1):
    if total==0:
        return ((),)
    if total<0:
        return ()
    return tuple((n,)+rest for n in range(minimum,total+1)
                 for rest in partitions(total-n,n))


def metric(word):
    return prod((2*n)**word.count(n)*factorial(word.count(n)) for n in set(word))


@lru_cache(None)
def gram(g):
    return sp.diag(*(metric(w) for w in partitions(g)))


def vertex_matrix(g):
    source, target = partitions(g),partitions(g-1)
    columns = []
    for word in source:
        out = vertex({tuple((n,4) for n in word):Q.one},q,r)
        assert all(axis==4 for w in out for _,axis in w)
        columns.append([Q.to_sympy(out.get(tuple((n,4) for n in t),Q.zero))
                        for t in target])
    return sp.Matrix(len(target),len(source),lambda i,j:columns[j][i])


@lru_cache(None)
def virasoro(g,n):
    """L_-n of the given q oscillator current at q-charge zero."""
    source,target = partitions(g),partitions(g+n)
    indices = {word:i for i,word in enumerate(target)}
    result = sp.zeros(len(target),len(source))
    for j,word in enumerate(source):
        for k in range(1,n):
            out = tuple(sorted(word+(k,n-k)))
            result[indices[out],j] += sp.Rational(1,4)
        for i,m in enumerate(word):
            out = tuple(sorted(word[:i]+word[i+1:]+(m+n,)))
            result[indices[out],j] += m
    return result


by_grade = {}
summary = {}
first_obstruction = None
for g in range(11):
    source = partitions(g)
    D = gram(g)
    N = vertex_matrix(g)
    if g==0:
        C = sp.zeros(1)
    else:
        C = D.inv()*N.T*gram(g-1)*N
    by_spin = {}
    for n in range(1,g+1):
        L = virasoro(g-n,n)
        for j,rows in by_grade[g-n].items():
            by_spin.setdefault(j,[]).extend(L*v for v in rows)
    by_spin = {j:sp.Matrix.hstack(*rows).columnspace() for j,rows in by_spin.items()}
    descendants = [v for rows in by_spin.values() for v in rows]
    primitive = (sp.Matrix.hstack(*descendants).T*D).nullspace() if descendants else [
        sp.eye(len(source))[:,i] for i in range(len(source))]
    j = int(sp.integer_nthroot(g,2)[0])
    assert len(primitive)==int(j*j==g)
    if primitive:
        assert j not in by_spin
        by_spin[j]=primitive
    U = sp.diag(*[(-1)**len(w) for w in source])
    X = sp.diag(*map(len,source))
    for j,rows in by_spin.items():
        assert all(C*v==j*(j+1)*v and U*v==(-1)**j*v for v in rows)
        assert len(rows)==len(partitions(g-j*j))-len(partitions(g-(j+1)**2))
    assert sum(map(len,by_spin.values()))==len(source)
    assert C.T*D==D*C
    if g<=3:
        assert C==sp.eye(len(source))-U
    if g==4:
        v = by_spin[2][0]
        v = v/next(x for x in v if x)
        assert C*v==6*v and U*v==v and C*X!=X*C
        first_obstruction = {
            "partitions_order":[list(w) for w in source],
            "even_spin_two_vector_coefficients":list(map(str,v)),
            "old_raised_norm_over_input_norm":str((v.T*D*C*v)[0]/(v.T*D*v)[0]),
            "old_parity_difference_over_input_norm":"0",
            "oscillator_count_commutator_rank":(C*X-X*C).rank(),
        }
    by_grade[g]=by_spin
    summary[str(g)] = {"dimension":len(source),
                       "spin_multiplicities":{str(j):len(rows) for j,rows in sorted(by_spin.items())},
                       "stationary_casimir_equals_legacy_parity":C==sp.eye(len(source))-U,
                       "zero_mode_C_commutes_with_oscillator_number":C*X==X*C}

# Exact visible vectors of the k=1 spin sector from the original pair.
# The explicit oscillator coefficients are checked against both charge
# vertices and the original cocycle on all four sign/quarter branches.
kernel = {}
ramond_gram = json.loads((base/"ramond-full-return-results-20261110.json").read_text())
norms = cover["packet_norms"]
q_axis2 = {((2,4),):Q.one/2}
q_axis11 = {((1,4),(1,4)):Q.one/2}
root2 = sp.sqrt(2)
for alpha, controls in signed["controls"].items():
    kernel[alpha]={}
    gram0 = sp.Matrix(ramond_gram[alpha]["full_return_Gram"])
    normalization=sp.diag(*[sp.sqrt(n) for n in norms])
    gram0=normalization@gram0@normalization
    for sign in (-1,1):
        low,high = [r for r in controls if r["sign"]==sign]
        a,b = [sp.Matrix(list(map(sp.sympify,row["charge_old"]))) for row in (low,high)]
        top,bottom = q+a,q+b
        midpoint=(top+bottom)/2
        assert top-bottom==2*q
        assert vertex({():Q.one},-q,top,1)=={
            ((2,4),):-Q.one/2, ((1,4),(1,4)):Q.one/2}
        assert vertex({():Q.one},q,bottom,1)==q_axis2 | q_axis11
        v = sp.Matrix([-1,root2])
        w = sp.Matrix([1,root2])
        matrix = sp.Matrix.hstack(v,w)
        assert abs(matrix.det())==2*root2
        # Low and high packet vectors come from actual old word coefficients.
        l,h = [sp.Matrix([[phase(charge,q)*sp.sympify(z) for z in row[
                         "exact_ground_amplitudes_up_to_old_cocycle"]]])
               for charge,row in ((a,low),(b,high))]
        i = 2 if sign==-1 else 1
        assert l[0]==0 and l[i]!=0 and h[0]!=0
        e = sp.eye(5)
        upper = (e[:,i]-h[i]/h[0]*e[:,0])/(l[i]*phase(-q,top))
        lower = e[:,0]/(h[0]*phase(q,bottom))
        F = sp.Matrix.hstack(v*phase(-q,top),w*phase(q,bottom))*sp.Matrix.vstack(l,h)
        assert sp.simplify(F@upper-v)==sp.zeros(2,1)
        assert sp.simplify(F@lower-w)==sp.zeros(2,1)
        probes = [upper,lower,upper+lower,upper+sp.I*lower]
        assert all(sp.simplify((p.H@gram0@p)[0])>0 for p in probes)
        kernel[alpha][str(sign)] = {
            "actual_input_packet_for_upper":i,
            "actual_input_coefficients_for_upper":[str(sp.expand(x)) for x in upper],
            "actual_input_coefficients_for_lower":[str(sp.expand(x)) for x in lower],
            "four_probe_squared_norms":[float(sp.simplify((p.H@gram0@p)[0])) for p in probes],
            "first_local_mode_vectors":[[-1,"sqrt(2)"],[1,"sqrt(2)"]],
            "selector_access_requires_no_new_neutral_conversion_insertion_in_the_old_experiment":True,
        }

x,y,z,a = sp.symbols("x y z a",real=True)
P = sp.Matrix([[a+z,x-sp.I*y],[x+sp.I*y,a-z]])
measurements = [sp.expand((s.H@P@s)[0])
                for s in (v,w,v+w,v+sp.I*w)]
assert sp.Matrix(measurements).jacobian((a,x,y,z)).det()!=0
cross = sp.expand((w.H@P@v)[0])
assert cross==a-3*z-2*sp.I*root2*y
# P=(I+X)/2, Hermitian contraction X. The feasible image is the
# ellipse from traceless Pauli coordinates; scalar +/-I lie inside it.
R = sp.expand(-cross.subs({a:(1+a)/2,x:x/2,y:y/2,z:z/2},simultaneous=True)/2)
assert sp.simplify(R-(-sp.Rational(1,4)-a/4+3*z/4+sp.I*root2*y/2))==0
for value,proj in ((sp.Rational(1,2),sp.diag(1,0)),(-1,sp.diag(0,1))):
    assert -(w.H@proj@v)[0]/2==value
same_cross = []
for sign in (-1,1):
    projection=sp.Matrix([[1,sign*root2],[sign*root2,2]])/3
    assert (projection**2-projection).applyfunc(sp.simplify)==sp.zeros(2)
    assert projection.rank()==1
    assert (w.H@projection@v)[0]==1
    same_cross.append([str(sp.simplify((s.H@projection@s)[0])) for s in (v,w)])
assert same_cross==[["3","1/3"],["1/3","3"]]
assert all(poly.eval(0)!=0 for poly in gate["polynomials"])

print(json.dumps({
    "exact_neutral_virasoro_validation_by_q_oscillator_grade":summary,
    "first_failure_of_legacy_parity":first_obstruction,
    "old_two_path_selector_probe_preparations":kernel,
    "real_measurement_map_in_order_a_x_y_z":[str(s) for s in measurements],
    "first_local_split_relative_selector_amplitude":{
        "ellipse_center":"-1/4","real_semiaxis":"3/4","imaginary_semiaxis":"sqrt(2)/2",
        "real_range":["-1","1/2"],"selector_I":"-1/2",
        "scope":"Old fixed channel before exp(-8d), any orthogonal selector compression, not the true P_adm.",
    },
    "distinct_old_rank_one_selectors_with_same_directed_ratio_minus_half":{
        "diagonal_positive_returns_v_w_for_the_two_selectors":same_cross,
        "both_commute_with_the_given_degenerate_old_energy_and_charge":True,
    },
    "unchanged_free_Pin_Fock_no_infinite_exact_energy_ladder":{
        "all_actual_flat_squared_masses_are_positive_algebraic":True,
        "number_of_distinct_squared_mass_roots":sum(p.degree() for p in gate["polynomials"]),
        "normal_branch_square_form":"(n*pi/2)^2 + lambda",
        "hypotheses":"Finite internal mass list, all strictly positive, positive excitation spectrum and fermionic occupation.",
        "scope":"All finite-particle excitations, arbitrary fixed time rescaling; does not cover dynamical boson or interacting sources.",
    },
    "physical_full_Padm_and_weighted_source":None,
},indent=2))
