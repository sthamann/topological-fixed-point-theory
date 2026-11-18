"""Original folded Higgs field as a convolution of two commuting root fields."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import runpy
import sys

repo = Path(sys.argv[1]).resolve()
root = repo / "experiments/theory-contracts"
with redirect_stdout(io.StringIO()):
    src = runpy.run_path(str(root / "chat-triangulation-20261108/check_tfpt_family_middle_ordered_cut_20261110.py"))
w = runpy.run_path(str(root / "universalraum-native-operations-ground-response-20260915/common.py"))
W = w["load_tensor"]()
_, boson_weights = w["weights"]()
sp, K, q, r, B, M, V = (src[k] for k in ("sp", "Q", "q", "r", "d1", "M", "V"))
old, combine, grade, vertex, inner = (src[k] for k in ("old", "combine", "grade", "vertex", "inner"))
metric = src["metric"]


def times(vector, scalar):
    return {word: a * scalar for word, a in vector.items() if a * scalar}


def add(*vectors):
    return combine((word, a) for v in vectors for word, a in v.items())


def cartan(vector, momentum, mode):
    if mode == 0:
        return times(vector, K.from_sympy(q.dot(momentum)))
    if mode < 0:
        return {tuple(sorted(word + ((-mode, 4),))): a for word, a in vector.items()}
    return combine((word[:j] + word[j+1:], mode * metric[4][i] * a)
                   for word, a in vector.items() for j, (n, i) in enumerate(word) if n == mode)


def norm(vector):
    return sp.simplify(inner(vector, vector))


def floor_degree(vector, momentum, charge):
    if not vector:
        return 0
    grades = {grade(word) for word in vector}
    assert len(grades) == 1
    return int(grades.pop() + momentum.dot(momentum)/2 - (momentum+charge).dot(momentum+charge)/2)


def descendant_parts(vector, momentum, total_mode, charge, cartan_field=False):
    # R_(total_mode-n) and the q-current J_n commute since q.r=0.
    # Positive modes stop once the intermediate oscillator degree is negative.
    n_min = total_mode - floor_degree(vector, momentum, r)
    n_max = floor_degree(vector, momentum, charge)
    result = {}
    for n in range(n_min, n_max+1):
        z = vertex(vector, r, momentum, n-total_mode)
        out = (cartan(z, momentum+r, n) if cartan_field
               else vertex(z, charge, momentum+r, -n))
        if out:
            result[n] = out
    return result


def equal(a, b):
    return not add(a, times(b, -1))


def double_weight(x):
    return tuple(int(t) for t in 2*M*x)


assert W.shape == (60, 2016)
bweights = {tuple(a) for a in boson_weights.tolist()}
assert double_weight(r) in bweights
assert double_weight(r+q) not in bweights and double_weight(B) not in bweights
assert max(abs(int(sum(a*b for a,b in zip(weight, double_weight(q)))//4)) for weight in bweights) == 1
assert q.dot(B) == -2 and q.dot(r) == 0
assert (M*q)[:5,0] == sp.zeros(5,1)

vacuum = {():K.one}
assert vertex(vacuum,r,sp.zeros(8,1),1) == vacuum
assert vertex(vacuum,r,r,1) == {}
# In the original canonical native boson algebra the same double creation
# from the vacuum has squared norm 2 by [b_r,b_r^*]=1.
result = {}
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    drift = sp.Matrix([j+alpha-4 for j in range(8)])
    inputs = [old["series_response"](f,[K.from_sympy(z) for z in old["cl"](V.T*drift/8)],4)
              for f in src["polys"]]
    packets = [old["vertex_zero"](f) for f in inputs]
    rows = []
    for i, v in enumerate(packets):
        pieces = descendant_parts(v, q, 0, -q)
        actual = vertex(v, B, q)
        assert equal(add(*pieces.values()), actual)
        assert equal(add(*descendant_parts(v, q, 0, sp.zeros(8,1), True).values()),src["middle"](v))
        mode_zero = pieces.get(0,{})
        # J_+J_- equals two on spin one at weight zero.
        # The old q family part occupies spins 0 and 1 on this band.
        raised = vertex(actual,q,r)
        spin_1 = times(vertex(raised,-q,r+q),K.from_sympy(sp.Rational(1,2)))
        assert equal(vertex(spin_1,q,r),raised)
        singlet = add(actual,times(spin_1,-1))
        assert not inner(singlet,spin_1)
        assert not inner(singlet,mode_zero)
        assert norm(vertex(mode_zero,q,r)) == 2*norm(mode_zero)
        assert all(norm(z)>=0 for z in (singlet,spin_1,mode_zero))
        modes = {}
        for m in (-2,-1,0,1,2):
            j_before = vertex(v,q,q,-m)
            early = vertex(j_before,B,2*q) if j_before else {}
            late = vertex(actual,q,r,-m)
            middle = add(*descendant_parts(v, q, m, sp.zeros(8,1), True).values())
            r_m = vertex(v,r,q,-m)
            # q,r commute; original cocycle product is constant at fixed momentum.
            assert equal(add(late,times(early,-1)),add(middle,times(r_m,m)))
            modes[str(m)] = {
                "early_norm2":str(norm(early)), "late_norm2":str(norm(late)),
                "middle_norm2":str(norm(middle)), "central_contact_norm2":str(norm(times(r_m,m))),
                "interference_2re":str(sp.simplify(2*m*sp.re(inner(middle,r_m)))),
                "difference_norm2":str(norm(add(late,times(early,-1)))),
            }
        # Exact comparison with a diagnostic high-pass filter of the old H.
        # The physical gap is a different operator. The common old weights
        # are removed by their eigenvalue factors before taking the difference.
        v_norm = norm(actual)
        plus_norm = sp.Rational(modes["-1"]["late_norm2"])
        minus_norm = norm(vertex(actual,-q,r,-1))
        assert plus_norm-minus_norm == v_norm
        assert drift.dot(q) == 0 and drift.dot(r) == 2*alpha-2
        energy = 32 + drift.dot(r)
        gaps = {"eigenvalues":list(map(str,(energy-8,energy,energy+8))),
                "ward_unmarked":str(v_norm),
                "after_old_high_pass_at_B_energy":str(plus_norm),
                "unrecovered_below_cutoff":str(minus_norm)}
        rows.append({"packet":i,"old_energy_cutoff_diagnostic":gaps,"convolution_indices":list(pieces),"B_norm2":str(norm(actual)),
                     "stationary_piece_norm2":str(norm(mode_zero)),"spin_0_norm2":str(norm(singlet)),
                     "spin_1_norm2":str(norm(spin_1)),"modes":modes})
    result[str(alpha)] = rows

print(json.dumps({"doubled_weights":{"q":double_weight(q),"r":double_weight(r),"B":double_weight(B)},
                  "first_native_r_vs_affine_r":{"one_particle_norm2":1,"native_two_particle_norm2":2,"affine_two_particle_norm2":0},"old_source":result,"physical_values":None},indent=2))
