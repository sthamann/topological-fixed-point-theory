"""Exact common folded charge/time test and local original R032 stress check."""
import ast
import importlib.util
import json
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp

from v1_e8_glue import e8_roots
from v1000_spin2_sk_witnesses import projector_matrices

repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[3]
sys.path.insert(0,str(repo))
from tfpt_explorer.current_block_geometry import _actual_current_source

contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
M = sp.Matrix(json.loads((repo / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json").read_text())["explicit_Weyl_matrix"])
syntax = ast.parse((repo/"verification/v634_st31_structure.py").read_text())
idx = next(ast.literal_eval(node.value) for node in syntax.body if isinstance(node, ast.Assign)
           and any(isinstance(t, ast.Name) and t.id=="SIGMA_IDX" for t in node.targets))
S = sp.eye(8)[idx,:]
assert M.T*M == sp.eye(8) and S**3 == sp.eye(8)

# Doubled coordinates of the original root set, in the existing target chart.
roots = sorted(tuple(int(x) for x in 2*r) for r in e8_roots())
assert set(roots)==set(_actual_current_source()["chevalley"].ridx)
ids = {r:i for i,r in enumerate(roots)}


def rotate(r):
    return tuple(r[i] for i in idx)


def disc(coords):
    if coords[0] % 2 == 0:
        return 2*(sum(x//2 for x in coords)%2)
    return 1+2*(sum((x-1)//2 for x in coords)%2)


grades = {}
for r in roots:
    grades[r] = disc(r[:5])
    assert grades[r] == disc(r[5:])
    assert rotate(r) in ids
counts = Counter((grades[r],grades[rotate(r)]) for r in roots)
grade_table = sp.Matrix(4,4,lambda i,j:counts[i,j])
assert sum(grade_table)==240
assert all((i-j)%2==0 for i,j in counts)
preserved_splits = []
for subset in combinations(range(8),5):
    target = tuple(i for i in range(8) if i not in subset)
    def class_in_split(root):
        value = disc(tuple(root[i] for i in subset))
        assert value == disc(tuple(root[i] for i in target))
        return value
    if all(class_in_split(r) == class_in_split(rotate(r)) for r in roots):
        preserved_splits.append(subset)
assert preserved_splits == [(0,2,4,6,7),(1,3,5,6,7)]
mixed = sp.simplify(sum(n*sp.I**(j-i) for (i,j),n in counts.items())/240)
phase = sp.SparseMatrix(240,240,{(i,i):sp.I**grades[r] for i,r in enumerate(roots)})
clock = sp.SparseMatrix(240,240,{(ids[rotate(r)],ids[r]):1 for r in roots})
assert clock.T*clock == sp.eye(240) and clock**3 == sp.eye(240)
neutral = sp.SparseMatrix(240,240,{(i,i):1 for i,r in enumerate(roots) if grades[r]==0})
covariant = clock*phase*clock.T
assert covariant*clock == clock*phase

seed = (2,0,0,2,0,0,0,0)
r1, r2 = rotate(seed), rotate(rotate(seed))
path_grades = [grades[r] for r in (seed,r1,r2)]
assert path_grades == [0,2,0]
e0 = sp.SparseMatrix(240,1,{(ids[seed],0):1})
u = sp.symbols("u", positive=True)  # exp(-t/32), t the original quotient time
cases = []
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    delta = sp.Matrix([(k+alpha)/8-sp.Rational(1,2) for k in range(8)])
    d = M*delta
    marked = tuple(2*x for x in M*(sp.eye(8)[:,0]+sp.eye(8)[:,2]))
    marked_orbit = (marked,rotate(marked),rotate(rotate(marked)))
    energy = {r:1+d.dot(sp.Matrix(r)/2) for r in roots}
    assert all((32*value).q==1 for value in energy.values())
    diagonal = sp.SparseMatrix(240,240,{(i,i):u**(32*energy[r]) for i,r in enumerate(roots)})
    norm = u**(64*energy[seed])
    word_first = (e0.T*phase.conjugate()*clock.T**2*diagonal*clock*phase*diagonal*clock*e0)[0]/norm
    word_last = (e0.T*phase.conjugate()*clock.T**2*phase*diagonal*clock*diagonal*clock*e0)[0]/norm
    exp = 32*(energy[r1]+energy[r2]-2*energy[seed])
    assert sp.simplify(word_first)==-u**exp and sp.simplify(word_last)==u**exp
    assert neutral*clock*e0 == sp.zeros(240,1)
    assert (e0.T*neutral*clock.T**2*diagonal*clock*phase*diagonal*clock*neutral*e0)[0]!=0
    moved = sp.SparseMatrix(240,240,{(i,i):u**(32*(1+(S*d).dot(sp.Matrix(r)/2)))
                                     for i,r in enumerate(roots)})
    assert moved*clock == clock*diagonal

    # The original folded Lhat_-1 response on a charged primary is its
    # shifted Cartan vector. Same d, energy and oscillator normalization.
    p = sp.Matrix(seed)/2
    vectors = [sp.Matrix(r)/2+d for r in (seed,r1,r2)]
    strengths = [v.dot(v) for v in vectors]
    assert all(strengths[j]-strengths[0]==2*(energy[r]-energy[seed])
               for j,r in enumerate((seed,r1,r2)))
    assert (S*p+S*d).dot(S*p+S*d)==strengths[0]
    # Check the ordered descendant insertion by its literal eight
    # one-oscillator coefficients, with unit oscillator Gram at level one.
    stress_inserted = sum(vectors[1][i]*vectors[1][i] for i in range(8))
    assert stress_inserted==strengths[1]
    assert energy[marked_orbit[1]]-energy[marked_orbit[0]]==sp.Rational(3,4)
    marked_g = [grades[r] for r in marked_orbit]
    marked_phases = [str(sp.I**(g-marked_g[0])) for g in marked_g]
    cases.append(dict(alpha=str(alpha), marked_r100_grades=marked_g, marked_r100_phases=marked_phases,
                      energies=[str(energy[r]) for r in (seed,r1,r2)],
                      first_phase_word=str(sp.simplify(word_first)),
                      last_phase_word=str(sp.simplify(word_last)),
                      stress_descendant_weights=list(map(str,strengths)),
                      stress_weight_change_under_first_step=str(strengths[1]-strengths[0])))


def circle_stress():
    spec = importlib.util.spec_from_file_location("quarter",contract/"check_quarter_reference.py")
    quarter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(quarter)
    DY,DM = quarter.higgs_matrix,quarter.mass_matrix
    charged = np.flatnonzero(np.abs(DM).sum(axis=0)==0)
    active = charged[np.abs(DY[np.ix_(charged,charged)]).sum(axis=0)>0]
    gamma = quarter.grading[np.ix_(active,active)]
    Y = -1j*gamma@DY[np.ix_(active,active)]
    plus = np.diag(gamma).real>0
    minus = ~plus
    assert plus.sum()==minus.sum()==42
    assert np.linalg.norm(gamma@Y+Y@gamma)<1e-14
    # The symmetric tangential metric insertion has internal coefficient gamma.
    # In the particle-hole matrix elements the momentum sign fixes the occupied grade.
    for sign_p in (-1,1):
        for sign_q in (-1,1):
            initial = minus if sign_p>0 else plus
            final = plus if sign_q>0 else minus
            yw = np.linalg.norm(Y[np.ix_(final,initial)])**2
            tw = np.linalg.norm(gamma[np.ix_(final,initial)])**2
            assert (yw>0)==(sign_p==sign_q)
            assert (tw>0)==(sign_p!=sign_q)
            assert yw*tw==0
    tags = quarter.native.P
    charges = np.array([(1 if i<len(tags) else -1) if tags[i%len(tags)][0] in ("u","d") else 0
                        for i in active])
    Z = np.diag(np.exp(2j*np.pi*charges/3))
    assert np.sum(charges[plus])%3==np.sum(charges[minus])%3==0
    assert np.linalg.norm(Z@Y-Y@Z)<1e-13 and np.linalg.norm(Z@gamma-gamma@Z)<1e-13

    # Smooth fixed support f on [-1,1], integral one, as in TFPTSTARTOVER20261109.
    nodes,weights = leggauss(96)
    f = weights*np.exp(-1/(1-nodes**2))
    f /= f.sum()
    radius, eps = 16384, .05
    s = np.arange(int(radius*eps))
    E = (s+1)/radius
    Fourier_squared = (np.cos(np.outer(E,nodes))@f)**2
    # Both opposite-sign directions, each of rank 42. n+m=s, shifted signs.
    signed_square_sum = (s+1)*(s*(s+2)/3 + .25)
    weight = (2*42/(2*np.pi*radius)**2)*signed_square_sum*Fourier_squared/(4*radius**2)
    mass_ratio = weight.sum()/eps**4
    inverse_ratio = (weight/E).sum()/eps**3
    expected_mass = 42/(96*np.pi**2)
    expected_inverse = 42/(72*np.pi**2)
    assert abs(mass_ratio/expected_mass-1)<.01
    assert abs(inverse_ratio/expected_inverse-1)<.01

    # One actual tangent perturbation; no new transverse directions inserted.
    tt = projector_matrices([.4,0,0])["TT"]
    tangent = np.array([1.,0.,0.,0.,0.,0.])
    assert np.max(np.abs(tt@tangent))==0
    assert np.linalg.matrix_rank(tt,tol=1e-12)==2
    return dict(tangent_mass_ratio=float(mass_ratio),tangent_mass_limit=float(expected_mass),
                tangent_inverse_ratio=float(inverse_ratio),tangent_inverse_limit=float(expected_inverse),
                mixed_local_Y_tangent_stress_zero=True,
                tensor_projection_of_circle_tangent=0,
                tensor_rank_of_separate_v1000_control=2)


print(json.dumps(dict(glue_table=grade_table.tolist(),roots_changing_fixed_grade=240-int(sp.trace(grade_table)),
                      mixed_fixed_grade_character=str(mixed),root_grade_path=path_grades,
                      grade_parity_preserved=True, separately_neutral_midcut_word=0,
                      alternative_coordinate_splits=preserved_splits,
                      folded_cases=cases,circle=circle_stress()),default=int,indent=2))
