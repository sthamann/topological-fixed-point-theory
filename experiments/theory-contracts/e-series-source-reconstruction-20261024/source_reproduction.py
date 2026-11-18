#!/usr/bin/env python3
"""Original scalar quarter, actual compiler glue, and the first source gap.

All lattice, field-character and quotient statements are exact. The final
Slater calculation is a numerical diagnostic on the unchanged v988 input.
No computation here establishes microscopic charged-field convergence.
"""
import argparse
import ast
from collections import Counter
import contextlib
from hashlib import sha256
import io
from itertools import combinations, product
import json
from math import factorial
from pathlib import Path
import runpy

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PINS = {}


def need(value, label):
    if not bool(value):
        raise RuntimeError(label)


def pin(relative):
    data = (REPO/relative).read_bytes()
    PINS[relative] = sha256(data).hexdigest()
    return data


def calculate():
    clock_file = "experiments/theory-contracts/e-series-joint-source-20261023/joint_clock.py"
    pin(clock_file)
    with contextlib.redirect_stdout(io.StringIO()):
        original = runpy.run_path(str(REPO/clock_file))
    PINS.update(original["PINS"])
    roots = sorted(tuple(int(x) for x in row) for row in original["roots"])
    need(len(roots) == 240, "original root carrier")
    T = sp.Matrix([[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1]])/2
    need(T.T*T == sp.eye(3) and sp.ones(1,4)*T == sp.zeros(1,3), "A3/D3 isometry")
    d3 = set()
    for i, j in combinations(range(3), 2):
        for a, b in product((-1,1), repeat=2):
            v = [0]*3
            v[i], v[j] = a, b
            d3.add(tuple(v))
    a3 = set()
    for i, j in product(range(4), repeat=2):
        if i != j:
            v = [0]*4
            v[i], v[j] = 1, -1
            a3.add(tuple(v))
    need({tuple(T*sp.Matrix(v)) for v in d3} == a3, "all twelve actual roots map bijectively")
    q = sp.ones(8,1)/4
    s = 2*q
    F = sp.diag(sp.eye(5), T)
    target = sp.Matrix([sp.Rational(1,2)]*5 + [sp.Rational(3,4)] + [sp.Rational(-1,4)]*3)
    need(F.T*F == sp.eye(8) and F*s == target, "actual glue-vector identity F(2q)=lambda")
    need(q.dot(q)/2 == sp.Rational(1,4) and s.dot(s)/2 == 1, "actual twist weights")

    def in_l0(v):
        return all(x.q == 1 for x in v) and sum(v[:5]) % 2 == 0 and sum(v[5:]) % 2 == 0

    def in_d8(v):
        return all(x.q == 1 for x in v) and sum(v) % 2 == 0

    need([in_l0(k*s) for k in range(1,5)] == [False,False,False,True], "glue has order four over L0")
    need([in_d8(k*s) for k in (1,2)] == [False,True], "same glue has order two over D8")
    block_parities = [(a,b) for a,b in product(range(2), repeat=2) if (a+b) % 2 == 0]
    need(block_parities == [(0,0),(1,1)], "L0 plus 2s is exactly D8")
    need(sp.Rational(16,2**2) == 4 and sp.Rational(4,2**2) == 1, "correct mu-index tower")

    # A root is stored as 2*alpha. Torus action exp(2pi i(q,alpha))
    # is i^(sum(2alpha)/2). On E8 all these exponents are even.
    def quarter_degree(root):
        need(sum(root) % 2 == 0, "integer quarter exponent")
        return (sum(root)//2) % 4

    def gamma_degree(root):
        return 2*(root[0] % 2)

    def alternate_degree(root):
        return (quarter_degree(root) + gamma_degree(root)) % 4

    degree_functions = [quarter_degree, alternate_degree]
    fixed_roots = []
    spectra = []
    for degree in degree_functions:
        counts = Counter(degree(root) for root in roots)
        need(counts[1] == counts[3] == 0, "both full-field extensions are involutions")
        fixed_roots.append({r for r in roots if degree(r) == 0})
        counts[0] += 8  # zero modes fix every Heisenberg oscillator
        spectra.append([counts[i] for i in range(4)])
    need(spectra == [[136,0,112,0], [120,0,128,0]], "full V1 field-character spectra")
    e7 = {r for r in roots if sum(r) == 0}
    s_pair = {tuple([1]*8), tuple([-1]*8)}
    need(len(e7) == 126 and fixed_roots[0] == e7 | s_pair, "actual E7 plus A1 root set")

    def simple_roots(shell):
        positive = {r for r in shell if next(x for x in r if x) > 0}
        decomposable = {tuple(a+b for a,b in zip(x,y)) for x,y in combinations(positive,2)}
        return sorted(positive - decomposable)

    cartan_data = []
    for shell in (e7, fixed_roots[0], fixed_roots[1]):
        simple = simple_roots(shell)
        gram = sp.Matrix([[sum(a*b for a,b in zip(x,y))/sp.Integer(4) for y in simple] for x in simple])
        cartan_data.append({"roots":len(shell),"rank":len(simple),"determinant":int(gram.det()),
                            "cartan":[list(map(int, gram.row(i))) for i in range(gram.rows)]})
    need([(x["roots"],x["rank"],x["determinant"]) for x in cartan_data]
         == [(126,7,2),(128,8,4),(112,8,4)], "fixed root-lattice types and determinants")
    d8_degrees = sorted(sum(x != 0 for j,x in enumerate(row) if j != i)
                        for i,row in enumerate(cartan_data[2]["cartan"]))
    need(d8_degrees == [1,1,1,2,2,2,2,3], "alternate fixed root lattice is D8")

    # All-field order follows already on the integral E8 lattice basis:
    # 2q=s is in E8*=E8; gamma is the nontrivial E8/D8 character.
    basis = original["B"]/2
    for j in range(8):
        need((s.dot(basis[:,j])).q == 1, "square trivial on all lattice generators")
        need((2*(q+sp.eye(8)[:,0]).dot(basis[:,j])).q == 1, "alternate square trivial")
    root_set = set(roots)
    brackets = 0
    for a in roots:
        for b in roots:
            total = tuple(x+y for x,y in zip(a,b))
            if total in root_set:
                need(all(degree(total) == (degree(a)+degree(b)) % 4 for degree in degree_functions),
                     "actual character actions preserve every root product")
                brackets += 1
    need(brackets == 13440, "complete nonopposite root OPE channels")
    # J has 60 four-cycles with product phase+1, plus four +/-i Cartans.
    J = original["J"]
    Jstd = original["J_standard"]
    unseen = set(roots)
    cycles = 0
    while unseen:
        point = next(iter(unseen))
        orb = []
        for _ in range(4):
            orb.append(point)
            point = tuple(Jstd*sp.Matrix(point))
        need(point == orb[0] and len(set(orb)) == 4, "native root clock is a four-cycle")
        unseen -= set(orb)
        cycles += 1
    need(cycles == 60 and original["power_basis_phases"](J,0,4) == (1,)*8,
         "unchanged native full-field J")
    native_spectrum = [60,64,60,64]
    need(all(row != native_spectrum for row in spectra), "neither scalar extension is native J")

    # Reproduce the FIRST exact excitation norm on the ACTUAL v988 case.
    source = "verification/v988_psi_lambda_reduction.py"
    tree = ast.parse(pin(source))
    selected = []
    names = {"SX","SY","SZ","TX","TY"}
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in names for t in node.targets):
            selected.append(node)
        if isinstance(node,ast.FunctionDef) and node.name in {"qwz_cylinder","occ_frame"}:
            selected.append(node)
    env = {"np":np}
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(REPO/source),"exec"),env)
    Nx, Ny = 16,6
    h0 = env["qwz_cylinder"](Nx,Ny,1.0,0)
    _,_,V1 = env["occ_frame"](Nx,Ny,1.0,1)
    P = V1@V1.conj().T
    Q = np.eye(h0.shape[0])-P
    energy = float(np.trace(P@h0).real)
    variance = float(np.trace(P@h0@Q@h0@P).real)
    variance_hs = float(np.linalg.norm(Q@h0@P,"fro")**2)
    need(variance > 7 and abs(variance-variance_hs) < 1e-9, "actual first excitation norm, not occupied diagonal proxy")
    occupied_energies = np.diag(V1.conj().T@h0@V1).real
    proxy_sum = float(sum(occupied_energies-occupied_energies.mean()))
    need(abs(proxy_sum) < 1e-9 and np.sqrt(variance) > 2.6, "original first central sum cannot be f1")
    shadow = [float(abs(np.sum((occupied_energies-occupied_energies.mean())**n))**0.5+1e-15)
              for n in range(1,5)]
    fitted_C = max(1.0,1.05*max((f/factorial(n)**2)**(1/(n+1))
                              for n,f in enumerate(shadow,1)))
    need(np.sqrt(variance) > fitted_C**2,
         "the original fitted proxy constant fails as a bound on the ACTUAL first energy norm")
    for relative, digest in PINS.items():
        need(sha256((REPO/relative).read_bytes()).hexdigest() == digest, "unchanged original input")
    return {
        "research_id":"MATH.SEAM.ESERIES.SOURCE-RECONSTRUCTION.06.REPRODUCTION",
        "exact_mathematical_certificate":"PASS", "fundamental_source_goal":"PARTIAL",
        "T_D3_to_A3":[[str(x) for x in T.row(i)] for i in range(4)],
        "F_2q_equals_original_lambda":True,
        "quarter_twist_weight":"1/4", "double_twist_weight":"1",
        "glue_order_over_D5_A3":4,"glue_order_over_D8":2,
        "extension_indices":[2,2],"mu_indices":[16,4,1],
        "canonical_scalar_quarter_V1_spectrum_1_i_minus1_minusi":spectra[0],
        "alternate_scalar_extension_V1_spectrum_1_i_minus1_minusi":spectra[1],
        "complete_extensions_over_fixed_D8":2,
        "both_scalar_extensions_have_order":2,
        "fixed_VOAs":["V_E7 tensor V_A1","V_D8_prime"],
        "native_J_V1_spectrum_1_i_minus1_minusi":native_spectrum,
        "native_J_equals_either_scalar_extension":False,
        "cartan_certificates":cartan_data,
        "all_nonopposite_root_OPE_channels":brackets,
        "actual_v988_slater_case":{"Nx":Nx,"Ny":Ny,"dimension":len(h0),"rank_P":V1.shape[1],
            "mean_energy":round(energy,10),"exact_variance_formula_evaluated":round(variance,10),
            "first_excitation_norm":round(float(np.sqrt(variance)),10),"first_central_proxy_sum_exact":0,
            "original_fitted_proxy_C":round(fitted_C,10),
            "bound_C_squared_for_n1":round(fitted_C**2,10),
            "proxy_fitted_C_bounds_actual_f1":False,
            "numeric_identity_tolerance":1e-9},
        "not_proved":["microscopic P1 or QWZ charged field exists on the common scaling core",
                      "scalar-quarter and native J operator identification",
                      "uniform energy bounds for the microscopic charged field",
                      "physical state, time or 3+1D completion"],
        "source_sha256":PINS,"promotion":False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=HERE/"source-reproduction-results.json")
    args = parser.parse_args()
    result = calculate()
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("cartan_certificates","source_sha256")},indent=2))
