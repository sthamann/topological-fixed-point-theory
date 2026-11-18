"""Evaluate the original source45 C6 resolvent, not a substitute full P1 source.

Exact symbolic identities, 70-digit original algebraic pole, and a read-only
capability call on the existing original-source consumer and diagnostic domain.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import mpmath as mp
import sympy as s

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
OUT = Path(__file__).resolve().parent
mp.mp.dps = 70
y, delta, tau = s.symbols('y delta tau', positive=True)
U = s.zeros(6)
for j in range(6):
    U[(j + 1) % 6, j] = 1
I = s.eye(6)
L = 2*I-U-U.T
c = (y-delta)**2+tau
A = (y*I-delta*U).T*(y*I-delta*U)+tau*I
assert s.simplify(A-c*I-y*delta*L) == s.zeros(6)
assert L.eigenvals() == {0:1, 1:2, 3:2, 4:1}
one = s.ones(6,1)
assert L*one == s.zeros(6,1)
assert (U-I).rank() == 5
# Check the displayed exact polynomial resolvent of original EM section005.
a = y*y+delta*delta+tau
b = y*delta
C = (U+U.T)/2
denom = (a*a-b*b)*(a*a-4*b*b)
num = (a**3-5*a*b*b)*I+2*b*(a*a-5*b*b)*C+4*a*b*b*C*C+8*b**3*C**3
assert s.expand((A*num-denom*I)) == s.zeros(6)

z = (mp.mpf(794)-7*mp.sqrt(9961))/2187
d = z**(mp.mpf(1)/6)
stationarity = sum(1/(v**6-z) for v in [mp.mpf(1), mp.mpf(1)/3, mp.mpf(2)/3])
assert abs(stationarity) < mp.mpf('1e-60')
rows = []
for label, yy in [('1',mp.mpf(1)),('1/3',mp.mpf(1)/3),('2/3',mp.mpf(2)/3)]:
    aa, bb = yy*yy+d*d, yy*d
    cc = (yy-d)**2
    Am = mp.matrix([[float(U[i,j]) for j in range(6)] for i in range(6)])
    Am = aa*mp.eye(6)-bb*(Am+Am.T)
    Ym = Am**-1
    Km = cc*Ym
    rr = (aa-mp.sqrt(aa*aa-4*bb*bb))/(2*bb)
    closed_row = [(rr**k+rr**(6-k))/(mp.sqrt(aa*aa-4*bb*bb)*(1-rr**6)) for k in range(6)]
    error = max(abs(Ym[0,k]-closed_row[k]) for k in range(6))
    assert error < mp.mpf('1e-60')
    assert min(Km) > 0
    assert max(abs(sum(Km[i,j] for j in range(6))-1) for i in range(6)) < mp.mpf('1e-60')
    rows.append(dict(y=label, c=mp.nstr(cc,35),
        Y_first_row=[mp.nstr(Ym[0,k],35) for k in range(6)],
        K_first_row=[mp.nstr(Km[0,k],35) for k in range(6)],
        K_spectrum=[dict(value=mp.nstr(cc/(cc+bb*ell),35),multiplicity=mult)
                    for ell,mult in [(0,1),(1,2),(3,2),(4,1)]],
        formula_error=mp.nstr(error,6)))

base = ROOT/'experiments/theory-contracts/original-paper-source-reconciliation-20261102/fieldsearch_20261106P'
spec = importlib.util.spec_from_file_location('original_source_evaluator',base/'source/evaluator.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
domain = json.loads((base/'domain/native-diagnostic-v1.json').read_text())
cap = mod.capability_report(domain)
contrib = mod.contribution_capability_report(domain)
assert not cap['numeric_bundle_ready'] and not cap['whole_source_bundle_supplied']
# Do not fabricate candidate/assembled bindings merely to call evaluate_source.
# The actual domain capability already decides whether it can evaluate Q0.
paths = [
 'documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex',
 'documents/archive--tfpt-45--source_extracts--03_em_flavor_source/sections/003-subsection-discrete-phase-operator-and-transport-grammar.tex',
 'documents/archive--tfpt-45--source_extracts--03_em_flavor_source/sections/004-source.tex',
 'documents/archive--tfpt-45--source_extracts--03_em_flavor_source/sections/005-source.tex',
 'documents/archive--tfpt-45--source_extracts--04_qft_source/sections/004-source.tex',
 'verification/v774_arf_spinor_compiler.py',
 str((base/'source/evaluator.py').relative_to(ROOT)),
 str((base/'domain/native-diagnostic-v1.json').relative_to(ROOT)),
]
result = dict(scope='Exact original local C6 and primitive selector; not whole P1 integral',
 delta=mp.nstr(d,55), z=mp.nstr(z,55), transport_stationarity_residual=mp.nstr(stationarity,6),
 exact_checks=dict(shifted_inverse_is_original_C6_Laplacian=True,
     original_polynomial_resolvent=True, centered_C6_no_local_zero=True,
     centered_eigenvalue_multiplicities=[2,2,1]),
 sectors=rows, original_source_capability=cap, original_fiber_capability=contrib,
 actual_domain_role=domain.get('input_role'), physical_source_identity=domain.get('physical_source_identity'),
 source_pins={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
(OUT/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(delta=result['delta'], exact_checks=result['exact_checks'],
 sectors=[dict(y=r['y'],K_first_row=r['K_first_row'],K_spectrum=r['K_spectrum']) for r in rows],
 whole_source_ready=cap['numeric_bundle_ready'],missing_operation=cap['missing_operations']),indent=2))
