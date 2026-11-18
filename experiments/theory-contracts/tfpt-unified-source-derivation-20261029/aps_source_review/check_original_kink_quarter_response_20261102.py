"""Original five-state Schur source, neutral Weyl line and quarter response.

No P1 field identification. The saved source states and CAR responses are read,
not reconstructed from a chosen or fitted source model.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import sys

import sympy as s

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / 'experiments/theory-contracts/joint-native-algebra-selection-20261016'
sys.path.insert(0, str(SOURCE))
import check_original_current_schur_kernel as original

hall_file = SOURCE / 'original_current_hall_kernel_results.json'
car_file = SOURCE / 'original_density_current_results.json'
hall = json.loads(hall_file.read_text())
car = json.loads(car_file.read_text())
expected = {
 '1/4': [Fraction(328635469,1048576), Fraction(2994963,4096), Fraction(189873,512), Fraction(174393,4096), Fraction(71417,65536)],
 '3/4': [Fraction(331641421,1048576), Fraction(3008787,4096), Fraction(194481,512), Fraction(215865,4096), Fraction(97529,65536)],
}
ratios = {}
for alpha in ('1/4','3/4'):
    case = hall['cases'][alpha]
    assert case['weighted_degree'] == 32
    state = {tuple(x['partition']): original.Q.from_sympy(s.sympify(x['coefficient']))
             for x in case['schur_coefficients'][0]}
    norm = original.Q.to_sympy(original.car.inner(state,state))
    assert norm == s.Rational(str(expected[alpha][0]))
    assert norm == s.sympify(case['normal_ordered_first_mode_jet'][0]['raw_gram'][0][0])
    values = []
    for index, mode in enumerate((8,16,24,32),1):
        lowered = original.lower(state, mode)
        value = original.Q.to_sympy(original.car.inner(lowered,lowered))
        assert value == s.Rational(str(expected[alpha][index]))
        # Previously saved, independently constructed original CAR occupation basis.
        assert value == s.sympify(car['cases'][alpha]['all_positive_density_mode_responses'][str(mode)]['raw_gram'][0][0])
        values.append(value/norm)
    ratios[alpha] = values
    print(alpha, 'norm:', norm, 'R_8,16,24,32:', *values)

delta = s.factor(ratios['1/4'][0]-ratios['3/4'][0])
assert delta == s.Rational(1141672523857920,108989133930161449) and delta > 0
# Verify the one-mode BCH exponent and the energy-order identity symbolically.
b, t, n = s.symbols('b t n', real=True, positive=True)
z = s.symbols('z', complex=True) # z=e^{int}, conjugate(z)=1/z
phi = -b/z
psi = b
comm = 8*n*(s.conjugate(phi)*psi - phi*s.conjugate(psi))
comm = comm.subs(s.conjugate(z),1/z)
d = phi + psi
gauss = -4*n*s.expand_complex(d*s.conjugate(d)).subs({s.re(z): (z+1/z)/2, s.im(z): (z-1/z)/(2*s.I)})
# Direct analytic identity: |1-z^-1|^2=2-z-z^-1.
exponent = -4*n*b*b*(2-z-1/z) - comm/2
assert s.simplify(exponent + 8*n*b*b*(1-z)) == 0
print('R8 difference:',delta)
print('Weyl BCH and CAR versus Hall:', 'PASS')
for path in (hall_file,car_file,SOURCE/'check_original_current_schur_kernel.py',SOURCE/'ORIGINAL-ORDERED-CURRENT-KERNEL-AND-SEAM-GRADING.txt'):
 print('sha256',hashlib.sha256(path.read_bytes()).hexdigest(),path.relative_to(ROOT))
