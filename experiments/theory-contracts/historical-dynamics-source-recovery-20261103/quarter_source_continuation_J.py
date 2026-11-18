"""Direct coefficient check of the existing TFPT v252 producer; no new model."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
ROOT = Path(os.environ['TFPT_REPO']) if os.environ.get('TFPT_REPO') else next(
    (p for p in Path(__file__).resolve().parents
     if (p / 'verification/v252_full_finite_triple.py').is_file()), None)
if ROOT is None:
    raise RuntimeError('Run inside the TFPT repository or set TFPT_REPO.')
sys.path.insert(0, str(ROOT / 'verification'))
sys.path.insert(0, str(ROOT))
import numpy as np
import v252_full_finite_triple as original
from tfpt_explorer.process import _code_basis

gamma = original.grading()
g = np.diag(gamma)
u = 1j * g
roots = [1, 1j, -1, -1j]

def conjugate_by_quarter(matrix):
    return u[:, None] * matrix * u.conjugate()[None, :]

phase = []
charge_rows = []
for name, chirality, color, generation in original.P:
    bl = F(-1) if name in ('nu', 'e') else F(1, 3)
    if chirality == 'L':
        y = F(-1, 2) if name in ('nu', 'e') else F(1, 6)
    else:
        y = {'nu': F(0), 'e': F(-1), 'u': F(2, 3), 'd': F(-1, 3)}[name]
    x = 5 * bl - 4 * y
    assert x.denominator == 1
    phase.append(roots[int(x) % 4])
    charge_rows.append([name, chirality, color, generation, str(y), str(bl), str(x)])
phase = np.asarray(phase + [z.conjugate() for z in phase])
assert np.array_equal(phase, u)
assert np.array_equal(u * u, -np.ones(original.NH))
K = original.Kswap()
assert np.array_equal(K * u.conjugate()[None, :], u[:, None] * K)

max_quarter_defect = 0.0
max_fixed_mr_formula_defect = 0.0
max_reference_transport_defect = 0.0
max_original_quaternion_orbit_defect = 0.0
basis_count = 0
fixed_mr_nonzero_defects = 0
uR = original.rep_PS(np.eye(2), np.diag([1j, -1j]), np.eye(3))
UR = uR @ K @ uR.conjugate() @ K
for slot in range(5):
    indices = list(combinations(range(3), 2)) + [(i, i) for i in range(3)] if slot == 4 else [(i, j) for i in range(3) for j in range(3)]
    for i, j in indices:
        for amplitude in (1, 1j):
            inputs = [np.zeros((3, 3), complex) for _ in range(5)]
            inputs[slot][i, j] = amplitude
            if slot == 4:
                inputs[slot][j, i] = amplitude
            D = original.full_D(*inputs)
            transformed = conjugate_by_quarter(D)
            target = original.full_D(*[-x for x in inputs])
            error = float(np.max(np.abs(transformed - target)))
            max_quarter_defect = max(max_quarter_defect, error)
            fixed_mr_target = original.full_D(*[-x for x in inputs[:4]], inputs[4])
            zero_y = [np.zeros((3, 3), complex) for _ in range(4)]
            majorana = original.full_D(*zero_y, inputs[4])
            remainder = transformed - fixed_mr_target
            e = float(np.max(np.abs(remainder + 2 * majorana)))
            max_fixed_mr_formula_defect = max(max_fixed_mr_formula_defect, e)
            fixed_mr_nonzero_defects += int(np.any(remainder != 0))
            assert error == 0 and e == 0
            # The actual R032 S reference is A gamma; T equals h_A on
            # the existing seam-even chart. The A term commutes with U.
            h_mass = -1j * g[:, None] * D
            h_mass_target = -1j * g[:, None] * target
            e_ref = float(np.max(np.abs(conjugate_by_quarter(h_mass) - h_mass_target)))
            max_reference_transport_defect = max(max_reference_transport_defect, e_ref)
            # Reuse the existing right-quaternion gauge path, not a new sigma.
            r_inputs = [-1j * inputs[0], 1j * inputs[1], 1j * inputs[2], -1j * inputs[3], -inputs[4]]
            e_ps = float(np.max(np.abs(UR @ D @ UR.conjugate().T - original.full_D(*r_inputs))))
            max_original_quaternion_orbit_defect = max(max_original_quaternion_orbit_defect, e_ps)
            assert e_ref == 0 and e_ps == 0
            basis_count += 1

V, words = _code_basis()
quartic_J = (1j ** 4) * V
assert np.array_equal(quartic_J, V)
native_five_J = V.T @ quartic_J
assert np.max(np.abs(native_five_J - np.eye(5))) < 1e-14
# The Spin(10) central generator acts as -I on the vector C5 chart.
# S I = -I S implies 2 S=0, so the complete intertwiner space is {0}.
native_five_to_spin_center_intertwiner_dimension = 0
spinor_weights = [5 - 2 * mask.bit_count() for mask in range(32) if mask.bit_count() % 2 == 0]
assert len(spinor_weights) == 16 and all(x % 4 == 1 for x in spinor_weights)
# Original determinant extension, PROOF sections3-4. This is reuse.
for e in (0, 1):
    for f in (0, 1):
        assert F(-e-f, 4) == F(-e*f, 2) - F(e+f-2*e*f, 4)

data = {
    'result_id': 'QUARTER20261104J',
    'producer': 'verification/v252_full_finite_triple.py',
    'producer_sha256': hashlib.sha256((ROOT / 'verification/v252_full_finite_triple.py').read_bytes()).hexdigest(),
    'full_internal_dimension': original.NH,
    'real_linear_parameter_basis': basis_count,
    'quarter_equals_i_original_grading': True,
    'quarter_commutes_with_actual_antiunitary_J': True,
    'quarter_squared': '-I on the one-particle coefficient space',
    'quarter_coefficient_defect': max_quarter_defect,
    'fixed_MR_remainder': '-2 D_MR',
    'fixed_MR_remainder_formula_defect': max_fixed_mr_formula_defect,
    'nonzero_MR_basis_defects': fixed_mr_nonzero_defects,
    'actual_R032_reference_quarter_transport_defect': max_reference_transport_defect,
    'existing_right_quaternion_orbit_defect': max_original_quaternion_orbit_defect,
    'existing_right_quaternion_is_not_U_X': not np.array_equal(UR, np.diag(u)),
    'original_holonomy_lift_reused': 'compiler-quartic-holonomy-20260919/PROOF.txt sections3-4',
    'original_reflection_lift_relation': 'L(r)^2=G^-1; L(r)^8=1, versus raw r^2=1',
    'joint_16_times_4_center_phase': 'i from the A3 partner for the section g=J',
    'native_quartic_clock_on_five': '+I5',
    'spin10_central_generator_on_vector_five': '-I5',
    'direct_native_vector_intertwiner_dimension': native_five_to_spin_center_intertwiner_dimension,
    'carrier_spinor_X_values': sorted(set(spinor_weights)),
    'physical_global_symmetry_proved': False,
    'complete_original_source_evaluated': False,
    'S_adm_H_tr_pos_and_common_spectrum_derived': False,
    'scope': 'Exact real-linear coefficient comparison on original full96 producer and actual native quartic code; no new field measure or physical identification.',
}
out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / 'quarter-source-results-J.json'
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(data, ensure_ascii=False, indent=2))
