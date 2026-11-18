"""Narrow type check on the EXISTING v252 constructor; no new physical source.

Does not call its run(), modify its ledger, or invent spacetime matrices.
The general sign proof and the C* state theorem are in PROOF.txt.
"""
from pathlib import Path
import json
import sys

import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "verification"))
import v252_full_finite_triple as original


def exact_gaussian(matrix):
    rows = []
    for row in matrix:
        values = []
        for z in row:
            z = complex(z)
            if not (z.real.is_integer() and z.imag.is_integer()):
                raise ValueError("Control inputs must give Gaussian integer entries")
            values.append(s.Integer(int(z.real)) + s.I * int(z.imag))
        rows.append(values)
    return s.Matrix(rows)


def main():
    # Deterministic Gaussian integer data exercise the original complex Yukawa
    # and symmetric Majorana construction. They are not a physical calibration.
    Y = np.array([[1+1j, 2-1j, 0], [0, 3+2j, -1j], [2, 0, 4-1j]])
    MR = np.array([[2, 1+1j, 0], [1+1j, 3, 2j], [0, 2j, 5]])
    D = exact_gaussian(original.full_D(Y, 2*Y, 3*Y, 4*Y, MR))
    U = exact_gaussian(original.Kswap())
    G = exact_gaussian(original.grading())
    B = U.conjugate().T * D
    plus = [i for i in range(original.NH) if G[i, i] == 1]
    Bplus = B.extract(plus, plus)
    UG = G * U
    BG = UG.conjugate().T * D
    checks = {
        "original_dimension_96": D.shape == (96, 96),
        "selfadjoint": D == D.conjugate().T,
        "J_square_plus": U * U.conjugate() == s.eye(96),
        "J_commutes_D": U * D.conjugate() == D * U,
        "J_anticommutes_grading": U * G.conjugate() == -G * U,
        "D_anticommutes_grading": D * G == -G * D,
        "finite_bilinear_symmetric": B.T == B,
        "chiral_bilinear_nonzero_but_symmetric": Bplus != s.zeros(48) and Bplus.T == Bplus,
        "grassmann_quadratic_coefficients_zero": Bplus - Bplus.T == s.zeros(48),
        "grading_modified_J_square_minus": UG * UG.conjugate() == -s.eye(96),
        "grading_modified_J_anticommutes_D": UG * D.conjugate() == -D * UG,
        "grading_modified_bilinear_still_symmetric": BG.T == BG,
    }
    # Standard real even KO signs. This is a type census, not a metric dimension
    # derivation: a nonzero chiral skew pairing needs both columns to pass.
    signs = {0: (1, 1, 1), 2: (-1, 1, -1), 4: (-1, 1, 1), 6: (1, 1, -1)}
    census = [dict(ko=k, epsilon=e, epsilon_prime=ep, epsilon_doubleprime=epp,
                   skew=e*ep == -1, preserves_chiral_pairing=epp == -1)
              for k, (e, ep, epp) in signs.items()]
    checks["only_KO2_has_both_required_signs"] = [
        r["ko"] for r in census if r["skew"] and r["preserves_chiral_pairing"]
    ] == [2]
    # A second output of that SAME constructor has simple Majorana eigenvalues.
    # Kramers' argument then excludes EVERY commuting square-minus antiunitary
    # on that output, not only the two explicitly tried J choices.
    args = [np.diag(x) for x in ([11,13,17], [19,23,29], [31,37,41],
                                 [1,2,3], [1,4,9])]
    DM = exact_gaussian(original.full_D(*args))
    lam = s.Symbol("lambda")
    characteristic = DM.charpoly(lam).as_expr()
    q = lam**4 - 3*lam**2 + 1
    indices = [original.IDX[("nu", chirality, 0, 0)] + offset
               for offset in (0, original.NP) for chirality in ("L", "R")]
    block = DM.extract(indices, indices)
    checks["actual_neutrino_block_has_q_characteristic"] = block.charpoly(lam).as_expr() == q
    checks["q_occurs_in_full_96_characteristic"] = s.rem(characteristic, q, lam) == 0
    checks["all_four_q_eigenvalues_are_simple_in_full_96"] = s.gcd(q, s.diff(characteristic, lam)) == 1
    Q = s.diag(*([1]*48 + [-1]*48))
    Ualt = Q * U
    args[-1] = np.zeros((3, 3))
    Dzero = exact_gaussian(original.full_D(*args))
    checks["no_majorana_control_allows_commuting_square_minus_J"] = (
        Ualt * Ualt.conjugate() == -s.eye(96)
        and Ualt * Dzero.conjugate() == Dzero * Ualt
        and Ualt * G.conjugate() == -G * Ualt)
    checks["same_J_fails_commutation_with_majorana_on"] = Ualt * DM.conjugate() != DM * Ualt
    for name, passed in checks.items():
        if not passed:
            raise RuntimeError(name)
    result = {
        "status": "PASS_SCOPED",
        "checks": checks,
        "chiral_dimension": len(plus),
        "chiral_nonzero_bilinear_entries": sum(x != 0 for x in Bplus),
        "ko_sign_census": census,
        "majorana_type_control": {
            "neutrino_block": str(block),
            "q": str(q),
            "gcd_q_and_full_characteristic_derivative": str(s.gcd(q, s.diff(characteristic, lam))),
            "scope": "An actual v252 constructor output, not the selected physical TFPT masses. A universal internal KO2 repair is excluded; no all-P1 or spacetime theorem.",
        },
        "scope": "Existing v252 finite constructor only; general identities proved in text. No original local P1 field action, spacetime domain, measure, reflection or full TFPT derived.",
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
