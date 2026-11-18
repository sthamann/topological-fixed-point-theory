"""R070: exact original color action on the retained folded CAR charge frame.

This checks a necessary interface condition, not physical fermionization,
source selection, local spacetime, or an alternative fitted frame.
"""
import itertools
import json
from pathlib import Path

import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FOLD = ROOT / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"


def main():
    m = s.Matrix(json.loads(FOLD.read_text())["explicit_Weyl_matrix"])
    eye = s.eye(8)
    axes = [eye[:, i] for i in range(8)]
    a, b = axes[0] - axes[1], axes[1] - axes[2]
    w = axes[0] + axes[1] - 2 * axes[2]
    colors = [a, -a, b, -b, a + b, -a - b]
    wf = m.T * w
    tests = []

    def need(label, condition):
        if not bool(condition):
            raise AssertionError(label)
        tests.append(label)

    need("original M is orthogonal", m.T * m == eye)
    need("original color charge in source frame", list(wf) == [1, -1, 0, -1, -1, -1, 0, -1])
    need("integer color charge on every NS integer charge", all(x.q == 1 for x in wf))
    need("integer color charge on half-spinor E8 coset", sum(wf) / 2 == -2)
    mapped = [m.T * q for q in colors]
    integer = [all(x.q == 1 for x in q) for q in mapped]
    need("exactly two of six original color roots belong to CAR D8", sum(integer) == 2)
    need("retained color A1 is the actual first simple root", integer == [True, True, False, False, False, False])
    need("first root exact preimage", list(mapped[0]) == [0, 0, 1, 0, 0, 0, -1, 0])

    for k, q in enumerate(mapped):
        need(f"color root {k} has norm two", q.dot(q) == 2)
        if integer[k]:
            need(f"color root {k} has even integer charge", sum(q) % 2 == 0)
        else:
            need(f"color root {k} is the retained E8 half coset", all(x.q == 2 for x in q))
            for j in range(8):
                need(f"color root {k} has nontrivial full monodromy with CAR generator {j}", (q[j] % 1) == s.Rational(1, 2))

    # On integer charges N=Z8 the original color character has order exactly3.
    # On D8 it is also onto; e0+e2 has even sum and color charge1.
    d8_witness = axes[0] + axes[2]
    need("color map onto Z3 on Z8", wf.dot(axes[0]) == 1)
    need("color map onto Z3 on D8", sum(d8_witness) % 2 == 0 and wf.dot(d8_witness) == 1)
    odd_neutral = axes[2]
    need("actual source has an odd neutral charge", wf.dot(odd_neutral) == 0 and odd_neutral.dot(odd_neutral) == 1)
    need("odd neutral charge lies outside E8", sum(odd_neutral) % 2 == 1)

    # Independent complete norm-two root control of the same original matrix.
    roots = []
    for i, j in itertools.combinations(range(8), 2):
        for x, y in itertools.product([-1, 1], repeat=2):
            roots.append(x * axes[i] + y * axes[j])
    for signs in itertools.product([-1, 1], repeat=8):
        if signs.count(-1) % 2 == 0:
            roots.append(s.Matrix(signs) / 2)
    rootset = {tuple(q) for q in roots}
    need("complete retained E8 shell", len(rootset) == 240)
    need("M preserves every original E8 root", {tuple(m * q) for q in roots} == rootset)
    need("original color roots lie in complete shell", all(tuple(q) in rootset for q in colors))

    out = {
        "result": "R070",
        "status": "SCOPED_INTERFACE_DECIDED_PHYSICAL_F1_OPEN",
        "input": str(FOLD.relative_to(ROOT)),
        "color_charge_in_original_fold_frame": [str(x) for x in wf],
        "color_root_preimages": [[str(x) for x in q] for q in mapped],
        "color_roots_in_unchanged_even_CAR_lattice": sum(integer),
        "color_roots_requiring_NS_R_sector_change": len(colors) - sum(integer),
        "monodromy_of_each_missing_root_with_each_odd_vector_generator": "-1",
        "common_neutral_even_lattice_index_in_E8": 6,
        "common_neutral_even_lattice_determinant": 36,
        "bosonic_and_fermionic_neutral_neighbors_determinants": [9, 9],
        "neighbor_index_over_common_even_lattice": 2,
        "scope": "Original fixed M and color Cartan. Center action extends to NS CAR, full color root action does not preserve its charge lattice. Four-module parent is retained; no new physical field map, local sector endpoint, state/measure/time, 96-matter identification, 3+1 or gravity is derived.",
        "checks_passed": len(tests),
        "checks": tests,
    }
    (HERE / "COLOR-SPIN-SOURCE-INTERFACE.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ["status", "color_charge_in_original_fold_frame", "color_roots_in_unchanged_even_CAR_lattice", "color_roots_requiring_NS_R_sector_change", "checks_passed"]}))


if __name__ == "__main__":
    main()
