"""Exact charge, centre and affine-grade typecheck for the native Gamma/Hodge bridge."""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORK = HERE.parent
REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")

GAMMA_CHECKER = WORK / "current" / "check_native_gamma_current.py"
GAMMA_RESULT = WORK / "current" / "native-gamma-current-results.json"
CURRENT_DESCENDANT = REPO / "experiments/theory-contracts/compiler-current-descendant-20260919/PROOF.txt"
ROOT_PHASE = REPO / "experiments/theory-contracts/compiler-root-phase-lift-20260919/PROOF.txt"
A4_RESPONSE = REPO / "experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/hodge_complement_20261028/E8-A4-PETERSEN-RESPONSE.txt"

PINS = {
    GAMMA_CHECKER: "e00ac829026c15d949672601bd5971187a767a23a25424abfbc5defc72be44b8",
    GAMMA_RESULT: "5b8ee6ccab21c3cdb02c55fe12d42359ce73a623f00336f7e999875e4d22ea4b",
    CURRENT_DESCENDANT: "ad45198c44b84365a862b709e8b7e49af708f2e1a58adff67dc55d7cf4fa459b",
    ROOT_PHASE: "d7465bd8ee3514d29bfcb4238b4772d7d2ee631301465682e50bcbfca5b4879f",
    A4_RESPONSE: "8311a51e2f9a0e5e9eb9307fc329cdce77e5d65a599da3fa5692b7460ec63a38",
}

checks = []


def need(value, label):
    if not bool(value):
        raise RuntimeError(label)
    checks.append(label)


def qstr(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def spectrum(values):
    return {qstr(key): count for key, count in sorted(Counter(values).items())}


def poly_from_roots(values):
    # Coefficients from highest degree to constant term.
    coeff = [Q(1)]
    for root in values:
        new = [Q(0)] * (len(coeff) + 1)
        for index, value in enumerate(coeff):
            new[index] += value
            new[index + 1] -= root * value
        coeff = new
    return [qstr(value) for value in coeff]


for path, expected in PINS.items():
    need(path.is_file(), f"source exists: {path.name}")
    need(sha256(path.read_bytes()).hexdigest() == expected, f"source pin: {path.name}")

gamma = json.loads(GAMMA_RESULT.read_text())
need(gamma["gamma_names"] == ["XI", "YX", "YY", "YZ", "ZI"], "actual five Gamma directions")
need(gamma["verdict"] == "NATIVE_HODGE_EQUALS_FIVE_ACTUAL_COLLECTIVE_PAULI_CURRENT_ACTIONS", "actual Gamma/Hodge operator identity")

descendant_text = CURRENT_DESCENDANT.read_text()
phase_text = ROOT_PHASE.read_text()
a4_text = A4_RESPONSE.read_text()
need("10_D5 tensor E(Bell10_A3)" in descendant_text and "conformal weight1/2+3/2=2" in descendant_text, "old Bell corner is E8 weight two")
need("J_0(X) E = E [X tensor I+I tensor X]|_Sym²" in descendant_text, "Gamma acts as an A3 horizontal zero mode")
need("The diagonal Z4 glue acts trivially: (-1)*(-1)=+1." in descendant_text, "old corner has glue classes two plus two")
need("Die bestehende P2-Einbettung zerlegt den D5-Vektor" in phase_text, "original P2 vector charge statement present")
need("Eingabe:     (10,bar5)" in a4_text and "Einsatz:     (5,10)" in a4_text and "Ausgabe:     (bar10,5)" in a4_text, "existing charged A4 by A4 current route present")
need("q_output + q_env = q_input" in a4_text, "existing mixed-current route has explicit charge balance")

# Original carrier hypercharge on the five D5 vector coordinates.
y5 = [Q(-1, 3)] * 3 + [Q(1, 2)] * 2
y_vector = y5 + [-value for value in y5]
y_wedge = [y5[i] + y5[j] for i in range(5) for j in range(i + 1, 5)]

need(spectrum(y_vector) == {"-1/2": 2, "-1/3": 3, "1/3": 3, "1/2": 2}, "D5 vector hypercharge spectrum")
need(spectrum(y_wedge) == {"-2/3": 3, "1/6": 6, "1": 1}, "spinor Lambda2(5) hypercharge spectrum")

shared_y = sorted(set(y_vector).intersection(y_wedge))
hom_y_dimension = sum(Counter(y_vector)[value] * Counter(y_wedge)[value] for value in shared_y)
need(shared_y == [] and hom_y_dimension == 0, "Hom_Y(vector10,Lambda2(5)) is zero")
need(poly_from_roots(y_vector) != poly_from_roots(y_wedge), "hypercharge characteristic polynomials differ")

# Spin(10) centre classes are recorded additively modulo four.
# vector10 is class 2, while the two half-spinors are classes 1 and 3.
old_d5_class = 2
old_a3_class = 2
spinor_classes = [1, 3]
paired_a3_classes = [3, 1]
need((old_d5_class + old_a3_class) % 4 == 0, "old vector by descendant corner is glue neutral")
need(all((d5 + a3) % 4 == 0 for d5, a3 in zip(spinor_classes, paired_a3_classes)), "mixed spinor by fundamental corners are glue neutral")
need(all(old_d5_class != value for value in spinor_classes), "D5 centre character differs for vector and spinor")

# Affine conformal weights in the original level-one embedding.
h_old_d5 = Q(1, 2)
h_old_a3_descendant = Q(3, 2)
h_spinor_d5 = Q(5, 8)
h_fundamental_a3 = Q(3, 8)
need(h_old_d5 + h_old_a3_descendant == 2, "old Bell corner has E8 weight two")
need(h_spinor_d5 + h_fundamental_a3 == 1, "mixed spinor corner has E8 weight one")
need(h_old_d5 + h_old_a3_descendant != h_spinor_d5 + h_fundamental_a3, "zero mode cannot identify the two affine grades")

result = {
    "verdict": "ABSTRACT_HODGE_CHANNEL_IDENTITY_ONLY_IN_OLD_BELL_CHART__NO_UNCHANGED_P2_CARRIER_IDENTIFICATION",
    "actual_gamma_directions": gamma["gamma_names"],
    "old_bell_corner": {
        "module": "(10_D5_vector, 10_A3_descendant)",
        "e8_conformal_weight": "2",
        "mu4_classes": [2, 2],
        "d5_hypercharge_spectrum": spectrum(y_vector),
    },
    "desired_p2_carrier_corner": {
        "module": "Lambda2(5)_SU5 inside 16_D5_spinor, paired with an A3 fundamental in E8",
        "e8_conformal_weight": "1",
        "mu4_class_options": [[1, 3], [3, 1]],
        "carrier_hypercharge_spectrum": spectrum(y_wedge),
    },
    "decisive_intertwiner_test": {
        "shared_hypercharges": [qstr(value) for value in shared_y],
        "dimension_Hom_Y": hom_y_dimension,
        "vector_characteristic_polynomial_coefficients": poly_from_roots(y_vector),
        "wedge_characteristic_polynomial_coefficients": poly_from_roots(y_wedge),
    },
    "meaning": "The five Gamma zero modes act in the A3 Bell descendant and leave the D5 vector charge, Spin(10) centre class and L0 grade unchanged. They therefore realize the Hodge channel only as an abstract operator identity in that old Bell chart. They do not furnish the unchanged physical P2 carrier readout.",
    "existing_positive_route": "The already verified A4 by A4 E8 bracket (5,10) x (10,bar5) -> (bar10,5) uses charged mixed currents, has explicit q_output+q_env=q_input on all 600 edges, and yields the Hodge carrier response after the stated inclusive family trace.",
    "scope": "This excludes only the neutral family-current identification with the unchanged D5/P2 charge dictionary. The positive A4 by A4 route is a different charged cut through the same E8 algebra, not an identification of the old D5 vector with Lambda2(5). Physical selection of its 50 directions remains open.",
    "checks": len(checks),
    "check_names": checks,
    "source_sha256": {str(path): expected for path, expected in PINS.items()},
}

output = json.dumps(result, indent=2, sort_keys=True) + "\n"
(HERE / "charge-typed-hodge-bridge-results.json").write_text(output)
print(output, end="")
