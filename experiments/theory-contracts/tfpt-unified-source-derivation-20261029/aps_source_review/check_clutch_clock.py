"""Independent exact finite check of Claude17, with the original marked D5+D3 split."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json

roots = []
for i, j in combinations(range(8), 2):
    for a, b in product((-2, 2), repeat=2):
        v = [0] * 8
        v[i], v[j] = a, b
        roots.append(tuple(v))
roots += [v for v in product((-1, 1), repeat=8) if v.count(-1) % 2 == 0]
if len(set(roots)) != 240:
    raise RuntimeError("Root reconstruction failed")

lk = (-1, -1, -1, 1, -1, 1, 1, 1)
ld = (0, 0, 0, 2, 0, 2, 0, 0)
native = (1, 1, 1, 1, 1, 0, 0, 0)

def ip(v, w):
    return Fraction(sum(a * b for a, b in zip(v, w)), 4)

def scalar_clock_spectrum(v):
    out = Counter({0: 8})  # Cartan, phase i^0
    for r in roots:
        k = ip(r, v)
        if k.denominator != 1:
            raise RuntimeError("Nonintegral root-cocharacter pairing")
        out[int(k) % 4] += 1
    return dict(sorted(out.items()))

def carrier_degrees(v, j):
    return tuple(Fraction(v[i] + v[j], 2) for i in range(5))

def allowed(v, j):
    d = carrier_degrees(v, j)
    return d[:3] == (0, 0, 0) and sorted(d[3:]) == [0, 1]

def kinds(seq):
    return dict(Counter("D8" if 0 in v else "half" for v in seq))

conjugate64 = [r for r in roots if 0 not in r and r[:5].count(-1) % 2 == 0]
actual64 = [r for r in roots if 0 not in r and r[:5].count(-1) % 2 == 1]
raw_phase = Counter({0: 8})
for r in roots:
    raw_phase[sum(a * b for a, b in zip(r, native)) % 4] += 1

out = {
    "source": "Independent full E8 enumeration, integer roots scaled by2; no physical input fitted",
    "lambda_kappa_is_root": lk in roots,
    "lambda_D8_is_root": ld in roots,
    "lambda_kappa_selected_carrier_j6": list(map(str, carrier_degrees(lk, 5))),
    "lambda_D8_carrier_j6": list(map(str, carrier_degrees(ld, 5))),
    "lambda_D8_carrier_j7": list(map(str, carrier_degrees(ld, 6))),
    "fixed_j6_clutch_roots": kinds([v for v in roots if allowed(v, 5)]),
    "some_j_clutch_roots": kinds([v for v in roots if any(allowed(v, j) for j in (5, 6, 7))]),
    "lambda_kappa_conjugate64_grades_used_by_Claude": dict(sorted(Counter(str(ip(r, lk)) for r in conjugate64).items())),
    "lambda_kappa_actual_exterior_even_64_grades": dict(sorted(Counter(str(ip(r, lk)) for r in actual64).items())),
    "chirality_source": "JOINT-LIFT-NATIVE-CLOCK-CONNECTION.txt:63-87: Lambda_even uses odd minus parity on D5 and paired family -sigma/2",
    "lambda_kappa_full_quarter_i_powers": scalar_clock_spectrum(lk),
    "all_root_full_quarter_spectrum_classes": len({tuple(scalar_clock_spectrum(v).items()) for v in roots}),
    "native_kappa_full_i_powers": dict(sorted(raw_phase.items())),
    "full_quarter_trace": scalar_clock_spectrum(lk)[0] - scalar_clock_spectrum(lk)[2],
    "native_kappa_trace": raw_phase[0] - raw_phase[2],
    "full_loop_exp2pii_lambda": "identity on all248 directions; all root pairings integral",
    "scope": "Root quarter != native kappa by full adjoint spectrum. Separate D5 endpoint kappa remains valid. Fixed marked j versus freely varying j counts are distinct.",
}
print(json.dumps(out, indent=2))
