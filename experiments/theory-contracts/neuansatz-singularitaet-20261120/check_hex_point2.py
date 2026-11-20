"""Zweite Spinhebung der Dreierwirkung (Ordnung 6) und Inhalt der 84 unter SU(5) x SU(4)."""
import json, itertools, collections
from fractions import Fraction as Fr
import check_hex_point as Hx, check_orbifold_zero_modes as Z
d6 = [Fr(0), Fr(1, 3), Fr(-1, 3), Fr(0), Fr(1, 3)]
hexa6 = Hx.analyse(d6 + Hx.a3, 6, (1, 5, 3, 2, 4))
print("hexagonal, Hebung Ordnung 6", hexa6["graddimensionen"], hexa6["fixierte_algebra"], "u1:", hexa6["u1"])
for k, lst in hexa6["grade"].items():
    print("  Grad", k, [(x["dims"], x["kubisch"], x["aus16"], x["aus16bar"]) for x in lst])
# Inhalt der 84 (Hebung Ordnung 3) nach Herkunft
h3 = Hx.d3 + Hx.a3
src = collections.Counter(t for t, v in Z.R if Z.grade(v, h3, 6) == 4)
print("84 nach Herkunft:", dict(src))
# Familienzahl in der 84: Zustaende aus (16,4) nach Familienindex und SM Art
fam = collections.Counter()
for t, v in Z.R:
    if Z.grade(v, h3, 6) == 4 and t in ("16x4", "16bx4b"):
        bar = t == "16bx4b"
        idx = tuple(v[5:]).index(max(v[5:])) if not bar else tuple(v[5:]).index(min(v[5:]))
        fam[(t, idx)] += 1
print("Verklebungsanteil der 84 nach Familienindex:", dict(sorted(fam.items())))
json.dump({"hexagonal_ordnung6": hexa6, "84_herkunft": dict(src), "84_familien": {str(k): n for k, n in fam.items()}},
          open("check_hex_point2.json", "w"), indent=1, ensure_ascii=False, default=str)
