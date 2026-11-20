"""Was beide singulaeren Punkte gemeinsam fest lassen (exakt)."""
import json, collections
from fractions import Fraction as Fr
import check_orbifold_zero_modes as Z
R = Z.R; T = {v: t for t, v in R}
h4 = Z.h
out = {}
for name, d3 in (("Hebung Ordnung 6", [Fr(0), Fr(1, 3), Fr(-1, 3), Fr(0), Fr(1, 3)]),
                 ("Hebung Ordnung 3", [Fr(0), Fr(1, 3), Fr(2, 3), Fr(0), Fr(1, 3)])):
    h3 = d3 + [Fr(1, 3), Fr(-1, 3), Fr(0), Fr(0)]
    both = [v for _, v in R if Z.grade(v, h4, 8) == 0 and Z.grade(v, h3, 6) == 0]
    comps = Z.components(both)
    types = sorted(Z.typ(c) for c in comps)
    src = [dict(collections.Counter(T[v] for v in c)) for c in comps]
    # alle Lagen des Dreierzyklus relativ zum Viererzyklus (welche Marke fest bleibt, Richtung)
    variants = collections.Counter()
    import itertools
    for perm in itertools.permutations(range(4)):
        a = [[Fr(1, 3), Fr(-1, 3), Fr(0), Fr(0)][perm[i]] for i in range(4)]
        for dperm in itertools.permutations(range(5)):
            dd = [d3[dperm[i]] for i in range(5)]
            b = [v for _, v in R if Z.grade(v, h4, 8) == 0 and Z.grade(v, dd + a, 6) == 0]
            variants["+".join(sorted(Z.typ(c) for c in Z.components(b))) or "nur Torus"] += 1
    out[name] = {"gemeinsam fixiert": types, "wurzeln": len(both), "herkunft": src, "alle Lagen": dict(variants)}
    print(name, types, len(both), src); print("   Lagen:", dict(variants))
json.dump(out, open("check_both_points.json", "w"), indent=1, ensure_ascii=False)
