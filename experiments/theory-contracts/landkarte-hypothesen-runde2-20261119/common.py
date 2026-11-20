"""Gemeinsame Helfer. Exakte oder hochgenaue Arithmetik, keine angepassten Eingaben."""
import json, sys

class Checks:
    def __init__(self, name): self.name, self.rows = name, []
    def ok(self, label, cond, detail=""):
        self.rows.append({"check": label, "pass": bool(cond), "detail": str(detail)})
        print(("PASS " if cond else "FAIL ") + label + ("  | " + str(detail) if detail != "" else ""))
    def done(self, extra=None):
        n = sum(r["pass"] for r in self.rows)
        print(f"{self.name}: {n}/{len(self.rows)} checks passed")
        json.dump({"name": self.name, "checks": self.rows, "data": extra or {}},
                  open(f"results_{self.name}.json", "w"), indent=1, ensure_ascii=False, default=str)
        sys.exit(0 if n == len(self.rows) else 1)

def e8():
    """E8 in Bourbaki Koordinaten: einfache Wurzeln, 240 Wurzeln, Koeffizienten jeder Wurzel in der einfachen Basis."""
    import itertools, sympy as sp
    h = sp.Rational(1, 2); e = sp.eye(8)
    al = [sp.Matrix([h, -h, -h, -h, -h, -h, -h, h]), e[:, 0] + e[:, 1], e[:, 1] - e[:, 0], e[:, 2] - e[:, 1], e[:, 3] - e[:, 2], e[:, 4] - e[:, 3], e[:, 5] - e[:, 4], e[:, 6] - e[:, 5]]
    roots = []
    for i, j in itertools.combinations(range(8), 2):
        for a in (1, -1):
            for b in (1, -1): v = [0] * 8; v[i] = a; v[j] = b; roots.append(sp.Matrix(v))
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0: roots.append(sp.Matrix([h * x for x in s]))
    Binv = sp.Matrix.hstack(*al).inv()
    coef = [tuple(int(x) for x in Binv * r) for r in roots]
    return al, roots, coef, Binv
