"""Shared helpers. Exact/high-precision arithmetic only; no fitted inputs."""
import json, sys
from mpmath import mp, mpf, pi, exp, log, findroot
mp.dps = 40
C3 = 1 / (8 * pi)
DTOP = 48 * C3**4
PHI_TREE = 1 / (6 * pi)
CODATA, SIG = mpf("137.035999177"), mpf("0.000000021")   # CODATA 2022, quoted in v2.2..v2.8
RB, RB_SIG = mpf("137.035999206"), mpf("0.000000011")    # Rb recoil 2020 (cited from memory)

def alpha_inv(phi_of_alpha):
    F = lambda a: a**3 - 2 * C3**3 * a**2 - mpf(4) / 5 * 41 * C3**6 * log(1 / phi_of_alpha(a))
    return 1 / findroot(F, 1 / mpf(137))

def z(x, ref=CODATA, s=SIG):
    return (x - ref) / s

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
