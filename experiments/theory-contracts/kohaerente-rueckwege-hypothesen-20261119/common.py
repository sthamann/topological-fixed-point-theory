"""Shared helpers. Exact or high-precision arithmetic; no fitted inputs."""
import json, sys
from itertools import product
from mpmath import mp, mpf, pi, exp, log, findroot
mp.dps = 40
C3 = 1 / (8 * pi)
DTOP = 48 * C3**4
PHI_TREE = 1 / (6 * pi)
CODATA, SIG = mpf("137.035999177"), mpf("0.000000021")   # CODATA 2022 (as in ORIGINIDEAS20261119)


def alpha_inv(damping):
    """Root of the fixed alpha equation with q = DTOP * damping(alpha)."""
    def phi(a):
        x = DTOP * damping(a)
        return PHI_TREE + x * (1 - x)**(-mpf(5) / 4)
    F = lambda a: a**3 - 2 * C3**3 * a**2 - mpf(4) / 5 * 41 * C3**6 * log(1 / phi(a))
    return 1 / findroot(F, 1 / mpf(137))


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


# ---- D8 half-spinor in the CAR chart of R113-FOUR-FAMILY-TRANSPORT-20261116 ----
# eight modes: 0,1,2 colour, 3,4 weak, 5,6,7 family; states = even occupation masks
STATES = tuple(s for s in range(256) if bin(s).count("1") % 2 == 0)
INDEX = {s: i for i, s in enumerate(STATES)}


def letter(state, mode, create):
    if (state >> mode) & 1 == create:
        return None
    return state ^ (1 << mode), (-1)**bin(state & ((1 << mode) - 1)).count("1")


def root_action(i, j, a, b):
    """Signed CAR two-letter word: letter b on mode j first, then letter a on mode i."""
    out = {}
    for s in STATES:
        f = letter(s, j, b)
        if f is None: continue
        g = letter(f[0], i, a)
        if g is not None: out[s] = (g[0], f[1] * g[1])
    return out


EVEN_ROOTS = tuple((i, j, a, b) for i in range(8) for j in range(i + 1, 8) for a, b in product((0, 1), repeat=2))
MIXED = tuple(r for r in EVEN_ROOTS if r[0] < 5 <= r[1])
occ = lambda s, modes: sum((s >> m) & 1 for m in modes)
KLF = lambda s: (occ(s, (0, 1, 2)), occ(s, (3, 4)), occ(s, (5, 6, 7)))
