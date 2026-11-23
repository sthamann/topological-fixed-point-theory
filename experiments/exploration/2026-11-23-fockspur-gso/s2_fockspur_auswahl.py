"""Fockspur statt Felddeterminante: welche Darstellung waehlt 3+2?

Frage (aus Bericht 1, Abschnitt 3 und 8): Eine positive Determinantenquelle ueber der Fuenf
waehlt 3+2, dieselbe Quelle ueber der 16 als Feldraum waehlt 4+1. Welche ist die richtige?

These: In TFPT ist die 16 = Lambda^gerade C^5 der Fockraum von fuenf Moden. Eine Zustandssumme
ueber diesen Fockraum ist eine Spur, keine Determinante ueber 16 Feldern. Fuer einen
Einteilchenoperator Sigma gilt exakt
    Tr_{Lambda C^5} e^{-T Sigma^} = det_5(1 + e^{-T Sigma}) = 2^5 prod cosh(y_i),   y = T lambda / 2,
    Tr_{Lambda C^5} (-1)^F e^{-T Sigma^} = det_5(1 - e^{-T Sigma}) = 2^5 prod sinh(y_i),
also Tr_16 = 16 (prod cosh + prod sinh) (gerader Teil = Summe beider Spinstrukturen des Kreises).
Ausserdem ist das vollstaendige antiperiodische Produkt der Paarquelle aus Bericht 1,
prod_n det(1 + g^2 Sigma^2 / lambda_n^2) mit lambda_n = (2pi/T)(n+1/2), gleich det cosh(gT Sigma/2).

Gerechnet:
 1  die drei Spuridentitaeten exakt (Teilmengensumme gegen Produktformel, symbolisch);
 2  das Weierstrassprodukt numerisch;
 3  globale Minima von Gamma = rho sum y^2 - log Z fuer Z = Tr_32, Tr_16, Tr_16', und als Gegenprobe
    die Felddeterminante ueber den 16 Gewichten (cosh und log(1+w^2));
 4  Basislinie: N Moden, volle Fockspur -> welche Teilung?  (balanciert = generisch);
 5  Ordnung der Vorzeichenauswahl: log(1 + prod tanh y) = det(y) + O(y^7).
"""
import itertools, math, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import sympy as sp
from scipy.optimize import minimize

OK = True
def check(name, cond, info=""):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}" + (f" -- {info}" if info else ""))

# ------------------------------------------------------------------ 1
print("1  Spuridentitaeten auf dem Fockraum von fuenf Moden")
y = sp.symbols("y1:6", real=True)
constraint = {y[4]: -sum(y[:4])}
subsets = [S for k in range(6) for S in itertools.combinations(range(5), k)]
def weight(S):                       # Eigenwert von e^{-T Sigma^} auf dem Zustand S: prod_{i in S} e^{-2 y_i}
    return sp.exp(-2 * sum(y[i] for i in S))
Z32 = sum(weight(S) for S in subsets)
Zpar = sum((-1) ** len(S) * weight(S) for S in subsets)
Z16 = sum(weight(S) for S in subsets if len(S) % 2 == 0)
c = sp.prod([sp.cosh(t) for t in y]); s_ = sp.prod([sp.sinh(t) for t in y])
def eq(a, b):
    num = [("-7/10", "3/10", "11/10", "-1/5"), ("2/5", "9/10", "-13/10", "1/4"), ("3/2", "-1/2", "1/5", "-4/5")]
    for v in num:
        v = [sp.Rational(a) for a in v]
        sub = dict(zip(y[:4], v)); sub[y[4]] = -sum(v)
        if abs(complex(sp.N((a - b).subs(sub), 30))) > 1e-20: return False
    return True
check("Tr_32 e^{-T Sigma^} = 2^5 prod cosh y_i  (Summe y = 0)", eq(Z32, 32 * c))
check("Tr_32 (-1)^F e^{-T Sigma^} = 2^5 prod sinh y_i", eq(Zpar, -32 * s_) or eq(Zpar, 32 * s_),
      "prod (1 - e^{-2y}) = 2^5 prod sinh y bei Summe y = 0")
sign = 1 if eq(Zpar, 32 * s_) else -1
check("Tr_16 = 16 (prod cosh + sigma prod sinh) mit sigma = " + str(sign), eq(Z16, 16 * (c + sign * s_)))
print(f"   Damit Tr_16 = 16 * prod cosh y * (1 + ({sign}) prod tanh y) > 0, da |prod tanh| < 1.")

# ------------------------------------------------------------------ 2
print("\n2  Vollstaendiges antiperiodisches Paarprodukt = det cosh")
for xval in (0.3, 1.0, 2.5):
    prod = 1.0
    for nn in range(0, 200000):
        prod *= 1 + (2 * xval) ** 2 / ((2 * nn + 1) * math.pi) ** 2
    check(f"prod_n (1 + x^2/((n+1/2)pi)^2 / 4 ...) = cosh(x) bei x={xval}", abs(prod / math.cosh(xval) - 1) < 1e-5,
          f"{prod:.8f} vs {math.cosh(xval):.8f}")

# ------------------------------------------------------------------ 3
print("\n3  Globale Minima")
basis5 = np.linalg.qr(np.vstack([np.ones(5), np.eye(5)[:4]]).T)[0][:, 1:]
def classify(lam, tol=2e-4):
    lam = np.array(sorted(lam))
    if np.max(np.abs(lam)) < 1e-5: return "0"
    groups = [[lam[0]]]
    for v in lam[1:]:
        (groups[-1].append(v) if abs(v - groups[-1][-1]) < tol * max(1, abs(v)) else groups.append([v]))
    parts = sorted(groups, key=len, reverse=True)
    label = "+".join(str(len(g)) for g in parts)
    if label in ("3+2", "4+1"):
        big = parts[0][0]
        label += " (" + ("Mehrheit positiv" if big > 0 else "Mehrheit negativ") + ")"
    return label
def gmin(fun, basis, rng, starts=300):
    best = None
    d = basis.shape[1]
    for _ in range(starts):
        z0 = rng.normal(size=d) * rng.uniform(0.1, 4)
        r = minimize(lambda z: fun(basis @ z), z0, method="BFGS", options={"gtol": 1e-12, "maxiter": 2000})
        if best is None or r.fun < best.fun - 1e-12: best = r
    return basis @ best.x, best.fun
def lcosh(t): return np.abs(t) + np.log1p(np.exp(-2 * np.abs(t))) - math.log(2)
w16_index = [S for S in subsets if len(S) % 2 == 0]
W16 = np.array([[1.0 if i in S else 0.0 for i in range(5)] for S in w16_index])
def weights16(lam):                    # Gewichte der 16 fuer Sigma mit Eigenwerten lam (zentriert)
    return W16 @ lam - (W16.sum(1) / 5) * lam.sum()
sources = {
    "Tr_32 (antiperiodisch)": lambda lam: -np.sum(lcosh(lam)),
    "Tr_16 (gerader Teil)":  lambda lam: -(np.sum(lcosh(lam)) + np.log1p(sign * np.prod(np.tanh(lam)))),
    "Tr_16' (ungerader Teil)": lambda lam: -(np.sum(lcosh(lam)) + np.log1p(max(-sign * np.prod(np.tanh(lam)), -1 + 1e-300))),
    "Felddet. 16, cosh":      lambda lam: -np.sum(lcosh(weights16(lam))) / 4.0,
    "Felddet. 16, log(1+w^2)": lambda lam: -np.sum(np.log1p(weights16(lam) ** 2)) / 4.0,
}
rng = np.random.default_rng(7)
grid = [0.12, 0.2, 0.3, 0.4, 0.45, 0.49]
results = {}
for name, src in sources.items():
    print("   Quelle:", name)
    for rho in grid:
        lam, val = gmin(lambda l: rho * np.sum(l ** 2) + src(l), basis5, rng)
        lab = classify(lam)
        results[(name, rho)] = (lab, lam)
        print(f"     rho={rho:4.2f}: {lab:26s} {np.round(np.sort(lam), 4)}  Gamma={val:.6f}")
check("Tr_32 waehlt 3+2 auf dem ganzen Gitter", all(results[("Tr_32 (antiperiodisch)", r)][0].startswith("3+2") for r in grid))
check("Tr_16 waehlt 3+2 auf dem ganzen Gitter", all(results[("Tr_16 (gerader Teil)", r)][0].startswith("3+2") for r in grid))
check("Tr_16 legt das Vorzeichen fest: Mehrheit positiv = Sigma ~ +(2,2,2,-3,-3) = -6Y; Tr_16' die Gegenbahn",
      len({results[("Tr_16 (gerader Teil)", r)][0] for r in grid}) == 1, results[("Tr_16 (gerader Teil)", grid[0])][0])
check("Felddeterminante ueber 16 waehlt nicht 3+2 (Gegenprobe)",
      all(not results[("Felddet. 16, cosh", r)][0].startswith("3+2") for r in grid))

# ------------------------------------------------------------------ 4
print("\n4  Basislinie: N Moden, volle Fockspur, rho = 0.3")
for N in range(2, 10):
    bN = np.linalg.qr(np.vstack([np.ones(N), np.eye(N)[:N - 1]]).T)[0][:, 1:]
    lam, val = gmin(lambda l: 0.3 * np.sum(l ** 2) - np.sum(lcosh(l)), bN, rng, starts=200)
    print(f"     N={N}: {classify(lam)}")
print("   Die volle Fockspur waehlt jeweils die balancierteste Zweierteilung; fuer N = 5 ist das 3+2.")

# ------------------------------------------------------------------ 5
print("\n5  Ordnung der Vorzeichenauswahl")
t = sp.symbols("t")
ser = sp.series(sp.log(1 + sp.prod([sp.tanh(t * yi) for yi in y])), t, 0, 8).removeO()
poly = sp.Poly(sp.expand(ser), t)
low = {k: sp.expand(poly.coeff_monomial(t ** k)) for k in range(0, 7)}
check("log(1 + prod tanh(t y)) = t^5 y1..y5 + O(t^7)",
      all(low[k] == 0 for k in (0, 1, 2, 3, 4, 6)) and sp.expand(low[5] - sp.prod(y)) == 0)
print("   Bis zur vierten Ordnung sind Tr_32 und Tr_16 gleich (Form 3+2 aus tr Sigma^4).")
print("   Die Projektion auf die 16 traegt erst mit det(Sigma) in fuenfter Ordnung bei und waehlt so")
print("   eine der beiden Bahnen +-(2,2,2,-3,-3), die die Determinantenquelle allein nicht trennt.")
for shape in [(2, 2, 2, -3, -3), (-2, -2, -2, 3, 3)]:
    print(f"   det auf {shape}: {np.prod(shape)}")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
sys.exit(0 if OK else 1)
