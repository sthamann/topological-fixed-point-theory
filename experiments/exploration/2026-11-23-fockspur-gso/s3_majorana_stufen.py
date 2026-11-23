"""Vorstufe: 2N reelle (Majorana) Moden, Quelle Sigma in so(2N), Spur ueber den Spinorfockraum.

Frage: Erzeugt dieselbe Fockspur-Quelle aus reellen Moden zuerst die komplexe Struktur J
(SO(2N) -> U(N)), und waehlt sie danach ohne weitere Quelle die Teilung 3+2?

Fuer Sigma in so(2N) mit Eigenwerten +-i mu_k gilt
    Tr_{2^N} e^{-T Sigma^} = prod_k 2 cosh(mu_k),   Tr_{S+} - Tr_{S-} = prod_k 2 sinh(mu_k)
(mu absorbiert T/2). Gamma = rho sum mu_k^2 - log Tr.  Es gibt keine Spurbedingung an die mu_k.

Gerechnet: globale Minima fuer N = 5 und N = 8; Hessian der so(10) Wirkung am J-Sattel in den
su(5) Richtungen; Basislinie, wann ueberhaupt eine Teilung entsteht.
"""
import math, sys, warnings
import numpy as np
from scipy.optimize import minimize, brentq
warnings.filterwarnings("ignore")
OK = True
def check(name, cond, info=""):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}" + (f" -- {info}" if info else ""))
def lcosh(t): return np.abs(t) + np.log1p(np.exp(-2 * np.abs(t))) - math.log(2)
rng = np.random.default_rng(3)
def gmin(fun, d, starts=300):
    best = None
    for _ in range(starts):
        r = minimize(fun, rng.normal(size=d) * rng.uniform(0.1, 3), method="BFGS", options={"gtol": 1e-12})
        if best is None or r.fun < best.fun - 1e-12: best = r
    return best.x, best.fun

print("1  Spinorfockspur von 2N Majoranamoden")
for N in (5, 8):
    for rho in (0.2, 0.35, 0.45):
        full = lambda m, rho=rho: rho * np.sum(m ** 2) - np.sum(lcosh(m))
        chiral = lambda m, rho=rho: rho * np.sum(m ** 2) - (np.sum(lcosh(m)) + np.log1p(np.prod(np.tanh(m))))
        mf, _ = gmin(full, N); mc, _ = gmin(chiral, N)
        mstar = brentq(lambda x: 2 * rho * x - math.tanh(x), 1e-6, 50)
        eq_full = np.allclose(np.abs(mf), mstar, atol=1e-5)
        eq_chi = np.allclose(np.abs(mc), np.abs(mc).mean(), atol=1e-5)
        neg = int(np.sum(mc < 0))
        print(f"   N={N} rho={rho}: |mu| voll {np.round(np.abs(mf),5)}  mu*={mstar:.5f}; chiral: |mu| = {np.abs(mc).mean():.5f}, Vorzeichen {np.sign(mc).astype(int)}")
        check(f"N={N}, rho={rho}: alle |mu_k| gleich (voll: = mu*; Sigma ~ J mit J^2 = -1, Rest U(N))", eq_full and eq_chi)
        check(f"N={N}, rho={rho}: chirale Spur waehlt gerade Anzahl negativer mu (eine Orientierungsklasse)", neg % 2 == 0)

print("\n2  Keine weitere Teilung am J-Sattel ohne zweite Quelle")
rho = 0.35
mstar = brentq(lambda x: 2 * rho * x - math.tanh(x), 1e-6, 50)
f2 = 2 * rho - 1 / math.cosh(mstar) ** 2
print(f"   Krümmung je Eigenwert f''(mu*) = 2 rho - sech^2(mu*) = {f2:.6f}")
check("Alle su(5) Richtungen (mu_k -> mu* + v y_k, Summe y = 0) sind stabil", f2 > 0)
print("   Eine so(10) Fockquelle bricht SO(10) -> U(5), aber nicht U(5) -> S(U(3) x U(2)).")
print("   3+2 braucht eine zweite, spurfreie su(5) Quelle (s2_fockspur_auswahl.py).")

print("\n3  Basislinie: Ursprung der Spurbedingung")
print("   Mit Spurbedingung (su(N)) entsteht die balancierte Teilung; ohne (so(2N)) entsteht J.")
print("   Die Zahl der Moden bleibt in beiden Stufen Eingabe: 10 = 2*5 reell beziehungsweise 5 komplex.")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
sys.exit(0 if OK else 1)
