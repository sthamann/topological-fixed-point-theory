"""Z10: Prüfung der zweiten externen Kritik (23. November 2026).

 1. Fehlende Gegenprobe zu Z7: Defekt mit einer windenden Komponente, aber ohne J (n+ = 1, n- = 0, kein
    Viertelfluss). Trägt er auch eine chirale Mode?
 2. Index als Invariante statt Zählung im Energiefenster: Windung von det M(phi) der Majoranamassen-
    matrix an jeder der acht Impulsecken (W = 0 an der physischen, W >= 2 r_w an den Verdopplern).
 3. Familien als chiraler Index der vier Marken: H^0(K(D)) für D = mu4 auf CP^1 (logarithmische
    Differentiale), daneben H^0(O(D)) (TFPTs Fünferträger). Wirkung der Vierteldrehung z -> i z.
 4. Twistgewicht h_R = N/16 von N Majorana-Weyl Fermionen: lokaler Strom mit h = 1 nur bei N = 16.
 5. Satzzahl: Gegenbeispiel zur Behauptung aus Z8, nur ein Satz trage die Familiencharaktere.
 6. Schleifenenergie E(R) = 2 pi sigma R + (h - c/24)/R in den Zuständen R-Grundzustand, Vakuum und im
    KMS-Zustand von H6.
"""
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import eigh

_src = Path(__file__).with_name("z7_ursprung_naht.py").read_text()
_ns = {"__file__": str(Path(__file__).with_name("z7_ursprung_naht.py"))}
exec(_src[:_src.index("L = 14")], _ns)
build = _ns["build"]
OK = []


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


def core_count(L, npl, nmi, fl, R=3.5, window=0.15, thr=0.6, dk=1e-3):
    H, (rad, ns) = build(L, 0.0, npl, nmi, fl)
    E, V = eigh(H, driver="evr")
    S = V[:, np.abs(E) < window]
    cm = np.tile(np.repeat(rad <= R, ns * 2), 2)
    cv, cvec = np.linalg.eigh(S.conj().T @ (cm[:, None] * S))
    Hp, _ = build(L, dk, npl, nmi, fl)
    Hm, _ = build(L, -dk, npl, nmi, fl)
    dH = (Hp - Hm) / (2 * dk)
    cs = S @ cvec[:, cv > thr]
    v = np.real(np.einsum("ij,ik,kj->j", cs.conj(), dH, cs))
    return int(np.sum(v > 0.2) - np.sum(v < -0.2)), np.round(v, 3), np.round(cv[cv > thr], 3)


print("1. Gegenprobe ohne J\n")
cases = {"mit J (1, 0, Viertelfluss)": (1, 0, 0.25), "ohne J (1, 0, kein Fluss)": (1, 0, 0.0),
         "gewöhnlich (1, 1)": (1, 1, 0.0), "nur J (0, 0, Viertelfluss)": (0, 0, 0.25)}
counts = {}
for name, (a, b, fl) in cases.items():
    n, v, w = core_count(14, a, b, fl)
    counts[name] = n
    print(f"   {name:28s}: netto {n}, Geschwindigkeiten {v}, Kerngewichte {w}")
check("ohne J ebenfalls eine Mode", counts["ohne J (1, 0, kein Fluss)"] == 1)
print("   Die Kritik hat recht: Die Modenzahl allein erzwingt J nicht.")

print("\n2. Index als Windung von det M(phi) an den acht Impulsecken (m = 0,6, r_w = 1)\n")
m, rw = 0.6, 1.0
corners = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]


def winding(n_plus, n_minus, W):
    ph = np.linspace(0, 2 * np.pi, 2001)
    det = (m * np.exp(1j * n_plus * ph)) * (m * np.exp(1j * n_minus * ph)) - W ** 2
    return int(round(np.sum(np.diff(np.unwrap(np.angle(det)))) / (2 * np.pi)))


for (npl, nmi) in ((1, 0), (1, 1), (0, 0), (2, 0)):
    per = [winding(npl, nmi, rw * 2 * sum(cn)) for cn in corners]
    chir = [(-1) ** sum(cn) for cn in corners]
    total = sum(w * c for w, c in zip(per, chir))
    print(f"   (n+, n-) = ({npl}, {nmi}): Windung je Ecke {per}, mit Eckenchiralität gewichtet: Index {total}")
check("Index = n+ + n-", all(sum(winding(a, b, rw * 2 * sum(cn)) * (-1) ** sum(cn) for cn in corners) == a + b
                             for a, b in ((1, 0), (1, 1), (0, 0), (2, 0))))
per0 = [winding(1, 0, 0.0) for cn in corners]
print(f"   ohne Wilson-Glied (W = 0 überall): Windung je Ecke {per0}, gewichtet "
      f"{sum(w * (-1) ** sum(cn) for w, cn in zip(per0, corners))}: Die Verdoppler heben sich auf.")
print("   Der Index hängt nur an der Windung an der physischen Ecke, nicht am Viertelfluss. Was J leistet, ist")
print("   etwas anderes: Mit J genügt ein einziges komplexes Massenfeld m e^{i theta} psi chi, dessen Phase um pi")
print("   windet; dann ist die Masse selbst zweiwertig (Phi -> -Phi), das sind die zwei Blätter auf Bosonischem,")
print("   und auf Fermionen wirkt mu4. Ohne J braucht die Konfiguration (1, 0) zwei unabhängige Majoranamassen")
print("   und hat keine zwei Blätter. J wird also durch die Blattstruktur von P1 gewählt, nicht durch die Zählung.")

print("\n3. Familien als chiraler Index der vier Marken D = mu4 auf CP^1\n")
z = sp.symbols("z")
basis_K = [sp.Integer(1), z, z ** 2]                 # omega = P(z) dz/(z^4 - 1), deg P <= 2
basis_O = [z ** k for k in range(5)]                 # f = Q(z)/(z^4 - 1), deg Q <= 4
# pullback under z -> i z: omega -> i P(i z) dz/(z^4 - 1); f -> Q(i z)/(z^4 - 1)
MK = sp.Matrix([[sp.Poly(sp.expand(sp.I * b.subs(z, sp.I * z)), z).coeff_monomial(c) for b in basis_K] for c in basis_K])
MO = sp.Matrix([[sp.Poly(sp.expand(b.subs(z, sp.I * z)), z).coeff_monomial(c) for b in basis_O] for c in basis_O])
evK = sorted([complex(e) for e in MK.eigenvals(multiple=True)], key=lambda x: (x.real, x.imag))
evO = [complex(e) for e in MO.eigenvals(multiple=True)]
half = (MO ** 2).eigenvals()
print(f"   h^0(K(D)) = h^0(O(2)) = {len(basis_K)}, h^1 = 0: Index 3 (Riemann-Roch). Vierteldrehung: {evK}")
print(f"   h^0(O(D)) = {len(basis_O)} (TFPTs Fünferträger). Vierteldrehung: {np.round(evO, 6)}")
print(f"   Halbdrehung auf dem Träger: Eigenwerte {half} -> Teilung 3 + 2")
check("Familien {-1, i, -i} aus dem Index", np.allclose(sorted(evK, key=lambda x: (x.real, x.imag)),
                                                       sorted([-1, 1j, -1j], key=lambda x: (x.real, x.imag))))
check("Träger 3 + 2", half == {1: 3, -1: 2})
dimL = [sp.binomial(5, k) for k in (0, 2, 4)]
print(f"   Lambda^gerade C^5 = {' + '.join(str(d) for d in dimL)} = {sum(dimL)} (Spin(10)-Spinor einer Generation)")
print(f"   Materie = H^0(K(D)) x Lambda^gerade H^0(O(D)): {3 * sum(dimL)} Komponenten")
check("3 x 16 = 48", 3 * sum(dimL) == 48)
print("   Spin-Dirac auf CP^1, gekoppelt an K^(1/2) x O(D) = O(3): ker D+ = H^0(O(2)) = 3, ker D- = 0.")
print("   Dieselben vier Marken geben also den Fünferträger mit 3 + 2 und einen chiralen Familienindex 3 mit")
print("   genau den Charakteren {-1, i, -i}, direkt aus z -> i z, ohne Spinstrukturen und ohne Wahl des Zentrums.")

print("\n4. Twistgewicht h_R = N/16\n")
for N in (8, 16, 24, 32, 48):
    h = sp.Rational(N, 16)
    print(f"   N = {N:2d}: c = {sp.Rational(N, 2)}, h_R = {h}" + ("  -> Gewicht 1: 128 Spinorströme, 120 + 128 = 248" if h == 1 else ""))
check("nur N = 16 gibt Gewicht 1", [N for N in (8, 16, 24, 32, 48) if sp.Rational(N, 16) == 1] == [16])
print("   Das wählt Blöcke aus 16 Majorana-Weyl Fermionen, also E8 je Block. Es entscheidet nicht zwischen einem")
print("   und drei Blöcken (E8 hoch 3 hat je Block wieder h_R = 1), wohl aber gegen eine gemeinsame 48er-Summe.")

print("\n5. Satzzahl: Gegenbeispiel aus der Kritik\n")
Mfam = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]], complex)      # Z1: g on the three family labels
big = np.kron(np.kron(Mfam, Mfam), Mfam)                           # three identical sets
diag_states = [np.kron(np.kron(np.eye(3)[i], np.eye(3)[i]), np.eye(3)[i]) for i in range(3)]
Dsub = np.array([[ds2.conj() @ big @ ds1 for ds1 in diag_states] for ds2 in diag_states])
ev_diag = np.sort_complex(np.round(np.linalg.eigvals(Dsub), 9))
print(f"   Unterraum 'alle Sätze tragen dasselbe Label' unter g: Eigenwerte {ev_diag}")
closed = np.allclose(big @ np.array(diag_states).T, np.array(diag_states).T @ Dsub)
print(f"   Der Unterraum ist invariant: {closed}")
check("Gegenbeispiel bestätigt", closed and np.allclose(ev_diag, np.sort_complex(np.array([-1, -1j, 1j]))))
print("   Die Kritik hat recht: Die Familiencharaktere schließen drei Sätze nicht aus. Die Behauptung aus Z8 ist")
print("   zu stark und wird zurückgenommen. Unterscheiden könnte das Niveau eines gemeinsamen E8 Stroms (1 gegen 3).")

print("\n6. Schleifenenergie E(R) = 2 pi sigma R + (h - c/24)/R, c = 8\n")
for name, h in (("Ramond-Grundzustand (h = 1)", 1.0), ("Vakuum (h = 0)", 0.0), ("KMS-Zustand von H6 (<L0> = c/24)", 8 / 24)):
    coef = h - 8 / 24
    print(f"   {name:34s}: Koeffizient von 1/R = {coef:+.4f}" +
          (f", Minimum bei R* = 1/sqrt(3 pi sigma)" if coef > 0 else ", kein stabilisierender Term"))
z4 = json.loads(Path(__file__).with_name("z4_ergebnis.json").read_text())
check("KMS-Zustand ohne 1/R-Term (aus Z4)", z4["checks"]["Energie null bei y = 1"])
print("   Im Zustand von H6 (Z4: <L0> = c/24 bei tau = i) verschwindet der Term; die Stabilisierung der Kritik")
print("   verlangt den Ramond-Grundzustand und gilt dann für eine andere Zustandswahl.")
print(f"\n   Schwelle der Kritik für ihr Zusatzmodell auf S^3: 2 - sqrt(8/3) = {2 - np.sqrt(8 / 3):.4f} (nur Arithmetik geprüft;")
print("   das Modell selbst liegt nicht vor).")

print("\nZusammenfassung")
for n, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
Path(__file__).with_name("z10_ergebnis.json").write_text(json.dumps({
    "checks": {n: c for n, c in OK}, "counts": counts, "family_index_eigs": [str(e) for e in evK],
    "carrier_half_turn": {str(k): v for k, v in half.items()}}, indent=1))
assert all(c for _, c in OK)
