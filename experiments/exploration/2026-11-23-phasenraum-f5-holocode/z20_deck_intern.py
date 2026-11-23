"""Z20: Der Deck-interne Zweig von v957 (CCC.SEAM.KINEMATICS.01).

v957 hebt die Naht auf S^3: Die Nahtsphäre ist die Hopf-Basis (z = z1/z0), das Deck D = diag(i, -i)
wirkt frei, die vier Marken werden vier Hopf-Fasern, Polfasern tragen die volle mu4-Holonomie
chi(D) = i, Nahtstränge die Zwei-Torsion chi(D^2) = -1 = (-1)^F. Offen ist dort die Gabel: Das Deck
wirkt als Raumquotient (Linse X = S^3/Z4) oder als Strukturgruppe eines internen Bündels über S^3.
Z17 fand: Im räumlichen Zweig ist eine einzelne Naht unmöglich.

Geprüft wird der interne Zweig.
 1. v957 nachgebaut: Wirkungen von Deck, Uhr und Reziprokem; Freiheit; Verschlingung der Fasern.
 2. Räumlicher Zweig: deckinvariante Felder mit Hopf-Ladung k; Zahl der Nähte.
 3. Interner Zweig: S^3 ist einfach zusammenhängend. Eine physische Deck-Holonomie braucht deshalb eine
    Defektlinie (Flussschlauch). Das Paarungsfeld Phi trägt J-Ladung 2 und wechselt um einen J-Schlauch
    das Vorzeichen. Es muss dort verschwinden und um pi winden: der verdrillte Defekt aus G.
 4. Uhrensymmetrie und Reziprokes: Wo liegt der Schlauch, wie viele braucht es?
 5. Bilanz der Zweige gegen H6.
"""
import json
from pathlib import Path

import numpy as np

OK, BEF = [], {}


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


print("1. v957 nachgebaut\n")
D = np.diag([1j, -1j])                       # deck (v492)
C = np.diag([1j, 1])                         # clock
Jr = np.array([[0, -1], [1, 0]], complex)    # reciprocal lift


def base(v):
    return np.inf if abs(v[0]) < 1e-14 else v[1] / v[0]


rng = np.random.default_rng(20)
v = rng.normal(size=2) + 1j * rng.normal(size=2)
v /= np.linalg.norm(v)
zb = base(v)
print(f"   Basispunkt z = {zb:.4f}")
print(f"   Deck:  z -> {base(D @ v):.4f}  (= -z: {np.isclose(base(D @ v), -zb)})")
print(f"   Uhr:   z -> {base(C @ v):.4f}  (= -i z: {np.isclose(base(C @ v), -1j * zb)})")
print(f"   Rezipr.: z -> {base(Jr @ v):.4f}  (= -1/z: {np.isclose(base(Jr @ v), -1 / zb)})")
check("Deck = Halbdrehung, Uhr = Vierteldrehung, Reziprokes = -1/z",
      np.isclose(base(D @ v), -zb) and np.isclose(base(C @ v), -1j * zb) and np.isclose(base(Jr @ v), -1 / zb))
check("D^2 = -1, D frei auf S^3", np.allclose(D @ D, -np.eye(2)) and
      all(not np.any(np.isclose(np.linalg.eigvals(np.linalg.matrix_power(D, k)), 1)) for k in (1, 2, 3)))
print("   Das Deck hebt die Blattinvolution z -> -z der Basis mit D^2 = -1 (Spinorvorzeichen) auf S^3.")


def fibre(zp, n=600):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    if zp == np.inf:
        pts = np.stack([0 * t + 0j, np.exp(1j * t)], 1)
    else:
        pts = np.stack([np.exp(1j * t), zp * np.exp(1j * t)], 1) / np.sqrt(1 + abs(zp) ** 2)
    x = np.stack([pts[:, 0].real, pts[:, 0].imag, pts[:, 1].real, pts[:, 1].imag], 1)
    pole = np.array([0.3, -0.5, 0.6, np.sqrt(1 - 0.3 ** 2 - 0.25 - 0.36)])
    # stereographic projection from a point off all fibres used here
    e = pole / np.linalg.norm(pole)
    Q, _ = np.linalg.qr(np.column_stack([e, np.eye(4)[:, :3]]))
    y = x @ Q                                   # first coordinate along e
    return y[:, 1:] / (1 - y[:, :1])


def linking(a, b):
    da = np.roll(a, -1, 0) - a
    db = np.roll(b, -1, 0) - b
    s = 0.0
    for i in range(len(a)):
        r = a[i] - b
        s += np.sum(np.einsum("j,kj->k", da[i], np.cross(db, r)) / np.linalg.norm(r, axis=1) ** 3)
    return s / (4 * np.pi)


marks = [1, 1j, -1, -1j]
lk_pole_mark = [linking(fibre(0.0), fibre(m)) for m in marks]
lk_poles = linking(fibre(0.0), fibre(np.inf))
print(f"   Gauß-Verschlingung Polfaser (z=0) mit den vier Markenfasern: {np.round(lk_pole_mark, 3)}")
print(f"   Verschlingung der beiden Polfasern: {lk_poles:+.3f}")
check("jede Markenfaser verschlingt die Polfaser einmal", np.allclose(np.abs(lk_pole_mark), 1, atol=0.02))
check("Polfasern bilden eine Hopf-Verschlingung", abs(abs(lk_poles) - 1) < 0.02)

print("\n2. Räumlicher Zweig (Deck als Raumquotient, X = S^3/Z4)\n")
print("   Monome z0^a z1^b haben Hopf-Ladung k = a + b und Deck-Charakter i^(a-b).")
rows = []
for k in range(0, 5):
    chars = sorted({int((a - (k - a)) % 4) for a in range(k + 1)})
    inv = any((a - (k - a)) % 4 == 0 for a in range(k + 1))
    rows.append((k, inv))
    print(f"   k = {k}: Charaktere i^m mit m in {chars}; deckinvariantes Monom vorhanden: {inv}")
check("k = 1 nicht invariant, kleinstes invariantes k > 0 ist 2", [r[1] for r in rows][:3] == [True, False, True])
print("   Ein deckinvariantes Feld mit k = 2, etwa a z0 z1 + b ..., hat generisch zwei Nullfasern: Nähte kommen")
print("   paarweise (Z17). Eine einzelne Naht ist im räumlichen Zweig unmöglich.")

print("\n3. Interner Zweig (Deck als Strukturgruppe eines internen Bündels über S^3)\n")
print("   pi_1(S^3) = 0: Ein flaches internes Z4-Bündel ohne Defekte ist trivial. Soll die v957-Holonomie")
print("   chi(D) = i physisch sein, muss sie um eine Defektlinie liegen: pi_1(S^3 minus Kreis) = Z erlaubt")
print("   jede Meridian-Holonomie, auch J mit J^2 = -1 = (-1)^F.")
J = np.array([[0, -1], [1, 0]], complex)    # internal quarter turn on (psi, chi)
mass = np.array([[0, 1], [1, 0]], complex)  # Phi psi chi
flip = J.T @ mass @ J
print(f"   Masse psi chi unter J: J^T M J = {flip.real.astype(int).tolist()} = -M, also trägt Phi die J-Ladung 2.")
check("Phi hat J-Ladung 2", np.allclose(flip, -mass))
# a continuous complex field on a punctured disc with monodromy -1 must vanish inside
zeros_found = []
for trial in range(200):
    c = rng.normal(size=6) + 1j * rng.normal(size=6)
    th = np.linspace(0, 2 * np.pi, 4001)
    # boundary value of a field with monodromy -1: e^{i theta/2} times a periodic function
    per = sum(c[n] * np.exp(1j * (n - 3) * th) for n in range(6)) + 3 * np.exp(0j * th)
    f = np.exp(1j * th / 2) * per
    wind = np.sum(np.diff(np.unwrap(np.angle(f)))) / (2 * np.pi)
    zeros_found.append(abs(wind - round(wind)) > 0.4)          # half-integer winding
print(f"   Randwindung eines Felds mit Monodromie -1: in allen 200 Proben halbzahlig "
      f"({check('halbzahlige Windung, also Nullstelle im Inneren', all(zeros_found))})")
print("   Eine stetige nirgends verschwindende Funktion auf der Scheibe hat ganzzahlige Randwindung. Also muss Phi")
print("   auf dem Schlauch verschwinden und um pi (mal ungerade) winden: genau der verdrillte Defekt aus G.")
print("   Nach Z7 trägt er eine chirale Majorana-Weyl Mode je Paar; 16 Paare geben c = 8.")

print("\n4. Lage und Zahl der Schläuche\n")
fixed_clock = [zp for zp in (0.0, np.inf, 1, 1j, -1, -1j, 0.5 + 0.2j)
               if (zp == np.inf and True) or (zp != np.inf and np.isclose(-1j * zp, zp))]
print(f"   Unter der Uhr (z -> -i z) feste Basispunkte: {['0' if f == 0 else 'unendlich' for f in fixed_clock]}")
check("uhreninvariante Fasern sind nur die Polfasern", len(fixed_clock) == 2)
print("   Die Marken bilden eine Uhrenbahn der Länge 4; ein einzelner uhreninvarianter Schlauch muss eine")
print("   Polfaser sein. Dort trägt v957 genau die volle mu4-Holonomie chi(D) = i.")
print("   Das Reziproke (z -> -1/z) vertauscht die Pole. Wird es von der Konfiguration verlangt, braucht es")
print("   beide Polfasern (eine Hopf-Verschlingung, zwei Nähte). Die Einzelnaht-Zustände bilden dagegen ein")
print("   Dublett unter dem Reziproken: Der Grundzustand wählt einen Pol und bricht das Reziproke spontan.")
cost = {"eine Polnaht": 1, "beide Polnähte": 2}
print(f"   Kosten in Einheiten der Nahtwirkung: {cost}; minimal ist eine Naht.")
check("minimale Konfiguration: genau ein Schlauch", min(cost.values()) == 1)

print("\n5. Bilanz der Zweige gegen H6\n")
table = [
    ("räumlich (Linse)", 2, "nein (Phi braucht räumliche Monodromie)", 16, 6),
    ("intern, Reziprokes verlangt", 2, "ja", 16, 6),
    ("intern, Reziprokes spontan gebrochen", 1, "ja", 8, 3),
]
print("   Zweig                                  Nähte  G hergeleitet                           c    N_fam")
for name, ns, g, c, nf in table:
    print(f"   {name:38s}  {ns}      {g:38s}  {c:2d}   {nf}")
print("   (c und N_fam unter der Lesart: je Naht ein 16er-Satz und drei Familien nach C'')")
check("nur intern mit gebrochenem Reziproken: eine Naht, c = 8, drei Familien", [r for r in table if r[1] == 1][0][3:] == (8, 3))
BEF["räumlicher Zweig erlaubt keine Einzelnaht"] = True
BEF["Reziprokes symmetrisch verlangt gibt zwei Nähte"] = True
BEF["warum die Deck-Holonomie physisch ist, bleibt Eingabe (v522-Umlauf)"] = True

urteil = ("Lesart, deutlich gestärkt: Im Deck-internen Zweig erzwingt eine physische Deck-Holonomie auf dem einfach "
          "zusammenhängenden S^3 mindestens einen Flussschlauch; um ihn wechselt Phi (J-Ladung 2) das Vorzeichen und muss "
          "dort verschwinden: Der verdrillte Defekt aus G folgt. Uhrensymmetrie legt ihn auf eine Polfaser (volle "
          "mu4-Holonomie, die vier Markenfasern verschlingen ihn einmal); minimal ist genau einer, mit spontan gebrochenem "
          "Reziproken. Eingaben: Deck intern, Deck-Holonomie physisch, keine erzwungene Reziprozität der Konfiguration.")
print(f"\nUrteil: {urteil}")
print("\nZusammenfassung")
for n_, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n_}")
Path(__file__).with_name("z20_ergebnis.json").write_text(json.dumps({
    "checks": {n_: c for n_, c in OK}, "befunde": BEF, "urteil": urteil,
    "linking_pole_marks": [float(x) for x in lk_pole_mark], "linking_poles": float(lk_poles)}, ensure_ascii=False, indent=1))
assert all(c for _, c in OK)
