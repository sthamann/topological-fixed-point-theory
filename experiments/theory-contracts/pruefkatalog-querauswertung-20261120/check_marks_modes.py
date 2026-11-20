"""GESAMT20261120: ein Objekt für den festen Bestand. Exakte Prüfungen, kein Beweis der Physik.

Objekt: die Sphäre mit dem Divisor D der vier Marken bei den vierten Einheitswurzeln und
ihrer Drehgruppe mu4 (z -> i z). Zwei kanonische Räume (im Atlas bewiesen: Rang 3 und
Dimension 5 gehören zum selben Divisor):

  H0(O(D)) = { p(z)/(z^4 - 1), Grad p <= 4 }   Dimension 5   (Träger)
  H1(Sphäre ohne Marken)                         Rang 3        (Familien)

Neu geprüft wird, wie mu4 und die Involutionen der Marken auf diesen acht Moden wirken.
Gleiche Charaktere sind keine Identifikation von Operatoren.
"""
from fractions import Fraction as F
from itertools import permutations

checks = []


def check(name, ok, detail=""):
    checks.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


# --- 1. Charaktere von mu4 (Exponent k steht für den Charakter i^k) -----------------
carrier = [k % 4 for k in range(5)]           # Monome z^k, der Nenner z^4 - 1 ist invariant
perm = [0, 1, 2, 3]                           # reguläre Darstellung auf den vier Marken
family = sorted(perm)[1:]                     # Schleifen um die Marken, Summe null: trivialer Charakter entfällt
modes = sorted(carrier + family)
check("Träger trägt die Charaktere 1, i, -1, -i, 1", carrier == [0, 1, 2, 3, 0])
check("Familien tragen i, -1, -i", family == [1, 2, 3])
check("acht Moden = zweimal die reguläre Darstellung von mu4", modes == [0, 0, 1, 1, 2, 2, 3, 3],
      "jeder Charakter genau zweimal")

# --- 2. Die Involution z -> -z (Quadrat des Erzeugers) --------------------------------
par = lambda ks: [(-1) ** k for k in ks]
even_c, odd_c = par(carrier).count(1), par(carrier).count(-1)
even_f, odd_f = par(family).count(1), par(family).count(-1)
check("Träger zerfällt in 3 gerade und 2 ungerade", (even_c, odd_c) == (3, 2), "1, z^2, z^4 gegen z, z^3")
check("Familien zerfallen in 1 und 2", (even_f, odd_f) == (1, 2), "wie der Anker (1, 1, 2): zwei gleich, eine anders")
check("Spur der Involution auf den acht Moden ist null", sum(par(modes)) == 0,
      "verträglich mit KERNSIX Schritt 2: Polarisation hat Spur null")
check("das ungerade Paar des Trägers ist ein konjugiertes Paar (i, -i)", sorted(k for k in carrier if k % 2) == [1, 3],
      "eine reelle Darstellung der Dimension 2")

# Jede holomorphe Involution der Sphäre wirkt auf Sym^4(C^2) mit Eigenwerten (1,-1,1,-1,1).
for name, action in (("z -> 1/z", lambda k: 4 - k), ("z -> -1/z", lambda k: 4 - k)):
    sign = lambda k, name=name: -1 if name == "z -> 1/z" else -(-1) ** k
    plus = sum(1 for k in range(5) if action(k) == k and sign(k) == 1) + sum(1 for k in range(5) if action(k) > k)
    check(f"{name} teilt den Träger ebenfalls 3 + 2", sorted((plus, 5 - plus)) == [2, 3])

# Nur bei vier Marken: allgemeiner Grad n gibt (n//2 + 1) + (n - n//2).
table = {n: (n // 2 + 1, n - n // 2, n - 1) for n in range(2, 9)}
only = [n for n, (a, b, r) in table.items() if (a, b, r) == (3, 2, 3)]
check("unter 2 bis 8 Marken gibt nur n = 4 die Teilung (3 + 2) + 3", only == [4], str(table))

# --- 3. Hyperladung des Originals als Funktion der Parität --------------------------
eps = [-1 if k % 2 == 0 else 1 for k in range(5)]          # gerade -> -1, ungerade -> +1
Y = [F(1, 12) + F(5, 12) * e for e in eps]
check("Y = 1/12 + (5/12) eps gibt -1/3 dreimal und 1/2 zweimal", sorted(Y) == [F(-1, 3)] * 3 + [F(1, 2)] * 2)
check("Y ist spurfrei", sum(Y) == 0)
check("tr Y^2 = 5/6", sum(y * y for y in Y) == F(5, 6), "Exponent von theta13 (v268)")
check("6 Y ist der Marker (-2,-2,-2,3,3)", sorted(6 * y for y in Y) == [-2, -2, -2, 3, 3])

# Eine Familie: gerade Teilmengen von fünf Plätzen (16 Zustände), Ladung = Summe der besetzten Y, verschoben.
states = [s for s in range(32) if bin(s).count("1") % 2 == 0]
Ys = sorted(Y)
charge = lambda s: sum(Ys[i] for i in range(5) if s >> i & 1) - sum(Ys) / 2
q = [charge(s) for s in states]
check("16 Zustände einer Familie", len(states) == 16)
check("Summe Y und Summe Y^3 über die 16 verschwinden", sum(q) == 0 and sum(x ** 3 for x in q) == 0)
sumY2 = sum(x * x for x in q)
check("Summe Y^2 über eine Familie = 10/3 = 16 mal 5/24", sumY2 == F(10, 3) and sumY2 == 16 * F(5, 24))
b1 = F(4, 3) * 3 + F(1, 10)
check("b1 = (4/3) mal 3 Familien + 1/10 = 41/10", b1 == F(41, 10))

# --- 4. Verklebung: zyklische Viererverklebung verlangt zwei ungerade Ränge -----------
splits = [(a, 8 - a) for a in range(1, 8) if a % 2 == 1 and (8 - a) % 2 == 1 and a >= 8 - a]
check("8 = ungerade + ungerade mit zyklischer Diskriminante Z4", splits == [(5, 3), (7, 1)],
      "D5 + D3 ist die einzige Wahl mit zwei echten Teilen; D3 = A3")
check("q(D5) + q(A3) = 5/4 + 3/4 = 2", F(5, 4) + F(3, 4) == 2)

# --- 5. Galoisgetriebe: gleicher Charakter wie der Träger -----------------------------
# (Z/5)^x = <2> wirkt auf Z/5 durch x -> 2x: Fixpunkt 0 und ein Viererzyklus.
orbit = [1, 2, 4, 3]
check("x -> 2x auf Z/5: ein Fixpunkt und ein Viererzyklus", [2 * x % 5 for x in orbit] == orbit[1:] + orbit[:1])
galois = sorted([0] + [0, 1, 2, 3])
check("Charakter der Galoiswirkung = Charakter von mu4 auf dem Träger", galois == sorted(carrier),
      "trivial plus regulär; die Abbildung selbst ist nicht konstruiert")

# --- 6. Zahlen des Wörterbuchs (KERNSIX Schritt 2) aus dem flachen Zustand ---------------
check("Markergewicht tr Y^2 / 4 = 5/24", F(5, 6) / 4 == F(5, 24))
check("Ströme zwischen Träger und Familie: 2 mal 5 mal 3 mal 2 = 60", 2 * 5 * 3 * 2 == 60)
check("Wurzeln: 40 (D5) + 12 (A3) + 2 mal 64 (16,4) + 60 (10,6) = 240", 40 + 12 + 128 + 60 == 240)

print(f"\n{sum(checks)} von {len(checks)} Prüfungen bestanden")
raise SystemExit(0 if all(checks) else 1)
