#!/usr/bin/env python3
"""A3 -- Anomaliefreiheit einer reflexionspositiven Theorie mit mu_4.

Prinzip: Die physische Theorie ist Rand einer invertierbaren, reflexionspositiven
Theorie; ihre Gesamtanomalie muss verschwinden. Mit einer Z_4, deren Quadrat die
Fermionparitaet ist (Skalar i auf den Fermionen), misst die eta Invariante des
getwisteten Diracoperators auf RP^(2k-1) die Anomalie eines Weylfermions in
d = 2k - 2 Dimensionen. Fixpunktformel (Atiyah, Patodi, Singer; Gilkey), exakt.

Ausfuehren: python3 a3_anomalie_eta.py
"""
from fractions import Fraction as Fr

import sympy as sp

OK = True


def check(name, cond, detail=""):
    global OK
    OK &= bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}{(' -- ' + detail) if detail else ''}")


def eta_rp(k, s):
    """eta auf S^(2k-1)/Z_2: Spinlift von -1 hat mu = i je komplexer Richtung,
    Twist mit Z_4 Ladung s (k ungerade: s ungerade) oder Vorzeichen (k gerade: s gerade)."""
    val = sp.Rational(1, 2) * sp.I ** s / (2 * sp.I) ** k
    val = sp.nsimplify(sp.expand(val))
    assert val.is_real
    return Fr(int(sp.numer(val)), int(sp.denom(val)))


print("1. eta(RP^(2k-1)) und Ordnung der Anomalie eines Weylfermions")
LIT = {2: "Z_8 (Pin-)", 4: "Z_16 (Pin+)", 6: "Z_16 (Pin-)", 8: "Z_2 x Z_32 (Pin+)"}
rows = []
for k in range(2, 6):
    s = 1 if k % 2 else 2
    e = eta_rp(k, s)
    rows.append((2 * k - 2, e, e.denominator))
    print(f"  d = {2 * k - 2}: eta = {e}, Ordnung {e.denominator}; Literatur Bordismus {LIT[2 * k - 2]}")
check("d = 4: eta = +-1/16, Ordnung 16", rows[1][1] in (Fr(1, 16), Fr(-1, 16)))
check("Ladung 3 (das konjugierte Fermion) traegt die entgegengesetzte Anomalie",
      eta_rp(3, 3) == -eta_rp(3, 1))
print("  Hinweis: fuer d >= 6 ist die Gruppe der Literatur um den Faktor 2 kleiner als die Ordnung der")
print("  komplexen eta Phase (reelle oder quaternionische Struktur der Spinoren). Mit mu_4 Twist")
print("  (k ungerade, d = 0, 4, 8, ...) tritt 16 nur bei d = 4 auf; ohne Twist auch bei d = 6.")

print("\n2. Kontrolle: Z_n Symmetrie in d = 4 auf Linsenraeumen S^5/Z_n (n ungerade, Spin)")


def eta_lens(n, s, k=3):
    tot = 0
    for j in range(1, n):
        mu = sp.exp(sp.I * sp.pi * sp.Rational(j * (n + 1), n))  # Spinlift der Ordnung n
        lam = sp.exp(2 * sp.I * sp.pi * sp.Rational(j, n))
        tot += lam ** s / (mu - 1 / mu) ** k
    val = sp.nsimplify(sp.N(-tot / n, 60), rational=True, tolerance=sp.Rational(1, 10 ** 40))
    assert abs(sp.N(-tot / n - val, 50)) < sp.Rational(1, 10 ** 40)
    return Fr(int(sp.numer(val)), int(sp.denom(val)))


orders = {}
for n in (3, 5, 7, 9):
    vals = [eta_lens(n, s) for s in range(1, n)]
    den = 1
    for v in vals:
        den = den * v.denominator // sp.igcd(den, v.denominator)
    orders[n] = den
    print(f"  n = {n}: eta(s) = {[str(v) for v in vals]}, gemeinsamer Nenner {den}")
check("kein ungerades n bis 9 erzeugt die Ordnung 16 (Basislinie fuer die 16)", 16 not in orders.values(),
      f"Nenner {orders}")

print("\n3. Eine Generation: mu_4 ist das Zentrum von Spin(10)")
# X = 5 - 2p auf Lambda^p C^5 (SU(5) x U(1)_X): 1 -> 5, 10 -> 1, 5bar -> -3
charges = {p: (5 - 2 * p) % 4 for p in (0, 2, 4)}
check("auf Lambda^gerade C^5 ist X = 5 - 2p immer 1 mod 4: die Z_4 wirkt als Skalar i", set(charges.values()) == {1})
e1 = eta_rp(3, 1)
dims = {0: 1, 2: 10, 4: 5}
check("16 Weylfermionen (mit nu_R): Anomalie 16 * eta ist ganzzahlig", (sum(dims.values()) * e1).denominator == 1,
      f"16 * ({e1}) = {sum(dims.values()) * e1}")
check("15 Weylfermionen (ohne nu_R): Anomalie bleibt", (15 * e1).denominator != 1, f"15 * ({e1}) = {15 * e1}")
print("  Reflexionspositivitaet + mu_4 + vier Dimensionen verlangen Vielfache von 16 Weylfermionen;")
print("  Lambda^gerade C^5 ist die kleinste Loesung. Die Familienzahl bleibt dabei frei.")

print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
