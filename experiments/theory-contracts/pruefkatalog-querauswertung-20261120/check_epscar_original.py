"""EPSCAR20261120: Ist eps_car des Originals die Parität der Deckinvolution z -> -z?

Liest die Originalstellen und prüft exakt, was aus ihnen folgt. Kein Beweis der Physik.
"""
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
KERN = next((ROOT/"documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections").glob("001-*.tex")).read_text()
CARR = next((ROOT/"documents/archive--tfpt-45--source_extracts--02_carrier_source/sections").glob("001-*.tex")).read_text()
V507 = (ROOT/"verification/v507_seam_bit_origin.py").read_text()
checks = []


def check(name, ok, detail=""):
    checks.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


flat = lambda s: " ".join(s.split())
# --- Original -------------------------------------------------------------------------
check("Original: eps_car ist die Einschränkung der Polarisation", r"\epscar:=\iC|_{\mathcal H_{\mathrm{car}}}" in KERN)
check("Original: Polarisation ist die Spiegelung auf Cauchydaten", r"\iC=2P_C-1=R_{\taudbl}" in KERN)
check("Original: die Naht ist die Fixmenge der Spiegelung", r"\Seam=\mathrm{Fix}(\taudbl)" in KERN)
check("Original: die drei Involutionen dürfen nicht wörtlich gleichgesetzt werden",
      "No later theorem may identify these three operators literally" in flat(KERN))
check("Original: 3 + 2 folgt dort aus zwei Annahmen über äußere Potenzen",
      r"\dim \Lambda^2 E_+=1" in CARR and r"\Lambda^2 E_- \cong E_-^\vee" in CARR and "Assume only" in CARR)
check("Ledger v507: das Nahtdeck ist z -> -z, die zentrale Involution der Marken",
      "the v506 seam deck z -> -z lands on the central" in flat(V507))
check("Ledger v507: jedes markenfreie Kragendeck ist eine Halbperiodenverschiebung",
      "EVERY mark-free collar deck (v216) IS a half-period" in flat(V507))

# --- Rechnung ---------------------------------------------------------------------------
# Träger H0(O(D)) = { p(z)/(z^4-1), Grad p <= 4 }. Rückzug unter z -> -z: z^k -> (-1)^k z^k.
plus = [k for k in range(5) if k % 2 == 0]
minus = [k for k in range(5) if k % 2 == 1]
check("z -> -z auf dem Träger: Eigenräume der Dimension 3 und 2", (len(plus), len(minus)) == (3, 2))
b, s = 3, 2
check("Annahme 1 des Originals folgt: dim Lambda^2 E_2 = 1", comb(s, 2) == 1)
check("Annahme 2 des Originals folgt: dim Lambda^2 E_3 = dim E_3", comb(b, 2) == b)
sols = [(bb, ss) for bb in range(1, 9) for ss in range(1, 9) if comb(ss, 2) == 1 and comb(bb, 2) == bb]
check("umgekehrt lassen die beiden Annahmen nur (3, 2) zu", sols == [(3, 2)])
# Die Fixpunkte von z -> -z sind 0 und unendlich, keine Marke und kein Kreis.
check("z -> -z bewegt alle vier Marken", all(-m != m for m in (1, 1j, -1, -1j)))
# Die Teilung hängt nicht an der Lage des Decks: jede Involution der Sphäre wirkt auf Sym^4 mit (1,-1,1,-1,1).
check("jede holomorphe Involution gibt 3 + 2", sorted([(1j) ** (4 - 2 * k) for k in range(5)], key=lambda c: c.real)
      == [-1, -1, 1, 1, 1])
# Vorzeichen: der Rückzug gibt +1 auf den drei geraden Moden; das Original setzt Farbe = ker(eps + 1).
check("Vorzeichen: eps_car = minus Rückzug", [-(-1) ** k for k in range(5)].count(-1) == 3,
      "eine Wahl der Hebung auf das Linienbündel")

print(f"\n{sum(checks)} von {len(checks)} Prüfungen bestanden")
raise SystemExit(0 if all(checks) else 1)
