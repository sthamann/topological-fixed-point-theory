"""Z172: das Gesetz als Eigenwertgleichung auf den drei Paarungen und die Netz-Suche über Verklebungen.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, C bedingt, L Lesart, N numerisch, O offen.

Die Suchen stehen in research/originlab/whole_solution_search.py (Universen Gesetzesformen und Verklebungen).

Teil A  Die Alpha-Gleichung als Eigenwertgleichung eines Transferoperators auf den drei Paarungen der vier Marken
Teil B  Gesetzesform-Suche: 1092 Transferoperatoren gegen den Messwert
Teil C  Netz-Suche: Naht mit Partnern, Raum der Geraden, reelle Strukturen
Teil D  Was Gesetz und Netz verbindet
Teil E  Bilanz und Urteil.
"""
import json
import os
import sys
import time
from itertools import permutations
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import whole_solution_search as ws                                         # noqa: E402

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 40
AINV = mp.mpf("137.0359992168407")
c3 = 1 / (8 * mp.pi)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


say("Teil A  Die Alpha-Gleichung als Eigenwertgleichung auf den drei Paarungen der vier Marken")
P = sp.Matrix(ws.PAIRINGS).T                                               # Spalten: {01|23}, {02|13}, {03|12}
check("A1 [E] die drei Paarungen sind eine orthogonale Basis der Richtungen mit Summe null",
      P.T * P == 4 * sp.eye(3) and sp.ones(1, 4) * P == sp.zeros(1, 3))
acts = {perm: sp.Matrix(ws.pairing_action(perm)) for perm in permutations(range(4))}
check("A2 [E] S4 wirkt auf den Paarungen durch 24 verschiedene vorzeichenbehaftete Permutationen; die Kleinsche "
      "Vierergruppe wirkt diagonal, der Quotient S3 vertauscht die Paarungen",
      len({tuple(m) for m in acts.values()}) == 24
      and all(m.applyfunc(abs).is_diagonal() == (sorted(p) == list(p) or all(p[p[i]] == i and p[i] != i
                                                                              for i in range(4)))
              for p, m in acts.items())
      and len({tuple(m.applyfunc(abs)) for m in acts.values()}) == 6)


def circ(row):
    return sp.Matrix(4, 4, lambda i, j: row[(j - i) % 4])


B, Cc, U = circ([2, -1, 0, -1]), circ([0, 1, 2, 1]), sp.diag(1, -1, 1, -1)
PI4 = (B + U * Cc * U) / 8
PI = sp.diag(0, 1, 0)
check("A3 [E] der Projektor Pi der Analyse 25 ist der Projektor auf die Paarung gegenüberliegender Marken {02|13}",
      PI4 == P * PI * P.T / 4)
C3cyc = acts[(0, 2, 3, 1)]
quarter = acts[(1, 2, 3, 0)]
check("A4 [E] der Dreierzyklus, der eine Marke festhält, vertauscht die drei Paarungen zyklisch ohne Vorzeichen",
      C3cyc == sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]]) and C3cyc ** 3 == sp.eye(3))
t, lam, cc = sp.symbols("t lambda c3", positive=True)
check("A5 [E] det(t - c3 Pi - lambda C) = t^2 (t - c3) - lambda^3, für beide Umlaufrichtungen",
      all(sp.expand((t * sp.eye(3) - cc * PI - lam * M).det() - (t ** 2 * (t - cc) - lam ** 3)) == 0
          for M in (C3cyc, C3cyc.T)))
check("A6 [E] die Vierteldrehung der Marken hält die Paarung {02|13} bis aufs Vorzeichen fest und vertauscht die "
      "anderen beiden; die Symmetrie der Naht (D4) enthält den Dreierzyklus nicht",
      quarter == sp.Matrix([[0, 0, 1], [0, -1, 0], [-1, 0, 0]])
      and len({tuple(acts[p].applyfunc(abs)) for p in permutations(range(4))
               if all((p[(i + 1) % 4] - p[i]) % 4 in (1, 3) for i in range(4))}) == 2)
N, b1, pw = 48, mp.mpf(41) / 10, mp.mpf(5) / 4
seam = lambda al: mp.mpf(4) / 3 * c3 + N * c3 ** 4 * mp.e ** (-2 * al) * (1 - N * c3 ** 4 * mp.e ** (-2 * al)) ** (-pw)
al = 1 / AINV
for _ in range(40):
    lam_n = (b1 * mp.log(1 / seam(al))) ** (mp.mpf(1) / 3)
    T = mp.matrix([[0, 0, lam_n], [lam_n, c3, 0], [0, lam_n, 0]])
    ev = mp.eig(T, left=False, right=False)
    tt = max(mp.re(x) for x in ev if abs(mp.im(x)) < mp.mpf("1e-30"))
    al = 2 * c3 ** 2 * tt
ZAHLEN["t"], ZAHLEN["lambda"] = mp.nstr(tt, 15), mp.nstr(lam_n, 15)
others = [x for x in ev if abs(mp.im(x)) > mp.mpf("1e-30")]
check("A7 [N] alpha = 2 c3^2 mal größter Eigenwert von T = c3 Pi + lambda C mit lambda^3 = b1 ln(1/phi_seam): "
      "alpha^-1 = 137,0359992168407; t = 2,30470345485719",
      abs(1 / al - AINV) < mp.mpf("1e-12") and abs(tt - mp.mpf("2.30470345485719")) < mp.mpf("1e-13"))
T5 = T ** 5
check("A8 [E/N] T ist nichtnegativ und primitiv (T^5 > 0): nach Perron-Frobenius genau ein positiver Eigenwert, "
      "einfach und betragsgrößter; die anderen beiden sind komplex",
      all(T5[i, j] > 0 for i in range(3) for j in range(3)) and len(others) == 2
      and all(abs(x) < tt for x in others))
z = 1 / tt
check("A9 [N] Erneuerungsform: 1 = c3 z + b1 ln(1/phi_seam) z^3 bei z = 1/t = c3/e^2; Eigenvektor mit Verhältnis "
      "lambda/t = (1 - c3/t)^(1/3) = 0,9942",
      abs(c3 * z + lam_n ** 3 * z ** 3 - 1) < mp.mpf("1e-30")
      and abs(lam_n / tt - (1 - c3 / tt) ** (mp.mpf(1) / 3)) < mp.mpf("1e-30") and abs(lam_n / tt - 0.9942) < 1e-4)
tau = mp.mpc("0.3", "1.1")
th8 = lambda k, x: mp.jtheta(k, 0, mp.e ** (1j * mp.pi * x)) ** 8
close = lambda a, b: abs(a / b - 1) < mp.mpf("1e-25")
s_ok = close(th8(3, -1 / tau), tau ** 4 * th8(3, tau)) and close(th8(2, -1 / tau), tau ** 4 * th8(4, tau)) \
    and close(th8(4, -1 / tau), tau ** 4 * th8(2, tau))
t_ok = close(th8(2, tau + 1), th8(2, tau)) and close(th8(3, tau + 1), th8(4, tau)) \
    and close(th8(4, tau + 1), th8(3, tau))
eta8 = lambda x: (mp.e ** (1j * mp.pi * x / 12) * mp.qp(mp.e ** (2j * mp.pi * x))) ** 8
check("A10 [E/N] dieselbe S3 wirkt als Modulgruppe auf den drei geraden Thetakanälen von E8: S hält theta_3^8 fest "
      "und vertauscht theta_2^8, theta_4^8; T hält theta_2^8 fest; ST ist ein Dreierzyklus. Die Phase von T auf "
      "theta^8/eta^8 ist exp(-2 pi i/3), das Umlaufgewicht von (ST)^3 ist eins",
      s_ok and t_ok and close(eta8(tau + 1), mp.e ** (2j * mp.pi / 3) * eta8(tau))
      and close(eta8(-1 / tau), tau ** 4 * eta8(tau)))
befund("Gesetz", "Eigenwertgleichung",
       "Die Alpha-Gleichung ist genau die Aussage: alpha/(2 c3^2) ist der Perron-Eigenwert des Transferoperators "
       "T = c3 Pi + lambda C auf den drei Paarungen der vier Marken. Pi ist die Paarung gegenüberliegender Marken, "
       "C der Dreierzyklus bei einer festgehaltenen Marke (4 = 3 + 1); auf den Thetakanälen von E8 sind das der von S "
       "festgehaltene Kanal und die Modultransformation ST. Die Eindeutigkeit der Wurzel folgt aus "
       "Perron-Frobenius. Nicht hergeleitet: das Umlaufgewicht lambda^3 = b1 ln(1/phi_seam), die Stärke c3 der "
       "Schleife und warum alpha der Eigenwert ist.")

say("\nTeil B  Gesetzesform-Suche: 1092 Transferoperatoren gegen den Messwert")
law = ws.run_law_form_search()
HIT = "c3 on opposite | 3-cycle | b1 ln(1/phi_seam) | alpha = 2 c3^2 t"
check("B1 [N] 1092 Formen (Diagonale, Stärke, Transport, Gewicht, Einheit): genau eine trifft den Messwert",
      law["objects"] == 1092 and law["hits"] == [HIT])
miss = list(law["nearest_miss"].items())
ZAHLEN["gesetz_naechste_fehlversuche"] = {k: float(v) for k, v in miss}
check("B2 [N] nächster Fehlversuch (phi_0 statt phi_seam) 24 000 Sigma, alle anderen über 100 000 Sigma",
      miss[0][0].count("phi_0") == 1 and abs(miss[0][1] - 23918) < 5 and all(abs(v) > 1e5 for _, v in miss[1:]))
check("B3 [E] Regelmaschine: 'alpha gemessen' allein wählt die Form eindeutig",
      any(r["selected"] == HIT and ["alpha_measured"] in r["minimal_rules"]
          for r in law["search"]["unique_selections"]))
mod = [o for o in ws.law_form_candidates()
       if o["label"] == "c3 on opposite | 3-cycle | 1 (modular) | alpha = 2 c3^2 t"][0]
ZAHLEN["alpha_inv_modulares_gewicht"] = mod["alpha_inverse"]
check("B4 [N negativ] mit dem nativen modularen Umlaufgewicht eins käme alpha^-1 = 311,6: die Gestalt des "
      "Transports ist modular (ST), sein Gewicht nicht", abs(mod["alpha_inverse"] - 311.6) < 0.1)
ratio = (mp.jtheta(3, 0, mp.e ** (-mp.pi)) / mp.jtheta(2, 0, mp.e ** (-mp.pi))) ** 8
t4 = mp.findroot(lambda x: x ** 2 * (x - c3) - b1 * mp.log(ratio), 1.8)
ZAHLEN["alpha_inv_kanalverhaeltnis"] = mp.nstr(1 / (2 * c3 ** 2 * t4), 6)
check("B5 [N negativ] das Verhältnis der Thetakanäle bei tau = i ist theta_3^8/theta_2^8 = 4; als Umlaufgewicht "
      "b1 ln 4 gäbe das alpha^-1 = 175,7: phi_seam ist kein Verhältnis der Kanalgewichte",
      abs(ratio - 4) < mp.mpf("1e-30") and abs(1 / (2 * c3 ** 2 * t4) - mp.mpf("175.656")) < 0.01)
befund("Gesetzesform", "eindeutig im Katalog",
       "Der Messwert legt die Form im Katalog fest; kein anderer Transferoperator aus denselben Bausteinen kommt "
       "näher als 24 000 Sigma. Das ist Auswahl durch die Messung, keine Herleitung.")

say("\nTeil C  Netz-Suche: Naht mit Partnern, Raum der Geraden, reelle Strukturen")
comp = ws.run_composition_search()
ZAHLEN["netz"] = {k: comp[k] for k in ("objects", "lorentzian", "euclidean", "real_slices",
                                       "signature_without_real_structure", "common_slice")}
n = sp.symbols("n")
check("C1 [E] 18 Verklebungen (0 bis 3 Partner, hermitesche Form (p,q) oder quaternionische Struktur)",
      comp["objects"] == 18)
check("C2 [E] die Nullgeraden einer hermiteschen Form bilden genau dann einen reellen Schnitt des Geradenraums, "
      "wenn 4(n - 3) = 2(n - 2), also n = 4: Naht plus genau ein Partner", sp.solve(4 * (n - 3) - 2 * (n - 2), n) == [4])
check("C3 [N] genau eine Verklebung gibt ein lorentzsches Netz: ein Partner mit Form (2,2), Klein-Form (2,4), also "
      "konforme Struktur 3+1 in vier Dimensionen", comp["lorentzian"] == ["1 partner | hermitian (2, 2)"])
check("C4 [N] Form (4,0): Klein-Form (6,0), keine reellen Geraden; Form (3,1): keine reelle Struktur",
      comp["signature_without_real_structure"] == {"(4,0)": (6, 0), "(3,1)": (0, 0)})
check("C5 [N] quaternionische Struktur mit einem Partner: euklidisches Netz S^4 (Klein-Form (1,5)); der gemeinsame "
      "Schnitt mit dem lorentzschen Netz hat Klein-Form (1,4), also eine dreidimensionale Sphäre (Zeitpunkt null)",
      comp["euclidean"] == ["1 partner | quaternionic"] and comp["common_slice"] == (1, 4))
befund("Netz", "ein Partner, Form (2,2)",
       "Im Katalog gibt nur die Naht mit genau einem Partner und ausgeglichener Form ein reelles lorentzsches Netz, "
       "und es ist vierdimensional. Das ist die Verklebung des Bestands (Tomita-Partner, CCR-Form (2,2), z79); "
       "die Modultheorie liefert genau einen Partner. Nicht geliefert: die Skala, der Einstein-Koeffizient und "
       "die Prämisse T' (I_pos). Die Firewall DIMENSION.UPLIFT.FIREWALL.01 bleibt bestehen.")

say("\nTeil D  Was Gesetz und Netz verbindet")
X = lambda x: N * c3 ** 4 * mp.e ** (-2 * x)
gl = 8 * mp.pi * X(AINV)
check("D1 [E] mit |Lambda| = 12 eps^2 in nativen Einheiten (z164) und der Lambda-Stufe: eps^2 G/a^2 = 4 c3^3 "
      "exp(-2/alpha)", abs(gl / 12 / (4 * c3 ** 3 * mp.e ** (-2 * AINV)) - 1) < mp.mpf("1e-35"))
lo, hi = mp.log10(gl / 12 * mp.mpf("0.2")), mp.log10(gl / 12 * 24)
ZAHLEN["log10_eps2"] = [mp.nstr(lo, 5), mp.nstr(hi, 5)]
check("D2 [N] das gibt eps^2 = 10^-123,3 bis 10^-121,2 für a^2/G = 0,2 bis 24, wie in z164",
      abs(lo + 123.3) < 0.1 and abs(hi + 121.25) < 0.1)
check("D3 [E] beide reinen Zahlen (Kopplung, eps^2 G/a^2) sind Werte derselben Funktion X bei alpha und 1/alpha",
      abs(gl / (8 * mp.pi) / X(AINV) - 1) < mp.mpf("1e-35") and X(1 / AINV) > X(AINV))
befund("Verbindung", "eine offene Frage",
       "Das Netz gibt Dimension und Signatur, aber keine Zahl. Seine freie Zahl und die Kopplung hängen beide am "
       "Gesetz. Offen bleibt genau: warum das Umlaufgewicht b1 ln(1/phi_seam) ist und warum X bei alpha und bei "
       "1/alpha gelesen wird.")

say("\nTeil E  Bilanz und Urteil")
bilanz = [
    ("E", "Gesetz", "Eigenwertgleichung eines Transferoperators auf den drei Paarungen; Eindeutigkeit aus "
     "Perron-Frobenius."),
    ("N", "Gesetzesform", "eine von 1092 Formen trifft; nächster Fehlversuch 24 000 Sigma."),
    ("E/N", "Netz", "nur ein Partner mit Form (2,2) gibt ein reelles lorentzsches Netz; vier Dimensionen folgen aus "
     "der Schnittbedingung; euklidischer Partner S^4 mit gemeinsamem Schnitt S^3."),
    ("O", "offen", "Umlaufgewicht b1 ln(1/phi_seam), Stärke c3, Auswertung bei alpha und 1/alpha; Skala und "
     "Einstein-Koeffizient; T' (I_pos); Ursprung der Naht (K19, K20)."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Beide Wege sind gerechnet. Das Gesetz hat eine native Operatorform: alpha/(2 c3^2) ist der "
          "Perron-Eigenwert des Transferoperators c3 Pi + lambda C auf den drei Paarungen der vier Marken; unter 1092 "
          "Formen trifft nur diese den Messwert. Hergeleitet ist es damit nicht, das Umlaufgewicht "
          "b1 ln(1/phi_seam) bleibt Eingabe. Das Netz: Nur die Naht mit genau einem Partner und Form (2,2) gibt ein "
          "reelles lorentzsches Netz, und es ist vierdimensional; das bestätigt die Verklebung des Bestands und "
          "macht die Vier zur Folge einer Schnittbedingung, liefert aber weder Skala noch Einstein-Koeffizient. "
          "Beide freien Zahlen hängen an derselben offenen Frage. Keine volle Lösung, Prüfbank unverändert, nichts "
          "befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z172_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z172.txt").write_text("\n".join(OUT) + "\n")
