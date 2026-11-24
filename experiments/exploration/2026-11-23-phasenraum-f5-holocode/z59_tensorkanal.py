"""Z59: Tensorkanal und Raumzeitrichtungen aus den Antworten des Originals (Fortsetzung der Rekonstruktion, z58).

Teil 1  Tensorkanal in 2D: Aus der originalen Kovarianz c(x) = 1/(1 + e^{4 pi x}) folgen die Virasoro-Antworten
        D_m = omega([L_m, L_-m]) = sum_x (x + m/2)^2 (c(x) - c(x+m)). Sie haben exakt die Form 2m(<L_0>_beta + h_alpha)
        + (c/12)(m^3 - m): die Zentralladung c = 1 je Quarter und das Twistgewicht h_alpha = (alpha - 1/2)^2/2 = 1/32 werden
        rekonstruiert, nicht gesetzt. c/2 ist die spurfreie (Spin-2-)Antwort <T T>.
Teil 2  Gravitation der Naht in 2D: Spurkanal (Fisher-Information, z58) und Tensorkanal (c) sind beide Antworten des
        Originals; die induzierte Wirkung ist die Polyakov-Wirkung mit Koeffizient c, die Weyl-Anomalie -c R/(24 pi).
        Reflexionspositivität in modularer Zeit (Gram-Matrix der rekonstruierten Stromkorrelation) geprüft.
Teil 3  Raumzeitrichtungen: Die Antworten des Originals sind Funktionen einer einzigen modularen Zeit und hängen nur von
        der Mode x ab (ein chiraler Zweig): Sie rekonstruieren genau eine positive lichtartige Translation. 3+1 braucht vier
        vertauschende; die vier Marken liefern Nullrichtungen (Tetraeder), aber keine Dynamik. Offen.
Teil 4  Bilanz.
"""
import json
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).parent
mp.mp.dps = 50
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


PI4 = 4*mp.pi
ALPHAS = (mp.mpf(1)/4, mp.mpf(3)/4)
R = 40


def occ(x):
    return 1/(1 + mp.exp(PI4*x))


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Tensorkanal in 2D -- Virasoro-Antworten des Originals")
tol = mp.mpf(10)**-35
ok_form, ok_c, ok_h = True, True, True
for a in ALPHAS:
    xs = [r + a for r in range(-R, R)]
    D = {m: mp.fsum((x + mp.mpf(m)/2)**2*(occ(x) - occ(x + m)) for x in xs) for m in range(1, 6)}
    L0 = mp.fsum(x*(occ(x) - (1 if x < 0 else 0)) for x in xs)
    h = (a - mp.mpf(1)/2)**2/2
    c_rec = 2*(D[2] - 2*D[1])
    A_rec = D[1]/2
    ok_c &= abs(c_rec - 1) < tol
    ok_h &= abs(A_rec - (L0 + h)) < tol and abs(h - mp.mpf(1)/32) < tol
    ok_form &= all(abs(D[m] - (2*m*A_rec + c_rec*(m**3 - m)/12)) < tol for m in range(1, 6))
    ZAHLEN[f"alpha={mp.nstr(a, 3)}"] = {"c": mp.nstr(c_rec, 20), "<L0>_beta": mp.nstr(L0, 12), "h": mp.nstr(h, 6)}
check("D_m = omega([L_m, L_-m]) hat für m = 1..5 exakt die Virasoro-Form 2m A + (c/12)(m^3 - m) (beide Quarter)", ok_form)
check("rekonstruierte Zentralladung c = 2(D_2 - 2 D_1) = 1 je Quarter (bis 1e-35)", ok_c)
check("A = D_1/2 = <L_0>_beta + h_alpha mit dem Twistgewicht h_alpha = (alpha - 1/2)^2/2 = 1/32 (beide Quarter)", ok_h)
befund("Die spurfreie Antwort <T T> = c/2 ist damit eine Zahl des Originals: c = 1 je Quarter, c = 2 für die Quarterquelle", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Gravitation der Naht in 2D -- Spur- und Tensorkanal aus denselben Antworten")
z58 = json.loads((HERE/"z58_ergebnis.json").read_text())
kappa2 = float(z58["zahlen"]["energiekumulanten_K"]["kappa2"])
check("Spurkanal (z58): Fisher-Information kappa_2 > 0; Tensorkanal: c > 0 -- beide positiv ohne Entwurf",
      kappa2 > 0 and ok_c)
taus = [mp.mpf(k)/20 for k in range(1, 10)]                        # modulare Zeiten in (0, 1/2), beta_mod = 1
for a in ALPHAS:
    xs = [r + a for r in range(-R, R)]

    def G(t):        # omega(J_1(-i t) J_-1) in modularer Zeit, aus den Fermionsummen
        return mp.fsum(occ(x)*(1 - occ(x + 1)) for x in xs)*mp.exp(-PI4*t) + mp.fsum(
            occ(x + 1)*(1 - occ(x)) for x in xs)*mp.exp(-PI4*(1 - t))
    gram = np.array([[float(G(ti + tj)) for tj in taus] for ti in taus])
    ev = np.linalg.eigvalsh(gram)
    if a == ALPHAS[0]:
        rp = ev.min() > -1e-12*ev.max()
    else:
        rp &= ev.min() > -1e-12*ev.max()
check("Reflexionspositivität in modularer Zeit: die Gram-Matrix G(t_i + t_j) der rekonstruierten Stromkorrelation ist "
      "positiv semidefinit (Osterwalder-Schrader in einer Richtung)", rp)
befund("Die induzierte 2D-Gravitation der Naht ist damit ganz rekonstruiert: Polyakov-Wirkung (c/96 pi) Int R (1/Box) R "
       "und Weyl-Anomalie -c R/(24 pi) mit c aus den Antworten (Quarterquelle c = 2, E8-Naht c = 8); in 2D hat sie keine "
       "propagierenden Gravitonen", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Raumzeitrichtungen -- wie viele Translationen rekonstruieren die Antworten?")
befund("Alle Antworten des Originals hängen nur von der Mode x und einer modularen Zeit ab (ein chiraler Zweig, z58): Sie "
       "rekonstruieren genau eine positive lichtartige Translation (Borchers: eine halbseitige modulare Inklusion)", True)
v = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])/np.sqrt(3)
n = np.hstack([np.ones((4, 1)), v])                                   # vier Nullvektoren (Tetraeder)
eta = np.diag([1, -1, -1, -1])
null_ok = all(abs(nk @ eta @ nk) < 1e-12 for nk in n)
rng = np.random.default_rng(59)
cone_ok = True
for _ in range(2000):
    p = rng.random(4)
    P = p @ n
    m2 = P @ eta @ P
    cone_ok &= m2 >= -1e-12 and P[0] > 0 and abs(m2 - (4/3)*(p.sum()**2 - (p**2).sum())) < 1e-10
rank = np.linalg.matrix_rank(n)
check("die vier Marken als Tetraeder-Nullrichtungen spannen R^{1,3} (Rang 4); positive Beiträge p_k längs n_k geben "
      "stets einen zukunftsgerichteten Impuls mit M^2 = (4/3)((Sum p)^2 - Sum p^2) >= 0", null_ok and rank == 4 and cone_ok)
befund("Das ist Kinematik, keine Rekonstruktion: Dass vier vertauschende positive Translationen aus den Antworten des "
       "Originals kommen, ist nicht gezeigt; die Quarterquelle allein trägt nur eine. 3+1 bleibt offen und braucht eine "
       "4D-Familie von Antworten (z. B. je Marke ein chirales Netz mit vertauschenden Translationen)", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Bilanz")
tafel = [
    ("Tensorkanal in 2D", "rekonstruiert: c = 1 je Quarter aus den Virasoro-Antworten, <TT> = c/2"),
    ("Spurkanal in 2D", "rekonstruiert: Fisher-Information (z58)"),
    ("2D-Gravitation der Naht", "vollständig: Polyakov-Wirkung mit c aus den Antworten"),
    ("Tensorkanal in 4D / t_n", "offen: braucht 4D-Antworten"),
    ("3+1", "offen: die Quelle trägt eine lichtartige Translation, gebraucht werden vier vertauschende"),
    ("Familiengewichte gamma", "offen: im Bestand keine Antwort des Originals, die gamma festlegt"),
]
for k_, v_ in tafel:
    say(f"   {k_}: {v_}")
ZAHLEN["bilanz"] = dict(tafel)
urteil = ("Der Tensorkanal ist in 2D aus den Antworten des Originals rekonstruiert: Die Virasoro-Antworten "
          "omega([L_m, L_-m]) der originalen Kovarianz haben exakt die Form 2m(<L_0> + h) + (c/12)(m^3 - m) mit c = 1 je "
          "Quarter und Twistgewicht h = 1/32. Zusammen mit der Fisher-Information (z58) ist die induzierte 2D-Gravitation der "
          "Naht vollständig festgelegt (Polyakov-Wirkung mit c), reflexionspositiv in modularer Zeit. Für 3+1 reicht die "
          "Quelle nicht: Ihre Antworten rekonstruieren genau eine lichtartige Translation; die vier Marken geben nur die "
          "Nullrichtungen (Kinematik). Offen bleiben der 4D-Tensorkanal, 3+1 und die Familiengewichte.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z59_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
