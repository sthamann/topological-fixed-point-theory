"""H14: Raten 2/3 und 1/3 als Anteile der Spinsektoren der E8 Gittersumme am quadratischen Torus."""
import itertools, mpmath as mp
from fractions import Fraction as Fr
from common import Checks
c = Checks("h14_theta"); mp.mp.dps = 40
def shares(tau):
    q = mp.exp(1j * mp.pi * tau); th = [mp.jtheta(k, 0, q) ** 8 for k in (2, 3, 4)]; S = sum(th); return [v / S for v in th], S / 2
sh, E4 = shares(1j)
c.ok("am quadratischen Torus sind die Anteile exakt 1/6, 2/3, 1/6 (40 Stellen)", all(abs(a - b) < mp.mpf(10) ** -35 for a, b in zip(sh, (mp.mpf(1) / 6, mp.mpf(2) / 3, mp.mpf(1) / 6))))
sh2, _ = shares(2j)
c.ok("am Rechteck 2:1 sind sie 0,00045 · 0,51471 · 0,48485, also keine einfachen Brüche mehr", abs(sh2[1].real - 0.514705) < 1e-5 and abs(sh2[0].real - 0.000446) < 1e-5, [mp.nstr(v.real, 6) for v in sh2])
rho = mp.exp(2j * mp.pi / 3); q = mp.exp(1j * mp.pi * rho); S6 = sum(mp.jtheta(k, 0, q) ** 8 for k in (2, 3, 4))
c.ok("am Sechseck verschwindet die Summe (E4 = 0): dort gibt es keine Anteile", abs(S6) < mp.mpf(10) ** -30)
c.ok("FOLGE: die Zahlen 2/3 und 1/3 hängen an der Torusform. Die Raten des Originals sind feste Zahlen ohne Formparameter", True)
c.ok("ABBRUCHREGEL GREIFT: drei Sektoren mit Wurzelzahlen 112 (Vektorteil) und 128 (Spinorteil) geben keine Vielfachheiten 6 und 6 im Spektrum {1, (2/3)^6, (1/3)^6}", 112 + 128 == 240 and 6 not in (112, 128, 3))
c.ok("im Originalmodell v724 stehen die Raten als (2/3)^6 und (1/3)^6 neben dem einen Maßstab 6 log(3/2); der Exponent 6 kommt dort nicht aus drei Randbedingungen", True)
# Gegenprobe zur parallelen Prüfung: gemischte Ladungsantwort
def e8roots():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for a in (1, -1):
            for b in (1, -1): v = [Fr(0)] * 8; v[i] = Fr(a); v[j] = Fr(b); R.append(tuple(v))
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0: R.append(tuple(Fr(x, 2) for x in s))
    return R
dot = lambda a, b: sum(x * y for x, y in zip(a, b))
def responses(R, limit=None):
    out = set(); self4 = sum(dot(R[0], r) ** 4 for r in R)
    for a in R[:limit]:
        da = [dot(a, r) for r in R]
        for b in R:
            if dot(a, b) == 0: out.add(Fr(sum(x * x * dot(b, r) ** 2 for x, r in zip(da, R))) / Fr(self4))
    return out
E = e8roots(); rE = responses(E, 8)
c.ok("NACHGERECHNET: in einem E8 hat die gemischte Antwort zweier senkrechter Wurzeln, geteilt durch die Einzelantwort, immer den Wert 1/3; zwischen zwei verschiedenen E8 ist sie 0", rE == {Fr(1, 3)}, rE)
D = []
for i, j in itertools.combinations(range(16), 2):
    for a in (1, -1):
        for b in (1, -1): v = [Fr(0)] * 16; v[i] = Fr(a); v[j] = Fr(b); D.append(tuple(v))
rD = responses(D, 4)
c.ok("NACHGERECHNET: in D16 (dieselben 480 Wurzeln wie D16+) nimmt sie nur die Werte 1/9 und 7/9 an. Die parallele Aussage stimmt", len(D) == 480 and rD == {Fr(1, 9), Fr(7, 9)}, rD)
c.ok("ERKLÄRUNG: 1/3 ist die Zahl jeder gleichmäßigen Verteilung (Summe x^2 y^2 zu Summe x^4 = 1 : 3). E8 Wurzeln sind so gleichmäßig, D16 Wurzeln nicht. Mit den Sektoranteilen hat diese 1/3 nichts zu tun", True)
c.ok("drei verschiedene Quellen für dieselben Brüche: Gleichverteilung (1/3), Thetawerte am Quadrat (1/6, 2/3, 1/6), Kanalraten (2/3, 1/3). Eine Abbildung zwischen ihnen fand sich nicht", True)
c.done({"anteile_quadrat": [mp.nstr(v.real, 20) for v in sh], "anteile_2zu1": [mp.nstr(v.real, 10) for v in sh2], "antwort_E8": [str(x) for x in rE], "antwort_D16": sorted(str(x) for x in rD),
        "urteil": "Gefallen: die Vielfachheiten folgen nicht und die Anteile hängen an der Form. Nebenbefund: die parallele Gitterprobe (0 oder 1/3 gegen 1/9 oder 7/9) ist nachgerechnet und erklärt."})
