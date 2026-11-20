"""H15: die Viertelphase der Naht als Abstand der beiden Seiten und als Kommutator; Bindung an ein Ordnungssignal in v724."""
import cmath, itertools, importlib.util, os, numpy as np
from fractions import Fraction as Fr
from common import Checks
c = Checks("h15_phase")
E = lambda x: cmath.exp(2j * cmath.pi * float(x))
def gauss(k): return sum(E(Fr(k * a * a, 8)) for a in range(4)) / 2            # Z4 mit q(a) = k a^2 / 8
def comm(k): return E(Fr(4 * k, 8) - 2 * Fr(k, 8))                              # b(1,1) = q(2) - 2 q(1)
c.ok("D5 hat q(1) = 5/8, A3 hat q(1) = 3/8: Gaußsummen exp(2 pi i 5/8) und exp(2 pi i 3/8), also Stufe 10 und 6 von sechzehn", abs(gauss(5) - E(Fr(5, 8))) < 1e-12 and abs(gauss(3) - E(Fr(3, 8))) < 1e-12)
c.ok("Produkt der beiden Phasen ist 1 (E8 ist phasenfrei), ihr Verhältnis ist i", abs(gauss(5) * gauss(3) - 1) < 1e-12 and abs(gauss(5) / gauss(3) - 1j) < 1e-12)
c.ok("der Kommutator von Schieber und Uhr ist auf der D5 Seite i, auf der A3 Seite -i", abs(comm(5) - 1j) < 1e-12 and abs(comm(3) + 1j) < 1e-12)
ok = all(abs(gauss(k) / gauss(-k % 8) - comm(k)) < 1e-12 for k in (1, 3, 5, 7))
c.ok("SATZ (alle vier möglichen Z4 Formen): sind zwei Seiten verklebbar (k und -k), dann ist das Verhältnis ihrer Phasen gleich dem Kommutator der ersten Seite. Die beiden i der Hypothese sind also zwingend dieselbe Zahl", ok)
c.ok("Umkehr der Orientierung (Seiten tauschen) gibt -i: die Orientierung der Naht ist das Vorzeichen dieser Phase", abs(gauss(3) / gauss(5) + 1j) < 1e-12)
# Bindung an das Originalmodell v724
root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
spec = importlib.util.spec_from_file_location("v724", os.path.join(root, "verification", "v724_phys_t3b_modular_flows.py")); v724 = importlib.util.module_from_spec(spec); import sys; sys.modules["v724"] = v724; spec.loader.exec_module(v724)
obj = v724.derive_internal_objects(); P = dict(zip(v724.CORNER_NAMES, obj.projectors)); M = obj.mother_covariances[float("inf")]
c.ok("v724 ist rein reell: Projektoren und Mutterkovarianz sind reelle symmetrische Matrizen", all(np.isrealobj(p) and np.allclose(p, p.T) for p in P.values()) and np.isrealobj(M) and np.allclose(M, M.T))
tr = lambda *ms: float(np.trace(np.linalg.multi_dot(ms)))
triples = list(itertools.combinations(v724.CORNER_NAMES, 3)); sgn = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}
anti = {t: sum(s * tr(M, P[t[p[0]]], P[t[p[1]]], P[t[p[2]]]) for p, s in sgn.items()) for t in triples}
c.ok("der umlaufabhängige Teil (vollständig antisymmetrisch in den drei Ecken) ist bei allen vier Tripeln exakt null. Das folgt schon aus der Reellheit: Umkehr der Reihenfolge ist Transponieren", all(abs(v) < 1e-12 for v in anti.values()), {"".join(t): f"{v:.1e}" for t, v in anti.items()})
mixed = {t: [tr(M, P[t[i]], P[t[j]], P[t[k]]) - tr(M, P[t[i]], P[t[k]], P[t[j]]) for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1))] for t in triples}
mags = {"".join(t): [round(x, 6) for x in v] for t, v in mixed.items()}
c.ok("ein Ordnungssignal gibt es trotzdem: Tr(M Pa Pb Pc) - Tr(M Pa Pc Pb) ist bei allen vier Tripeln von null verschieden (deckt sich mit der parallelen Aussage)", all(max(abs(x) for x in v) > 1e-6 for v in mixed.values()), mags)
big = sorted(max(abs(x) for x in v) for v in mixed.values())
c.ok("ABBRUCHREGEL GREIFT: die vier Signale haben verschiedene Beträge, und ihr Vorzeichen folgt keiner Umlaufrichtung (die vier Ecken liegen auf einer Zeitgeraden t = -2, -1, 0, 1, nicht auf einem Kreis)", big[-1] / big[0] > 1.05, [round(x, 6) for x in big])
c.ok("FOLGE: das Ordnungssignal von v724 ist nicht die Viertelphase. In einem reellen Modell kann i nicht als Spur dreier Projektoren erscheinen", True)
c.ok("GRENZE: die parallele Sitzung hat ihr Signal nicht im Repository hinterlegt. Geprüft ist die hier definierte Größe, nicht notwendig dieselbe", True)
c.done({"antisymmetrisch": {"".join(t): v for t, v in anti.items()}, "gemischt": mags, "urteil": "Kern exakt (die zwei i sind zwingend dieselbe Zahl), Bindung an das Ordnungssignal von v724 gefallen."})
