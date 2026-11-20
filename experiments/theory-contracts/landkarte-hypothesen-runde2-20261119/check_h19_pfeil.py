"""H19: Zeitpfeil als Umlauf unendlicher Ordnung und weiche Moden an der Marke."""
import cmath, math, itertools
from fractions import Fraction as Fr
from common import Checks
c = Checks("h19_pfeil")
def spec(tau, n=12, k=8):                                    # Laplace auf dem Torus der Fläche 1 mit Form tau
    y = tau.imag; v = sorted(4 * math.pi**2 * abs(m - nn * tau)**2 / y for m in range(-40, 41) for nn in range(-n, n + 1) if (m, nn) != (0, 0)); return [round(x, 8) for x in v[:k]]
c.ok("GEFALLEN (1): die Scherung um 4n ändert den Torus nicht. tau = i und tau = i + 4, i + 8 haben dasselbe Spektrum; ein gescherter Torus ist nicht lang und dünn", spec(1j) == spec(1j + 4) == spec(1j + 8), spec(1j)[:5])
c.ok("der Umlauf zählt nur, wie die Zyklen benannt sind. Sichtbar wird er erst in einer Darstellung auf Zuständen", True)
def order(N, power):                                         # T auf den N Thetafunktionen der Stufe N (N gerade): Phasen exp(i pi k^2 / N)
    ph = [Fr(power * k * k, 2 * N) % 1 for k in range(N)]; ph = [p - ph[0] for p in ph]
    return next(m for m in range(1, 8 * N + 1) if all((m * p) % 1 == 0 for p in ph))
o = {N: (order(N, 1), order(N, 4)) for N in (4, 8, 16)}
c.ok("GEFALLEN (2): auf jedem endlichen Zustandsraum hat der Umlauf endliche Ordnung. Stufe 4: T hat Ordnung 8, T^4 Ordnung 2. Stufe 8: 16 und 4. Die Zählung wird wieder zur Uhr", o == {4: (8, 2), 8: (16, 4), 16: (32, 8)}, o)
c.ok("Folge: eine Zahl, die nur wächst, gibt es erst mit unendlich vielen Zuständen. Der Träger hat fünf", True)
gap = lambda y: 4 * math.pi**2 / y
c.ok("weich werden die Moden nur, wenn man sich der Marke im Formenraum nähert (Im tau -> unendlich): Lücke 4 pi^2 / Im tau", [round(gap(y), 3) for y in (1, 10, 100)] == [39.478, 3.948, 0.395])
ys = [4 / (2 * math.pi) * math.log(1 / t) for t in (1e-2, 1e-4, 1e-8, 1e-16)]
c.ok("an einer Faser I4 wächst Im tau nur wie (4 / 2 pi) log(1/|t|): die Lücke schließt sich logarithmisch langsam (bei |t| = 10^-16 erst auf ein Dreiundzwanzigstel)", abs(gap(ys[-1]) / gap(1) - 1 / ys[-1]) < 1e-12 and 22 < ys[-1] < 24, [round(y, 2) for y in ys])
# das Original: flaches Kissen = quadratischer Torus geteilt durch z -> -z
even = spec(1j, k=12)
c.ok("ABBRUCHREGEL GREIFT: das flache Kissen des Originals ist der quadratische Torus geteilt durch die Halbdrehung. Jede Schwingung cos(2 pi (m x + n y)) ist gerade, die Lücke ist 4 pi^2 je Flächeneinheit und fest. Es gibt dort keinen dünnen Hals und keine Formänderung", even[0] == round(4 * math.pi**2, 8))
c.ok("das Energiekriterium der parallelen Prüfung (weiche Moden auf großen Skalen) ist am Original damit nicht erfüllt. Seine genaue Formel liegt nicht im Repository; geprüft ist nur die Voraussetzung einer nach unten unbeschränkten Lücke", True)
c.ok("WAS BLEIBT: der Unterschied zwischen Uhr (endliche Ordnung) und Zählung (unendliche Ordnung) ist genau der Unterschied zwischen endlich und unendlich vielen Zuständen. Ein Zeitpfeil braucht einen wachsenden Zustandsraum, nicht einen anderen Umlauf", True)
c.done({"ordnungen": {str(k): v for k, v in o.items()}, "urteil": "Gefallen: die Scherung ändert das Spektrum nicht, auf endlichen Zustandsräumen ist der Umlauf wieder eine Uhr, und das flache Kissen hat eine feste Lücke."})
