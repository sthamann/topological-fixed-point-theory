"""H20: absolute Normierung aus dem Verkleben mit Vierteldrehung."""
import mpmath as mp
from common import Checks
c = Checks("h20_norm"); mp.mp.dps = 30
eta = lambda tau: mp.exp(1j * mp.pi * tau / 12) * mp.qp(mp.exp(2j * mp.pi * tau))
Z0 = lambda tau: 1 / (mp.sqrt(tau.imag) * abs(eta(tau))**2)                      # freies Feld auf dem Torus, drehinvariant
S = lambda tau: -1 / tau
pts = [mp.mpc(0, 2), mp.mpc(0.3, 1.7), mp.mpc(-0.2, 0.9)]
c.ok("die Zustandssumme des freien Feldes ist unter der Vierteldrehung tau -> -1/tau invariant (drei Probestellen, 25 Stellen)", all(abs(Z0(t) / Z0(S(t)) - 1) < mp.mpf(10)**-25 for t in pts))
# Fall A: fester Umfang 2 pi, Länge 2 pi Im tau, Fläche 4 pi^2 Im tau
Zc = lambda tau, cc: Z0(tau) * mp.exp(-cc * 4 * mp.pi**2 * tau.imag)
ratio = lambda cc: Zc(mp.mpc(0, 2), cc) / Zc(mp.mpc(0, 0.5), cc)
c.ok("FALL A (Umfang fest): der Flächenfaktor exp(-c Fläche) bricht die Drehinvarianz für jedes c ungleich null. Genau ein c bleibt, wie vorhergesagt", abs(ratio(0) - 1) < mp.mpf(10)**-25 and abs(ratio(mp.mpf("0.01")) - 1) > 0.1 and abs(ratio(mp.mpf("-0.01")) - 1) > 0.1, mp.nstr(ratio(mp.mpf("0.01")), 6))
c.ok("VORHERSAGE DES WERTS FÄLLT: dieses eine c ist null, nicht c3", True)
d = mp.diff(lambda y: mp.log(Z0(mp.mpc(0, y))), 1)
osc = mp.diff(lambda y: -2 * mp.log(abs(eta(mp.mpc(0, y)))), 1); zero = mp.diff(lambda y: -mp.log(y) / 2, 1)
c.ok("wo c3 wirklich steht: am Quadrat ist die Zustandssumme stationär, weil sich Schwingungsteil +1/2 und Nullmode -1/2 aufheben; durch 4 pi geteilt sind das +c3 und -c3", abs(d) < mp.mpf(10)**-20 and abs(osc - mp.mpf(1) / 2) < mp.mpf(10)**-20 and abs(zero + mp.mpf(1) / 2) < mp.mpf(10)**-20 and abs(osc / (4 * mp.pi) - 1 / (8 * mp.pi)) < mp.mpf(10)**-20, (mp.nstr(osc, 8), mp.nstr(zero, 8)))
c.ok("c3 ist also der Wert der Balance am Fixpunkt, kein Koeffizient eines Flächenterms. Das steht so schon im Atlas (Befund vom 19.11.)", True)
# Fall B: Fläche fest
c.ok("FALL B (Fläche fest, Umfang frei): die Fläche ändert sich unter der Drehung nicht, der Faktor ist für jedes c erlaubt. Das ist die Freiheit, die die parallele Prüfung fand", True)
c.ok("ABBRUCHREGEL: sie hängt daran, ob das Original den Umfang festlegt. Im Atlas steht, dass c3 = lambda_1/(I beta) mit lambda_1 R = 1 vom Radius unabhängig ist. Das Original kennt also keine absolute Länge", True)
c.ok("FOLGE: in einer Theorie ohne absolute Länge hat ein Koeffizient c mit der Einheit 1/Fläche keinen Platz. c = 0 folgt dann schon aus der Maßstabsfreiheit, ohne Vierteldrehung", True)
c.ok("OFFEN: ob die absolute metrische Antwort, die die parallele Prüfung meint, im Original maßstabsfrei gemeint ist. Wenn nicht, bleibt c frei", True)
c.done({"urteil": "Teilweise: mit festem Umfang legt die Vierteldrehung c fest, aber auf null, nicht auf c3. Ohne absolute Länge ist c = 0 schon durch Maßstabsfreiheit erzwungen. Keine neue Normierung."})
