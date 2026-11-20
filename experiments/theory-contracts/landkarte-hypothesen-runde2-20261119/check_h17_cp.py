"""H17: messbare CP Verletzung als ungerader Teil unter Umkehr der Orientierung. Gerechnet mit der Originalformel aus v88 / v561."""
import os, sys, mpmath as mp
from common import Checks
root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")); sys.path.insert(0, os.path.join(root, "verification"))
from tfpt_constants import phi0
c = Checks("h17_cp"); mp.mp.dps = 30; deg = 180 / mp.pi
lam = mp.sqrt(phi0 * (1 - phi0)); s12, s23, s13 = lam, phi0 / (1 + lam), lam**3 / 3
c12, c23, c13 = (mp.sqrt(1 - s**2) for s in (s12, s23, s13)); K = s12 * c12 * s23 * c23 * s13 * c13**2; eps = 3 * lam**2
Jp, Jm = K * mp.sin(mp.pi / 3 + eps), K * mp.sin(-mp.pi / 3 + eps)
c.ok("NACHGERECHNET: J+ = 3,3270e-5 (Phase pi/3 + 3 lambda^2) und J- = -2,7896e-5 (Phase -pi/3 + 3 lambda^2). Die berichteten Zahlen stimmen", abs(Jp - mp.mpf("3.32702e-5")) < 1e-9 and abs(Jm + mp.mpf("2.7896e-5")) < 1e-9, (mp.nstr(Jp, 6), mp.nstr(Jm, 6)))
c.ok("Ursache des ungleichen Betrags: bei der Umkehr kippt nur pi/3, der Zusatz 3 lambda^2 behält sein Vorzeichen. Die Mischungswinkel bleiben gleich", abs(Jm / Jp + mp.sin(mp.pi / 3 - eps) / mp.sin(mp.pi / 3 + eps)) < 1e-25)
odd, even = (Jp - Jm) / 2, (Jp + Jm) / 2
c.ok("ungerader Teil = K sin(pi/3) cos(3 lambda^2) = 3,0583e-5, gerader Teil = K cos(pi/3) sin(3 lambda^2) = 0,2687e-5", abs(odd - K * mp.sin(mp.pi / 3) * mp.cos(eps)) < 1e-25 and abs(even - K * mp.cos(mp.pi / 3) * mp.sin(eps)) < 1e-25, (mp.nstr(odd, 6), mp.nstr(even, 6)))
JP, JE, GP, GE = mp.mpf("3.08e-5"), mp.mpf("0.14e-5"), mp.mpf("65.7"), mp.mpf("3.0")          # eingefrorene Vergleichswerte aus v88
pJ_full, pJ_odd = (Jp - JP) / JE, (odd - JP) / JE
c.ok("gegen den im Repository eingefrorenen Messwert J = 3,08(14)e-5: voll +1,76 sigma, ungerader Teil -0,16 sigma", abs(pJ_full - 1.764) < 0.01 and abs(pJ_odd + 0.155) < 0.01, (mp.nstr(pJ_full, 4), mp.nstr(pJ_odd, 4)))
d_full = (mp.pi / 3 + eps) * deg; d_odd = mp.asin(mp.sin(mp.pi / 3) * mp.cos(eps)) * deg
pg_full, pg_odd, pg_60 = (d_full - GP) / GE, (d_odd - GP) / GE, (60 - GP) / GE
c.ok("VORHER FESTGELEGTER ZWEITTEST, die Phase selbst gegen gamma = 65,7(3,0) Grad: voll 68,65 Grad (+0,98 sigma), ungerade Lesart 58,89 Grad (-2,27 sigma) oder 60 Grad (-1,90 sigma)", abs(pg_full - 0.985) < 0.01 and abs(pg_odd + 2.27) < 0.01 and abs(pg_60 + 1.9) < 0.01, (mp.nstr(d_odd, 5), mp.nstr(pg_odd, 4)))
chi_full, chi_odd = pJ_full**2 + pg_full**2, pJ_odd**2 + pg_odd**2
c.ok("ABBRUCHREGEL GREIFT: beide Größen zusammen passen mit der ungeraden Lesart schlechter (5,2) als mit der vollen (4,1). Der Gewinn bei J wird bei der Phase mehr als verloren", chi_odd > chi_full, (mp.nstr(chi_full, 3), mp.nstr(chi_odd, 3)))
Js = Jp * mp.mpf("0.0408") / s23
c.ok("das Original erklärt den Überschuss bei J schon anders: s23 ist dort bei der Quellskala 0,0434 statt 0,0408 bei M_Z. Mit dem Wert bei M_Z wird J = 3,13e-5 (+0,34 sigma), ohne Eingriff in die Phase", abs((Js - JP) / JE - 0.34) < 0.03, mp.nstr(Js, 5))
def vtd(d): return abs(s12 * s23 - c12 * c23 * s13 * mp.e**(1j * d))
c.ok("der gerade Teil lässt sich nicht wegtransformieren: die zwei Orientierungen haben verschiedene Beträge |V_td| (CP gerade, messbar). Sie sind zwei verschiedene Theorien, nicht Spiegelbilder", abs(vtd(mp.pi / 3 + eps) / vtd(-mp.pi / 3 + eps) - 1) > 0.01, (mp.nstr(vtd(mp.pi / 3 + eps), 6), mp.nstr(vtd(-mp.pi / 3 + eps), 6)))
dl = 4 * mp.pi / 3
c.ok("bei den Leptonen ist die Phase 4 pi/3 ohne Zusatz: die Umkehr kippt sie genau, der gerade Teil ist null. Dort sagt die Hypothese nichts Neues", abs(mp.sin(dl) + mp.sin(-dl)) < 1e-25)
c.ok("EINORDNUNG: die Nähe von 3,058e-5 zum Messwert fiel nachträglich auf und hält dem zweiten Test nicht stand", True)
c.done({"J_plus": mp.nstr(Jp, 8), "J_minus": mp.nstr(Jm, 8), "ungerade": mp.nstr(odd, 8), "gerade": mp.nstr(even, 8), "zug_J": [mp.nstr(pJ_full, 4), mp.nstr(pJ_odd, 4)], "zug_gamma": [mp.nstr(pg_full, 4), mp.nstr(pg_odd, 4)],
        "urteil": "Gefallen. Die berichteten Zahlen sind nachgerechnet; die ungerade Lesart trifft J, verfehlt aber die Phase, und der gerade Teil ist messbar."})
