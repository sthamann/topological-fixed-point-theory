"""Vorrechnungen zu den zehn Hypothesen der zweiten Runde (19.11.2026). Kleine exakte Proben, keine Prüfung der Hypothesen.
Zahlen aus der parallelen Sitzung (J+, J-, OPH Quotienten) sind dort berichtet und hier nur weiterverrechnet."""
import sympy as sp, mpmath as mp, numpy as np, json
out = {}
# 1. Nachbarknoten im erweiterten E8 Diagramm: Marke 4 und Marke 5 streichen
e = sp.eye(8); h = sp.Rational(1, 2)
al = [sp.Matrix([h, -h, -h, -h, -h, -h, -h, h]), e[:, 0] + e[:, 1], e[:, 1] - e[:, 0], e[:, 2] - e[:, 1], e[:, 3] - e[:, 2], e[:, 4] - e[:, 3], e[:, 5] - e[:, 4], e[:, 6] - e[:, 5]]
marks = [2, 3, 4, 6, 5, 4, 3, 2]; theta = sum((m * a for m, a in zip(marks, al)), sp.zeros(8, 1)); a0 = -theta
keep = [a0, al[7], al[6], al[3], al[2], al[0], al[1]]                       # ohne alpha6 (Marke 4) und alpha5 (Marke 5)
G = sp.Matrix(7, 7, lambda i, j: keep[i].dot(keep[j]))
B = sp.Matrix.hstack(*al); v = sp.Matrix.vstack(*[k.T for k in keep]).nullspace()[0]
c = B.solve(v); c = c * sp.ilcm(*[sp.fraction(x)[1] for x in c]); c = c / sp.igcd(*list(c)); v = B * c
out["nachbarknoten"] = {"theta_norm": int(theta.dot(theta)), "det_rest": int(G.det()), "rang": G.rank(), "senkrechter_gittervektor_norm": int(v.dot(v)), "koeffizienten": [int(x) for x in c]}
# 2. elliptische Fläche mit vollem 2 Torsionsgitter und Punkt der Ordnung 4: y^2 = x (x + t^2)(x + 1)
t, x = sp.symbols("t x"); f = sp.expand(x * (x + t**2) * (x + 1)); disc = sp.factor(sp.discriminant(f, x)); c4 = sp.expand(16 * ((t**2 + 1)**2 - 3 * t**2))
out["flaeche_stufe4"] = {"diskriminante": str(disc), "fasern_bei": "0 (I4), 1 und -1 (je I2), unendlich (I4, Grad 12 - 8)", "c4_grad": int(sp.degree(c4, t)), "c4_bei_0": int(c4.subs(t, 0)),
                         "doppelverhaeltnis_0_unendlich_1_-1": -1, "punkt_ordnung_4": bool(sp.simplify(f.subs(x, t) - (t * (t + 1))**2) == 0), "kusp_breiten_summe": 4 + 4 + 2 + 2}
# 3. Geschlecht der Doppelüberlagerung über n Marken
out["geschlecht"] = {n: (n - 2) // 2 for n in (2, 4, 6, 8, 12)}
# 4. Anteile der drei Sektoren der E8 Gittersumme am quadratischen Torus
mp.mp.dps = 30; q = mp.exp(-mp.pi); th = [mp.jtheta(k, 0, q) ** 8 for k in (2, 3, 4)]; S = sum(th)
out["theta_anteile_tau_i"] = [mp.nstr(v / S, 15) for v in th]
# 5. Phase zwischen den Kitaev Stufen von D5 und A3
out["phasen"] = {"exp(2 pi i (10-6)/16)": str(sp.exp(2 * sp.pi * sp.I * sp.Rational(10 - 6, 16))), "exp(2 pi i (10+6)/16)": str(sp.exp(2 * sp.pi * sp.I * sp.Rational(16, 16)))}
# 6. quadratische Form der sechs Paare auf vier Alternativen
M = sp.ones(4, 4) - sp.eye(4); ev = M.eigenvals()
out["paarform"] = {"eigenwerte": {str(k): int(m) for k, m in ev.items()}, "signatur": "(1,3)", "det": int(M.det()), "auf_ebene_summe_null": "negativ definit (-1 dreifach)"}
# 7. berichtete Jarlskog Werte beider Orientierungen (parallele Sitzung)
Jp, Jm = 3.3270e-5, -2.7896e-5
out["jarlskog"] = {"ungerader_teil": (Jp - Jm) / 2, "gerader_teil": (Jp + Jm) / 2, "verhaeltnis_gerade_zu_ungerade": (Jp + Jm) / (Jp - Jm)}
# 8. berichtete OPH Quotienten
out["oph_quotient"] = {N: 8 - 96 / N for N in (240, 960, 3840)}
out["oph_quotient"]["berichtet"] = [7.6, 7.9, 7.975]; out["oph_quotient"]["tfpt"] = 1.125
# 9. E2 am Quadrat
E2 = 1 - 24 * mp.nsum(lambda n: n * mp.exp(-2 * mp.pi * n) / (1 - mp.exp(-2 * mp.pi * n)), [1, mp.inf])
out["E2_stern_bei_i"] = mp.nstr(E2 - 3 / mp.pi, 5)
print(json.dumps(out, indent=1, ensure_ascii=False, default=str)); json.dump(out, open("vorrechnung_runde2.json", "w"), indent=1, ensure_ascii=False, default=str)
