"""KATALOGQUER20261120: numerische Beobachtung am Reaktorwinkel, kein Beweis.

Frage: v328 hält fest, dass der aus sin^2 theta13 rückgerechnete Seed
e^{5/6} * sin^2 theta13_exp = 0.0505 etwa 5 % unter phi0 liegt. Ist dieser Wert ein
schon vorhandenes TFPT Objekt?

Befund: 0.0505 = phi0 (1 - phi0) = lambda_Y^2 = lambda_C^2 innerhalb 0.12 sigma.
Das ist eine nachträgliche Beobachtung an einem Datenpunkt. Sie leitet nichts her
und ändert keine eingefrorene Vorhersage (REG.FREEZE.01 bleibt).

Messwerte wie im Register (v328, Scorecard M075, M093): NuFIT 6.0 NO.
"""
from math import exp, pi, sqrt

PHI0 = 1 / (6 * pi) + 3 / (256 * pi**4)
GAMMA = 5 / 6  # tr_E Y^2, v268
S13_OBS, S13_ERR = 0.02195, 0.00058  # NuFIT 6.0 NO, wie v328
LAMC_OBS, LAMC_ERR = 0.2243, 0.0008  # Kaonmittel, wie M093

checks = []


def check(name, ok, detail):
    checks.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def pull(pred):
    return (pred - S13_OBS) / S13_ERR


lam2 = PHI0 * (1 - PHI0)
frozen = PHI0 * exp(-GAMMA)
ladder = lam2 * exp(-GAMMA)
seed_back = S13_OBS * exp(GAMMA)
seed_err = S13_ERR * exp(GAMMA)

check("eingefrorene Vorhersage reproduziert", abs(frozen - 0.023108) < 1e-6,
      f"phi0 e^-5/6 = {frozen:.6f}, Zug {pull(frozen):+.2f} sigma")
check("rückgerechneter Seed wie v328", abs(seed_back - 0.0505) < 1e-4,
      f"{seed_back:.6f} +- {seed_err:.6f}")
check("rückgerechneter Seed = lambda_C^2", abs(seed_back - lam2) < 0.5 * seed_err,
      f"lambda_C^2 = {lam2:.6f}, Abstand {(seed_back - lam2) / seed_err:+.2f} sigma; "
      f"gegen phi0 {(seed_back - PHI0) / seed_err:+.2f} sigma")
check("Leiterlesart trifft den Messwert", abs(pull(ladder)) < 0.5,
      f"lambda_C^2 e^-5/6 = {ladder:.6f}, Zug {pull(ladder):+.2f} sigma")

# Seedfreie Form: Verhältnis zweier Messgrößen gegen eine reine Zahl.
ratio_obs = sqrt(S13_OBS) / LAMC_OBS
ratio_err = ratio_obs * sqrt((0.5 * S13_ERR / S13_OBS) ** 2 + (LAMC_ERR / LAMC_OBS) ** 2)
ratio_pred = exp(-GAMMA / 2)
check("seedfrei: sin theta13 / lambda_C = e^-5/12", abs(ratio_obs - ratio_pred) < ratio_err,
      f"Messung {ratio_obs:.4f} +- {ratio_err:.4f}, Zahl {ratio_pred:.4f}; "
      f"eingefrorene Lesart verlangt {ratio_pred / sqrt(1 - PHI0):.4f}")

# Ehrliche Reichweite: die Daten bestimmen nur den Koeffizienten erster Ordnung.
k_fit = (1 - seed_back / PHI0) / PHI0
k_err = seed_err / PHI0**2
check("Koeffizient k in phi0 (1 - k phi0)", abs(k_fit - 1) < k_err and abs(k_fit) > 1.9 * k_err,
      f"k = {k_fit:.2f} +- {k_err:.2f}; k = 1 verträglich, k = 0 bei {k_fit / k_err:.1f} sigma")

# Kontrolle gegen Beliebigkeit: benachbarte Basen bei festem Exponenten 5/6.
neighbours = {
    "phi0 (eingefroren)": PHI0,
    "1/(6 pi) (Baumwert)": 1 / (6 * pi),
    "phi0 (1 - phi0) = lambda_C^2": lam2,
    "phi0 / (1 + phi0)": PHI0 / (1 + PHI0),
    "phi0 (1 - 1/(4 pi)) = Omega_b": PHI0 * (1 - 1 / (4 * pi)),
    "phi0 (1 - 2 phi0)": PHI0 * (1 - 2 * PHI0),
    "phi0 (1 - phi0/2)": PHI0 * (1 - PHI0 / 2),
    "phi0 (1 - phi0)^2": PHI0 * (1 - PHI0) ** 2,
    "phi0 (1 + phi0)": PHI0 * (1 + PHI0),
}
print("\nNachbarbasen bei Exponent 5/6 (Zug in sigma):")
hits = 0
for name, base in neighbours.items():
    p = pull(base * exp(-GAMMA))
    hits += abs(p) < 0.5
    print(f"  {name:42s} {base * exp(-GAMMA):.6f}  {p:+.2f}")
check("Kontrolle ist nicht trennscharf zwischen Formen gleicher erster Ordnung", hits >= 2,
      f"{hits} von {len(neighbours)} Nachbarn unter 0.5 sigma: phi0(1-phi0) und phi0/(1+phi0) "
      "sind mit heutigen Daten nicht unterscheidbar")

# Wirkung auf die gemeinsame Seedprüfung (M066: chi2 = 4.51 bei 4 Freiheitsgraden).
chi2_new = 4.51 - pull(frozen) ** 2 + pull(ladder) ** 2
check("Seed chi2 ohne den theta13 Riss", chi2_new < 1.0, f"4.51 -> {chi2_new:.2f}")

print(f"\n{sum(checks)} von {len(checks)} Prüfungen bestanden")
raise SystemExit(0 if all(checks) else 1)
