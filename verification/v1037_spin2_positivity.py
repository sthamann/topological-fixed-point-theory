"""v1037 -- positivity note for the entire-form-factor spin-2 rows (GRAV.SPIN2.UNITARITY.01,
GRAV.KMS.HESSIAN.01, QGAMB.IDG.01, GRAV.AMPLITUDE.01, GRAV.LOOP.FINITE.01, GRAV.NONLOCAL.01,
GRAV.SPIN2.EMERGENCE.01): pole freedom is not spectral positivity.

Question.  v304/v370/v380/v386/v389/v400 dress the spin-2 propagator with the entire seam KMS
form factor, P(p^2) = e^{-p^2/M^2}/p^2 = 1/(p^2 a), a(u) = e^u.  v386 item 4 reads this as
"the propagator's spectral density is a single positive-residue pole at p^2 = 0".  Is that a
positive (Kaellen-Lehmann, reflection-positive) spin-2 propagator?

Method.  A positive two-point function of a field with a massless pole has the Kaellen-Lehmann
form C(Q^2) = Z/Q^2 + Int_0^oo drho(s)/(Q^2 + s), rho >= 0 (Euclidean Q^2 > 0).  Everything
below follows from that representation; the controls are explicit positive and ghost measures.

  [E] 1. KL INEQUALITY.  Q^2 C(Q^2) = Z + Int Q^2 drho/(Q^2 + s) >= Z for every Q^2 > 0, and
        C(Q^2) - Z/Q^2 >= 0.  Checked on 50 random positive measures (point masses and a
        continuum), all Q^2 in [1e-4, 1e4].
  [E] 2. A SINGLE POLE IS EXACTLY Z/Q^2.  If the spectral density is only the pole at 0 (rho = 0),
        C = Z/Q^2.  The entire form factor gives C - Z/Q^2 = (e^{-u} - 1)/Q^2 < 0: it is not "a
        single positive-residue pole" in the spectral sense (v386 item 4 is a pole statement).
  [X] 3. ENTIRE FORM FACTOR VIOLATES KL.  Q^2 P = Z e^{-Q^2/M^2} < Z for every Q^2 > 0 (and every
        a(u) > 1 on u > 0 does the same), so e^{-p^2/M^2}/p^2 has no positive spectral representation.
  [E] 4. INVERSE PROPAGATOR SIGN.  KL with Int drho/s finite gives C^{-1} = Q^2/Z - (Int drho/s) Q^4/Z^2
        + O(Q^6): the Q^4 coefficient is <= 0.  The entire factor (and its IR projection a = 1 + u,
        the R + R^2 Stelle sign) gives + Q^4/(M^2 Z): the wrong sign for positivity.
  [E] 5. POSITIVE CONTROL.  C = Z/Q^2 + r/(Q^2 + M^2), r > 0, satisfies 1 and 4 and has
        C^{-1} = Q^2/Z - r Q^4/(M^2 Z^2) + O(Q^6) (symbolic).
  [E] 6. GHOST CONTROL.  The Stelle truncation 1/(Q^2 (1 + Q^2/M^2)) = 1/Q^2 - 1/(Q^2 + M^2) has a
        negative weight and violates 1 as well: ghosts and entire factors fail the same test, so
        removing the ghost pole does not restore positivity.
  [E] 7. TIMELIKE GROWTH.  At p^2 = -s (s > 0) the entire propagator grows like e^{s/M^2}/s, while a KL
        function is bounded there by its spectral integral away from the support; e^{-p^2/M^2}/p^2 has no
        cut at all, so there is no spectral density to be positive.
  [O] 8. SCOPE.  The rows' pole statements (one pole at p^2 = 0, residue 1, no Stelle ghost) and the
        power counting stand; what fails is the reading as a physical positive (KL / reflection-positive)
        spin-2 propagator and the inference "single pole => Im A >= 0 => unitarity" beyond the pole.
        A positive spin-2 sector needs a positive spectral measure (e.g. Z/Q^2 + Int drho/(Q^2+s)).

Verdict.  Read-only on the ledger: this script checks the mathematics behind the note; it changes no
number of v304/v370/v380/v386/v389/v400.  Python only (numpy, sympy).
"""
import csv
import os

import numpy as np
import sympy as sp

from tfpt_constants import check, summary, reset

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = ("GRAV.SPIN2.UNITARITY.01", "GRAV.KMS.HESSIAN.01", "QGAMB.IDG.01", "GRAV.AMPLITUDE.01",
        "GRAV.LOOP.FINITE.01", "GRAV.NONLOCAL.01", "GRAV.SPIN2.EMERGENCE.01")


def run():
    reset()
    rng = np.random.default_rng(1037)
    Q2 = np.geomspace(1e-4, 1e4, 801)
    Z = 1.0

    # 1. KL inequality on random positive measures
    ok1 = True
    for _ in range(50):
        s = rng.uniform(1e-3, 1e3, size=6)
        w = rng.uniform(0, 5, size=6)
        cont = np.trapezoid(np.exp(-np.linspace(0, 40, 400)) / (Q2[:, None] + np.linspace(1e-3, 40, 400)[None, :]),
                            np.linspace(1e-3, 40, 400), axis=1)
        C = Z / Q2 + np.sum(w[:, None] / (Q2[None, :] + s[:, None]), axis=0) + cont
        ok1 &= bool(np.all(Q2 * C >= Z * (1 - 1e-14)) and np.all(C - Z / Q2 >= -1e-14))
    check("KL INEQUALITY [E]: C = Z/Q^2 + Int drho/(Q^2+s), rho >= 0, gives Q^2 C >= Z and C - Z/Q^2 >= 0 "
          "(50 random positive measures, Q^2 in [1e-4, 1e4])", ok1)

    # 2./3. entire form factor
    M2 = 1.0
    P = np.exp(-Q2 / M2) / Q2
    check("SINGLE POLE = Z/Q^2 [E]: with rho = 0 the KL form is exactly Z/Q^2; the entire factor gives "
          "C - Z/Q^2 = (e^{-u} - 1)/Q^2 < 0 everywhere, so it is not 'a single positive-residue pole' spectrally",
          bool(np.all(P - Z / Q2 < 0)))
    a_list = {"e^u": np.exp, "1+u+u^2": lambda u: 1 + u + u * u, "cosh(u)+u": lambda u: np.cosh(u) + u}
    with np.errstate(over="ignore"):            # a(u) = inf for large u means Q^2 P = 1/a = 0 < Z
        viol = {k: bool(np.all(1 / a(Q2 / M2) < Z)) for k, a in a_list.items()}
    check("ENTIRE FORM FACTOR VIOLATES KL [X]: Q^2 P = e^{-Q^2/M^2} < 1 = Z for every Q^2 > 0; likewise for "
          "any a(u) > 1 on u > 0 (%s)" % viol, all(viol.values()))

    # 4./5. inverse propagator expansion (symbolic)
    q2, z, r, m2, k = sp.symbols("Q2 Z r M2 k", positive=True)
    rho_s = sp.symbols("I1", positive=True)          # I1 = Int drho/s
    C_gen = z / q2 + rho_s - k * q2                  # small-Q^2 expansion of the KL integral: Int drho/(s) - Q^2 Int drho/s^2 ...
    inv_gen = sp.series(1 / C_gen, q2, 0, 3).removeO()
    c4_gen = sp.simplify(inv_gen.coeff(q2, 2))
    inv_ent = sp.series(q2 * sp.exp(q2 / m2) / z, q2, 0, 3).removeO()
    c4_ent = sp.simplify(inv_ent.coeff(q2, 2))
    inv_pos = sp.series(1 / (z / q2 + r / (q2 + m2)), q2, 0, 3).removeO()
    check("INVERSE PROPAGATOR SIGN [E]: KL gives C^{-1} = Q^2/Z - I1 Q^4/Z^2 + ... (Q^4 coefficient %s <= 0); the "
          "entire factor gives +Q^4/(M^2 Z) (coefficient %s > 0), the Stelle/R+R^2 sign" % (c4_gen, c4_ent),
          sp.simplify(c4_gen + rho_s / z**2) == 0 and sp.simplify(c4_ent - 1 / (m2 * z)) == 0)
    check("POSITIVE CONTROL [E]: C = Z/Q^2 + r/(Q^2+M^2) has C^{-1} = Q^2/Z - r Q^4/(M^2 Z^2) + O(Q^6) and "
          "satisfies Q^2 C >= Z",
          sp.simplify(inv_pos - (q2 / z - r * q2**2 / (m2 * z**2))) == 0
          and bool(np.all(Q2 * (1 / Q2 + 0.7 / (Q2 + 2.0)) >= 1)))

    # 6. ghost control
    stelle = 1 / (Q2 * (1 + Q2 / M2))
    check("GHOST CONTROL [E]: the Stelle truncation 1/Q^2 - 1/(Q^2+M^2) (negative weight) gives Q^2 C = "
          "M^2/(Q^2+M^2) < Z too: removing the ghost pole does not restore positivity",
          bool(np.all(Q2 * stelle < Z)) and np.allclose(stelle, 1 / Q2 - 1 / (Q2 + M2)))

    # 7. timelike growth
    s_t = np.array([1.0, 5.0, 10.0, 20.0])
    grow = np.abs(np.exp(s_t / M2) / (-s_t))
    kl_t = np.abs(1 / (-s_t) + 0.7 / (-s_t + 50.0))  # positive KL control evaluated below its support (s < 50)
    check("TIMELIKE GROWTH [E]: at p^2 = -s the entire propagator grows like e^{s/M^2}/s (%s), a KL control "
          "stays O(1/s) below its support (%s); the entire factor has no cut, hence no spectral density"
          % (np.round(grow, 2).tolist(), np.round(kl_t, 4).tolist()),
          bool(np.all(np.diff(grow) > 0) and grow[-1] > 1e7 and np.all(kl_t < 1.1)))

    # 8. ledger rows (read-only) and scope
    with open(os.path.join(HERE, "status_ledger.csv"), newline="", encoding="utf-8") as fh:
        led = {row["claim_id"]: row for row in csv.DictReader(fh)}
    present = [rid for rid in ROWS if rid in led]
    amp = led.get("GRAV.AMPLITUDE.01", {}).get("claim", "")
    check("LEDGER ROWS (read-only) [E]: all seven rows exist and rest on the entire form factor; "
          "GRAV.AMPLITUDE.01 names P(p^2)=e^{-p^2/M^2}/p^2",
          len(present) == 7 and "e^{-p^2/M^2}/p^2" in amp
          and all(("entire" in led[rid]["claim"].lower() or "untruncated" in led[rid]["claim"].lower()
                   or "infinite-derivative" in led[rid]["claim"].lower()) for rid in present))
    check("SCOPE [O]: the pole statements (one pole at p^2=0, residue 1, no Stelle ghost) and the power "
          "counting stand; the reading as a positive KL / reflection-positive spin-2 propagator and the "
          "inference 'single pole => Im A >= 0 => unitarity' beyond the pole do not; a positive spin-2 sector "
          "needs a positive spectral measure", True)

    return summary("v1037 spin-2 positivity: entire form factor is pole-free but not Kaellen-Lehmann positive")


if __name__ == "__main__":
    raise SystemExit(1 if run() else 0)
