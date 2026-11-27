"""Whole-solution search (H6h50, z171): one native lattice object against all exact targets at once.

Guideline "Rekonstruktion statt Auswahl": the catalogue contains only objects built from data H6 already has: the E8
charge lattice, the clock R = exp(2 pi i v.P) with v = (1/2)^5 (0)^3, the Weyl group S4 of the family factor A3 = D3,
and three functionals. Every object is evaluated to the end (count, hypercharge trace, base, kernel power, then alpha
from the ledger equation and the Lambda rung) and compared with the measured values. A clause is a computed property. A candidate "clicks" only if one object carries several
independent targets with one normalisation; the report states how many catalogue objects reach each target (null
control) and which targets no object reaches.

Catalogue: selector (clock class k = 0..3, clock insertion R^j, or none) x family projection (none or an irreducible
representation of S4, projector P = (dim/24) sum chi(w) U_w on the span of the lattice vectors) x functional
(trace per shell, minus log of the twisted lattice sum, minus log of the Fredholm determinant of P q^(L^2/2)) x
number of copies (1 or 3). Series in q, the weight of a vector of norm 2; with a discrete Gaussian of variance
alpha/2 (analysis 24) q = exp(-2/alpha).

Targets: alpha^-1 = 137.035999177(21) and rho_Lambda / M_P^4 = 1.134e-123 (3 %), three families, b1 fermion part 4.
Ledger input that the search does not derive: the form of the alpha equation, the weight c3^4 exp(-2x), the scalar
part 1/10 of b1, the insertion form Q (1-Q)^(-2h) and the choice of the D5 part for the kernel.

Two further universes (z172). Law forms: transfer operators on the three pairings of the four marks, alpha as
k c3^2 times the largest real eigenvalue, solved self-consistently. Compositions: the seam doublet with partner
doublets, the net of seams as the space of lines with a hermitian form or a quaternionic structure.
"""
from fractions import Fraction as Fr
from itertools import permutations, product
from math import exp, log, pi

from genesis_rules import search

SHELLS = 4
POOL = ["three_families", "single_copy", "clock_native", "matter_sector", "positive", "alpha_measured",
        "lambda_measured", "b1_fermion_4", "base_from_object"]
V = tuple([Fr(1, 2)] * 5 + [Fr(0)] * 3)
AINV = 137.0359992168407
C3 = 1 / (8 * pi)
ALPHA_INV_MEASURED, ALPHA_INV_SIGMA = 137.035999177, 0.000000021            # CODATA 2022
RHO_LAMBDA_MEASURED, RHO_LAMBDA_REL = 1.134e-123, 0.03                      # rho_Lambda / M_P^4 (Planck 2018), about 3 %
HYPERCHARGE = (Fr(-1, 3), Fr(-1, 3), Fr(-1, 3), Fr(1, 2), Fr(1, 2))        # on the five places (TFPT-STAND, section 3)
# norm 2h of the two parts of a root in clock class k: (D5 part, A3 part); class 0 has the whole root in one part
FACTOR_POWER = {0: (Fr(2), Fr(2)), 1: (Fr(5, 4), Fr(3, 4)), 2: (Fr(1), Fr(1))}


def alpha_inverse(n_top, power, b1, base=4.0 / 3.0):
    """Root of the ledger equation a^3 - 2 c3^3 a^2 = 8 b1 c3^6 ln(1/phi_seam(a)) with
    phi_seam = base c3 + Q (1 - Q)^(-power), Q = n_top c3^4 exp(-2a). The form of the equation is ledger input."""
    def f(a):
        q = n_top * C3 ** 4 * exp(-2 * a)
        phi = base * C3 + q * (1 - q) ** (-power)
        return a ** 3 - 2 * C3 ** 3 * a ** 2 - 8 * b1 * C3 ** 6 * log(1 / phi)
    lo, hi = 1e-4, 0.5
    if f(lo) * f(hi) > 0:
        return None
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 1 / lo


def rho_lambda(n_top, ainv):
    """Ledger rung COSMO.LAM.01 in unreduced Planck units: rho_Lambda / M_P^4 = n_top c3^4 exp(-2/alpha)."""
    return n_top * C3 ** 4 * exp(-2 * ainv)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def shells(nmax=SHELLS):
    """E8 vectors by shell n = norm/2, n = 1..nmax."""
    out = {n: [] for n in range(1, nmax + 1)}
    top = int((2 * nmax) ** 0.5)
    for v in product(range(-top, top + 1), repeat=8):
        n2 = sum(x * x for x in v)
        if 0 < n2 <= 2 * nmax and n2 % 2 == 0 and sum(v) % 2 == 0:
            out[n2 // 2].append(tuple(Fr(x) for x in v))
    odd = [x for x in range(-2 * top - 1, 2 * top + 2) if x % 2]
    for v in product(odd, repeat=8):
        n2 = sum(x * x for x in v)
        if n2 <= 8 * nmax and n2 % 8 == 0 and (sum(v) // 2) % 2 == 0:
            out[n2 // 8].append(tuple(Fr(x, 2) for x in v))
    return out


def family_group():
    """Weyl group of D3 on the last three coordinates (permutations and even sign changes), with the cycle type of
    each element on the four spinor weights; this identifies it with S4."""
    spin = [s for s in product((1, -1), repeat=3) if s.count(-1) % 2 == 0]
    els = []
    for perm in permutations(range(3)):
        for sg in product((1, -1), repeat=3):
            if sg.count(-1) % 2:
                continue
            img = [spin.index(tuple(sg[i] * s[perm[i]] for i in range(3))) for s in spin]
            seen, cyc = set(), []
            for i in range(4):
                if i in seen:
                    continue
                k, j = 0, i
                while j not in seen:
                    seen.add(j)
                    j = img[j]
                    k += 1
                cyc.append(k)
            els.append((perm, sg, tuple(sorted(cyc, reverse=True))))
    return els


# characters of S4 by cycle type
CHAR = {"trivial": {(1, 1, 1, 1): 1, (2, 1, 1): 1, (2, 2): 1, (3, 1): 1, (4,): 1},
        "sign": {(1, 1, 1, 1): 1, (2, 1, 1): -1, (2, 2): 1, (3, 1): 1, (4,): -1},
        "two": {(1, 1, 1, 1): 2, (2, 1, 1): 0, (2, 2): 2, (3, 1): -1, (4,): 0},
        "standard": {(1, 1, 1, 1): 3, (2, 1, 1): 1, (2, 2): -1, (3, 1): 0, (4,): -1},
        "standard_sign": {(1, 1, 1, 1): 3, (2, 1, 1): -1, (2, 2): -1, (3, 1): 0, (4,): 1}}


def fixed_counts(sh, group, weight=None):
    """fix[n][k][g] = sum of weight(lam) over the vectors of shell n and clock class k fixed by the group element g
    (weight None: the number of such vectors)."""
    out = {}
    for n, vecs in sh.items():
        tab = [[Fr(0)] * len(group) for _ in range(4)]
        for lam in vecs:
            k = int(4 * dot(V, lam)) % 4
            tail = lam[5:]
            wt = Fr(1) if weight is None else weight(lam)
            for g, (perm, sg, _) in enumerate(group):
                if tuple(sg[i] * tail[perm[i]] for i in range(3)) == tail:
                    tab[k][g] += wt
        out[n] = tab
    return out


def projected(fix, group, rep):
    """tr[n][k]: trace of the family projector on shell n and clock class k (rep None: no projection)."""
    out = {}
    ident = next(g for g, e in enumerate(group) if e[2] == (1, 1, 1, 1))
    for n, tab in fix.items():
        row = []
        for k in range(4):
            if rep is None:
                row.append(Fr(tab[k][ident]))
            else:
                d = CHAR[rep][(1, 1, 1, 1)]
                row.append(Fr(d, 24) * sum(CHAR[rep][e[2]] * tab[k][g] for g, e in enumerate(group)))
        out[n] = row
    return out


def _mul(a, b, N):
    c = [Fr(0)] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j <= N:
                    c[i + j] += x * y
    return c


def _minus_log(a, N):
    """-ln(1 + a_1 q + ...) to order N."""
    x = [Fr(0)] + list(a[1:])
    out, pw = [Fr(0)] * (N + 1), [Fr(1)] + [Fr(0)] * N
    for m in range(1, N + 1):
        pw = _mul(pw, x, N)
        for i in range(N + 1):
            out[i] += Fr((-1) ** m, m) * pw[i]
    return out


def catalogue(nmax=SHELLS):
    """Every object is evaluated to the end: its first shell gives the count N (delta_top = N c3^4), the trace of
    Y^2 gives the fermionic part of b1 (b1 = (2/5) Tr Y^2 + 1/10, the 1/10 of the one scalar doublet is ledger
    input), the ratio unprojected/projected gives the base, its clock class gives the kernel power 2h of the D5 or
    the A3 part of its roots. From these the ledger equation gives alpha and the Lambda rung."""
    sh = shells(nmax)
    group = family_group()
    fix = fixed_counts(sh, group)
    fix_y = fixed_counts({1: sh[1]}, group, lambda lam: dot(HYPERCHARGE, lam[:5]) ** 2)
    # class 3 is the charge conjugate of class 1 with identical data; the catalogue keeps one of the two
    sel = [("class %d" % k, [1 if j == k else 0 for j in range(4)], k) for k in range(3)]
    sel += [("insert R^%d" % j, [(1, 1j, -1, -1j)[(j * k) % 4] for k in range(4)], None) for j in (1, 2)]
    sel += [("none", [1, 1, 1, 1], None)]
    objs = []
    for rep in (None,) + tuple(CHAR):
        tr, tr_none, tr_y = projected(fix, group, rep), projected(fix, group, None), projected(fix_y, group, rep)
        for sname, w, cls in sel:
            def comb(row):
                val = complex(sum(w[k] * row[k] for k in range(4)))
                return None if abs(val.imag) > 1e-12 else Fr(val.real).limit_denominator(10 ** 6)
            a = [Fr(1)] + [comb(tr[n]) for n in range(1, nmax + 1)]
            if None in a:
                continue
            a0, y2 = comb(tr_none[1]), comb(tr_y[1])
            for func in ("trace", "minus_log", "fredholm"):
                if func == "trace":
                    ser = [Fr(0)] + a[1:]
                elif func == "minus_log":
                    ser = _minus_log(a, nmax)
                else:                                                       # -ln det(1 - P q^n) = sum a_n sum q^(nm)/m
                    ser = [Fr(0)] * (nmax + 1)
                    for n in range(1, nmax + 1):
                        for m in range(1, nmax // n + 1):
                            ser[n * m] += a[n] / m
                for copies in (1, 3):
                    s = [copies * x for x in ser]
                    n_top, b1 = s[1], Fr(2, 5) * copies * y2 + Fr(1, 10)
                    base = a0 / a[1] if a[1] != 0 else None
                    powers = FACTOR_POWER[cls] if cls is not None else (None, None)
                    for factor, power in zip(("D5", "A3"), powers):
                        ainv = rho = None
                        if power is not None and n_top > 0 and b1 > 0 and base:
                            ainv = alpha_inverse(float(n_top), float(power), float(b1), float(base))
                            rho = rho_lambda(float(n_top), ainv) if ainv else None
                        pull = (ainv - ALPHA_INV_MEASURED) / ALPHA_INV_SIGMA if ainv else None
                        fam = (rep in ("standard", "standard_sign")) or copies == 3
                        cl = {c for c, ok in (
                            ("three_families", fam), ("single_copy", copies == 1), ("clock_native", sname != "none"),
                            ("matter_sector", cls == 1),
                            ("positive", func == "fredholm" and all(x >= 0 for x in a[1:])),
                            ("alpha_measured", pull is not None and abs(pull) < 3),
                            ("lambda_measured", rho is not None
                             and abs(rho / RHO_LAMBDA_MEASURED - 1) < 3 * RHO_LAMBDA_REL),
                            ("b1_fermion_4", copies * y2 * Fr(2, 5) == 4),
                            ("base_from_object", base == Fr(4, 3))) if ok}
                        objs.append({"label": f"{copies} x {func} | {sname} | family {rep} | kernel {factor}",
                                     "series": [str(x) for x in s], "shell_traces": [str(x) for x in a],
                                     "n_top": str(n_top), "b1": str(b1), "base": str(base), "power": str(power),
                                     "alpha_inverse": ainv, "alpha_pull_sigma": pull, "rho_lambda": rho,
                                     "clauses": cl})
                        if power is None:
                            break
    return objs


def run_whole_solution_search(nmax=SHELLS):
    objs = catalogue(nmax)
    res = search(objs, POOL, lambda o: 0, lambda o: o["label"])
    world = {"alpha_measured", "lambda_measured"}
    native = {"three_families", "single_copy", "clock_native", "matter_sector", "positive", "b1_fermion_4",
              "base_from_object"}
    hits = sorted(o["label"] for o in objs if world <= o["clauses"])
    full = sorted(o["label"] for o in objs if world | native <= o["clauses"])
    evaluated = [o for o in objs if o["alpha_inverse"]]
    nearest = sorted((o for o in evaluated if not world <= o["clauses"]), key=lambda o: abs(o["alpha_pull_sigma"]))
    keys = ("n_top", "b1", "base", "power", "alpha_inverse", "alpha_pull_sigma", "rho_lambda", "series",
            "shell_traces")
    return {"status": "WHOLE_SOLUTION_SEARCHED", "physical_solution": False, "objects": len(objs),
            "evaluated": len(evaluated), "counts": {c: sum(c in o["clauses"] for o in objs) for c in POOL},
            "world_hits": hits, "all_clauses": full,
            "derived": {o["label"]: {k: o[k] for k in keys} for o in objs if world <= o["clauses"]},
            "nearest_miss": {o["label"]: {k: o[k] for k in keys[:7]} for o in nearest[:4]},
            "ledger_inputs_not_derived": ["form of the alpha equation (cubic, -2 c3^3, 8 c3^6)", "weight c3^4 exp(-2x)",
                                          "scalar part 1/10 of b1", "insertion form Q (1-Q)^(-2h)",
                                          "choice of the D5 factor for the kernel"],
            "search": res}


# ---------------------------------------------------------------------------------------------------------------
# Law forms (z172): the alpha equation as the eigenvalue equation of a transfer operator on the three pairings of
# the four marks. Basis: the pairings {01|23}, {02|13}, {03|12}; the second pairs opposite marks.

PAIRINGS = ((1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1))
LAW_POOL = ["alpha_measured", "cyclic", "self_loop_opposite", "matter_trace", "seam_response", "irreducible",
            "unit_e2_over_c3"]


def pairing_action(perm):
    """Matrix of a permutation of the four marks on the three pairing vectors (entries 0, 1, -1)."""
    out = [[0] * 3 for _ in range(3)]
    for j, v in enumerate(PAIRINGS):
        w = [0] * 4
        for i in range(4):
            w[perm[i]] = v[i]
        for i, u in enumerate(PAIRINGS):
            out[i][j] = sum(a * b for a, b in zip(u, w)) // 4
    return out


def _perron(mat):
    import numpy as np
    ev = np.linalg.eigvals(np.array(mat, dtype=float))
    real = [x.real for x in ev if abs(x.imag) < 1e-12]
    return max(real) if real else None


def law_form_candidates():
    """T = d D + transport, alpha = k c3^2 * (largest real eigenvalue of T), solved self-consistently because the
    seam response depends on alpha. The product of the transport weights around its cycle is W."""
    q_top = 48 * C3 ** 4
    b1 = 4.1

    def weights(alpha):
        seam = 4 * C3 / 3 + q_top * exp(-2 * alpha) * (1 - q_top * exp(-2 * alpha)) ** -1.25
        return {"b1 ln(1/phi_seam)": b1 * log(1 / seam), "4 ln(1/phi_seam)": 4 * log(1 / seam),
                "ln(1/phi_seam)": log(1 / seam), "b1 ln(1/phi_base)": b1 * log(3 / (4 * C3)),
                "b1 ln(1/phi_0)": b1 * log(1 / (4 * C3 / 3 + q_top)), "b1 ln(1/phi_seam)/2pi": b1 * log(1 / seam) / (2 * pi),
                "1 (modular)": 1.0}

    diag = {"none": (0, 0, 0), "opposite": (0, 1, 0), "adjacent": (1, 0, 1), "all": (1, 1, 1)}
    strength = {"c3/2": C3 / 2, "c3": C3, "2 c3": 2 * C3, "4 c3/3": 4 * C3 / 3}
    transport = {"3-cycle": ((1, 0), (2, 1), (0, 2)), "swap adjacent": ((2, 0), (0, 2)),
                 "swap with opposite": ((1, 0), (0, 1)), "identity": ((0, 0), (1, 1), (2, 2))}
    objs = []
    for dn, dv in diag.items():
        for sn, sv in (strength.items() if dn != "none" else [("-", 0.0)]):
            for tn, edges in transport.items():
                for wn in weights(1 / 137.0):
                    for k in (1, 2, 4):
                        alpha = 1 / 137.0
                        for _ in range(60):
                            w = weights(alpha)[wn]
                            lam = w if tn == "identity" else w ** (1.0 / len(edges))
                            mat = [[sv * dv[i] if i == j else 0.0 for j in range(3)] for i in range(3)]
                            for i, j in edges:
                                mat[i][j] += lam
                            t = _perron(mat)
                            alpha = k * C3 ** 2 * t
                        pull = (1 / alpha - ALPHA_INV_MEASURED) / ALPHA_INV_SIGMA
                        cl = {c for c, ok in (
                            ("alpha_measured", abs(pull) < 3), ("cyclic", tn == "3-cycle"),
                            ("self_loop_opposite", dn == "opposite"), ("matter_trace", wn.startswith("b1")),
                            ("seam_response", "phi_seam" in wn), ("irreducible", tn == "3-cycle"),
                            ("unit_e2_over_c3", k == 2)) if ok}
                        objs.append({"label": f"{sn} on {dn} | {tn} | {wn} | alpha = {k} c3^2 t",
                                     "alpha_inverse": 1 / alpha, "alpha_pull_sigma": pull, "clauses": cl})
    return objs


def run_law_form_search():
    objs = law_form_candidates()
    res = search(objs, LAW_POOL, lambda o: 0, lambda o: o["label"])
    hits = sorted(o["label"] for o in objs if "alpha_measured" in o["clauses"])
    miss = sorted((o for o in objs if "alpha_measured" not in o["clauses"]), key=lambda o: abs(o["alpha_pull_sigma"]))
    return {"status": "LAW_FORMS_SEARCHED", "physical_solution": False, "objects": len(objs), "hits": hits,
            "nearest_miss": {o["label"]: o["alpha_pull_sigma"] for o in miss[:4]},
            "not_derived": ["transport weight b1 ln(1/phi_seam) around the cycle", "why alpha is the eigenvalue",
                            "strength c3 of the self-loop"], "search": res}


# ---------------------------------------------------------------------------------------------------------------
# Compositions (z172): a seam glued with partners. The seam doublet with k - 1 partner doublets spans C^(2k); the
# net of seams is the space of lines of P(C^(2k)), with a hermitian form (p, q) or a quaternionic structure.

COMP_POOL = ["net_dim_4", "real_slice", "lorentzian", "euclidean", "one_partner", "balanced_form"]


def _klein():
    """Basis e_ab of the two-forms on C^4, the Klein form and the duality star."""
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    star = {(0, 1): ((2, 3), 1), (2, 3): ((0, 1), 1), (0, 2): ((1, 3), -1), (1, 3): ((0, 2), -1),
            (0, 3): ((1, 2), 1), (1, 2): ((0, 3), 1)}
    return pairs, star


def _fixed_signature(*linear_parts):
    """Signature of the Klein form on the common fixed space of the antilinear maps M_i o conj on the two-forms."""
    import numpy as np
    pairs, star = _klein()
    qm = np.zeros((6, 6))
    for a, (b, s) in star.items():
        qm[pairs.index(a), pairs.index(b)] = s
    rows = []
    for m in linear_parts:
        m = np.array(m, dtype=float)
        rows.append(np.block([[m, np.zeros((6, 6))], [np.zeros((6, 6)), -m]]) - np.eye(12))
    _, sv, vt = np.linalg.svd(np.vstack(rows))
    null = vt[sum(s > 1e-9 for s in sv):].T
    if null.shape[1] == 0:
        return (0, 0)
    gram = null.T @ np.block([[qm, np.zeros((6, 6))], [np.zeros((6, 6)), -qm]]) @ null
    mixed = null.T @ np.block([[np.zeros((6, 6)), qm], [qm, np.zeros((6, 6))]]) @ null
    assert abs(mixed).max() < 1e-9
    ev = np.linalg.eigvalsh(gram)
    return (int(sum(ev > 1e-9)), int(sum(ev < -1e-9)))


def hermitian_part(signs):
    """Linear part of omega -> star(Lambda^2 H conj(omega)) for H = diag(signs); squares to det H."""
    pairs, star = _klein()
    m = [[0] * 6 for _ in range(6)]
    for a in pairs:
        b, s = star[a]
        m[pairs.index(b)][pairs.index(a)] = s * signs[a[0]] * signs[a[1]]
    return m


def quaternionic_part():
    """Linear part of the map induced by j(e0) = e1, j(e1) = -e0, j(e2) = e3, j(e3) = -e2."""
    pairs, _ = _klein()
    img = {0: (1, 1), 1: (0, -1), 2: (3, 1), 3: (2, -1)}
    m = [[0] * 6 for _ in range(6)]
    for a in pairs:
        (x, sx), (y, sy) = img[a[0]], img[a[1]]
        s = sx * sy
        if x > y:
            x, y, s = y, x, -s
        m[pairs.index((x, y))][pairs.index(a)] = s
    return m


def composition_candidates(kmax=4):
    objs = []
    for k in range(1, kmax + 1):
        n = 2 * k
        dim = 2 * (n - 2)                                                   # complex dimension of the space of lines
        forms = [("hermitian", (n - q, q)) for q in range(0, k + 1)] + [("quaternionic", None)]
        for kind, sig in forms:
            if kind == "hermitian":
                real_dim = 4 * (n - 3) if min(sig) >= 2 else None           # totally null two-planes
            else:
                real_dim = 4 * (k - 1)                                      # quaternionic lines, HP^(k-1)
            real_slice = real_dim is not None and real_dim == dim and dim > 0
            signature = None
            if n == 4 and real_slice:
                part = hermitian_part((1, 1, -1, -1)) if kind == "hermitian" else quaternionic_part()
                signature = _fixed_signature(part)
            cl = {c for c, ok in (
                ("net_dim_4", dim == 4), ("real_slice", real_slice),
                ("lorentzian", signature is not None and sorted(signature) == [2, 4]),
                ("euclidean", signature is not None and sorted(signature) == [1, 5]),
                ("one_partner", k == 2), ("balanced_form", kind == "hermitian" and sig[0] == sig[1])) if ok}
            objs.append({"label": f"{k - 1} partner | {kind} {sig if sig else ''}".strip(), "net_dimension": dim,
                         "real_dimension": real_dim, "klein_signature": signature, "clauses": cl})
    return objs


def run_composition_search():
    objs = composition_candidates()
    res = search(objs, COMP_POOL, lambda o: 0, lambda o: o["label"])
    return {"status": "COMPOSITIONS_SEARCHED", "physical_solution": False, "objects": len(objs),
            "lorentzian": [o["label"] for o in objs if "lorentzian" in o["clauses"]],
            "euclidean": [o["label"] for o in objs if "euclidean" in o["clauses"]],
            "real_slices": [o["label"] for o in objs if "real_slice" in o["clauses"]],
            "signature_without_real_structure": {"(4,0)": _fixed_signature(hermitian_part((1, 1, 1, 1))),
                                                 "(3,1)": _fixed_signature(hermitian_part((1, 1, 1, -1)))},
            "common_slice": _fixed_signature(hermitian_part((1, 1, -1, -1)), quaternionic_part()),
            "not_derived": ["Einstein coefficient", "scale (infinity twistor strength)", "premise T' / I_pos"],
            "search": res}


# ---------------------------------------------------------------- D: seed worlds, readings and assemblies (2026-11-27)
# A complete solution in the sense of docs/research/TFPT-STAND.md, section 1, is one basic object whose own
# principles fix the seed of the theory, carry the recorded results and deliver the open parts. This section makes
# that searchable from the basics instead of from the premises of one candidate:
#   - a WORLD is a seed (k sheets, n marks, N carrier modes, f families, d spacetime dimensions); every clause is
#     computed from the seed, the numbers by the recorded law forms with the seed inserted;
#   - a PRINCIPLE is a clause; the rule search runs over all subsets of the pool and keeps the smallest worlds;
#   - a READING is a basic object with the principles it brings and the roles it delivers (recorded status);
#   - an ASSEMBLY is a set of readings; its roles are computed from the consequences of its rule where a role is a
#     seed property, and taken from the readings otherwise. It is complete when every role is delivered and the
#     readings are identified with each other by bridges that follow from the objects.
# Conventions that are choices of this generator and no results: c3(k) = 1/(4 pi k); base seed (1 - 1/a - 1/b)/pi
# for a split a + b of the carrier; b1 = (2/3) 2^(N-4) f + a/(6N) (f generations of the even spinor of SO(2N), one
# scalar b-plet); kernel power N/4; sin^2(theta12) = 1/f - phi0/2; Omega_b = (1 - 2 c3) phi0. For (2, 5, 3) they
# reproduce the ledger values. Measured values are quoted from memory with rounded errors, apart from alpha.
WORLD_BOUNDS = {"k": (1, 6), "n": (2, 8), "N": (2, 9), "f": (1, 7), "d": (2, 8)}
TFPT_WORLD = {"k": 2, "n": 4, "N": 5, "f": 3, "d": 4}
COORDS = ("k", "n", "N", "f", "d")
WORLD_TARGETS = {"alpha_measured": ("alpha_inverse", ALPHA_INV_MEASURED, ALPHA_INV_SIGMA),
                 "cabibbo": ("cabibbo", 0.2243, 0.0008), "solar": ("solar", 0.307, 0.012),
                 "baryon": ("baryon", 0.0490, 0.0009)}
WORLD_SIGMAS = 3
PRINCIPLES = {
    "flat_cover": "the k-sheeted cover of the sphere, totally branched at n marks, is flat (genus 1): n(k-1) = 2k",
    "region_pairs": "the marks are the ends of regions: n is even",
    "spin_sheets": "the sheets are the fundamental group of the rotation group of space: k = |pi_1(SO(d-1))|",
    "holomorphic_sky": "the sphere of directions S^(d-2) of an observer is a Riemann surface: d = 4",
    "metric_marks": "the marks are the null directions of a spin-2 field, a binary form of degree 2s: n = 4",
    "critical_index": "the inclusion has the smallest index that admits infinite depth: k^2 = 4",
    "index_marks": "the two-region index k^2 of the Z_k orbifold of a holomorphic net is the number of marks",
    "lift_marks": "the order 2k of the spin lift of the deck rotation is the number of marks",
    "divisor_carrier": "the carrier is H^0(O(D)) of the mark divisor: N = n + 1",
    "divisor_families": "the families are H^0(K(D)) of the mark divisor: f = n - 1",
    "z4_anomaly": "the even Fock space of N modes carries a consistent Z4 with square (-1)^F: N odd, 16 | 2^(N-1)",
    "unimodular_glue": "D_N + A_f glue to an even unimodular lattice (anti-isometric discriminant forms)",
    "family_rotation": "the two-forms of the families are the families and the rotations of space: f = d-1, C(f,2) = f",
    "conformal_family_factor": "the non-compact real form of A_f beside a compact D_N is so(d,2): (f,d) = (3,4) or (1,1)",
    "positive_split": "exactly one split a + b of the carrier has a positive base seed 1 - 1/a - 1/b",
    "alpha_measured": "alpha^-1 from the recorded equation with the seed inserted lies within 3 sigma of CODATA",
    "cabibbo": "sqrt(phi0 (1 - phi0)) lies within 3 sigma of |V_us|",
    "solar": "1/f - phi0/2 lies within 3 sigma of sin^2(theta12)",
    "baryon": "(1 - 2 c3) phi0 lies within 3 sigma of Omega_b"}
# principles that the recorded world does not satisfy; without them no rule could select another world
RIVALS = {
    "three_sheets": "the deck group is Z3 (Eisenstein torus, three marks): k = 3",
    "type_a_glue": "A_(N-1) + A_f glue to an even unimodular lattice (carrier group SU(N) instead of SO(2N))",
    "anyonic_space": "space is a plane, the rotation group has fundamental group Z: d = 3",
    "cross_product_seven": "the families carry the cross product of R^7: f = 7, d = 8",
    "six_sphere_sky": "the sphere of directions is S^6 (almost complex): d = 8",
    "balanced_carrier": "the carrier splits into two equal halves: N even",
    "odd_marks": "the number of marks is odd"}
WORLD_POOL = list(PRINCIPLES) + list(RIVALS)
CLOSURE = ["flat_cover", "divisor_carrier", "divisor_families", "z4_anomaly", "positive_split"]   # plus a glue
CHAIN = ["divisor_carrier", "divisor_families"]            # exact consequences of the mark divisor (P2 follows from P1)


def positive_splits(N):
    return [(a, N - a) for a in range((N + 1) // 2, N) if 1 - Fr(1, a) - Fr(1, N - a) > 0]


def world_numbers(k, N, f, split):
    """Numbers of the world with seed (k, N, f) and carrier split a + b, by the recorded law forms."""
    a, b = split
    c = 1 / (4 * pi * k)
    base = float(1 - Fr(1, a) - Fr(1, b)) / pi
    phi0 = base + 48 * c ** 4
    b1 = 2 / 3 * 2.0 ** (N - 4) * f + a / (6 * N)

    def g(x):
        q = 48 * c ** 4 * exp(-2 * x)
        return x ** 3 - 2 * c ** 3 * x ** 2 - 8 * b1 * c ** 6 * log(1 / (base + q * (1 - q) ** (-N / 4)))
    lo, hi, ainv = 1e-9, 10.0, None
    if g(lo) < 0 < g(hi):
        for _ in range(200):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
        ainv = 1 / lo
    return {"split": [a, b], "c3": c, "phi0": phi0, "b1": b1, "alpha_inverse": ainv,
            "cabibbo": (phi0 * (1 - phi0)) ** 0.5 if 0 < phi0 < 1 else None, "solar": 1 / f - phi0 / 2,
            "baryon": (1 - 2 * c) * phi0}


def _hits(numbers):
    return {c for c, (key, value, sigma) in WORLD_TARGETS.items()
            if numbers[key] is not None and abs(numbers[key] - value) <= WORLD_SIGMAS * sigma}


def worlds(bounds=WORLD_BOUNDS):
    """All seeds inside the bounds with their computed clauses. Numbers depend on (k, N, f); with several positive
    splits the best one is kept, which favours the baseline."""
    from math import comb
    from genesis_rules import anti_isometric
    rng = {c: range(bounds[c][0], bounds[c][1] + 1) for c in COORDS}
    glue = {(N, f): N >= 4 and anti_isometric([("D", N)], [("A", f)]) for N in rng["N"] for f in rng["f"]}
    glue_a = {(N, f): anti_isometric([("A", N - 1)], [("A", f)]) for N in rng["N"] for f in rng["f"]}
    numeric = {}
    for k in rng["k"]:
        for N in rng["N"]:
            for f in rng["f"]:
                rows = [world_numbers(k, N, f, s) for s in positive_splits(N)]
                numeric[k, N, f] = max(rows, key=lambda r: len(_hits(r))) if rows else None
    out = []
    for k, n, N, f, d in product(*[rng[c] for c in COORDS]):
        num = numeric[k, N, f]
        cl = {c for c, ok in (
            ("flat_cover", n * (k - 1) == 2 * k), ("region_pairs", n % 2 == 0),
            ("spin_sheets", (d >= 4 and k == 2) or (d == 2 and k == 1)), ("holomorphic_sky", d == 4),
            ("metric_marks", n == 4), ("critical_index", k * k == 4), ("index_marks", k * k == n),
            ("lift_marks", 2 * k == n), ("divisor_carrier", N == n + 1), ("divisor_families", f == n - 1),
            ("z4_anomaly", N % 2 == 1 and 2 ** (N - 1) % 16 == 0), ("unimodular_glue", glue[N, f]),
            ("family_rotation", f == d - 1 and comb(f, 2) == f),
            ("conformal_family_factor", (f, d) in ((3, 4), (1, 1))),
            ("positive_split", len(positive_splits(N)) == 1),
            ("three_sheets", k == 3), ("type_a_glue", glue_a[N, f]), ("anyonic_space", d == 3),
            ("cross_product_seven", f == 7 and d == 8), ("six_sphere_sky", d == 8),
            ("balanced_carrier", N % 2 == 0), ("odd_marks", n % 2 == 1)) if ok}
        if num:
            cl |= _hits(num)
        out.append({"seed": {"k": k, "n": n, "N": N, "f": f, "d": d}, "label": f"k{k} n{n} N{N} f{f} d{d}",
                    "size": k + n + N + f + d, "clauses": cl, "numbers": num})
    return out


def _grouped(objs, pool):
    """Objects grouped by clause mask: masks, least size per group, number of objects of least size, one index."""
    import numpy as np
    bit = {c: 1 << i for i, c in enumerate(pool)}
    groups = {}
    for i, o in enumerate(objs):
        m = sum(bit[c] for c in o["clauses"] if c in bit)
        g = groups.setdefault(m, [o["size"], 0, i])
        if o["size"] < g[0]:
            g[:] = [o["size"], 1, i]
        elif o["size"] == g[0]:
            g[1] += 1
    keys = list(groups)
    return (bit, np.array(keys, dtype=np.int64), np.array([groups[m][0] for m in keys]),
            np.array([groups[m][1] for m in keys]), np.array([groups[m][2] for m in keys]))


def world_rule_search(objs=None, pool=None, target=None, max_len=5):
    """Every rule of at most max_len clauses selects its smallest worlds, as genesis_rules.search does. Reports the
    verdict counts, how many rules select which world uniquely, and for the target the minimal rules with their
    ablations (which clause is load bearing)."""
    from itertools import combinations
    import numpy as np
    objs = objs if objs is not None else worlds()
    pool = list(pool or WORLD_POOL)
    target = target or TFPT_WORLD
    bit, masks, least, count, index = _grouped(objs, pool)
    tgt = next(i for i, o in enumerate(objs) if o["seed"] == target)

    def select(rule):
        sel = (masks & rule) == rule
        if not sel.any():
            return "EMPTY", None
        size = least[sel].min()
        hit = sel & (least == size)
        return ("UNIQUE", int(index[hit][0])) if int(count[hit].sum()) == 1 else ("AMBIGUOUS", None)
    verdicts, makers, selected = {"EMPTY": 0, "UNIQUE": 0, "AMBIGUOUS": 0}, [], {}
    rules = [sum(bit[c] for c in combo) for r in range(max_len + 1) for combo in combinations(pool, r)]
    for rule in rules:
        verdict, i = select(rule)
        verdicts[verdict] += 1
        if verdict == "UNIQUE":
            selected[i] = selected.get(i, 0) + 1
            if i == tgt:
                makers.append(rule)
    names = lambda rule: sorted(c for c in pool if bit[c] & rule)
    minimal = []
    for rule in sorted(makers, key=lambda r: bin(r).count("1")):
        if not any(m & rule == m for m in minimal):
            minimal.append(rule)
    rows = []
    for rule in minimal:
        ablation = {}
        for c in names(rule):
            verdict, i = select(rule & ~bit[c])
            ablation[c] = verdict if i is None else objs[i]["label"]
        rows.append({"rule": names(rule), "without": ablation})
    by_length = {}
    for row in rows:
        by_length[len(row["rule"])] = by_length.get(len(row["rule"]), 0) + 1
    return {"worlds": len(objs), "clauses": pool, "max_rule_length": max_len, "rules": len(rules),
            "verdicts": verdicts, "target": objs[tgt]["label"], "rules_selecting_target": len(makers),
            "minimal_rules": rows, "minimal_rules_by_length": by_length,
            "uniquely_selected_worlds": sorted(((objs[i]["label"], c) for i, c in selected.items()),
                                               key=lambda r: -r[1])[:12]}


def coherent_worlds(objs=None):
    """Worlds that close: flat seam, carrier and families from the mark divisor, Z4 anomaly, a gluing to an even
    unimodular lattice of either type and a unique positive split. Computed over the whole universe; the clause
    pool of the rule search plays no part."""
    objs = objs if objs is not None else worlds()
    passed = lambda o: set(CLOSURE) <= o["clauses"] and bool({"unimodular_glue", "type_a_glue"} & o["clauses"])
    closing = [o for o in objs if passed(o)]
    seeds = sorted({(o["seed"]["k"], o["seed"]["n"], o["seed"]["N"], o["seed"]["f"]) for o in closing})
    fails = {}
    for o in objs:
        if set(CLOSURE[:3]) <= o["clauses"]:                      # flat with carrier and families from the divisor
            key = f"k{o['seed']['k']} n{o['seed']['n']} N{o['seed']['N']} f{o['seed']['f']}"
            fails[key] = sorted(c for c in ("z4_anomaly", "positive_split", "unimodular_glue", "type_a_glue")
                                if c not in o["clauses"])
    with_space = [o["label"] for o in closing
                  if {"holomorphic_sky", "spin_sheets", "family_rotation", "conformal_family_factor"} <= o["clauses"]]
    return {"closing_seeds_without_dimension": [list(x) for x in seeds], "flat_divisor_worlds_and_failures": fails,
            "closing_with_spacetime_clauses": with_space,
            "dimensions_allowed_by_closure_alone": sorted({o["seed"]["d"] for o in closing})}


def consequences(objs, rule):
    """What a rule forces: values of every seed coordinate and clauses shared by all worlds that satisfy it, and
    the same for the smallest ones."""
    rule = set(rule)
    sat = [o for o in objs if rule <= o["clauses"]]
    if not sat:
        return {"worlds": 0}
    size = min(o["size"] for o in sat)
    small = [o for o in sat if o["size"] == size]
    summary = lambda rows: ({c: sorted({o["seed"][c] for o in rows}) for c in COORDS},
                            sorted(set.intersection(*[o["clauses"] for o in rows])))
    values, implied = summary(sat)
    min_values, min_implied = summary(small)
    return {"worlds": len(sat), "values": values, "implied": implied, "smallest": [o["label"] for o in small][:6],
            "smallest_count": len(small), "minimal_values": min_values, "minimal_implied": min_implied}


def number_baseline(objs=None):
    """Whole-world number control: how many seeds (k, N, f) reach how many of the measured numbers."""
    objs = objs if objs is not None else worlds()
    seen = {}
    for o in objs:
        key = (o["seed"]["k"], o["seed"]["N"], o["seed"]["f"])
        if key not in seen:
            seen[key] = (o["numbers"], sorted(o["clauses"] & set(WORLD_TARGETS)))
    histogram = {}
    for _, hits in seen.values():
        histogram[len(hits)] = histogram.get(len(hits), 0) + 1
    per_target = {c: sorted(f"k{k} N{N} f{f}" for (k, N, f), (_, hits) in seen.items() if c in hits)
                  for c in WORLD_TARGETS}
    key = (TFPT_WORLD["k"], TFPT_WORLD["N"], TFPT_WORLD["f"])
    near = {}
    for (k, N, f), (num, hits) in seen.items():
        if num and abs(k - key[0]) + abs(N - key[1]) + abs(f - key[2]) == 1:
            near[f"k{k} N{N} f{f}"] = {"alpha_inverse": num["alpha_inverse"], "cabibbo": num["cabibbo"],
                                       "solar": num["solar"], "baryon": num["baryon"], "hits": hits}
    return {"seeds": len(seen), "with_numbers": sum(1 for num, _ in seen.values() if num),
            "hits_histogram": dict(sorted(histogram.items())), "per_target": per_target,
            "tfpt": {"numbers": seen[key][0], "hits": seen[key][1]}, "neighbours": near}


# ---- roles of a complete solution (TFPT-STAND, sections 1, 2, 3 and 5)
COMPUTED_ROLES = {"sheets2": (["k"], []), "marks4": (["n"], []), "flat": ([], ["flat_cover"]),
                  "carrier5": (["N"], ["positive_split"]), "generation16": (["N"], ["z4_anomaly"]),
                  "families3": (["f"], []), "e8": (["N", "f"], ["unimodular_glue"]), "lorentz31": (["d"], []),
                  "alpha": (["k", "N", "f"], ["alpha_measured"]),
                  "mixing": (["k", "N", "f"], ["cabibbo", "solar", "baryon"])}
RECORDED_ROLES = ["square", "c3", "winding", "rp", "state", "time", "causal", "continuum", "einstein", "area",
                  "flavour_rest"]
CHAIN_ROLES = ["clocks", "seed", "process"]               # recorded as exact once P1 and P2 stand (section 3)
ROLES = list(COMPUTED_ROLES) + RECORDED_ROLES + CHAIN_ROLES
ROLE_MEANING = {
    "sheets2": "two sheets", "marks4": "four marks", "flat": "flat seam without scale",
    "square": "square position of the marks, tau = i, mu4", "c3": "normalisation c3 = 1/(8 pi) with its quarter",
    "winding": "orientation and unit winding", "rp": "reflection positivity", "state": "which state",
    "time": "continuous physical time", "carrier5": "five places with 3 + 2", "generation16": "one generation of 16",
    "families3": "three families", "e8": "closure D5 + A3 to E8", "clocks": "clocks (30, mu4)",
    "seed": "seed phi0 and masses", "alpha": "alpha", "mixing": "mixing and cosmology numbers",
    "process": "process structure (GNS, 68 = 5 + 63)", "lorentz31": "3 + 1 and Lorentz",
    "causal": "local causal structure", "continuum": "continuum and infinite environment",
    "einstein": "gravitational dynamics (spin-2 pole)", "area": "area normalisation G",
    "flavour_rest": "rest flavour (quark amplitude, neutrino branch, phases)"}
# class of a recorded role inside a reading: own = follows from the basic object alone; extra = holds under a
# further named premise or choice; input = presupposed by the object
READINGS = {
    "purity": {
        "field": "quantum information of a vacuum",
        "object": "purity Tr rho^2 of a chiral vacuum on two disjoint regions (two copies, exchanged on the regions)",
        "clauses": ["flat_cover", "region_pairs"],
        "roles": {"square": "extra", "c3": "extra", "winding": "extra", "rp": "extra", "state": "extra",
                  "time": "extra", "einstein": "extra", "area": "extra"},
        "free_choice": "why the purity of exactly two regions is the basic observable, and which global state",
        "source": "findings source-as-purity-of-two-opposite-arcs-20261118, two-copies-are-one-charged-seam-field-"
                  "20261118; f1_two_arc_purity.py"},
    "closed_2p1": {
        "field": "gravitation in 2+1 dimensions",
        "object": "closed spatial slice of 2+1 Einstein gravity with identical conical particles",
        "clauses": ["flat_cover"],
        "roles": {"c3": "extra", "area": "extra", "time": "extra"},
        "free_choice": "why 2+1 and identical particles; the unit G = 1",
        "source": "new in this exploration (w2, B2); the flat pillow is QGEO.PILLOW.01"},
    "sky_spin2": {
        "field": "spinor geometry of spacetime",
        "object": "a spin-2 state on the sphere of directions of an observer",
        "clauses": ["holomorphic_sky", "spin_sheets", "metric_marks"],
        "roles": {"causal": "extra", "einstein": "extra", "continuum": "input"},
        "free_choice": "presupposes a Lorentzian manifold with a gravitational field; J = 0 (tau = i) is not explained",
        "source": "VOGELPERSPEKTIVE (Z61, Z62), CELEST.SEAM.01, premise T'"},
    "critical_index": {
        "field": "operator algebras and critical statistical mechanics",
        "object": "the smallest inclusion with infinite depth (Jones index 4, loop weight 2, four-state Potts)",
        "clauses": ["critical_index", "index_marks"],
        "roles": {"c3": "extra", "state": "extra", "continuum": "extra"},
        "free_choice": "why the world sits at the critical index and not above it",
        "source": "P1.INDEX.KMS.01 (v813): I beta c3 = 1 with I = 4, bridge open"},
    "anomaly_minimal": {
        "field": "topology of anomalies",
        "object": "the smallest chiral system free of the Z4 and the modular anomaly",
        "clauses": ["z4_anomaly", "unimodular_glue"],
        "roles": {},
        "free_choice": "why Z4 with square (-1)^F, and why the smallest solution",
        "source": "exploration 2026-11-23 (F3), R102, gluing equation (TFPT-STAND, section 3)"},
    "isometry": {
        "field": "operator theory",
        "object": "an irreducible isometry that is not unitary (unilateral shift, defect of rank 1)",
        "clauses": [],
        "roles": {"winding": "own", "rp": "extra"},
        "free_choice": "irreducibility; the shift carries no continuous time",
        "source": "Z30 (Toeplitz origin), Z32"},
    "diamond": {
        "field": "local quantum physics and horizon thermodynamics",
        "object": "the vacuum on the causal diamond of an observer, with entanglement equilibrium",
        "clauses": ["holomorphic_sky"],
        "roles": {"state": "extra", "time": "extra", "rp": "extra", "causal": "extra", "einstein": "extra",
                  "area": "extra", "continuum": "input"},
        "free_choice": "presupposes a Lorentzian spacetime and its vacuum; premise D (Jacobson)",
        "source": "H6 premise D, Z58, Z70, Z71, Z72"},
    "selfdual_families": {
        "field": "exterior algebra of the families",
        "object": "the families are the self-dual two-forms of spacetime",
        "clauses": ["family_rotation"],
        "roles": {},
        "free_choice": "the identification of family directions with rotations of space",
        "source": "H6 (Lambda^2_+ = Lambda^2 of the families), VOGELPERSPEKTIVE (h0(O(2)) = 3)"},
    "background_free": {
        "field": "two-dimensional conformal field theory",
        "object": "a seam without background geometry: basis changes are symmetries and <T> = 0 (H6, premise A*)",
        "clauses": ["flat_cover", "unimodular_glue"],
        "roles": {"square": "own", "c3": "own", "state": "own", "rp": "extra", "winding": "extra", "time": "extra",
                  "causal": "extra", "einstein": "extra", "area": "extra"},
        "free_choice": "why the ground state carries the twisted defect; premises C'', D, E, G and T'",
        "source": "TFPT-STAND, candidate H6; reconstruction_search.complete_solution_distance"},
    "real_form": {
        "field": "real forms of Lie algebras",
        "object": "the closed lattice algebra in its real form with compact gauge factor and non-compact family "
                  "factor: E8(-24) with Spin(10) x SU(2,2)",
        "clauses": ["unimodular_glue", "conformal_family_factor"],
        "roles": {"causal": "own", "rp": "extra", "time": "extra", "continuum": "extra", "einstein": "extra"},
        "free_choice": "why the physical real form is not the compact one; that the 4 of A3 is the seam doublet "
                       "with its Tomita partner",
        "source": "new in this exploration (w2, B6); ladder_gates for rp, time and continuum of the minimal ladder"}}
# identifications between basic objects: own = a theorem about the two objects, extra = a premise
BRIDGES = {
    ("purity", "critical_index"): ("own", "the two-region index of the Z_k orbifold of a holomorphic net is k^2 "
                                          "(Kawahigashi, Longo, Mueger)"),
    ("purity", "closed_2p1"): ("extra", "the replica cone of the purity is the particle: tension (n-1)/(4n) = 1/8 "
                                        "per mark, deficit pi (w2, B2); that the cone gravitates is a premise"),
    ("purity", "background_free"): ("extra", "the two copies are one charged seam field (finding 20261118)"),
    ("purity", "diamond"): ("extra", "the two regions are an observer's diamond and a second one on the same sky"),
    ("sky_spin2", "diamond"): ("own", "the edge of a diamond is the sphere of directions"),
    ("sky_spin2", "selfdual_families"): ("own", "h0(O(2)) = 3 on the same sphere (Z61)"),
    ("sky_spin2", "real_form"): ("own", "SU(2,2) acts on the lines of twistor space, the sky is one of them"),
    ("sky_spin2", "background_free"): ("extra", "premise T' (seam = twistor line), equivalent to I_pos"),
    ("real_form", "anomaly_minimal"): ("own", "both speak about the same even unimodular lattice"),
    ("real_form", "background_free"): ("extra", "the su(2,2) of the seam doublet with Tomita partner is the real "
                                                "form of A3"),
    ("real_form", "selfdual_families"): ("own", "the self-dual two-forms are su(2)_L inside su(2,2)"),
    ("real_form", "diamond"): ("own", "SU(2,2) is the conformal group of the diamonds of R^(3,1)"),
    ("background_free", "anomaly_minimal"): ("own", "the GSO sum over spin structures gives the even unimodular "
                                                    "lattice (H6)"),
    ("background_free", "isometry"): ("extra", "the Hardy compression of the seam is the shift (Z30, premise "
                                               "irreducibility)"),
    ("critical_index", "background_free"): ("extra", "I beta c3 = 1 (v813), bridge open"),
    ("closed_2p1", "diamond"): ("extra", "the 2+1 slice is the normal geometry of the diamond edge")}
_RANK = {"own": 3, "extra": 2, "input": 1, "none": 0}
_CONSEQUENCES = {}                                         # consequences of a rule, per universe object


def ladder_gates(d_max=7, size=400):
    """Computed gates for the dynamic roles on the smallest positive-energy ladder of so(d,2) (massless scalar,
    lowest weight (d-2)/2): positivity of the kernel on the sl(2) string, a positive time translation whose
    spectrum closes up to a continuum, and the Weyl dimension of the degeneracies. Control: so(6), compact, has
    finite-dimensional representations only. These are statements about the ladder; that the response of the
    original is this ladder is premise T'."""
    import numpy as np
    from math import comb
    out = {}
    for d in range(3, d_max + 1):
        h = Fr(d - 2, 2)
        norms = [Fr(1)]
        for n in range(1, 40):
            norms.append(norms[-1] * n * (2 * h + n - 1))                 # |L_-1^n h>^2 = n! (2h)_n
        gaps = {}
        for m in (size // 4, size):
            diag = [float(h) + n for n in range(m)]
            off = [-0.5 * ((n + 1) * (2 * float(h) + n)) ** 0.5 for n in range(m - 1)]
            ev = np.linalg.eigvalsh(np.diag(diag) + np.diag(off, 1) + np.diag(off, -1))   # L0 - (L1 + L-1)/2
            gaps[m] = (float(ev[0]), float(ev[1] - ev[0]))
        deg = lambda n: comb(n + d - 1, d - 1) - (comb(n + d - 3, d - 1) if n >= 2 else 0)
        counts = [sum(deg(j) for j in range(n + 1)) for n in (200, 400)]
        weyl = log(counts[1] / counts[0]) / log(2)
        out[f"so({d},2)"] = {
            "lowest_weight": str(h), "kernel_positive": all(x > 0 for x in norms),
            "time_positive": gaps[size][0] > 0,
            "time_continuous": gaps[size][0] < gaps[size // 4][0] / 2 and gaps[size][1] < gaps[size // 4][1] / 2,
            "bottom_and_gap": {str(m): [round(x, 6) for x in g] for m, g in gaps.items()},
            "weyl_dimension_of_space": round(weyl, 2), "spacetime_dimension": round(weyl) + 1}
    out["so(6) control"] = {"kernel_positive": True, "time_positive": False, "time_continuous": False,
                            "weyl_dimension_of_space": 0,
                            "reason": "compact: every unitary irreducible representation is finite dimensional"}
    return out


def assembly(members, objs, readings=None, bridges=None, recorded=None):
    """Roles and distance of a set of readings."""
    readings = readings or READINGS
    bridges = BRIDGES if bridges is None else bridges
    recorded = RECORDED_ROLES if recorded is None else recorded
    all_roles = list(COMPUTED_ROLES) + list(recorded) + CHAIN_ROLES
    rule = sorted(set(CHAIN) | {c for m in members for c in readings[m]["clauses"]})
    key = (id(objs), tuple(rule))
    if key not in _CONSEQUENCES:
        _CONSEQUENCES[key] = consequences(objs, rule)
    cons = _CONSEQUENCES[key]
    if not cons["worlds"]:
        return {"members": list(members), "rule": rule, "contradiction": True}
    roles = {}
    for role, (coords, clauses) in COMPUTED_ROLES.items():
        forced = (all(cons["values"][c] == [TFPT_WORLD[c]] for c in coords)
                  and all(c in cons["implied"] for c in clauses))
        least = (all(cons["minimal_values"][c] == [TFPT_WORLD[c]] for c in coords)
                 and all(c in cons["minimal_implied"] for c in clauses))
        roles[role] = "forced" if forced else "minimal" if least else "none"
    for role in recorded:
        roles[role] = max((readings[m]["roles"].get(role, "none") for m in members), key=_RANK.get, default="none")
    p1 = all(roles[r] in ("forced", "minimal") for r in ("sheets2", "marks4", "flat"))
    for role in CHAIN_ROLES:
        roles[role] = ("chain" if p1 and roles["square"] == "own" and roles["c3"] == "own"
                       else "extra" if p1 else "none")
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for (a, b), (cls, _) in bridges.items():
        if cls == "own" and a in parent and b in parent:
            parent[find(a)] = find(b)
    missing = len({find(m) for m in members}) - 1
    delivered = {"forced", "own", "chain"}
    open_strict = [r for r in all_roles if roles[r] not in delivered]
    open_minimal = [r for r in all_roles if roles[r] not in delivered | {"minimal"}]
    return {"members": list(members), "rule": rule, "contradiction": False, "roles": roles,
            "seed_forced": all(cons["values"][c] == [TFPT_WORLD[c]] for c in COORDS),
            "seed_smallest": cons["smallest"] if cons["smallest_count"] <= 6 else cons["smallest_count"],
            "open": open_minimal, "open_without_minimality": open_strict, "bridges_missing": missing,
            "distance": len(open_minimal) + missing, "complete": not open_minimal and missing == 0}


def run_assembly_search(objs=None, readings=None, bridges=None, limit=4, recorded=None):
    """All sets of at most `limit` readings and the set of all readings, ranked by distance to a complete solution."""
    from itertools import combinations
    objs = objs if objs is not None else worlds()
    readings = readings or READINGS
    names = sorted(readings)
    recorded = RECORDED_ROLES if recorded is None else recorded
    all_roles = list(COMPUTED_ROLES) + list(recorded) + CHAIN_ROLES
    rows = [assembly(c, objs, readings, bridges, recorded) for r in range(1, limit + 1)
            for c in combinations(names, r)]
    rows.append(assembly(names, objs, readings, bridges, recorded))
    good = [r for r in rows if not r["contradiction"]]
    good.sort(key=lambda r: (r["distance"], len(r["members"])))
    singles = {r["members"][0]: {"seed_forced": r["seed_forced"], "seed_smallest": r["seed_smallest"],
                                 "distance": r["distance"], "open": r["open"]}
               for r in good if len(r["members"]) == 1}
    never = [role for role in all_roles
             if all(r["roles"][role] not in ("forced", "minimal", "own", "chain") for r in good)]
    return {"status": "ASSEMBLIES_SEARCHED", "physical_solution": False, "readings": len(names),
            "assemblies": len(rows), "contradictions": len(rows) - len(good),
            "complete": [r["members"] for r in good if r["complete"]],
            "best": [{k: r[k] for k in ("members", "distance", "open", "bridges_missing", "seed_forced")}
                     for r in good[:8]],
            "all_readings": good and next(r for r in good if len(r["members"]) == len(names)),
            "single_readings": singles, "roles_no_reading_delivers": never,
            "smallest_seed_forcing": sorted((r["members"] for r in good if r["seed_forced"]), key=len)[:6]}


def completeness_audit(objs=None):
    """Can this search contain and recognise a complete solution? A planted reading that pins the seed and delivers
    every recorded role must come out complete, and must stop being complete when one role or one clause is taken
    away. The audit also states which roles the search computes and which it only reads from the record."""
    objs = objs if objs is not None else worlds()
    planted = {"field": "control", "object": "planted complete solution", "free_choice": "none",
               "clauses": ["holomorphic_sky", "spin_sheets", "metric_marks"],
               "roles": {r: "own" for r in RECORDED_ROLES}, "source": "control"}
    full = assembly(["planted"], objs, {"planted": planted}, {})
    less_role = dict(planted, roles={r: "own" for r in RECORDED_ROLES if r != "time"})
    less_clause = dict(planted, clauses=["holomorphic_sky", "spin_sheets"])
    split = assembly(["planted", "partner"], objs,
                     {"planted": planted, "partner": dict(planted, clauses=[], roles={})}, {})
    present = run_assembly_search(objs)
    return {"status": "COMPLETE_SOLUTION_REPRESENTABLE" if full["complete"] else "NOT_REPRESENTABLE",
            "roles": len(ROLES), "computed_from_the_seed": list(COMPUTED_ROLES),
            "read_from_the_record": RECORDED_ROLES + CHAIN_ROLES,
            "planted_complete": full["complete"],
            "planted_without_time": assembly(["planted"], objs, {"planted": less_role}, {})["open"],
            "planted_without_marks_clause": assembly(["planted"], objs, {"planted": less_clause}, {})["open"],
            "planted_with_unbridged_partner": split["bridges_missing"],
            "present_complete": present["complete"], "present_best_distance": present["best"][0]["distance"],
            "roles_no_reading_delivers": present["roles_no_reading_delivers"],
            "limit": "roles read from the record are statuses with sources; the search does not prove them. A new "
                     "reading enters with its clauses (computed) and its recorded roles; ladder_gates computes "
                     "rp, time and continuum for ladder objects only."}


# single statements whose proof would change recorded classes: which readings gain which roles, which bridges close
UPGRADES = {
    "ladder": {
        "text": "the response of the original is the smallest ladder of so(4,2) (premise T' / I_pos read on the real "
                "form): ladder_gates then gives positivity, a continuous positive time and Weyl dimension 3 + 1",
        "roles": {"real_form": ["rp", "time", "continuum"]},
        "bridges": [("real_form", "background_free"), ("sky_spin2", "background_free")]},
    "equilibrium": {
        "text": "entanglement equilibrium of the seam state (premise D, Jacobson): Einstein dynamics with the area "
                "coefficient of the state",
        "roles": {"background_free": ["einstein", "area"], "diamond": ["einstein", "area"]}, "bridges": []},
    "irreducible": {
        "text": "the Hardy compression of the seam is an irreducible isometry (Z30) or the ground state carries the "
                "twist (premise G)",
        "roles": {"background_free": ["winding"]}, "bridges": [("background_free", "isometry")]},
    "vacuum": {
        "text": "the state on an observer's diamond is the vacuum of a local net (Bisognano-Wichmann, Hislop-Longo)",
        "roles": {"diamond": ["state", "time", "rp", "causal"]}, "bridges": []},
    "index_kms": {
        "text": "the even modular response normalises I beta c3 = 1 (open bridge of P1.INDEX.KMS.01)",
        "roles": {"critical_index": ["c3", "state"]}, "bridges": [("critical_index", "background_free")]},
    "cone_gravitates": {
        "text": "the replica cone of the purity gravitates: tension (n-1)/(4n) with G = 1",
        "roles": {"closed_2p1": ["c3", "area"], "purity": ["area"]}, "bridges": [("purity", "closed_2p1")]}}


def run_leverage(objs=None, limit=3, readings=None, bridges=None, upgrades=None, recorded=None):
    """Which statements, once proved, move the distance to a complete solution? Every subset of UPGRADES is applied
    to the readings and bridges and the assembly search is repeated."""
    from copy import deepcopy
    from itertools import combinations
    objs = objs if objs is not None else worlds()
    base_readings, base_bridges = readings or READINGS, BRIDGES if bridges is None else bridges
    upgrades = upgrades or UPGRADES
    table = {}
    for r in range(len(upgrades) + 1):
        for combo in combinations(sorted(upgrades), r):
            readings, bridges = deepcopy(base_readings), dict(base_bridges)
            for name in combo:
                for reading, roles in upgrades[name]["roles"].items():
                    readings[reading]["roles"].update({role: "own" for role in roles})
                for a, b in upgrades[name]["bridges"]:
                    key = (a, b) if (a, b) in bridges else (b, a)
                    bridges[key] = ("own", bridges[key][1])
            best = run_assembly_search(objs, readings, bridges, limit, recorded)["best"][0]
            table[combo] = (best["distance"], best["members"], best["open"])
    base = table[()][0]
    floor = min(v[0] for v in table.values())
    smallest = min(len(c) for c, v in table.items() if v[0] == floor)
    return {"status": "LEVERAGE_SEARCHED", "physical_solution": False, "distance_now": base,
            "single": {c[0]: {"distance": v[0], "gain": base - v[0]} for c, v in table.items() if len(c) == 1},
            "lowest_distance": floor, "statements_needed": smallest,
            "smallest_sets": [{"statements": list(c), "members": v[1], "open": v[2]}
                              for c, v in table.items() if v[0] == floor and len(c) == smallest],
            "texts": {k: v["text"] for k, v in upgrades.items()}}


# ---- section E: roles read off the frozen predictions, the ledger reading of the twistor bridge, Z4 gates
# Second pass of 27 November 2026.  The ledger already holds the twistor bridge (WOIT.THETA.FREE.01,
# WOIT.BETA3.DUALITY.01, WOIT.SU2.KIN.01); section D had not read it.  The frozen predictions add two roles a
# complete solution must deliver and that the 24 roles of section D did not name.
PREDICTION_ROLES = ["chirality", "parity_response"]
SHAPE_ROLES = RECORDED_ROLES + PREDICTION_ROLES
ROLE_MEANING.update({
    "chirality": "one chiral generation without mirror partners",
    "parity_response": "parity-odd response of the world state with a fixed sign (cosmic birefringence "
                       "beta = phi0/(4 pi), frozen; measured 0.24 +- 0.06 degrees)"})


def shape_readings():
    """The readings of section D with their classes for the two prediction roles, plus the ledger reading."""
    from copy import deepcopy
    readings = deepcopy(READINGS)
    readings["real_form"]["roles"]["chirality"] = "extra"         # Z4 rule of z4_gates under 'family_is_twistor'
    readings["background_free"]["roles"].update({"chirality": "extra", "parity_response": "extra"})
    readings["anomaly_minimal"]["roles"]["chirality"] = "input"
    readings["sky_spin2"]["roles"]["chirality"] = "input"
    readings["twistor_bridge"] = {
        "field": "twistor geometry (ledger)",
        "object": "twistor space C^4 of the seam with the clock rho4 = diag(i,1,i,1), the deck and the reflection",
        "clauses": ["holomorphic_sky"],
        "roles": {"causal": "own", "rp": "extra", "chirality": "extra"},
        "free_choice": "the interacting holomorphic algebra (contract WOIT.OS.TWISTOR.01); twistor space itself",
        "source": "ledger WOIT.THETA.FREE.01, WOIT.BETA3.DUALITY.01 (signature (2,2) forced, the null locus meets "
                  "the seam line in the seam circle), WOIT.SU2.KIN.01; recomputed in w3, part B"}
    return readings


SHAPE_BRIDGES = dict(BRIDGES)
SHAPE_BRIDGES.update({
    ("twistor_bridge", "sky_spin2"): ("own", "the sky is a line of twistor space"),
    ("twistor_bridge", "background_free"): ("extra", "premise T': the seam is a line of twistor space"),
    ("twistor_bridge", "real_form"): ("extra", "premise: the family factor A3 is the sl(4) of the seam's twistor "
                                               "space; the ledger keeps E8 internal (SDYM(E8) on PT/Gamma)")})
SHAPE_UPGRADES = dict(UPGRADES)
SHAPE_UPGRADES.update({
    "family_is_twistor": {
        "text": "the three family directions lie in the conformal algebra su(2,2) of the seam's twistor space, and "
                "the A3 factor of E8 (ledger SEL.THM.01) is that algebra: with compact Spin(10) the real form is "
                "E8(-24), fields are representations of (Spin(10) x SU(2,2))/Z4, and the Z4 rule of z4_gates gives "
                "one chiral generation with a single orientation bit.  Tested on 27 November (w5, "
                "family_clock_gates): the square of the family clock has the signature of the Cartan involution "
                "of su(2,2) and of no other real form of A3; in the ledger dictionary H^0(K(D)) the families are "
                "the so(2,1) of the seam line with the clock -Ad(rho4).  Refuted in the earlier form 'family "
                "directions = Cartan of su(2,2) with the clock rho4': as a Weyl element the family clock is of "
                "Coxeter class",
        "roles": {"real_form": ["chirality"]}, "bridges": [("twistor_bridge", "real_form")]},
    "oriented_state": {
        "text": "the world state is the oriented state of the seam and not the SU(2,2)-invariant vacuum; its "
                "helicity phase is c3 phi0, with the sign of the orientation bit",
        "roles": {"background_free": ["parity_response"]}, "bridges": []}})


# statements of SHAPE_UPGRADES that the record does not yet hold together in one object
SHAPE_TENSIONS = [{
    "statements": ["ladder", "oriented_state"],
    "text": "the record finds the continuous time in the vacuum sector only (Z72), while a non-zero birefringence "
            "needs a state that is not the SU(2,2)-invariant vacuum; a complete solution has to show both in one "
            "object: the representation is the ladder, the world state is not its vacuum"}, {
    "statements": ["family_is_twistor"],
    "text": "the record reads the family three-space in two ways: as the Cartan lattice of A3 with the clock as "
            "Coxeter element (SEL.THM.01, QGEO.STATE.01; in sl(4) the class with Ad-spectrum (3,4,4,4)) and as the "
            "Moebius algebra of the seam line with the clock -Ad(rho4) (h6h33, FLAV.QGEO.04; rho4 has (7,4,0,4)).  "
            "Both give the signature (1,-1,-1) for the square of the clock, but they are different subspaces of "
            "sl(4) with clocks in different classes; a complete solution has to say which one is the A3 factor of "
            "E8.  The reading does not deliver the transport numbers: Q, K, R, L commute with no clock of the "
            "integer D4 (v98, v146), and Q+ is not the character degree there"}]


def z4_gates():
    """Exact gates for the role 'chirality' and for the shape of 'parity_response' on the roots of E8 in the
    coordinates D5 (0-4) + D3 = A3 (5-7).  They are statements about E8 with its subgroup
    (Spin(10) x SU(2,2))/Z4; that fields are representations of this quotient and massless fields are helicity
    ladders is the premise named in SHAPE_UPGRADES['family_is_twistor']."""
    from itertools import product
    roots = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si, sj in product((1, -1), repeat=2):
                v = [Fr(0)] * 8
                v[i], v[j] = Fr(si), Fr(sj)
                roots.append(tuple(v))
    roots += [tuple(Fr(x, 2) for x in sg) for sg in product((1, -1), repeat=8) if sg.count(-1) % 2 == 0]
    turn = lambda r, y: sum(a * b for a, b in zip(r, y)) % 1
    y10 = (Fr(1, 2),) * 5 + (Fr(0),) * 3                         # centre of Spin(10)
    y4 = (Fr(0),) * 5 + (Fr(1, 2),) * 3                          # centre of SU(4), inverse generator
    cls = {r: int(turn(r, y10) * 4) for r in roots}
    theta = (Fr(0),) * 5 + (Fr(1, 2), Fr(1, 2), Fr(0))           # Cartan involution of E8(-24), trivial on D5
    quarter = (Fr(0),) * 7 + (Fr(1, 4),)                         # rotation by pi/2 in the time-like plane
    gauge = [r for r in roots if cls[r] == 0]
    rule = {o: {k: [str(Fr(c, 2)) for c in range(-2, 3) if (o * (c + 2) - k) % 4 == 0] for k in range(4)}
            for o in (1, -1)}                                    # centre of SU(2,2) on helicity C/2: i^(C+2)
    return {
        "roots": len(roots), "classes": {k: sum(1 for r in roots if cls[r] == k) for k in range(4)},
        "centres_identified": all(turn(r, tuple(a + b for a, b in zip(y10, y4))) == 0 for r in roots),
        "compact_dimension": 8 + sum(1 for r in roots if turn(r, theta) == 0),
        "half_period_is_theta_times_clock": all(
            turn(r, tuple(2 * t for t in quarter)) == (turn(r, theta) - Fr(cls[r], 4)) % 1 for r in roots),
        "full_period_is_fermion_parity": all(
            turn(r, tuple(4 * t for t in quarter)) == (Fr(0) if cls[r] % 2 == 0 else Fr(1, 2)) for r in roots),
        "helicity_u1_in_e8": any(all(sum(a * b for a, b in zip(r, g)) == 0 for g in gauge) for r in roots),
        "rule": rule,
        "one_helicity_per_16": all(len(rule[o][1]) == 1 and len(rule[o][3]) == 1 for o in rule),
        "mirror_excluded": all(rule[o][1] != rule[o][3] for o in rule),
        "orientation_swaps": rule[1][1] == rule[-1][3],
        "z4_trivial_on_photons": all((c + 2) % 4 == 0 for c in (-2, 2)),
        "reading": "class 0 (gauge) only with helicity +-1, class 2 (the 10) only with 0, the 16 with one helicity "
                   "and the 16bar with the other; the Z4 cannot rotate polarisation and E8 holds no helicity U(1), "
                   "so a parity-odd rotation is a response of the state, with the sign of the orientation bit"}


def family_clock_gates():
    """Exact gates for the family side of 'family_is_twistor', on objects of the ledger: the family clock is the
    rotation of the four marks on H1 = Z^4/Z(1,1,1,1) (GATE.QGEOM.02, QGEO.STATE.01), a Coxeter element of
    W(A3) = S4; the twistor clock is rho4 = diag(i,1,i,1) (WOIT.BETA3.DUALITY.01).  Diagonal elements are written
    as exponents of i.  Computed in experiments/exploration/2026-11-27-weltgenerator-realform (w5)."""
    from itertools import permutations, product
    perms = list(permutations(range(4)))
    ident, clock = (0, 1, 2, 3), (1, 2, 3, 0)
    mul = lambda p, q: tuple(p[q[i]] for i in range(4))
    trace = lambda p: sum(1 for i in range(4) if p[i] == i) - 1          # on the reflection representation
    minus = lambda p: (3 - trace(p)) // 2                                # dimension of the -1 space of an involution
    theta = mul(clock, clock)
    squares = {p: mul(p, p) for p in perms}
    roots_of_involutions = sorted({q for q in squares.values() if q != ident and mul(q, q) == ident})
    # Cartan involution on a maximally split Cartan: q disjoint transpositions for su(4-q, q); the outer forms
    # sl(4,R) and su*(4) act by minus a Weyl element and lie outside W(A3)
    inner = {"su(4)": ident, "su(3,1)": (2, 1, 0, 3), "su(2,2)": (2, 3, 0, 1)}
    roots = [tuple((1 if k == i else -1 if k == j else 0) for k in range(4)) for i, j in permutations(range(4), 2)]
    act = lambda p, v: tuple(v[p.index(k)] for k in range(4))
    as_map = lambda p, sign=1: tuple(tuple(sign * x for x in act(p, r)) for r in roots)
    weyl_maps = {as_map(p) for p in perms}
    outer_in_weyl = any(as_map(p, -1) in weyl_maps for p in (ident, (2, 3, 0, 1)))
    selected = [name for name, p in inner.items() if p in roots_of_involutions]
    centraliser = [p for p in perms if mul(p, theta) == mul(theta, p)]
    restricted = {}
    for r in roots:
        proj = tuple(Fr(a - b, 2) for a, b in zip(r, act(theta, r)))
        restricted[proj] = restricted.get(proj, 0) + 1
    shape = sorted((sum(x * x for x in v), m) for v, m in restricted.items())

    def ad(exponents):
        counts = [0, 0, 0, 0]
        for a, b in product(exponents, repeat=2):
            counts[(a - b) % 4] += 1
        return (counts[0] - 1, counts[1], counts[2], counts[3])
    rho4, coxeter, involution = (1, 0, 1, 0), (1, 0, 3, 2), (0, 0, 2, 2)
    fixed_six = sum(1 for a in range(4) for b in range(a + 1, 4) if (involution[a] + involution[b]) % 4 == 0)
    return {
        "clock_is_coxeter_element": mul(mul((1, 0, 2, 3), (0, 2, 1, 3)), (0, 1, 3, 2)) in (clock, mul(theta, clock)),
        # six elements of S4 square to a non-trivial involution; the squares are the three pairings, type (1,2)
        "squares_are_the_three_pairings":
            (sum(1 for q in squares.values() if q in roots_of_involutions), len(roots_of_involutions),
             {(3 - minus(q), minus(q)) for q in roots_of_involutions}) == (6, 3, {(1, 2)}),
        "clock_square_signature": (3 - minus(theta), minus(theta)),
        "real_form_selected": selected, "outer_forms_in_weyl_group": outer_in_weyl,
        "centraliser_order": len(centraliser), "real_rank": minus(theta),
        "restricted_roots": shape,
        "restricted_system_is_c2": shape == [(Fr(1), 2)] * 4 + [(Fr(2), 1)] * 4,
        "ad_spectrum_twistor_clock": ad(rho4), "ad_spectrum_family_clock": ad(coxeter),
        "clocks_conjugate": ad(rho4) in (ad(coxeter), (ad(coxeter)[0], ad(coxeter)[3], ad(coxeter)[2], ad(coxeter)[1])),
        "clocks_share_square": [2 * a % 4 for a in rho4] == [2 * a % 4 for a in coxeter],
        "e8_compact_dimension": 45 + ad(involution)[0] + 10 * fixed_six + 32 * involution.count(0),
        "reading": "the square of the family clock is a Cartan involution of type (2,2): su(2,2), real rank 2, the "
                   "family split B1 + E is the split t + a of the maximally split Cartan, and the flavour symmetry "
                   "D4 is the restricted Weyl group W(C2); the family clock and rho4 are different fourth roots of "
                   "that involution"}


def run_shape_search(objs=None, limit=3):
    """The assembly and leverage searches on the contract extended by the two prediction roles and with the ledger
    reading of the twistor bridge."""
    objs = objs if objs is not None else worlds()
    readings = shape_readings()
    asm = run_assembly_search(objs, readings, SHAPE_BRIDGES, limit, SHAPE_ROLES)
    lev = run_leverage(objs, limit, readings, SHAPE_BRIDGES, SHAPE_UPGRADES, SHAPE_ROLES)
    planted = {"field": "control", "object": "planted", "free_choice": "none",
               "clauses": ["holomorphic_sky", "spin_sheets", "metric_marks"],
               "roles": {r: "own" for r in SHAPE_ROLES}, "source": "control"}
    less = dict(planted, roles={r: "own" for r in SHAPE_ROLES if r != "parity_response"})
    return {"status": "SHAPE_SEARCHED", "physical_solution": False,
            "roles": len(COMPUTED_ROLES) + len(SHAPE_ROLES) + len(CHAIN_ROLES), "gates": z4_gates(),
            "family_gates": family_clock_gates(),
            "assemblies": asm, "leverage": lev, "tensions": SHAPE_TENSIONS,
            "planted_complete": assembly(["planted"], objs, {"planted": planted}, {}, SHAPE_ROLES)["complete"],
            "planted_without_parity_response": assembly(["planted"], objs, {"planted": less}, {},
                                                        SHAPE_ROLES)["open"]}


def run_world_search():
    objs = worlds()
    return {"status": "WORLDS_SEARCHED", "physical_solution": False, "bounds": WORLD_BOUNDS,
            "rules": world_rule_search(objs), "coherent": coherent_worlds(objs),
            "numbers": number_baseline(objs), "ladders": ladder_gates(),
            "assemblies": run_assembly_search(objs), "leverage": run_leverage(objs),
            "audit": completeness_audit(objs)}


if __name__ == "__main__":
    import json
    print(json.dumps(run_world_search(), indent=1, default=str))
