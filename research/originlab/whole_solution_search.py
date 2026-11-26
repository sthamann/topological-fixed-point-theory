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
