"""Erzeuger für die Szenariosuche. Jede Funktion gen_*() liefert Objekte {id, fam, name, feats}.
feats ist eine Menge von Tripeln (art, wert, wofür): art = 'zahl' | 'bruch' | 'typ'.
Neue Erzeuger: eine Funktion gen_*() hier oder in einer Datei gen_*.py ergänzen. Kein Erzeuger kennt den Fingerabdruck."""
import itertools, math
from fractions import Fraction as Fr
from functools import reduce

def Z(v, w): return ("zahl", int(v), w)
def B(v, w):
    v = Fr(v); return ("zahl", int(v), w) if v.denominator == 1 else ("bruch", f"{v.numerator}/{v.denominator}", w)
def T(v, w): return ("typ", str(v), w)
def obj(id_, fam, name, feats): return {"id": id_, "fam": fam, "name": name, "feats": {f for f in feats if f is not None}}

# ---------- ADE Bausteine ----------
ROOTS = lambda t, n: {"A": n * (n + 1), "D": 2 * n * (n - 1), "E": {6: 72, 7: 126, 8: 240}.get(n, 0)}[t]
DET = lambda t, n: {"A": n + 1, "D": 4, "E": {6: 3, 7: 2, 8: 1}.get(n, 0)}[t]
COX = lambda t, n: {"A": n + 1, "D": 2 * n - 2, "E": {6: 12, 7: 18, 8: 30}.get(n, 0)}[t]
def graph(t, n, ext=False):
    if t == "A":
        e = [(i, i + 1) for i in range(1, n)]
        if ext: e += [(0, 1)] + ([(0, n)] if n > 1 else [])
    elif t == "D":
        e = [(i, i + 1) for i in range(1, n - 2)] + [(n - 2, n - 1), (n - 2, n)]
        if ext: e += [(0, 2)]
    else:
        chain = [1, 3, 4, 5, 6, 7, 8][: n - 1]; e = list(zip(chain, chain[1:])) + [(2, 4)]
        if ext: e += [(0, {6: 2, 7: 1, 8: 8}[n])]
    return set(range(0 if ext else 1, n + 1)), e
def classify(nodes, edges):
    out = []; nodes = set(nodes); adj = {v: set() for v in nodes}
    for a, b in edges:
        if a in nodes and b in nodes: adj[a].add(b); adj[b].add(a)
    seen = set()
    for s in nodes:
        if s in seen: continue
        comp = {s}; stack = [s]
        while stack:
            for w in adj[stack.pop()]:
                if w not in comp: comp.add(w); stack.append(w)
        seen |= comp; m = len(comp); br = [v for v in comp if len(adj[v]) == 3]
        if not br: out.append(("A", m)); continue
        arms = []
        for w in adj[br[0]]:
            k, prev, cur = 1, br[0], w
            while True:
                nxt = [x for x in adj[cur] if x != prev]
                if not nxt: break
                prev, cur = cur, nxt[0]; k += 1
            arms.append(k)
        arms.sort(); out.append(("D", m) if arms[:2] == [1, 1] else ("E", m))
    return tuple(sorted(out))
def subsystems(parent):
    start = (parent,); seen = {start}; todo = [start]
    while todo:
        cur = todo.pop()
        for i, (t, n) in enumerate(cur):
            rest = cur[:i] + cur[i + 1:]
            for ext in (False, True):
                nodes, edges = graph(t, n, ext)
                for d in nodes:
                    new = tuple(sorted(rest + classify(nodes - {d}, edges)))
                    if new and new not in seen: seen.add(new); todo.append(new)
    return seen
def tname(c): return " + ".join(f"{t}{n}" for t, n in sorted(c, key=lambda x: (-x[1], x[0])))
def gen_teilgitter():
    for parent in (("E", 8), ("E", 7), ("E", 6)):
        prank = parent[1]
        for c in subsystems(parent):
            rank = sum(n for _, n in c)
            if rank < prank - 1: continue
            roots = sum(ROOTS(t, n) for t, n in c); det = reduce(lambda a, b: a * b, (DET(t, n) for t, n in c), 1)
            f = [T(f"{t}{n}", "Baustein") for t, n in c] + [T(f"E{prank}", "Elterngitter"), Z(roots, "Wurzeln"), Z(ROOTS(*parent) - roots, "Wurzeln außerhalb"), Z(rank, "Rang"), Z(roots + rank, "Dimension"), Z(len(c), "Zahl der Bausteine")]
            for t, n in c: f += [Z(n, f"Rang {t}{n}"), Z(ROOTS(t, n), f"Wurzeln {t}{n}"), Z(ROOTS(t, n) + n, f"Dimension {t}{n}"), Z(COX(t, n), f"Coxeterzahl {t}{n}"), Z(DET(t, n), f"Determinante {t}{n}")]
            for t, n in c:
                if t == "A": f += [B(Fr(k * (n + 1 - k), n + 1), f"Klebenorm A{n} Klasse {k}") for k in range(1, n // 2 + 2) if k <= n]
                if t == "D": f += [B(Fr(n, 4), f"Klebenorm D{n} Spinor"), B(Fr(n, 8), f"Spinfeldgewicht D{n}")]
                if t == "A" and n % 2 == 1: f.append(B(Fr(n, 8), f"Gewicht A{n} Mittelklasse/2"))
            if rank == prank:
                idx = math.isqrt(det // DET(*parent))
                if idx * idx * DET(*parent) == det:
                    f.append(Z(idx, "Index im Elterngitter"))
                    if len(c) == 2 and idx > 1 and {DET(t, n) for t, n in c} == {idx}: f.append(T(f"Z{idx}", "zyklische Klebung zweier Bausteine"))
            else: f.append(Z(1, "freie Richtung"))
            yield obj(f"gitter:{tname(c)}<E{prank}", "Teilgitter", f"{tname(c)} in E{prank}", f)

# ---------- Spiegelungsgruppen (Grade nach Shephard und Todd) ----------
DEG = {"G4": [4, 6], "G5": [6, 12], "G6": [4, 12], "G7": [12, 12], "G8": [8, 12], "G9": [8, 24], "G10": [12, 24], "G11": [24, 24], "G12": [6, 8], "G13": [8, 12], "G14": [6, 24], "G15": [12, 24],
       "G16": [20, 30], "G17": [20, 60], "G18": [30, 60], "G19": [60, 60], "G20": [12, 30], "G21": [12, 60], "G22": [12, 20], "G23 (H3)": [2, 6, 10], "G24": [4, 6, 14], "G25": [6, 9, 12], "G26": [6, 12, 18],
       "G27": [6, 12, 30], "G28 (F4)": [2, 6, 8, 12], "G29": [4, 8, 12, 20], "G30 (H4)": [2, 12, 20, 30], "G31": [8, 12, 20, 24], "G32": [12, 18, 24, 30], "G33": [4, 6, 10, 12, 18], "G34": [6, 12, 18, 24, 30, 42],
       "G35 (E6)": [2, 5, 6, 8, 9, 12], "G36 (E7)": [2, 6, 8, 10, 12, 14, 18], "G37 (E8)": [2, 8, 12, 14, 18, 20, 24, 30]}
def gen_spiegelungsgruppen():
    table = dict(DEG)
    for n in range(1, 9): table[f"W(A{n})"] = list(range(2, n + 2))
    for n in range(2, 9): table[f"W(B{n})"] = list(range(2, 2 * n + 1, 2))
    for n in range(4, 9): table[f"W(D{n})"] = sorted(list(range(2, 2 * n - 1, 2)) + [n])
    for m in (5, 6, 8, 10, 12): table[f"I2({m})"] = [2, m]
    for name, d in table.items():
        order = reduce(lambda a, b: a * b, d, 1)
        f = [Z(order, "Ordnung"), Z(len(d), "Rang"), Z(sum(x - 1 for x in d), "Zahl der Spiegelungen"), Z(reduce(math.gcd, d), "Ordnung des Zentrums"), Z(max(d), "höchster Grad")] + [Z(x, "Grad") for x in d]
        yield obj(f"spiegel:{name}", "Spiegelungsgruppe", f"Spiegelungsgruppe {name}", f)

# ---------- endliche Drehgruppen, McKay ----------
def gen_mckay():
    data = {"binär tetraedrisch (E6)": (24, [1, 1, 1, 2, 2, 2, 3], (2, 3, 3)), "binär oktaedrisch (E7)": (48, [1, 1, 2, 2, 2, 3, 3, 4], (2, 3, 4)), "binär ikosaedrisch (E8)": (120, [1, 2, 2, 3, 3, 4, 4, 5, 6], (2, 3, 5))}
    for n in range(2, 13): data[f"zyklisch Z{n} (A{n-1})"] = (n, [1] * n, (n, n))
    for n in range(2, 9): data[f"binär diedrisch Dic{n} (D{n+2})"] = (4 * n, [1, 1, 1, 1] + [2] * (n - 1), (2, 2, n))
    for name, (order, marks, sig) in data.items():
        f = [Z(order, "Ordnung"), Z(len(marks), "Zahl der Darstellungen"), Z(sum(marks), "Coxeterzahl (Summe der Marken)"), Z(max(marks), "größte Marke")] + [Z(m, "Marke") for m in set(marks) if m > 1] + [Z(s, "Drehordnung") for s in sig]
        if order % 2 == 0 and "zyklisch" not in name: f.append(Z(order // 2, "Ordnung der Drehgruppe"))
        yield obj(f"mckay:{name}", "Drehgruppe", name, f)

# ---------- Linsenräume ----------
def gen_linsen():
    for p in range(2, 17):
        for k in range(0, p // 2 + 1):
            eta = Fr(p * p - 1, 6 * p) - Fr(k * (p - k), p)
            f = [Z(p, "Ordnung des Schnitts"), B(abs(eta), "Randphase eta") if eta else Z(0, "Randphase verschwindet"), B(Fr(k * (p - k), p), "Klebenorm der Windung") if k else None, B(Fr(p * p - 1, 6 * p), "eta ohne Windung")]
            if k: f.append(Z(k, "Windung"))
            yield obj(f"linse:S3/Z{p},k={k}", "Linsenraum", f"Linsenraum S³/Z{p} mit Windung {k}", f)

# ---------- Modulkurven ----------
def _pf(n): return {p for p in range(2, n + 1) if n % p == 0 and all(p % q for q in range(2, int(p**.5) + 1))}
def _phi(n): return sum(1 for a in range(1, n + 1) if math.gcd(a, n) == 1)
def gen_modulkurven():
    for N in range(1, 26):
        mu = Fr(N) * reduce(lambda a, b: a * b, (1 + Fr(1, p) for p in _pf(N)), 1); cusps = sum(_phi(math.gcd(d, N // d)) for d in range(1, N + 1) if N % d == 0)
        def leg(a, p): return 0 if a % p == 0 else (1 if pow(a % p, (p - 1) // 2, p) == 1 else -1)
        nu2 = 0 if N % 4 == 0 else reduce(lambda a, b: a * b, (1 + (leg(-1, p) if p > 2 else 0) for p in _pf(N)), 1)
        nu3 = 0 if N % 9 == 0 else reduce(lambda a, b: a * b, (1 + (leg(-3, p) if p > 3 else (0 if p == 3 else -1)) for p in _pf(N)), 1)
        g = 1 + mu / 12 - Fr(nu2, 4) - Fr(nu3, 3) - Fr(cusps, 2)
        yield obj(f"modul:X0({N})", "Modulkurve", f"Modulkurve X₀({N})", [Z(N, "Stufe"), Z(mu, "Index"), Z(cusps, "Spitzen"), Z(g, "Geschlecht"), Z(nu2, "Punkte der Ordnung 2"), Z(nu3, "Punkte der Ordnung 3")])
        if N >= 3:
            mu = Fr(N**3, 2) * reduce(lambda a, b: a * b, (1 - Fr(1, p * p) for p in _pf(N)), 1); cusps = mu / N; g = 1 + mu / 12 - cusps / 2
            yield obj(f"modul:X({N})", "Modulkurve", f"Modulkurve X({N})", [Z(N, "Stufe"), Z(mu, "Index = Ordnung der Symmetriegruppe"), Z(cusps, "Spitzen"), Z(g, "Geschlecht"), Z(N, "Breite jeder Spitze")])
        if N >= 5:
            mu = Fr(N * N, 2) * reduce(lambda a, b: a * b, (1 - Fr(1, p * p) for p in _pf(N)), 1); cusps = Fr(sum(_phi(d) * _phi(N // d) for d in range(1, N + 1) if N % d == 0), 2); g = 1 + mu / 12 - cusps / 2
            yield obj(f"modul:X1({N})", "Modulkurve", f"Modulkurve X₁({N})", [Z(N, "Stufe"), Z(mu, "Index"), Z(cusps, "Spitzen"), Z(g, "Geschlecht"), T(f"Z{N}", "Torsion der Kurvenfamilie")])

# ---------- elliptische Flächen mit endlich vielen Schnitten (Miranda und Persson) ----------
FIB = {"I": lambda n: ("A", n - 1) if n > 1 else None, "I*": lambda n: ("D", n + 4), "II*": lambda n: ("E", 8), "III*": lambda n: ("E", 7), "IV*": lambda n: ("E", 6), "III": lambda n: ("A", 1), "IV": lambda n: ("A", 2), "II": lambda n: None}
EUL = {"I": lambda n: n, "I*": lambda n: n + 6, "II*": lambda n: 10, "III*": lambda n: 9, "IV*": lambda n: 8, "III": lambda n: 3, "IV": lambda n: 4, "II": lambda n: 2}
SURF = {"X22": ([("II", 0), ("II*", 0)], [1]), "X33": ([("III", 0), ("III*", 0)], [2]), "X44": ([("IV", 0), ("IV*", 0)], [3]), "X11": ([("I*", 0), ("I*", 0)], [2, 2]), "X211": ([("II*", 0), ("I", 1), ("I", 1)], [1]),
        "X411": ([("I*", 4), ("I", 1), ("I", 1)], [2]), "X9111": ([("I", 9), ("I", 1), ("I", 1), ("I", 1)], [3]), "X8211": ([("I", 8), ("I", 2), ("I", 1), ("I", 1)], [4]), "X6321": ([("I", 6), ("I", 3), ("I", 2), ("I", 1)], [6]),
        "X5511": ([("I", 5), ("I", 5), ("I", 1), ("I", 1)], [5]), "X4422": ([("I", 4), ("I", 4), ("I", 2), ("I", 2)], [4, 2]), "X3333": ([("I", 3)] * 4, [3, 3]), "X321": ([("III*", 0), ("I", 2), ("I", 1)], [2]),
        "X431": ([("IV*", 0), ("I", 3), ("I", 1)], [3]), "X222": ([("I*", 2), ("I", 2), ("I", 2)], [2, 2]), "X141": ([("I*", 1), ("I", 4), ("I", 1)], [4])}
def gen_flaechen():
    for name, (fibs, tors) in SURF.items():
        lat = [FIB[t](n) for t, n in fibs]; lat = [x for x in lat if x]; tor = reduce(lambda a, b: a * b, tors, 1)
        f = [Z(len(fibs), "Ausnahmefasern"), Z(tor, "Ordnung der Torsion"), Z(sum(n for _, n in lat), "Rang der Fasergitter"), Z(sum(EUL[t](n) for t, n in fibs), "Eulerzahl")] + [T(f"{t}{n}", "Fasergitter") for t, n in lat]
        f += [Z(n, "Breite einer Faser") for t, n in fibs if t == "I" and n > 1]
        if len(tors) == 1 and tor > 1: f.append(T(f"Z{tor}", "zyklische Torsion"))
        if sum(n for _, n in lat) == 8: f.append(T("E8", "Fasergitter füllen E8"))
        yield obj(f"flaeche:{name}", "Elliptische Fläche", f"elliptische Fläche {name} (" + ", ".join(t.replace("I", f"I{n}", 1) if t in ("I", "I*") else t for t, n in fibs) + ")", f)

# ---------- Codes ----------
def gen_codes():
    for m in range(2, 8):
        for r in range(0, m):
            n = 2**m; k = sum(math.comb(m, i) for i in range(r + 1)); d = 2**(m - r)
            f = [Z(n, "Länge"), Z(k, "Dimension"), Z(d, "Abstand"), B(Fr(k, n), "Rate"), Z(m, "Zahl der Variablen")]
            if 2 * r == m - 1: f += [T("selbstdual", "Code"), Z(m, "selbstdual bei m")]
            yield obj(f"code:RM({r},{m})", "Code", f"Reed Muller Code RM({r},{m}) [{n},{k},{d}]", f)
    for name, n, k, d in (("Hamming [7,4,3]", 7, 4, 3), ("erweiterter Hamming [8,4,4]", 8, 4, 4), ("Golay [23,12,7]", 23, 12, 7), ("erweiterter Golay [24,12,8]", 24, 12, 8), ("Hexacode [6,3,4] über F4", 6, 3, 4), ("ternärer Golay [12,6,6]", 12, 6, 6)):
        yield obj(f"code:{name}", "Code", name, [Z(n, "Länge"), Z(k, "Dimension"), Z(d, "Abstand"), B(Fr(k, n), "Rate")] + ([T("selbstdual", "Code")] if 2 * k == n else []))
    for name, n, k, d in (("[[5,1,3]]", 5, 1, 3), ("[[7,1,3]] Steane", 7, 1, 3), ("[[9,1,3]] Shor", 9, 1, 3), ("[[4,2,2]]", 4, 2, 2), ("[[8,3,3]]", 8, 3, 3), ("[[15,1,3]]", 15, 1, 3), ("[[6,4,2]]", 6, 4, 2), ("[[16,10,3]]", 16, 10, 3), ("[[8,0,4]]", 8, 0, 4)):
        yield obj(f"qcode:{name}", "Quantencode", f"Quantencode {name}", [Z(n, "Qubits"), Z(k, "logische Qubits") if k else None, Z(d, "Abstand"), Z(n - k, "Stabilisatoren"), Z(2**k, "Codedimension") if k else None, B(Fr(k, n), "Rate") if k else None])

# ---------- Qubits, Stabilisatorzustände, Clifford ----------
def gen_qudits():
    for d in (2, 3, 4, 5):
        for n in range(1, 6 if d == 2 else 3):
            D = d**n; f = [Z(D, "Dimension"), Z(D * D - 1, "Paulirichtungen"), Z(n, "Zahl der Systeme"), Z(d, "Stufen je System")]
            if d in (2, 3, 5):
                st = D * reduce(lambda a, b: a * b, (d**k + 1 for k in range(1, n + 1)), 1); f.append(Z(st, "Stabilisatorzustände"))
                cl = d**(n * n + 2 * n) * reduce(lambda a, b: a * b, (d**(2 * j) - 1 for j in range(1, n + 1)), 1); f.append(Z(cl, "Ordnung der Cliffordgruppe bis auf Phasen"))
                f.append(Z(D + 1, "unverzerrte Basen")); f.append(Z(D * (D + 1), "Strahlen in unverzerrten Basen"))
            for kk in (2, 3, 4): f.append(Z(math.comb(D + kk - 1, kk), f"symmetrischer Teil von {kk} Kopien"))
            yield obj(f"qudit:{n}x{d}", "Qubits", f"{n} System(e) mit je {d} Stufen", f)

# ---------- Anyonen und Randphasen ----------
def gen_anyonen():
    for nu in range(1, 17):
        yield obj(f"kitaev:{nu}", "Kitaev Stufe", f"Kitaev Stufe {nu} von sechzehn", [Z(nu, "Stufe"), B(Fr(nu, 2), "zentrale Ladung"), B(Fr(nu, 16), "Gewicht des Wirbels"), Z(4 if nu % 2 == 0 else 3, "Zahl der Sektoren"), T("Z4" if nu % 4 == 2 else ("Z2xZ2" if nu % 4 == 0 else "Ising"), "Verschmelzung")])
    for k in range(1, 13):
        f = [Z(k, "Stufe"), B(Fr(3 * k, k + 2), "zentrale Ladung"), Z(k + 1, "Zahl der Sektoren"), Z(k + 2, "Nenner")] + [B(Fr(j * (j + 2), 4 * (k + 2)), "Gewicht") for j in range(1, k + 1)]
        yield obj(f"su2:{k}", "Anyonen", f"SU(2) auf Stufe {k}", f)
    for n in range(2, 13):
        for k in range(1, 2 * n, 2 if n % 2 == 0 else 1):
            if n % 2 == 1 and k % 2: continue
            hs = {Fr(k * a * a, 2 * n) % 1 for a in range(1, n)}
            yield obj(f"zn:Z{n},k={k}", "Klebegruppe", f"Klebegruppe Z{n} mit Form {k}a²/{2*n}", [Z(n, "Ordnung"), T(f"Z{n}", "Klebegruppe")] + [B(h, "Gewicht") for h in hs if h])

# ---------- flache und runde Kegelflächen ----------
def gen_kegel():
    for sig in itertools.combinations_with_replacement(range(2, 9), 3):
        chi = 2 - sum(1 - Fr(1, m) for m in sig)
        if chi > 0: yield obj(f"kegel:{sig}", "Kegelfläche", f"runde Kugel mit Kegelpunkten {sig}", [Z(2 / chi, "Ordnung der Drehgruppe"), Z(4 / chi, "Ordnung der Doppelgruppe")] + [Z(m, "Kegelordnung") for m in sig])
        if chi == 0: yield obj(f"kegel:{sig}", "Kegelfläche", f"flache Kugel mit Kegelpunkten {sig}", [Z(3, "Kegelpunkte"), Z(max(sig), "höchste Kegelordnung"), T("flach", "Geometrie")] + [Z(m, "Kegelordnung") for m in sig])
    yield obj("kegel:(2,2,2,2)", "Kegelfläche", "flache Kugel mit vier Halbdrehpunkten", [Z(4, "Kegelpunkte"), Z(2, "Kegelordnung"), T("flach", "Geometrie"), Z(1, "freier Parameter der Lage"), T("Z2", "Drehung")])
