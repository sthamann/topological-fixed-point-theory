"""Szenariosuche S42. Sucht nicht über handgeschriebene Hypothesen, sondern über erzeugte Objekte.
Aufruf: python3 -B suche.py [--fakes 400] [--seed 1]
Liest fingerprint.json, sammelt alle gen_*() aus generators.py und gen_*.py, schreibt results.json und ERGEBNIS.md."""
import json, sys, math, glob, inspect, importlib, os, random, re, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, HERE)
ARG = lambda k, d: type(d)(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
NF, SEED, RARE = ARG("--fakes", 400), ARG("--seed", 1), 12

def universe():
    objs = {}
    for path in ["generators.py"] + sorted(os.path.basename(p) for p in glob.glob(os.path.join(HERE, "gen_*.py"))):
        mod = importlib.import_module(path[:-3])
        for n, f in inspect.getmembers(mod, inspect.isfunction):
            if n.startswith("gen_") and f.__module__ == mod.__name__:
                for o in f(): objs[o["id"]] = o
    return [objs[k] for k in sorted(objs)]

def key(a, w): return (a, int(w) if a == "zahl" else str(w))

class Suche:
    def __init__(self, objs):
        self.objs = objs; self.N = len(objs); self.have = collections.defaultdict(set)
        for i, o in enumerate(objs):
            for a, w, _ in sorted(o["feats"], key=str): self.have[(a, w)].add(i)
        self.pool = {a: sorted({k[1] for k in self.have if k[0] == a}, key=str) for a in ("zahl", "bruch", "typ")}
        # Brücken: seltene gemeinsame Merkmale, unabhängig vom Fingerabdruck
        self.C = np.zeros((self.N, self.N), bool); self.bridge = collections.defaultdict(list)
        for k, s in self.have.items():
            if 2 <= len(s) <= RARE and not (k[0] == "zahl" and k[1] < 7):
                for i in s:
                    for j in s:
                        if i < j and objs[i]["fam"] != objs[j]["fam"]: self.C[i, j] = self.C[j, i] = True; self.bridge[(i, j)].append(k)
    def matrix(self, keys):
        M = np.zeros((self.N, len(keys)), bool)
        for c, k in enumerate(keys):
            for i in self.have.get(k, ()): M[i, c] = True
        n = M.sum(0); w = np.log2(self.N / (1.0 + n)); return M, w
    def best(self, keys, full=False):
        M, w = self.matrix(keys); s1 = M @ w
        P = (M[:, None, :] | M[None, :, :]) @ w; P = np.where(self.C, P, -1.0); iu = np.triu_indices(self.N, 1); pv = P[iu]
        order = np.argsort(-pv)[:150]; pairs = [(int(iu[0][t]), int(iu[1][t]), float(pv[t])) for t in order if pv[t] > 0]
        tri = {}
        for i, j, _ in pairs[:60]:
            m = M[i] | M[j]; cand = np.where(self.C[i] | self.C[j])[0]
            sc = (m[None, :] | M[cand]) @ w
            for t in np.argsort(-sc)[:8]:
                k = int(cand[t])
                if k not in (i, j): tri[tuple(sorted((i, j, k)))] = float(sc[t])
        tr = sorted(tri.items(), key=lambda x: -x[1])
        if not full: return float(s1.max()), (pairs[0][2] if pairs else 0.0), (tr[0][1] if tr else 0.0), s1
        return M, w, s1, pairs, tr
    def fake(self, items, rng):
        out, used = [], set()
        for a, v in items:
            pool = self.pool[a]
            if a == "zahl": near = [x for x in pool if v / 2 <= x <= 2 * v and x != v] or pool
            else: near = [x for x in pool if x != v]
            for _ in range(50):
                x = rng.choice(near)
                if (a, x) not in used: break
            used.add((a, x)); out.append((a, x))
        return out

def run(S, items, names, label, rng):
    keys = [key(a, w) for a, w in items]; M, w, s1, pairs, tr = S.best(keys, full=True)
    null = [S.best(S.fake(keys, rng)) for _ in range(NF)]
    n1 = np.array([x[0] for x in null]); n2 = np.array([x[1] for x in null]); n3 = np.array([x[2] for x in null]); nobj = np.array([x[3] for x in null])
    p = lambda arr, v: float(((arr >= v - 1e-9).sum() + 1) / (len(arr) + 1))
    total = float(w.sum()); cov = lambda m: [names[c] for c in range(len(keys)) if m[c]]; mis = lambda m: [names[c] for c in range(len(keys)) if not m[c]]
    why = lambda i, m=None: sorted({f"{names[c]}: {wf}" for c in range(len(keys)) if M[i, c] for a, v, wf in S.objs[i]["feats"] if (a, v) == keys[c]})
    singles = []
    for i in np.argsort(-s1)[:25]:
        i = int(i)
        if s1[i] <= 0: break
        singles.append({"id": S.objs[i]["id"], "name": S.objs[i]["name"], "fam": S.objs[i]["fam"], "punkte": round(float(s1[i]), 2), "anteil": round(float(s1[i]) / total, 3), "p_einzeln": p(nobj[:, i], s1[i]), "p_gesamt": p(n1, s1[i]), "traegt": why(i), "merkmale": len({(a, v) for a, v, _ in S.objs[i]["feats"]})})
    def scen(idx, sc, narr):
        m = np.zeros(len(keys), bool)
        for i in idx: m |= M[i]
        br = sorted({f"{k[1]}" + (" (" + ", ".join(sorted({wf for i in idx for a, v, wf in S.objs[i]['feats'] if (a, v) == k})) + ")") for (a, b), ks in S.bridge.items() if a in idx and b in idx for k in ks})
        return {"objekte": [S.objs[i]["id"] for i in idx], "namen": [S.objs[i]["name"] for i in idx], "familien": sorted({S.objs[i]["fam"] for i in idx}), "punkte": round(sc, 2), "anteil": round(sc / total, 3), "p_gesamt": p(narr, sc), "bruecken": br, "traegt": cov(m), "fehlt": mis(m), "beitrag": {S.objs[i]["id"]: why(i) for i in idx}}
    never = [names[c] for c in range(len(keys)) if not M[:, c].any()]
    return {"lauf": label, "merkmale": len(keys), "gesamtgewicht": round(total, 2), "fakes": NF, "null_median": {"einzeln": round(float(np.median(n1)), 2), "paar": round(float(np.median(n2)), 2), "tripel": round(float(np.median(n3)), 2)},
            "gewichte": {names[c]: round(float(w[c]), 2) for c in range(len(keys))}, "niemand_traegt": never, "einzeln": singles, "paare": [scen((i, j), sc, n2) for i, j, sc in pairs[:15]], "tripel": [scen(idx, sc, n3) for idx, sc in tr[:15]]}

def known_text():
    txt = ""
    for pat in ("tfpt_landkarte/content-*.js", "docs/research/fundamental-plan/RESULTS.md", "docs/research/fundamental-plan/SEARCH-METHODS.md"):
        for f in glob.glob(os.path.join(ROOT, pat)): txt += open(f, encoding="utf-8", errors="ignore").read()
    return re.sub(r"\s+", " ", txt)
def tokens(o):
    s = o["id"].split(":", 1)[1]; t = [s.split("<")[0].strip(), s.split(",")[0]]
    m = re.search(r"\(([^)]*)\)", o["name"]); t += [m.group(1)] if m and len(m.group(1)) > 3 else []
    return [x for x in dict.fromkeys(t) if len(x) > 2]

def satz(sc, S, byid):
    a = " und ".join(sc["namen"]); br = "; ".join(sc["bruecken"][:3]) or "kein gemeinsames seltenes Merkmal"
    return {"annahme": f"Die Bausteine der Theorie sind die Merkmale von {a}. Die Objekte hängen über {br} zusammen.",
            "pruefung": "Finde eine Abbildung zwischen den Objekten, die das gemeinsame Merkmal wirklich identifiziert (nicht nur dieselbe Zahl), und leite daraus eines der fehlenden Merkmale her: " + (", ".join(sc["fehlt"][:6]) or "keines fehlt") + ".",
            "faellt": "Fällt, wenn die gemeinsame Zahl in beiden Objekten etwas Verschiedenes zählt und keine strukturerhaltende Abbildung existiert, oder wenn die Abbildung ein getragenes Merkmal verschiebt."}

if __name__ == "__main__":
    rng = random.Random(SEED); fp = json.load(open(os.path.join(HERE, "fingerprint.json"), encoding="utf-8")); objs = universe(); S = Suche(objs); byid = {o["id"]: o for o in objs}
    items = [(x["art"], x["wert"]) for x in fp["items"]]; names = [f'{x["wert"]} ({x["name"]})' for x in fp["items"]]
    out = {"objekte": len(objs), "familien": dict(collections.Counter(o["fam"] for o in objs)), "laeufe": [run(S, items, names, "alles", rng)]}
    anchor = byid.get(fp.get("anker", ""))
    if anchor:
        ak = {(a, v) for a, v, _ in anchor["feats"]}; rest = [(it, nm) for it, nm in zip(items, names) if key(*it) not in ak]
        out["anker"] = {"id": anchor["id"], "traegt": [nm for it, nm in zip(items, names) if key(*it) in ak]}
        out["laeufe"].append(run(S, [r[0] for r in rest], [r[1] for r in rest], "Rest ohne Anker", rng))
    txt = known_text()
    for L in out["laeufe"]:
        for e in L["einzeln"]: e["im_text"] = [t for t in tokens(byid[e["id"]]) if t in txt]
        for sc in L["paare"] + L["tripel"]:
            sc["im_text"] = {i: [t for t in tokens(byid[i]) if t in txt] for i in sc["objekte"]}; sc["neu"] = [i for i in sc["objekte"] if not sc["im_text"][i]]; sc["hypothese"] = satz(sc, S, byid)
    json.dump(out, open(os.path.join(HERE, "results.json"), "w"), ensure_ascii=False, indent=1)
    md = [f"# Szenariosuche, automatisch erzeugt\n\n{len(objs)} Objekte aus {len(out['familien'])} Familien, {NF} Zufallsfingerabdrücke je Lauf, Startwert {SEED}.\n"]
    for L in out["laeufe"]:
        md += [f"\n## Lauf: {L['lauf']} ({L['merkmale']} Merkmale, Gewicht {L['gesamtgewicht']})\n", f"Zufall erreicht im Mittel: einzeln {L['null_median']['einzeln']}, Paar {L['null_median']['paar']}, Tripel {L['null_median']['tripel']}.\n", "Kein Objekt trägt: " + (", ".join(L["niemand_traegt"]) or "nichts") + "\n", "\n### Einzelne Objekte\n\n| Objekt | Punkte | Anteil | p einzeln | p gesamt | im Text |\n|---|---|---|---|---|---|"]
        md += [f"| {e['name']} | {e['punkte']} | {e['anteil']:.0%} | {e['p_einzeln']:.3f} | {e['p_gesamt']:.3f} | {'ja' if e['im_text'] else 'nein'} |" for e in L["einzeln"][:15]]
        for kind in ("paare", "tripel"):
            md += [f"\n### Beste {kind.capitalize()}\n\n| Objekte | Punkte | Anteil | p gesamt | Brücke | neu |\n|---|---|---|---|---|---|"]
            md += [f"| {' + '.join(sc['namen'])} | {sc['punkte']} | {sc['anteil']:.0%} | {sc['p_gesamt']:.3f} | {'; '.join(sc['bruecken'][:2])} | {len(sc['neu'])} |" for sc in L[kind][:10]]
    open(os.path.join(HERE, "ERGEBNIS.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")
    for L in out["laeufe"]:
        print("==", L["lauf"], L["null_median"], "| niemand:", L["niemand_traegt"])
        for e in L["einzeln"][:10]: print(f"  {e['punkte']:6.1f} p1={e['p_einzeln']:.3f} pg={e['p_gesamt']:.3f} {e['name']}")
        for sc in L["paare"][:6]: print(f"  P {sc['punkte']:6.1f} pg={sc['p_gesamt']:.3f} {' + '.join(sc['namen'])} | {sc['bruecken'][:2]}")
        for sc in L["tripel"][:6]: print(f"  T {sc['punkte']:6.1f} pg={sc['p_gesamt']:.3f} {' + '.join(sc['namen'])}")
