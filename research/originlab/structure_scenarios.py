"""Generate structure scenarios without a model and test them against the reviewed Atlas record.

The fingerprint is computed from the Atlas: every integer and fraction that
established findings of several independent problems state. Generated objects
(structure_objects.py) carry computed features. An object, pair or triple is
scored by the rarity-weighted fingerprint it carries and compared with random
fingerprints drawn from values the object pool itself carries, so the search
over the whole pool is part of every p-value. A layer decomposition removes
what the best object carries and repeats on the remainder. A withheld test
selects on one part of the fingerprint and predicts the other part.

Equal numbers are no identification. Every scenario is a declared hypothesis
for atlas_scenarios.assess and for source review; nothing is promoted.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import argparse
import inspect
import json
import math
import random
import re

import numpy as np

import structure_objects
from atlas_scenarios import PROGRAM, knowledge

SCOPE = ("numerical coincidence search over generated objects and reviewed Atlas statements; "
         "a carried number is recorded evidence for a comparison and proves no identification")
# Bookkeeping of the record is no structure: identifiers, dates, decimals, error sizes, counts of
# checks, digits of agreement, moduli of exact arithmetic and line references.
NOISE = re.compile(r"\b\d{4}-\d{2}-\d{2}\b|\b20\d{6}\b|\b[A-Za-z_]+\d+[A-Za-z_\d]*\b|\b[0-9a-f]{8,}\b"
                   r"|\d+(?:[.,]\d+)?e[−\-+]?\d+|\d+[.,]\d+|(?<![\d/])\.\d+|\d+\s*%"
                   r"|\b\d+(?:/\d+)?\s+(?:(?:exakte|gezielte|native|unabhängige|symbolische|numerische|direkte)\s+){0,2}"
                   r"(?:Kontrollen|Vergleiche|Dezimalstellen|Energiekontrollen|Tests)\b"
                   r"|\bmod(?:ulo)?\s+\d+|\bbei\s+65537|:\d+\s*[–-]\s*\d+|\bZeilen?\s*\d+(?:\s*und\s*\d+)?")
FRACTION = re.compile(r"(?<![\d/^])(?<!\d[.,])(\d{1,3})/(\d{1,3})(?![\d/]|[.,]\d)")
INTEGER = re.compile(r"(?<![\d^_/\-])(?<!\d[.,])\b(\d{1,7})\b(?![\d/]|[.,]\d)")
RARE = 12


def stated_numbers(text):
    """Integers and reduced fractions a statement names; identifiers, dates and decimals are dropped."""
    text = NOISE.sub(" ", text)
    fractions = {Fraction(int(a), int(b)) for a, b in FRACTION.findall(text) if int(b)}
    integers = {int(x) for x in INTEGER.findall(FRACTION.sub(" ", text))}
    integers |= {int(f) for f in fractions if f.denominator == 1}
    return ({("int", n) for n in integers if n > 2} |
            {("frac", f"{f.numerator}/{f.denominator}") for f in fractions if f.denominator > 1})


def atlas_fingerprint(kb, min_problems=3):
    """Numbers that established findings of at least min_problems problems state, with their premises."""
    problems, premises = defaultdict(set), defaultdict(Counter)
    for finding in kb["findings"]:
        if finding["status"] != "established_scoped":
            continue
        text = " ".join((finding["title"]["de"], finding["statement"]["de"]))
        for key in stated_numbers(text):
            problems[key].add(finding["problem"])
            premises[key].update(finding["premise_ids"])
    rows = [{"kind": kind, "value": value, "problems": len(problems[(kind, value)]),
             "premise": premises[(kind, value)].most_common(1)[0][0] if premises[(kind, value)] else None}
            for kind, value in problems if len(problems[(kind, value)]) >= min_problems]
    return sorted(rows, key=lambda row: (-row["problems"], row["kind"], str(row["value"])))


def universe(ratios=False, bound=64):
    """All generated objects in a fixed order; optionally with ratios of two of their own integers."""
    objects = {}
    for name, function in inspect.getmembers(structure_objects, inspect.isfunction):
        if name.startswith("gen_"):
            for row in function():
                objects[row["id"]] = {**row, "feats": set(row["feats"])}
    if ratios:
        for row in objects.values():
            integers = sorted({(value, role) for kind, value, role in row["feats"] if kind == "int" and 2 <= value <= 4096})
            first = {}
            for value, role in integers:
                first.setdefault(value, role)
            for a in first:
                for b in first:
                    ratio = Fraction(a, b)
                    if a != b and ratio.denominator > 1 and max(ratio.numerator, ratio.denominator) <= bound:
                        row["feats"].add(("frac", f"{ratio.numerator}/{ratio.denominator}",
                                          f"ratio {first[a]} : {first[b]}"))
    # An object with declared parts also enters as its closure: itself with everything it contains.
    for row in list(objects.values()):
        names, todo = [], list(row.get("parts", ()))
        while todo:   # parts of parts belong to the closure as well
            name = todo.pop(0)
            if name in objects and name not in names and name != row["id"]:
                names.append(name)
                todo.extend(objects[name].get("parts", ()))
        parts = [objects[name] for name in names]
        if parts:
            feats = set(row["feats"]).union(*(part["feats"] for part in parts))
            objects[f"closure:{row['id']}"] = {"id": f"closure:{row['id']}", "fam": "Closure", "feats": feats,
                                               "name": f"{row['name']} with its parts",
                                               "parts": [part["id"] for part in parts]}
    return [objects[name] for name in sorted(objects)]


class Search:
    def __init__(self, objects):
        self.objects, self.n = objects, len(objects)
        self.have = defaultdict(set)
        for index, row in enumerate(objects):
            for kind, value, _ in sorted(row["feats"], key=str):
                self.have[(kind, value)].add(index)
        self.pool = {kind: sorted({key[1] for key in self.have if key[0] == kind}, key=str)
                     for kind in ("int", "frac", "type")}
        # A bridge is a rare feature shared by objects of different families.
        self.linked = np.zeros((self.n, self.n), bool)
        self.bridge = defaultdict(list)
        for key, members in self.have.items():
            if 2 <= len(members) <= RARE and not (key[0] == "int" and key[1] < 7):
                for i in members:
                    for j in members:
                        if i < j and objects[i]["fam"] != objects[j]["fam"]:
                            self.linked[i, j] = self.linked[j, i] = True
                            self.bridge[(i, j)].append(key)

    def weight(self, key):
        return math.log2(self.n/(1.0+len(self.have.get(key, ()))))

    def matrix(self, keys):
        carried = np.zeros((self.n, len(keys)), bool)
        for column, key in enumerate(keys):
            for index in self.have.get(key, ()):
                carried[index, column] = True
        return carried, np.array([self.weight(key) for key in keys])

    def fake(self, keys, rng):
        """A random fingerprint of values the pool carries, of the same kind and similar size."""
        rows, used = [], set()
        for kind, value in keys:
            pool = self.pool[kind]
            near = ([x for x in pool if value/2 <= x <= 2*value and x != value] if kind == "int"
                    else [x for x in pool if x != value]) or pool
            for _ in range(50):
                draw = rng.choice(near)
                if (kind, draw) not in used:
                    break
            used.add((kind, draw))
            rows.append((kind, draw))
        return rows

    def best(self, keys, full=False):
        carried, weights = self.matrix(keys)
        single = carried@weights
        pair = np.where(self.linked, (carried[:, None, :] | carried[None, :, :])@weights, -1.0)
        upper = np.triu_indices(self.n, 1)
        values = pair[upper]
        order = np.argsort(-values, kind="stable")[:150]
        pairs = [(int(upper[0][t]), int(upper[1][t]), float(values[t])) for t in order if values[t] > 0]
        triples = {}
        for i, j, _ in pairs[:60]:
            joined = carried[i] | carried[j]
            candidates = np.where(self.linked[i] | self.linked[j])[0]
            scores = (joined[None, :] | carried[candidates])@weights
            for t in np.argsort(-scores, kind="stable")[:8]:
                k = int(candidates[t])
                if k not in (i, j):
                    triples[tuple(sorted((i, j, k)))] = float(scores[t])
        triples = sorted(triples.items(), key=lambda row: (-row[1], row[0]))
        if not full:
            return float(single.max()), pairs[0][2] if pairs else 0.0, triples[0][1] if triples else 0.0
        return carried, weights, single, pairs, triples


def p_value(null, value):
    return (sum(1 for x in null if x >= value-1e-9)+1)/(len(null)+1)


def layers(search, keys, rng, fakes, level=0.05, limit=8):
    """Remove what the best object carries and repeat while the best object beats the search-wide null."""
    rows, rest = [], list(keys)
    while rest and len(rows) < limit:
        carried, weights = search.matrix(rest)
        single = carried@weights
        index = int(np.argmax(single))
        null = [search.best_single(search.fake(rest, rng)) for _ in range(fakes)]
        p = p_value(null, float(single[index]))
        row = {"object": search.objects[index]["id"], "name": search.objects[index]["name"],
               "score": round(float(single[index]), 2), "share_of_rest": round(float(single[index]/weights.sum()), 3),
               "null_median": round(float(np.median(null)), 2), "p_search_wide": round(p, 4),
               "carries": [value for column, (_, value) in enumerate(rest) if carried[index, column]],
               "significant": p <= level}
        rows.append(row)
        if p > level:
            break
        rest = [key for column, key in enumerate(rest) if not carried[index, column]]
    explained = [row for row in rows if row["significant"]]
    return rows, rest if rows and rows[-1]["significant"] else rest


def best_single(self, keys):
    carried, weights = self.matrix(keys)
    return float((carried@weights).max())


Search.best_single = best_single


def withheld(search, keys, rng, splits, fakes):
    """Select the best object and pair on two thirds of the fingerprint, predict the withheld third."""
    carried, weights = search.matrix(keys)
    pair_all = None
    wins = {"single": [0, 0], "pair": [0, 0]}
    chosen = {"single": Counter(), "pair": Counter()}
    columns = list(range(len(keys)))
    for _ in range(splits):
        rng.shuffle(columns)
        cut = (2*len(columns))//3
        fit, held = sorted(columns[:cut]), sorted(columns[cut:])
        single = carried[:, fit]@weights[fit]
        i = int(np.argmax(single))
        pair = np.where(search.linked, (carried[:, None, fit] | carried[None, :, fit])@weights[fit], -1.0)
        a, b = map(int, np.unravel_index(int(np.argmax(pair)), pair.shape))
        picks = {"single": (i,), "pair": tuple(sorted((a, b)))}
        held_keys = [keys[c] for c in held]
        for name, members in picks.items():
            chosen[name][members] += 1
            real = float((carried[list(members)][:, held].any(0))@weights[held])
            for _ in range(fakes):
                score = sum(search.weight(key) for key in search.fake(held_keys, rng)
                            if search.have.get(key, set()) & set(members))
                wins[name][0] += score >= real-1e-9
                wins[name][1] += 1
    rows = {}
    for name in wins:
        members, count = chosen[name].most_common(1)[0]
        rows[name] = {"p_withheld": round((wins[name][0]+1)/(wins[name][1]+1), 4),
                      "most_selected": [search.objects[m]["id"] for m in members],
                      "selected_share": round(count/splits, 3), "distinct_selections": len(chosen[name])}
    return rows


def margin(carried, weights, rng, draws=2000):
    """Resample the fingerprint: how often does the best object stay ahead of the next one of another kind?"""
    single = carried@weights
    order = [int(i) for i in np.argsort(-single, kind="stable")[:2]]
    if len(order) < 2:
        return None
    columns, wins = list(range(len(weights))), 0
    for _ in range(draws):
        sample = [rng.choice(columns) for _ in columns]
        wins += float(carried[order[0], sample]@weights[sample]) > float(carried[order[1], sample]@weights[sample])
    return {"first": order[0], "second": order[1], "first_ahead_share": round(wins/draws, 3)}


def describe(search, keys, members, carried, weights, score, null, premise):
    mask = carried[list(members)].any(0)
    roles = {search.objects[m]["id"]: sorted({f"{value}: {role}" for column, (kind, value) in enumerate(keys)
                                             if carried[m, column] for k, v, role in search.objects[m]["feats"]
                                             if (k, v) == (kind, value)}) for m in members}
    bridges = sorted({f"{key[1]} ({', '.join(sorted({role for m in members for k, v, role in search.objects[m]['feats'] if (k, v) == key}))})"
                      for (i, j), found in search.bridge.items() if i in members and j in members for key in found})
    return {"objects": [search.objects[m]["id"] for m in members], "names": [search.objects[m]["name"] for m in members],
            "score": round(score, 2), "share": round(score/float(weights.sum()), 3), "p_search_wide": round(p_value(null, score), 4),
            "bridges": bridges, "carries": [keys[c][1] for c in range(len(keys)) if mask[c]],
            "missing": [keys[c][1] for c in range(len(keys)) if not mask[c]], "roles": roles,
            "adopts": sorted({premise[keys[c]] for c in range(len(keys)) if mask[c] and premise.get(keys[c])})}


def proposal(row, revision):
    """A scenario for atlas_scenarios.assess: the objects enter as one declared hypothesis."""
    name = sha256(json.dumps([row["objects"], revision]).encode()).hexdigest()[:12]
    names = " and ".join(row["names"])
    link = "; ".join(row["bridges"][:3]) or "no shared rare feature"
    return {"id": f"structure-{name}",
            "statement": (f"The numbers {', '.join(map(str, row['carries'][:12]))} of the record are features of {names}. "
                          f"Recorded link between the objects: {link}. Search-wide p = {row['p_search_wide']}."),
            "adopts": row["adopts"],
            "new_premises": [{"id": f"hypothesis:structure-{name}",
                              "statement": f"One source realises {names} together.",
                              "rule": "Give a structure-preserving map between the objects that identifies the shared feature, "
                                      "and derive one missing number from it. The hypothesis falls when the shared number "
                                      "counts different things in the two objects and no such map exists.",
                              "predicts": f"missing so far: {', '.join(map(str, row['missing'][:12]))}"}]}


def run(root, fakes=300, splits=60, seed=1, min_problems=3, ratios=False, limit=8):
    root = Path(root)
    if not (root/PROGRAM).is_file():
        return {"status": "ATLAS_UNAVAILABLE", "physical_solution": False}
    kb = knowledge(root)
    rng = random.Random(seed)
    fingerprint = atlas_fingerprint(kb, min_problems)
    search = Search(universe(ratios))
    keys = [(row["kind"], row["value"]) for row in fingerprint]
    premise = {(row["kind"], row["value"]): row["premise"] for row in fingerprint}
    carried, weights, single, pairs, triples = search.best(keys, full=True)
    null = [search.best(search.fake(keys, rng)) for _ in range(fakes)]
    nulls = [[row[k] for row in null] for k in range(3)]
    ranked = lambda members, score, k: describe(search, keys, members, carried, weights, score, nulls[k], premise)
    singles = [ranked((int(i),), float(single[i]), 0) for i in np.argsort(-single, kind="stable")[:limit] if single[i] > 0]
    pair_rows = [ranked((i, j), score, 1) for i, j, score in pairs[:limit]]
    triple_rows = [ranked(members, score, 2) for members, score in triples[:limit]]
    layer_rows, rest = layers(search, keys, rng, fakes)
    significant = [row for row in layer_rows if row["significant"]]
    scenarios = [row for row in (*singles[:2], *pair_rows[:3], *triple_rows[:3]) if row["adopts"]]
    ahead = margin(carried, weights, rng)
    if ahead:
        ahead = {"first": search.objects[ahead["first"]]["id"], "second": search.objects[ahead["second"]]["id"],
                 "first_ahead_share": ahead["first_ahead_share"]}
    return {
        "status": "STRUCTURE_SCENARIOS_RANKED", "physical_solution": False, "scope": SCOPE,
        "atlas_revision": kb["revision"], "seed": seed, "fakes": fakes, "ratios": ratios,
        "pool": {"objects": search.n, "families": dict(sorted(Counter(row["fam"] for row in search.objects).items()))},
        "fingerprint": {"items": len(keys), "min_problems": min_problems, "weight": round(float(weights.sum()), 2),
                        "carried_by_no_object": [keys[c][1] for c in range(len(keys)) if not carried[:, c].any()],
                        "rows": fingerprint},
        "null_median": {"single": round(float(np.median(nulls[0])), 2), "pair": round(float(np.median(nulls[1])), 2),
                        "triple": round(float(np.median(nulls[2])), 2)},
        "singles": singles, "pairs": pair_rows, "triples": triple_rows,
        "margin_over_runner_up": ahead,
        "layers": layer_rows, "significant_layers": len(significant),
        "unexplained_after_layers": [value for _, value in rest],
        "withheld": withheld(search, keys, rng, splits, max(10, fakes//10)),
        "proposals": [proposal(row, kb["revision"]) for row in scenarios],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--fakes", type=int, default=300)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--min-problems", type=int, default=3)
    parser.add_argument("--ratios", action="store_true", help="add ratios of two integers of the same object")
    args = parser.parse_args()
    print(json.dumps(run(args.repo, args.fakes, seed=args.seed, min_problems=args.min_problems, ratios=args.ratios),
                     ensure_ascii=False, indent=1))
