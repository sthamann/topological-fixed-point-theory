"""KI1 (reduced) - planning with interventions in a switch world.

World: 16 cells with symbols 0..3; active four-cell rules (five families, combined) must sum
to zero in the group.  The visible cells are SWITCHES (a nearly minimal set that fixes the
world), the hidden cells are LAMPS.  An action sets one switch to another symbol; the world
answers with the new consistent lamp state, or rejects the action when no consistent state
exists.  Goal: make one lamp show a given symbol.

Every model is used as a world model inside the SAME planner: for each of 8 candidate actions
the model predicts the lamps, the planner takes the action with the highest predicted
probability of the goal.  Success is judged by the exact world (brute-force enumeration).
A two-step variant has 64 candidate action pairs and two goals at once.

Models (same inputs, same targets, same training budget):
  transformer   global attention, told which rule families are active
  plain         relational net with learned local updates
  gen_same      the same net with the ordinary exact factor message of the TRUE group rule
  gen_other     the same net with the factor message of the other group (misinformed)
  tfpt          the same net with the quartic-code projector message
  rand          the same net with a random rank-5 projector
Groups: new worlds with trained rule combinations, new rule combinations, two-step plans.
Not built from the package: the warehouse world with keys and doors, graphs of 24-48 nodes,
horizons 16-64, calibrated sets of possible states.
"""
import json, sys, collections, itertools
import numpy as np, torch
import worlds as W, models as M

NTRAIN, NTEST, STEPS, SEEDS, GROUPS, HMAX = 1000, 1000, 2500, [0, 1, 2], ["xor", "z4"], 7
MODELS = ["transformer", "plain", "gen_same", "gen_other", "tfpt", "rand"]
if "--quick" in sys.argv:
    NTRAIN, NTEST, STEPS, SEEDS = 200, 100, 150, [0]
CAND = np.array(list(itertools.product(range(4), repeat=HMAX)), dtype=np.int64)


def switch_puzzle(combo, group, rng):
    quads = W.FAM[list(combo)].reshape(-1, 4)
    w = W.sample_world(combo, group, rng)
    masked = []
    for c in rng.permutation(16):
        if len(masked) == HMAX:
            break
        if W.n_solutions(w, masked + [int(c)], quads, group) == 1:
            masked.append(int(c))
    known = np.ones(16, bool); known[masked] = False
    return w, known, quads, np.array(masked)


def completions(w, masked, quads, group):
    """all hidden assignments and their rule sums for the current switches."""
    cand = CAND[: 4 ** len(masked), HMAX - len(masked):] if len(masked) < HMAX else CAND
    cand = np.unique(cand, axis=0)
    Wm = np.repeat(w[None], len(cand), 0); Wm[:, masked] = cand
    return Wm, W.gsum(Wm[:, quads], group)


def outcome(Wm, S, quads, group, w, changes):
    """exact answer of the world to setting switches; None if rejected or ambiguous."""
    S2 = S.copy()
    for c, u in changes:
        inq = (quads == c).any(1)
        S2[:, inq] = (S2[:, inq] ^ (w[c] ^ u)) if group == "xor" else (S2[:, inq] + (u - w[c])) % 4
    ok = np.flatnonzero((S2 == 0).all(1))
    if len(ok) != 1:
        return None
    out = Wm[ok[0]].copy()
    for c, u in changes:
        out[c] = u
    return out


def make_set(n, combos, group, seed, plan=True):
    rng = np.random.default_rng(seed)
    D = collections.defaultdict(list)
    while len(D["x"]) < n:
        combo = combos[rng.integers(len(combos))]
        w, known, quads, masked = switch_puzzle(combo, group, rng)
        rule = np.zeros(5, bool); rule[list(combo)] = True
        if not plan:
            D["x"].append(w); D["known"].append(known); D["rule"].append(rule); continue
        Wm, S = completions(w, masked, quads, group)
        acts = [(int(c), u) for c in np.flatnonzero(known) for u in range(4) if u != w[c]]
        res = {a: outcome(Wm, S, quads, group, w, [a]) for a in acts}
        good = [(a, t) for a, r in res.items() if r is not None for t in masked if r[t] != w[t]]
        if not good:
            continue
        a_star, t = good[rng.integers(len(good))]
        v = res[a_star][t]
        others = [a for a in acts if a != a_star]
        cand = [a_star] + [others[i] for i in rng.permutation(len(others))[:7]]
        cand = [cand[i] for i in rng.permutation(8)]
        cx = np.repeat(w[None], 8, 0)
        for i, (c, u) in enumerate(cand):
            cx[i, c] = u
        succ = np.array([res[a] is not None and res[a][t] == v for a in cand])
        # two-step: pairs of the 8 candidates on different switches, two goals
        pairs = [(a, b) for a in cand for b in cand]
        res2 = [outcome(Wm, S, quads, group, w, [a, b]) if a[0] != b[0] else None for a, b in pairs]
        ok2 = [i for i, r in enumerate(res2) if r is not None and (r[masked] != w[masked]).sum() >= 2]
        if not ok2:
            continue
        i2 = ok2[rng.integers(len(ok2))]
        ch = masked[res2[i2][masked] != w[masked]]
        t2 = rng.choice(ch, 2, replace=False)
        v2 = res2[i2][t2]
        cx2 = np.repeat(w[None], 64, 0)
        for i, (a, b) in enumerate(pairs):
            cx2[i, a[0]] = a[1]; cx2[i, b[0]] = b[1]
        succ2 = np.array([r is not None and (r[t2] == v2).all() for r in res2])
        for k, val in (("x", w), ("known", known), ("rule", rule), ("cx", cx), ("t", t), ("v", v), ("succ", succ),
                       ("cx2", cx2), ("t2", t2), ("v2", v2), ("succ2", succ2),
                       ("admissible", np.array([res[a] is not None for a in cand]))):
            D[k].append(val)
    out = {k: np.array(v) for k, v in D.items()}
    out["active"] = np.repeat(out["rule"], 4, axis=1)
    return out


def data(group):
    mk = lambda name, n, combos, seed, plan: W.cached(f"ki1p_{group}_{name}_{n}", lambda: make_set(n, combos, group, seed, plan))
    return dict(train=mk("train", NTRAIN, W.TRAIN_COMBOS, 1, False),
                new_worlds=mk("id", NTEST, W.TRAIN_COMBOS, 2, True),
                new_rule_combinations=mk("ood", NTEST, W.OOD_COMBOS, 3, True))


@torch.no_grad()
def logp(net, kw_fn, x, known, d_idx, quads):
    out = []
    for s in range(0, len(x), 512):
        o = net(torch.as_tensor(x[s:s + 512]), torch.as_tensor(known[s:s + 512]), quads, **kw_fn(d_idx[s:s + 512]))[-1]
        out.append(torch.log_softmax(o, -1))
    return torch.cat(out).numpy()


def job(a):
    model, group, seed = a
    torch.set_num_threads(2); torch.manual_seed(seed)
    d = data(group)
    quads = torch.as_tensor(W.QUADS16)
    other = "z4" if group == "xor" else "xor"
    if model == "transformer":
        net = M.Transformer()
        M.train(net, d["train"], quads, 2 * STEPS, lr=5e-4, seed=seed, use_rule=True)
    else:
        mode = {"gen_same": f"gen_{group}", "gen_other": f"gen_{other}"}.get(model, model)
        net = M.RelNet(mode, seed=seed)
        M.train(net, d["train"], quads, STEPS, seed=seed, use_active=True)
    net.eval()
    out = dict(model=model, group=group, seed=seed)
    for split in ("new_worlds", "new_rule_combinations"):
        t = d[split]; n = len(t["x"])
        kw = (lambda i: dict(rule=torch.as_tensor(t["rule"][i]))) if model == "transformer" else \
             (lambda i: dict(active=torch.as_tensor(t["active"][i]), steps=8))
        # one step: 8 candidates
        idx = np.repeat(np.arange(n), 8)
        lp = logp(net, kw, t["cx"].reshape(n * 8, 16), t["known"][idx], idx, quads).reshape(n, 8, 16, 4)
        score = lp[np.arange(n)[:, None], np.arange(8)[None], t["t"][:, None], t["v"][:, None]]
        pick = score.argmax(1)
        out[split] = t["succ"][np.arange(n), pick].tolist()
        # rule violations in the predicted states of admissible candidates
        pred = np.where(t["known"][idx], t["cx"].reshape(n * 8, 16), lp.reshape(n * 8, 16, 4).argmax(-1))
        viol = (W.gsum(pred[:, W.QUADS16], group) != 0) & t["active"][idx]
        out[split + "_violation_rate"] = float(viol.any(1)[t["admissible"].reshape(-1)].mean())
        # two steps: 64 candidate pairs, two goals
        idx = np.repeat(np.arange(n), 64)
        lp = logp(net, kw, t["cx2"].reshape(n * 64, 16), t["known"][idx], idx, quads).reshape(n, 64, 16, 4)
        ar = np.arange(n)[:, None]
        score = sum(lp[ar, np.arange(64)[None], t["t2"][:, k][:, None], t["v2"][:, k][:, None]] for k in (0, 1))
        out[split + "_two_step"] = t["succ2"][np.arange(n), score.argmax(1)].tolist()
    return out


if __name__ == "__main__":
    for g in GROUPS:
        d = data(g)
    res = M.run_parallel(job, [(m, g, s) for g in GROUPS for m in MODELS for s in SEEDS], nproc=8)
    rng = np.random.default_rng(0)
    table = {}
    for g in GROUPS:
        d = data(g)
        for split in ("new_worlds", "new_rule_combinations"):
            for suffix, sk in (("", "succ"), ("_two_step", "succ2")):
                key = split + suffix
                per = {m: np.mean([r[key] for r in res if r["model"] == m and r["group"] == g], 0) for m in MODELS}
                row = {m: round(float(v.mean()), 4) for m, v in per.items()}
                row["random choice"] = round(float(d[split][sk].mean()), 4)
                row["exact search (upper bound)"] = round(float(d[split][sk].any(1).mean()), 4)
                for a, b in (("tfpt", "gen_same"), ("tfpt", "plain"), ("gen_same", "plain"), ("tfpt", "rand")):
                    diff = per[a] - per[b]
                    bs = [diff[rng.integers(len(diff), size=len(diff))].mean() for _ in range(2000)]
                    row[f"{a} minus {b} (paired, 95% interval)"] = [round(float(diff.mean()), 4),
                                                                   round(float(np.quantile(bs, .025)), 4), round(float(np.quantile(bs, .975)), 4)]
                table[f"{g} | {key}"] = row
            table[f"{g} | {split} | rule violations in predicted states"] = {
                m: round(float(np.mean([r[split + "_violation_rate"] for r in res if r["model"] == m and r["group"] == g])), 4) for m in MODELS}
        table[f"{g} | share of candidate actions the world accepts"] = round(float(d["new_worlds"]["admissible"].mean()), 3)
        table[f"{g} | hidden lamps per world (mean)"] = round(float((~d["new_worlds"]["known"]).sum(1).mean()), 2)
    json.dump(dict(setup=dict(train_worlds=NTRAIN, test_tasks_per_group=NTEST, seeds=len(SEEDS), steps=STEPS), table=table),
              open("out/ki1_planning.json", "w"), indent=1)
    print(json.dumps(table, indent=1))
