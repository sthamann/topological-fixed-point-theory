"""Reproduce the minimal finite reversible process and operation split, using only integers."""
from collections import Counter
from itertools import permutations, product


def reduced_latin_squares(n):
    """Enumerate all reduced Latin squares. Identity has label 0."""
    labels = tuple(range(n))
    possible_rows = [
        [(first,) + rest for rest in permutations(x for x in labels if x != first)]
        for first in labels
    ]
    used = [{x} for x in labels]
    rows = [labels]

    def extend():
        if len(rows) == n:
            yield tuple(rows)
            return
        for row in possible_rows[len(rows)]:
            if any(value in column for value, column in zip(row, used)):
                continue
            rows.append(row)
            for column, value in zip(used, row):
                column.add(value)
            yield from extend()
            for column, value in zip(used, row):
                column.remove(value)
            rows.pop()

    yield from extend()


def canonical_table(t):
    labels = range(len(t))
    relabelings = ((0,) + rest for rest in permutations(labels[1:]))

    def relabel(p):
        inv = tuple(p.index(i) for i in labels)
        return tuple(inv[t[p[i]][p[j]]] for i in labels for j in labels)

    return min(map(relabel, relabelings))


def classify(n):
    latin = 0
    labeled = 0
    groups = {}
    labels = range(n)
    for t in reduced_latin_squares(n):
        latin += 1
        if not all(t[t[a][b]][c] == t[a][t[b][c]] for a, b, c in product(labels, repeat=3)):
            continue
        labeled += 1
        groups[canonical_table(t)] = t
    return latin, labeled, list(groups.values())


def audit():
    all_models = []
    print("order  reduced Latin squares  labelled groups  group types  noncommuting types")
    for n in range(1, 7):
        latin, labelled, groups = classify(n)
        nonabelian = [t for t in groups if any(t[i][j] != t[j][i] for i,j in product(range(n),repeat=2))]
        all_models.extend(groups)
        print(n, latin, labelled, len(groups), len(nonabelian))
    assert len(nonabelian) == 1
    t, = nonabelian
    ops = range(1, len(t))
    order_two = {g for g in ops if t[g][g] == 0}
    assert (len(order_two), len(ops) - len(order_two)) == (3, 2)

    aut = []
    for tail in permutations(ops):
        p = (0,) + tail
        if all(p[t[i][j]] == t[p[i]][p[j]] for i,j in product(range(len(t)),repeat=2)):
            aut.append(p)
    orbits = {frozenset(p[g] for p in aut) for g in ops}
    assert orbits == {frozenset(order_two), frozenset(set(ops)-order_two)}
    print("aut size / nonidentity orbit sizes:", len(aut), sorted(map(len,orbits),reverse=True))

    def hits(table):
        size = len(table)
        return size == 6 and sum(table[i][i] == 0 for i in range(1,size)) == 3

    print("3 and 2, order 6:", sum(hits(table) for table in groups), "/", len(groups))
    print("3 and 2, all types through 6:", sum(map(hits,all_models)), "/", len(all_models))

    # Full syntactic family with one product, one equality and three atom positions.
    realized = Counter()
    expressions = {}
    for u,v,w in product("xe", repeat=3):
        def atom(letter,g):
            return g if letter == "x" else 0
        subset = frozenset(g for g in ops if t[atom(u,g)][atom(v,g)] == atom(w,g))
        expressions[u+v+"="+w] = subset
        realized[len(subset)] += 1
    print("true-set sizes among eight one-product equalities:", dict(sorted(realized.items())))
    print("predicates selecting exactly the involutions:", [k for k,v in expressions.items() if v==order_two])
    print("distinct semantic predicates:", len(set(expressions.values())))
    assert realized == {0:4, 3:1, 5:3}


if __name__ == "__main__":
    audit()
