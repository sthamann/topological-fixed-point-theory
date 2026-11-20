"""Binary associativity, its first overlap, and contextual consistency. Python 3."""
from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb


@lru_cache(None)
def trees(leaves):
    if len(leaves) == 1:
        return (leaves[0],)
    return tuple((left, right) for k in range(1, len(leaves))
                 for left in trees(leaves[:k]) for right in trees(leaves[k:]))


def rewrites(tree, permutation, target_shape="right"):
    """Apply ((x y) z) -> p0 (p1 p2), at every possible site."""
    if isinstance(tree, str):
        return
    left, right = tree
    if isinstance(left, tuple):
        pieces = (*left, right)
        x, y, z = (pieces[i] for i in permutation)
        yield (x, (y, z)) if target_shape == "right" else ((x, y), z)
    for target in rewrites(left, permutation, target_shape):
        yield (target, right)
    for target in rewrites(right, permutation, target_shape):
        yield (left, target)


def degree_pairs(vertices, edges):
    indeg = Counter(v for _, v in edges)
    outdeg = Counter(u for u, _ in edges)
    return {v: (indeg[v], outdeg[v]) for v in vertices}


def compact(tree):
    return tree if isinstance(tree, str) else '(' + ''.join(map(compact, tree)) + ')'


def shape_graph(permutation, leaves, seed, target_shape="right"):
    universe = [t for order in permutations(leaves) for t in trees(order)]
    all_edges = {(u, v) for u in universe for v in rewrites(u, permutation, target_shape)}
    neighbors = defaultdict(set)
    for u, v in all_edges:
        neighbors[u].add(v)
        neighbors[v].add(u)
    visited = {seed}
    queue = deque([seed])
    while queue:
        for v in neighbors[queue.popleft()] - visited:
            queue.append(v)
            visited.add(v)
    return visited, {(u, v) for u, v in all_edges if u in visited}


def substitute(tree, values):
    if isinstance(tree, str):
        return values[ord(tree) - ord('a')]
    return tuple(substitute(child, values) for child in tree)


def audit():
    canonical = (0, 1, 2)
    for arity in range(1, 9):
        labels = tuple('abcdefgh'[:arity])
        critical = [t for t in trees(labels) if len(set(rewrites(t, canonical))) > 1]
        print('arity', arity, 'bracketings', len(trees(labels)), 'branching', len(critical))
        if critical:
            break
    assert len(critical) == 1 and arity == 4
    seed = critical[0]
    signature = Counter({(0, 2): 1, (1, 1): 3, (2, 0): 1})
    n_hits = 0
    for p in permutations(range(3)):
        vertices, edges = shape_graph(p, labels, seed)
        counts = Counter(degree_pairs(vertices, edges).values())
        hit = counts == signature and len(edges) == 5
        n_hits += hit
        print('permutation', ''.join(str(i + 1) for i in p), 'component', len(vertices),
              'arrows', len(edges), 'exact_signature', hit)
    print('three-variable right-bracketed control hits', n_hits, '/ 6')
    # Other possible three-variable linear right sides, modulo naming the
    # variables on the left. The unchanged tautology is omitted.
    extra = 0
    for p in permutations(range(3)):
        if p == canonical:
            continue
        verts, arrows = shape_graph(p, labels, seed, 'left')
        counts = Counter(degree_pairs(verts, arrows).values())
        extra += counts == signature and len(arrows) == 5
    print('all different three-variable binary right sides:', n_hits + extra, '/ 11')
    assert extra == 0

    vertices, edges = shape_graph(canonical, labels, seed)
    roles = {v: 'transit' if pair == (1, 1) else 'extremal'
             for v, pair in degree_pairs(vertices, edges).items()}
    assert Counter(roles.values()) == {'transit': 3, 'extremal': 2}
    for v in sorted(vertices, key=compact):
        print(' ', compact(v), roles[v], '->', ','.join(sorted(compact(w) for u,w in edges if u == v)))

    orientation_hits = exact_arm_hits = 0
    for signs in product(range(2), repeat=5):
        cycle_edges = [(i, (i+1) % 5) if s else ((i+1) % 5, i)
                       for i, s in enumerate(signs)]
        degrees = degree_pairs(range(5), cycle_edges)
        if Counter(degrees.values()) == signature:
            orientation_hits += 1
            s, = (v for v, d in degrees.items() if d == (0, 2))
            t, = (v for v, d in degrees.items() if d == (2, 0))
            exact_arm_hits += (s-t) % 5 in (2, 3)
    print('cycle directions: degree split', orientation_hits, '/ 32; and arm lengths 2,3',
          exact_arm_hits, '/ 32')

    # Each arity-five pentagon fixes one interval of length two or four.
    five = list('abcde')
    embeddings = [('first four', lambda v: (v, 'e')),
                  ('last four', lambda v: ('a', substitute(v, five[1:])))]
    for i in range(4):
        expanded = five[:i] + [tuple(five[i:i+2])] + five[i+2:]
        embeddings.append(('fused ' + ''.join(five[i:i+2]),
                           lambda v, data=expanded: substitute(v, data)))
    face_roles = defaultdict(list)
    for name, embed in embeddings:
        for t, role in roles.items():
            face_roles[embed(t)].append((name, role))
    assert set(face_roles) == set(trees(tuple(five)))
    conflicts = {t: rs for t, rs in face_roles.items() if len(set(r for _, r in rs)) > 1}
    print('arity 5: pentagons', len(embeddings), 'vertices', len(face_roles),
          'vertices with inconsistent inherited roles', len(conflicts))
    for t, rs in sorted(conflicts.items(), key=lambda item:compact(item[0])):
        print(' ', compact(t), rs)
    # All local named coordinate 2-planes; changing the readout cannot make
    # context independence work. Each singleton support is counted once.
    local = sorted(vertices, key=compact)
    images = [{t: embed(t) for t in local} for _, embed in embeddings]
    conflicts_by_plane = []
    for two in combinations(local, 2):
        assigned = defaultdict(set)
        for face in images:
            for t, u in face.items():
                assigned[u].add(t in two)
        conflicts_by_plane.append(sum(len(x) > 1 for x in assigned.values()))
    print('all ten local two-sets: conflict counts', sorted(conflicts_by_plane))

    # A sharp compression of an effect forces invariance of the whole embedded
    # face. Thus an off-diagonal entry between unequal face membership masks
    # vanishes. All arithmetic here is exact combinatorics.
    mask = {v: tuple(v in face.values() for face in images) for v in face_roles}
    surviving = [(a, b) for a, b in combinations(local, 2)
                 if all(mask[face[a]] == mask[face[b]] for face in images)]
    print('possible off-diagonal local matrix entries (unordered)', len(surviving), '/ 10')

    images_of = {t: {face[t] for face in images} for t in local}
    connected = {local[0]}
    while True:
        grown = connected | {t for t in local if any(images_of[t] & images_of[u]
                                                     for u in connected)}
        if grown == connected:
            break
        connected = grown
    print('size of forced diagonal equality class', len(connected), '/ 5')
    print('smallest basis-state lift retaining the original tree and its local role',
          sum(len({role for _, role in assignments}) for assignments in face_roles.values()))
    # For arbitrary arity, every embedded copy through a tree picks three
    # connected internal vertices. Below we shrink each of the four hanging
    # pieces to a leaf and recover the original four-leaf bracketing.
    def cuts(tree, n):
        if n == 1:
            yield '*'
        elif not isinstance(tree, str):
            for i in range(1, n):
                for left in cuts(tree[0], i):
                    for right in cuts(tree[1], n-i):
                        yield (left, right)
    def forget_labels(tree):
        return '*' if isinstance(tree, str) else tuple(map(forget_labels, tree))
    shape_roles = {forget_labels(t): role for t, role in roles.items()}
    @lru_cache(None)
    def available_roles(tree):
        if isinstance(tree, str):
            return frozenset()
        local_roles = {shape_roles[s] for s in cuts(tree, 4)}
        return frozenset(local_roles | available_roles(tree[0]) | available_roles(tree[1]))
    for n in range(4, 10):
        all_trees = trees(tuple('abcdefghi'[:n]))
        multiplicities = [len(available_roles(t)) for t in all_trees]
        catalan = comb(2*(n-1), n-1) // n
        assert min(multiplicities) >= 1 and len(all_trees) == catalan
        assert sum(multiplicities) == 2*catalan-(n+1)
        assert Counter(next(iter(available_roles(t))) for t in all_trees
                       if len(available_roles(t)) == 1) == {'extremal': 2, 'transit': n-1}
        print('arity', n, 'trees', len(all_trees), 'mixed', multiplicities.count(2),
              'minimal tree/role states', sum(multiplicities))
    assert not surviving and len(connected) == 5
    assert min(conflicts_by_plane) == 7 and max(conflicts_by_plane) == 9
    assert len(conflicts) == 8 and n_hits == 1
    assert (orientation_hits, exact_arm_hits) == (20, 10)


if __name__ == '__main__':
    audit()
