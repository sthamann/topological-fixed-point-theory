"""Exact two-clock comparison in the original Gaussian E8 quartic coordinates."""
from collections import Counter
from fractions import Fraction
from itertools import product, permutations


def compose(left, right):
    return tuple(left[i] for i in right)


def inverse(perm):
    return tuple(perm.index(i) for i in range(len(perm)))


def conjugate(left, right):
    return compose(compose(left, right), inverse(left))


def perm_power(value, exponent):
    result = tuple(range(len(value)))
    for _ in range(exponent):
        result = compose(result, value)
    return result


def commutator(left, right):
    return compose(conjugate(left, right), inverse(right))


def closure(generators):
    seen = {tuple(range(len(generators[0])))}
    queue = list(seen)
    for element in queue:
        for generator in generators:
            candidate = compose(element, generator)
            if candidate not in seen:
                seen.add(candidate)
                queue.append(candidate)
    return seen


identity6 = tuple(range(6))
quarter = (1, 2, 3, 0, 4, 5)
swap = (0, 1, 2, 3, 5, 4)
groups = {frozenset(perm_power(p, k) for k in range(5))
          for p in permutations(identity6) if p != identity6 and perm_power(p, 5) == identity6}
normalized = [h for h in groups if {conjugate(quarter, p) for p in h} == h]
good = sorted((h for h in normalized if all(conjugate(quarter, p) == perm_power(p, 2)
                                            for p in h)), key=min)
assert (len(groups), len(normalized), len(good)) == (36, 4, 2)
a = min(good[0] - {identity6})
b = conjugate(swap, a)
assert b in good[1]
assert [i for i in identity6 if a[i] == i] == [4]
assert [i for i in identity6 if b[i] == i] == [5]
assert len(closure((a, b))) == 360
commutators = [commutator(perm_power(a, i), perm_power(b, j)) for i, j in product(range(1, 5), repeat=2)]
assert all(c[4] != 4 and c[5] != 5 for c in commutators)
assert sum(perm_power(c, 2) == identity6 for c in commutators) == 8
assert sum(perm_power(c, 5) == identity6 for c in commutators) == 8
print("S6: 36 subgroups of order five, 4 normalized by G, 2 squared by G; fixed labels 4 and 5.")
print("The two subgroups generate A6 (order 360); all 16 generator commutators move both fixed labels.")
print("Eight commutators have order 2 and eight have order 5.")

# RM(1,3) and the exact quartic vectors of the original source.
points = list(product((0, 1), repeat=3))
supports = [tuple(i for i, x in enumerate(points) if (sum(a*b for a,b in zip(linear, x))+offset)%2)
            for linear in points if any(linear) for offset in (0, 1)]
roots8 = []
for i, sign in product(range(8), (-1, 1)):
    roots8.append(tuple(2*sign*(j == i) for j in range(8)))
for support in supports:
    for signs in product((-1, 1), repeat=4):
        root = [0]*8
        for j, sign in zip(support, signs):
            root[j] = sign
        roots8.append(tuple(root))
roots = {tuple(complex(x, y) for x, y in zip(r[::2], r[1::2])) for r in roots8}
rays = {min((tuple(z*unit for z in r) for unit in (1, -1, 1j, -1j)),
            key=lambda r: tuple((z.real, z.imag) for z in r)) for r in roots}
assert len(roots) == 240 and len(rays) == 60

W = ((2,2,-1,-1,-1,-1),(0,0,-1,-1,1,1),(0,0,-1,1,-1,1),(0,0,1,-1,-1,1),(1,-1,0,0,0,0))
gram = (4,12,12,12,24)
assert all(sum(d*row[i]*row[j] for d,row in zip(gram,W)) == (40 if i == j else -8)
           for i,j in product(range(6),repeat=2))


def polynomial_values(z):
    t = (sum(x**4 for x in z),
         6*(z[0]**2*z[1]**2+z[2]**2*z[3]**2),
         6*(z[0]**2*z[2]**2+z[1]**2*z[3]**2),
         6*(z[0]**2*z[3]**2+z[1]**2*z[2]**2),
         24*z[0]*z[1]*z[2]*z[3])
    return tuple(sum(x*w[j] for x,w in zip(t,W)) for j in range(6))


def apply(matrix, vector):
    return tuple(sum(a*b for a,b in zip(row, vector)) for row in matrix)


def matmul(left, right):
    return tuple(tuple(sum(a*b for a,b in zip(row, col)) for col in zip(*right)) for row in left)


def adjoint(matrix):
    return tuple(tuple(z.conjugate() for z in row) for row in zip(*matrix))


identity4 = tuple(tuple(complex(i == j) for j in range(4)) for i in range(4))


def mat_power(matrix, exponent):
    result = identity4
    for _ in range(exponent):
        result = matmul(result, matrix)
    return result


def mat_conj(left, right):
    return matmul(matmul(left, right), adjoint(left))


def action(matrix, probes):
    probe = (1, 2, 5, 11)
    originals = polynomial_values(probe)
    assert len(set(originals)) == 6
    transform = tuple(zip(*matrix))
    p = tuple(originals.index(value) for value in polynomial_values(apply(transform, probe)))
    for probe in probes:
        old = polynomial_values(probe)
        assert polynomial_values(apply(transform, probe)) == tuple(old[i] for i in p)
    return p


# Newton interpolation on the 70 integer nodes with sum <= 4 determines a quartic.
# All evaluations on these dyadic matrices are exact in binary64: integer
# numerators are bounded by 2**53 and the denominator is a power of two.
# Each matrix entry is in (1/2) Z[i]: each native isometry sends 2 e_j to a root.
# This property is preserved under all products of the verified root isometries.
probes = [z for z in product(range(5), repeat=4) if sum(z) <= 4]
assert len(probes) == 70
transpositions = {}
for ray in sorted(rays, key=lambda r: tuple((z.real, z.imag) for z in r)):
    reflection = tuple(tuple(identity4[i][j] - ray[i]*ray[j].conjugate()/2
                             for j in range(4)) for i in range(4))
    assert {apply(reflection, root) for root in roots} == roots
    p = action(reflection, probes)
    assert sum(i != j for i,j in enumerate(p)) == 2 and perm_power(p,2)==identity6
    transpositions.setdefault(p, reflection)
assert len(transpositions) == 15


def lift(perm):
    lifts = {identity6: identity4}
    queue = [identity6]
    for start in queue:
        for p, reflection in transpositions.items():
            nxt = compose(start, p)
            if nxt not in lifts:
                lifts[nxt] = matmul(lifts[start], reflection)
                if nxt == perm:
                    result = lifts[nxt]
                    assert action(result, probes) == perm
                    assert {apply(result, root) for root in roots} == roots
                    return result
                queue.append(nxt)
    return lifts[perm]


# Normalize the chosen representatives to exact order five; 16 is 1 mod 5.
ua, ub = (mat_power(lift(p), 16) for p in (a, b))
for perm, matrix in ((a, ua),(b,ub)):
    assert mat_power(matrix,5)==identity4
    assert action(matrix,probes)==perm
u_comm = matmul(matmul(matmul(ua,ub),adjoint(ua)),adjoint(ub))
assert action(u_comm,probes)==commutator(a,b)
print("All 60 native reflections preserve all 240 E8 roots; all 15 exact transpositions recovered.")
print("Both clock representatives lift to exact order-five native Gaussian isometries.")
print("Fixed same-vacuum normalized quartic amplitudes on marker 4: C_A=1, C_B=-1/5, commutator=-1/5.")
assert all((Fraction(1) if p[4]==4 else Fraction(-1,5))**2 == want
           for p,want in ((a,1),(b,Fraction(1,25)),(commutator(a,b),Fraction(1,25))))


paulis = (((1,0),(0,1)),((0,1),(1,0)),((0,-1j),(1j,0)),((1,0),(0,-1)))
kernel = {tuple(tuple(phase * p[i//2][j//2]*q[i%2][j%2] for j in range(4)) for i in range(4))
          for p,q,phase in product(paulis,paulis,(1,-1,1j,-1j))}
assert len(kernel)==64
qa = lift(quarter)
as_five = [u for k in kernel if mat_power(u:=matmul(k,ua),5)==identity4]
bs_five = [u for k in kernel if mat_power(u:=matmul(k,ub),5)==identity4]
pair_counts = Counter()
for k in kernel:
    g = matmul(k,qa)
    aa = sum(mat_conj(g,u)==mat_power(u,2) for u in as_five)
    bb = sum(mat_conj(g,u)==mat_power(u,2) for u in bs_five)
    pair_counts[(aa,bb)] += 1
print("All native Gaussian lifts by the 64-element kernel: order-five candidates per clock",
      len(as_five),len(bs_five))
print("Number of lifts G grouped by exact compatible lifts of each clock:",dict(sorted(pair_counts.items())))
central = {tuple(tuple(phase*z for z in row) for row in identity4) for phase in (1,-1,1j,-1j)}
criteria = Counter()
for k in kernel:
    g = matmul(k,qa)
    g4 = mat_power(g,4)
    ca = [u for u in as_five if mat_conj(g,u)==mat_power(u,2)]
    cb = [u for u in bs_five if mat_conj(g,u)==mat_power(u,2)]
    assert (bool(ca) and bool(cb)) == (g4 in central)
    criteria[(g4[0][0] if g4 in central else "noncentral",len(ca),len(cb))] += 1
    if ca and cb:
        witness = (g,ca[0],cb[0])
assert pair_counts == {(0,0):32,(2,2):32}
print("Criterion for both exact squares: G^4 is scalar. Counts:",dict(criteria))
g, ua, ub = witness
assert action(g,probes)==quarter and action(ua,probes)==a and action(ub,probes)==b
assert mat_conj(g,ua)==mat_power(ua,2) and mat_conj(g,ub)==mat_power(ub,2)
assert all({apply(u,r) for r in roots}==roots for u in (g,ua,ub))

r0 = lift(swap)
r_lifts = [u for k in kernel if mat_power(u:=matmul(k,r0),2)==identity4]
exchange = Counter()
for k in kernel:
    g = matmul(k,qa)
    ca = [u for u in as_five if mat_conj(g,u)==mat_power(u,2)]
    cb = [u for u in bs_five if mat_conj(g,u)==mat_power(u,2)]
    if ca and cb:
        exchange[sum(matmul(r,g)==matmul(g,r) and mat_conj(r,ua) in cb
                     for r,ua in product(r_lifts,ca))] += 1
print("Involutive lifts of the marker swap:",len(r_lifts))
print("Simultaneous exact swapping witnesses per successful quarter lift:",dict(exchange))


def conjugate_entries(matrix):
    return tuple(tuple(z.conjugate() for z in row) for row in matrix)


assert {tuple(z.conjugate() for z in r) for r in roots} == roots
mirror_types = Counter()
exchange_mod_phase = Counter()
for k in kernel:
    g = matmul(k,qa)
    if mat_power(g,4) != identity4:
        continue
    ca = [u for u in as_five if mat_conj(g,u)==mat_power(u,2)]
    cb = [u for u in bs_five if mat_conj(g,u)==mat_power(u,2)]
    count_by_type = Counter()
    phase_by_type = Counter()
    for s in kernel:
        r = matmul(s,r0)
        for kind, change in (("linear", lambda m:m),("antilinear",conjugate_entries)):
            square = matmul(r,change(r))
            turn = mat_conj(r,change(g))
            if square==identity4 and turn==g:
                count_by_type[kind] += sum(mat_conj(r,change(u)) in cb for u in ca)
            if square in central and matmul(turn,adjoint(g)) in central:
                phase_by_type[kind] += sum(mat_conj(r,change(u)) in cb for u in ca)
    mirror_types[(count_by_type["linear"],count_by_type["antilinear"])] += 1
    exchange_mod_phase[(phase_by_type["linear"],phase_by_type["antilinear"])] += 1
assert mirror_types == {(0, 8): 32} and exchange_mod_phase == {(0, 16): 32}
print("Exact involutive exchange counts (linear, antilinear) for the 32 common lifts:",dict(mirror_types))
print("Projective exchange counts (linear, antilinear):",dict(exchange_mod_phase))
