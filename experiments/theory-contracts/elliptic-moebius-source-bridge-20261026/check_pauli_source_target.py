#!/usr/bin/env python3
"""Exact finite check for the source Jac(C)[2] versus native Pauli target.

All complex entries used below are Gaussian integers.  The F2 part is exhaustive:
all 35 two-planes and all 2^16 four by four binary matrices are inspected.
"""

from itertools import combinations, permutations, product


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print("PASS", label)


def matmul(a, b, mod=None):
    out = tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(len(b)))
              for j in range(len(b[0])))
        for i in range(len(a))
    )
    if mod is None:
        return out
    return tuple(tuple(x % mod for x in row) for row in out)


def dagger(a):
    return tuple(tuple(a[j][i].conjugate() for j in range(len(a)))
                 for i in range(len(a[0])))


def kron(a, b):
    return tuple(tuple(a[i // len(b)][j // len(b[0])] *
                       b[i % len(b)][j % len(b[0])]
                       for j in range(len(a[0]) * len(b[0])))
                 for i in range(len(a) * len(b)))


def identity(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def f2_apply(a, v):
    return tuple(sum(a[i][j] * v[j] for j in range(len(v))) % 2
                 for i in range(len(a)))


def f2_inverse(a):
    n = len(a)
    aug = [list(a[i]) + [int(i == j) for j in range(n)] for i in range(n)]
    for col in range(n):
        pivot = next(row for row in range(col, n) if aug[row][col])
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(n):
            if row != col and aug[row][col]:
                aug[row] = [x ^ y for x, y in zip(aug[row], aug[col])]
    return tuple(tuple(row[n:]) for row in aug)


def from_columns(*columns):
    return tuple(tuple(column[row] for column in columns)
                 for row in range(len(columns[0])))


def block_diag(a, b):
    za = tuple(0 for _ in range(len(a)))
    zb = tuple(0 for _ in range(len(b)))
    return tuple(tuple(a[i][j] if j < len(a) else 0
                       for j in range(len(a) + len(b)))
                 for i in range(len(a))) + tuple(
        tuple(0 if j < len(a) else b[i][j - len(a)]
              for j in range(len(a) + len(b)))
        for i in range(len(b))
    )


def symplectic(v, w):
    return (v[0] * w[1] + v[1] * w[0] +
            v[2] * w[3] + v[3] * w[2]) % 2


def polynomial_mul(p, q):
    out = 0
    i = 0
    while p:
        if p & 1:
            out ^= q << i
        p >>= 1
        i += 1
    return out


def characteristic_polynomial_f2(a):
    # Polynomials are bit strings: bit k is the coefficient of x^k.
    n = len(a)
    det = 0
    for perm in permutations(range(n)):
        term = 1
        for i in range(n):
            entry = (0b10 if i == perm[i] else 0) ^ a[i][perm[i]]
            term = polynomial_mul(term, entry)
        # In characteristic two the permutation sign is invisible.
        det ^= term
    return det


def main():
    # Frozen v783 chart: z'=(z2,z0,z1,z3).
    sigma = (
        (0, 0, 1, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 0, 1),
    )
    check("frozen sigma has trace 1 and order 3",
          sum(sigma[i][i] for i in range(4)) == 1 and
          matmul(matmul(sigma, sigma), sigma) == identity(4))

    i = 1j
    one_qubit = {
        "I": ((1, 0), (0, 1)),
        "X": ((0, 1), (1, 0)),
        "Z": ((1, 0), (0, -1)),
        "Y": ((0, -i), (i, 0)),
    }
    paulis = {a + b: kron(one_qubit[a], one_qubit[b])
              for a in "IXZY" for b in "IXZY"}
    one_name = {(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}
    bit_name = {v: one_name[v[:2]] + one_name[v[2:]]
                for v in product((0, 1), repeat=4)}

    generator_bits = ((1, 0, 0, 0), (0, 1, 0, 0),
                      (0, 0, 1, 0), (0, 0, 0, 1))
    images = []
    signed_images = []
    for bits in generator_bits:
        image = matmul(matmul(sigma, paulis[bit_name[bits]]), dagger(sigma))
        hit = next((candidate, sign) for candidate in product((0, 1), repeat=4)
                   for sign in (1, -1)
                   if image == tuple(tuple(sign * x for x in row)
                                     for row in paulis[bit_name[candidate]]))
        images.append(hit[0])
        signed_images.append((bit_name[bits], sign_name :=
                              ("+" if hit[1] == 1 else "-"),
                              bit_name[hit[0]]))

    m = from_columns(*images)
    expected_m = (
        (0, 0, 1, 0),
        (0, 1, 0, 1),
        (1, 0, 1, 0),
        (0, 1, 0, 0),
    )
    check("exact Pauli generator conjugation",
          signed_images == [("XI", "+", "IX"), ("ZI", "-", "ZZ"),
                            ("IX", "+", "XX"), ("IZ", "+", "ZI")])
    check("frozen Pauli symplectic matrix", m == expected_m)

    i4 = identity(4)
    m2 = matmul(m, m, 2)
    minpoly_zero = tuple(tuple(m2[r][c] ^ m[r][c] ^ i4[r][c]
                               for c in range(4)) for r in range(4))
    fixed_rank_full = f2_inverse(tuple(tuple(m[r][c] ^ i4[r][c]
                                             for c in range(4))
                                      for r in range(4)))
    check("minimal polynomial x^2+x+1 and no nonzero fixed Pauli",
          minpoly_zero == ((0, 0, 0, 0),) * 4 and fixed_rank_full is not None)
    check("characteristic polynomial (x^2+x+1)^2=x^4+x^2+1",
          characteristic_polynomial_f2(m) == 0b10101)

    vectors = tuple(product((0, 1), repeat=4))
    nonzero = tuple(v for v in vectors if any(v))
    planes = set()
    for a, b in combinations(nonzero, 2):
        planes.add(frozenset(((0, 0, 0, 0), a, b,
                              tuple(x ^ y for x, y in zip(a, b)))))
    invariant = []
    for plane in planes:
        if {f2_apply(m, v) for v in plane} != set(plane):
            continue
        nz = tuple(v for v in plane if any(v))
        kind = "symplectic" if any(symplectic(a, b)
                                    for a, b in combinations(nz, 2)) else "lagrangian"
        invariant.append((kind, frozenset(bit_name[v] for v in nz), plane))

    expected_lagrangians = {
        frozenset(("IX", "XI", "XX")),
        frozenset(("IZ", "ZI", "ZZ")),
        frozenset(("XZ", "YY", "ZX")),
    }
    expected_symplectic = {
        frozenset(("ZY", "XY", "YI")),
        frozenset(("IY", "YX", "YZ")),
    }
    check("exactly three invariant contexts and two invariant symplectic planes",
          {names for kind, names, _ in invariant if kind == "lagrangian"} == expected_lagrangians and
          {names for kind, names, _ in invariant if kind == "symplectic"} == expected_symplectic)
    symplectic_planes = [plane for kind, _, plane in invariant if kind == "symplectic"]
    check("the two symplectic planes are orthogonal complements",
          len(symplectic_planes) == 2 and
          all(symplectic(v, w) == 0
              for v in symplectic_planes[0] for w in symplectic_planes[1]))

    c_source = ((1, 1), (1, 0))
    q_source = ((1, 0), (1, 1))
    # Ordered bases (ZY,YI) and (IY,YZ).
    p = from_columns((0, 1, 1, 1), (1, 1, 0, 0),
                     (0, 0, 1, 1), (1, 1, 0, 1))
    c_block = block_diag(c_source, c_source)
    q_block = block_diag(q_source, q_source)
    pinv = f2_inverse(p)
    check("source cubic block equals native sigma in the displayed bases",
          matmul(matmul(pinv, m, 2), p, 2) == c_block)
    q_native = matmul(matmul(p, q_block, 2), pinv, 2)
    check("displayed native Q is symplectic, involutive, and reverses sigma",
          matmul(q_native, q_native, 2) == i4 and
          matmul(matmul(q_native, m, 2), q_native, 2) == m2)

    omega = ((0, 1, 0, 0), (1, 0, 0, 0),
             (0, 0, 0, 1), (0, 0, 1, 0))
    sp4 = []
    for entries in product((0, 1), repeat=16):
        a = tuple(tuple(entries[4 * r + c] for c in range(4)) for r in range(4))
        try:
            f2_inverse(a)
        except StopIteration:
            continue
        if matmul(matmul(tuple(zip(*a)), omega, 2), a, 2) == omega:
            sp4.append(a)
    check("full Sp(4,2) enumeration", len(sp4) == 720)
    c_intertwiners = [a for a in sp4 if matmul(a, c_block, 2) == matmul(m, a, 2)]
    cq_intertwiners = [a for a in c_intertwiners
                       if matmul(a, q_block, 2) == matmul(q_native, a, 2)]
    q_candidates = [a for a in sp4
                    if matmul(a, a, 2) == i4
                    and matmul(matmul(a, m, 2), a, 2) == m2
                    and all({f2_apply(a, v) for v in plane} == set(plane)
                            for plane in symplectic_planes)]
    check("18 c-equivariant symplectic identifications", len(c_intertwiners) == 18)
    check("9 plane-preserving native Q candidates", len(q_candidates) == 9)
    check("2 simultaneous (c,Q) identifications after Q is fixed",
          len(cq_intertwiners) == 2)

    code_sigma = sigma
    code_fixed = tuple(tuple(code_sigma[r][c] ^ i4[r][c] for c in range(4))
                       for r in range(4))
    code_fixed_vectors = [v for v in vectors if f2_apply(code_fixed, v) == (0, 0, 0, 0)]
    pauli_fixed_vectors = [v for v in vectors if f2_apply(m, v) == v]
    check("code fixed dimension 2 versus Pauli fixed dimension 0",
          len(code_fixed_vectors) == 4 and len(pauli_fixed_vectors) == 1)

    print("sigma Pauli action:", signed_images)
    print("M_pauli =", m)
    print("Q_native =", q_native)
    print("invariant planes =", sorted((kind, sorted(names))
                                       for kind, names, _ in invariant))
    print("counts: Sp4=720, c-intertwiners=18, Q-candidates=9, (c,Q)-intertwiners=2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
