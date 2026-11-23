"""W2: Die TFPT Uhren- und Ladungsoperatoren als Clifford-Gruppe eines Quantensystems der
Dimension 5.

Prüft:
  a) Die Weil-Darstellung von SL(2,5) enthält die Dilatationen x -> 2x, x -> 3x bis auf Vorzeichen.
  b) Verschiebung X (x -> x+1, TFPT C5) und Dilatation G (x -> 2x) erfüllen G X G^-1 = X^2 und
     erzeugen AGL(1,5) = F20 (TFPT: G C5 G^-1 = C5^2, G^4 = I).
  c) Die sechs gegenseitig unverzerrten Basen (MUBs) von C^5; die unitäre Clifford-Gruppe
     permutiert sie als PSL(2,5) = A5 (Ordnung 60, Stabilisator 10); auch die antiunitäre
     Konjugation K ergibt nichts Neues. TFPT hat auf sechs Fasern PGL(2,5) mit Stabilisator F20;
     das ist mit einem einzelnen Fünfersystem nicht erreichbar (negatives Ergebnis).
  d) Y = (I - 5 P)/12 mit der Parität P = S^2 kommutiert mit SL(2,5), nicht mit der Verschiebung.
"""
import itertools
import numpy as np

p = 5
z = np.exp(2j * np.pi / p)
X = np.zeros((p, p)); [X.__setitem__(((x + 1) % p, x), 1) for x in range(p)]
Zc = np.diag([z ** x for x in range(p)])
F = np.array([[z ** ((x * y) % p) for y in range(p)] for x in range(p)]) / np.sqrt(p)
T = np.diag([z ** ((3 * x * x) % p) for x in range(p)])          # 3 = 1/2 mod 5


def dil(a):
    M = np.zeros((p, p))
    for x in range(p):
        M[(a * x) % p, x] = 1
    return M


G = dil(2)
P = dil(p - 1)                                                     # parity x -> -x


def key(m):
    r = np.round(m, 7)
    return tuple(np.concatenate([r.real.ravel(), r.imag.ravel()]) + 0.0)


def generate(gens, cap=200000):
    e = np.eye(p, dtype=complex)
    seen, frontier = {key(e): e}, [e]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                m = g @ h
                k = key(m)
                if k not in seen:
                    seen[k] = m
                    nxt.append(m)
        frontier = nxt
        if len(seen) > cap:
            raise RuntimeError("too large")
    return list(seen.values())


# a) Weil image of SL(2,5) (scalar on F chosen as in w1: the group must have order 120)
for k in range(40):
    S = np.exp(2j * np.pi * k / 40) * F
    W = generate([S, T])
    if len(W) == 120:
        break
in_W = {name: [s for s in (+1, -1) if any(np.allclose(w, s * M) for w in W)]
        for name, M in (("x->2x", dil(2)), ("x->3x", dil(3)), ("x->4x (Parität)", P))}
print(f"a) |Weil-Bild| = {len(W)}; Dilatationen im Bild mit Vorzeichen: {in_W}")

# b) clock relation
print(f"b) G X G^-1 = X^2: {np.allclose(G @ X @ np.linalg.inv(G), X @ X)};  G^4 = I: "
      f"{np.allclose(np.linalg.matrix_power(G, 4), np.eye(p))};  |<X, G>| = {len(generate([X, G]))}")

# c) mutually unbiased bases and their permutation by the Clifford group
bases = [np.eye(p)]
for k in range(p):
    w, v = np.linalg.eig(X @ np.linalg.matrix_power(Zc, k))
    bases.append(v / np.linalg.norm(v, axis=0))
unbiased = all(np.allclose(np.abs(bases[i].conj().T @ bases[j]) ** 2, 1 / p, atol=1e-9)
               for i, j in itertools.combinations(range(6), 2))


def which_basis(U, B):
    img = U @ B
    for j, C in enumerate(bases):
        overlap = np.abs(C.conj().T @ img) ** 2
        if np.allclose(np.sort(overlap.max(axis=0)), 1, atol=1e-7):
            return j
    raise RuntimeError("no basis")


def as_perm(U):
    return tuple(which_basis(U, B) for B in bases)


def compose(a, b):
    return tuple(b[i] for i in a)


def closure(gens):
    e = tuple(range(6))
    seen, frontier = {e}, [e]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                m = compose(g, h)
                if m not in seen:
                    seen.add(m)
                    nxt.append(m)
        frontier = nxt
    return seen


gens6 = {"X (C5)": as_perm(X), "Z": as_perm(Zc), "S (Fourier)": as_perm(S), "T (Chirp)": as_perm(T),
         "G (x->2x)": as_perm(G)}
sl = closure([gens6["S (Fourier)"], gens6["T (Chirp)"]])
full = closure(list(gens6.values()))
stab = [g for g in full if g[0] == 0]
print(f"c) sechs MUBs gegenseitig unverzerrt: {unbiased}")
print(f"   Wirkung auf den sechs Basen: SL(2,5) -> Gruppe der Ordnung {len(sl)}; "
      f"mit X, Z, G -> Ordnung {len(full)}, transitiv: {len({g[0] for g in full}) == 6}, "
      f"Punktstabilisator Ordnung {len(stab)}")
print(f"   Erzeuger als Permutationen der sechs Basen: {gens6}")

# c2) the antiunitary complex conjugation K (Tomita/orientation type) on the six bases
def as_perm_conj():
    out = []
    for B in bases:
        img = B.conj()
        for j, C in enumerate(bases):
            if np.allclose(np.sort((np.abs(C.conj().T @ img) ** 2).max(axis=0)), 1, atol=1e-7):
                out.append(j)
                break
    return tuple(out)


K = as_perm_conj()
extended = closure(list(gens6.values()) + [K])
stab_ext = [g for g in extended if g[0] == 0]
print(f"   mit antiunitärer Konjugation K = {K}: Ordnung {len(extended)}, "
      f"Punktstabilisator Ordnung {len(stab_ext)}")

# d) hypercharge
Y = (np.eye(p) - 5 * P) / 12
print(f"d) Y Eigenwerte {sorted(set(np.round(np.linalg.eigvalsh(Y), 9)))}; "
      f"[Y, SL(2,5)] = 0: {all(np.allclose(Y @ w, w @ Y) for w in W)}; "
      f"[Y, X] = 0: {np.allclose(Y @ X, X @ Y)}; [Y, G] = 0: {np.allclose(Y @ G, G @ Y)}")
assert len(W) == 120 and unbiased and len(full) == 60 and len(stab) == 10
# Negative result: for p = 1 mod 4 the conjugation already acts inside PSL(2,5); the six-basis
# action stays A5. TFPT's six-fibre PGL(2,5) needs a similitude with non-square multiplier.
assert len(extended) == 60 and len(stab_ext) == 10
