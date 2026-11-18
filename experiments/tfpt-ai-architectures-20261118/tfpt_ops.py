"""Exact TFPT quartic-code objects, as defined in
_newest2/TFPT_Gesamtdokumentation_Code_Quartik_Rekursion_2026-09-26.md (sections 4-6).

  H = (C^4)^{x4}  (dim 256),  S = Sym^4 projector (rank 35),
  Q_P = (1/16) sum_A A^{x4} over the 16 two-qubit Paulis (rank 16),
  P = S Q_P (rank 5),  M4 = (1/60) sum_l (|psi_l><psi_l|)^{x4},  40 M4 = S + P,
  H_mom = I - 20 M4  (energies 0, 1/2, 1).

Everything here is a finite exact construction; nothing is a physical claim.
"""
import itertools
import numpy as np

_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.diag([1, -1]).astype(complex)


def paulis2():
    return np.array([np.kron(a, b) for a in (_I, _X, _Y, _Z) for b in (_I, _X, _Y, _Z)])


def stabilizer_states():
    """The 60 two-qubit stabilizer states as the Clifford orbit of |00>."""
    H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    Sg = np.diag([1, 1j]).astype(complex)
    CN = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
    gens = [np.kron(H, _I), np.kron(_I, H), np.kron(Sg, _I), np.kron(_I, Sg), CN]

    def canon(v):
        k = np.argmax(np.abs(v) > 1e-9)
        v = v * np.exp(-1j * np.angle(v[k]))
        return v / np.linalg.norm(v)

    start = canon(np.array([1, 0, 0, 0], dtype=complex))
    seen = {tuple(np.round(start, 6)): start}
    frontier = [start]
    while frontier:
        nxt = []
        for v in frontier:
            for g in gens:
                w = canon(g @ v)
                key = tuple(np.round(w, 6))
                if key not in seen:
                    seen[key] = w
                    nxt.append(w)
        frontier = nxt
    return np.array(list(seen.values()))


def kron4(a, b, c, d):
    return np.kron(np.kron(a, b), np.kron(c, d))


def sym_projector():
    E = np.eye(256).reshape(4, 4, 4, 4, 256)
    S = np.zeros((256, 256))
    for p in itertools.permutations(range(4)):
        S += E.transpose(p + (4,)).reshape(256, 256)
    return S / 24


def pauli_projector():
    return sum(kron4(A, A, A, A) for A in paulis2()).real / 16


def code_basis():
    """Integer 256x5 indicator matrix B of the five supports and isometry V = B G^{-1/2}."""
    B = np.zeros((256, 5))
    for idx, w in enumerate(itertools.product(range(4), repeat=4)):
        vals = sorted(set(w))
        counts = sorted(w.count(v) for v in vals)
        if len(vals) == 1:
            B[idx, 0] = 1
        elif len(vals) == 2 and counts == [2, 2]:
            B[idx, vals[0] ^ vals[1]] = 1
        elif len(vals) == 4:
            B[idx, 4] = 1
    G = B.T @ B
    return B, B / np.sqrt(np.diag(G))


def build(verify=True):
    stab = stabilizer_states()
    S, Q = sym_projector(), pauli_projector()
    P = S @ Q
    M4 = sum(kron4(*(np.outer(v, v.conj()),) * 4) for v in stab) / 60
    B, V = code_basis()
    Hmom = np.eye(256) - 20 * M4.real
    if verify:
        rk = lambda A: int(np.linalg.matrix_rank(A, tol=1e-9))
        assert len(stab) == 60
        assert (rk(S), rk(Q), rk(P)) == (35, 16, 5)
        assert np.allclose(S @ Q, Q @ S) and np.allclose(P @ P, P)
        assert np.allclose(40 * M4, S + P), "40 M4 = S + P failed"
        assert np.allclose(V @ V.T, P), "V V^T = P failed"
        ev = np.round(np.linalg.eigvalsh(Hmom), 9)
        assert sorted(set(ev)) == [0.0, 0.5, 1.0]
        assert [int((ev == e).sum()) for e in (0.0, 0.5, 1.0)] == [5, 30, 221]
    return dict(stab=stab, S=S, Q=Q, P=P, M4=M4.real, B=B, V=V, Hmom=Hmom, pauli=paulis2())


# ---------------------------------------------------------------- torch side
def torch_ops(kind="tfpt", seed=0):
    """Tensors used by the networks.  kind='rand' swaps V for a random rank-5 isometry
    (the specificity control); everything else stays identical."""
    import torch
    o = build(verify=False)
    V = o["V"]
    if kind == "rand":
        rng = np.random.default_rng(1000 + seed)
        V = np.linalg.qr(rng.standard_normal((256, 5)))[0]
    return dict(
        V4=torch.tensor(V.reshape(4, 4, 4, 4, 5), dtype=torch.float32),
        stab=torch.tensor(o["stab"], dtype=torch.complex64),
        pauli=torch.tensor(o["pauli"], dtype=torch.complex64),
        had=torch.tensor([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]], dtype=torch.float32),
    )


def code_feat(V4, z):
    """z: (...,4,4) real amplitudes of a four-register product state -> V^T psi, (...,5)."""
    import torch
    return torch.einsum("abcdk,...a,...b,...c,...d->...k", V4, z[..., 0, :], z[..., 1, :], z[..., 2, :], z[..., 3, :])


def code_messages(V4, z):
    """Parameter-free transition operator: d<P>/dz_i (up to 2), shape (...,4 slots,4)."""
    import torch
    a, b, c, d = (z[..., i, :] for i in range(4))
    g = torch.stack([
        torch.einsum("abcdk,...b,...c,...d->...ak", V4, b, c, d),
        torch.einsum("abcdk,...a,...c,...d->...bk", V4, a, c, d),
        torch.einsum("abcdk,...a,...b,...d->...ck", V4, a, b, d),
        torch.einsum("abcdk,...a,...b,...c->...dk", V4, a, b, c)], dim=-3)
    f = torch.einsum("...ak,...a->...k", g[..., 0, :, :], a)
    return (g * f[..., None, None, :]).sum(-1), f


def exp_P(V4, z):
    return (code_feat(V4, z) ** 2).sum(-1)


def exp_Q(pauli, z):
    """<Q_P> on a product state; z (...,4,4) real or complex."""
    import torch
    zc = z.to(torch.complex64)
    e = torch.einsum("...ia,Aab,...ib->...iA", zc.conj(), pauli, zc).real
    return e.prod(-2).sum(-1) / 16


def exp_20M4(stab, z):
    """<20 M4> = (1/3) sum_l prod_i |<psi_l|z_i>|^2 ; H_mom = 1 - this."""
    import torch
    ov = torch.einsum("la,...ia->...il", stab.conj(), z.to(torch.complex64)).abs() ** 2
    return ov.prod(-2).sum(-1) / 3


def prob_xor0(had, p):
    """Classical shadow of 4 Q_P: Pr[x1^x2^x3^x4 = 0] for independent p (...,4,4)."""
    return (p @ had).prod(-2).sum(-1) / 4


def prob_alldiff(p):
    """Classical shadow of 24 |c4><c4|: Pr[all four symbols distinct]."""
    out = 0
    for pi in itertools.permutations(range(4)):
        out = out + p[..., 0, pi[0]] * p[..., 1, pi[1]] * p[..., 2, pi[2]] * p[..., 3, pi[3]]
    return out


# ---------------------------------------------------------------- E8 lattice
def e8_roots():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = np.zeros(8); v[i], v[j] = si, sj; R.append(v)
    for s in itertools.product((0.5, -0.5), repeat=8):
        if sum(x < 0 for x in s) % 2 == 0:
            R.append(np.array(s))
    return np.array(R)


def _d8(x):
    f = np.rint(x)
    odd = (f.sum(-1).astype(np.int64) % 2) != 0
    err = x - f
    k = np.abs(err).argmax(-1)
    g = f.copy()
    rows = np.arange(x.shape[0])
    g[rows, k] += np.where(err[rows, k] >= 0, 1.0, -1.0)
    return np.where(odd[:, None], g, f)


def nearest(x, lattice):
    """Nearest lattice point for x (n,8); lattices have unit covolume after the stated scale."""
    if lattice == "Z8":
        return np.rint(x)
    if lattice == "D8":               # covolume 2 -> scale by 2^(-1/8)
        s = 2 ** (-1 / 8)
        return _d8(x / s) * s
    if lattice == "E8":               # covolume 1
        a, b = _d8(x), _d8(x - 0.5) + 0.5
        pick = ((x - a) ** 2).sum(-1) <= ((x - b) ** 2).sum(-1)
        return np.where(pick[:, None], a, b)
    raise ValueError(lattice)


if __name__ == "__main__":
    o = build()
    print("ok: 60 states, ranks 35/16/5, 40 M4 = S + P, spec(H_mom) = {0:5, 1/2:30, 1:221}")
    R = e8_roots()
    print("E8 roots", len(R), "norm^2", set(np.round((R ** 2).sum(1), 9)))
