"""W: Weil-(Oszillator-)Darstellung von SL(2, F_p) auf C^p.

Baut die Darstellung aus der diskreten Fouriertransformation S und dem Chirp T auf,
erzeugt die Gruppe vollständig (Breitensuche), zerlegt nach dem zentralen Element -I in
einen bosonischen Teil (-I wirkt als +1) und einen fermionischen Teil (-I wirkt als -1) und
prüft für p = 3, 5, 7, 11, 13:
  - Gruppenordnung p(p^2 - 1) (Darstellung linear, nicht nur projektiv),
  - Dimensionen und Treue beider Teile, SU(2)- bzw. SO(3)-Wertigkeit,
  - ob der bosonische Teil Sym^2 des fermionischen ist,
  - Kommutant auf C^p und die einzige spurfreie Größe darin (Hyperladung),
  - Spektrum der Fourier-Vierteldrehung S.
Numerisch mit Toleranz 1e-9; Gruppenordnungen und Multiplizitäten sind ganzzahlig.
"""
import itertools
import numpy as np

TOL = 1e-9


def key(m):
    r = np.round(m, 7)
    return tuple(np.concatenate([r.real.ravel(), r.imag.ravel()]) + 0.0)


def generate(gens, limit=60000):
    n = gens[0].shape[0]
    start = np.eye(n, dtype=complex)
    seen = {key(start): start}
    frontier = [start]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                m = g @ h
                k = key(m)
                if k not in seen:
                    seen[k] = m
                    nxt.append(m)
                    if len(seen) > limit:
                        raise RuntimeError("group too large")
        frontier = nxt
    return list(seen.values())


def weil(p, a=1):
    """Weil representation for the additive character x -> exp(2 pi i a x / p)."""
    z = np.exp(2j * np.pi / p)
    F = np.array([[z ** ((a * x * y) % p) for y in range(p)] for x in range(p)]) / np.sqrt(p)
    half = (p + 1) // 2                      # 1/2 mod p
    T = np.diag([z ** ((a * half * x * x) % p) for x in range(p)])
    target = p * (p * p - 1)
    # The DFT realises the Weyl element only up to a scalar; pick the scalar that makes
    # the generated group exactly SL(2, p).
    for k in range(8 * p):
        lam = np.exp(2j * np.pi * k / (8 * p))
        S = lam * F
        try:
            G = generate([S, T], limit=3 * target)
        except RuntimeError:
            continue
        if len(G) == target:
            return S, T, G
    raise RuntimeError(f"no linear Weil representation found for p={p}")


def central_minus_one(G, gens, n):
    for g in G:
        if np.allclose(g, np.eye(n), atol=TOL):
            continue
        if all(np.allclose(g @ h, h @ g, atol=TOL) for h in gens) and np.allclose(g @ g, np.eye(n), atol=TOL):
            return g
    raise RuntimeError("no central involution")


def block_basis(Z, sign):
    w, v = np.linalg.eigh((Z + Z.conj().T) / 2)
    return v[:, np.abs(w - sign) < 1e-6]


def restrict(g, V):
    return V.conj().T @ g @ V


def galois5(v):
    """Galois conjugate sqrt5 -> -sqrt5 of a real number of the form (m + n sqrt5)/2."""
    for m in range(-12, 13):
        for k in range(-12, 13):
            if abs(v - (m + k * np.sqrt(5)) / 2) < 1e-7:
                return (m - k * np.sqrt(5)) / 2
    raise ValueError(v)


def analyse(p, a=1):
    S, T, G = weil(p, a)
    n = p
    Z = central_minus_one(G, [S, T], n)
    parity = np.zeros((p, p))
    for x in range(p):
        parity[(-x) % p, x] = 1
    z_is_parity = np.allclose(Z, parity, atol=TOL)
    Vb, Vf = block_basis(Z, +1), block_basis(Z, -1)
    out = {"p": p, "order": len(G), "W(-I)": "+Parität" if z_is_parity else "-Parität",
           "dim_bos": Vb.shape[1], "dim_ferm": Vf.shape[1]}
    # faithfulness and determinants
    for name, V in (("bos", Vb), ("ferm", Vf)):
        if V.shape[1] == 0:
            out[f"kernel_{name}"] = len(G)
            continue
        imgs = [restrict(g, V) for g in G]
        out[f"kernel_{name}"] = sum(np.allclose(m, np.eye(V.shape[1]), atol=1e-7) for m in imgs)
        dets = [np.linalg.det(m) for m in imgs]
        out[f"det1_{name}"] = all(abs(d - 1) < 1e-7 for d in dets)
        chars = [np.trace(m) for m in imgs]
        out[f"chars_{name}"] = chars
        # Frobenius-Schur indicator: (1/|G|) sum chi(g^2)
        out[f"fs_{name}"] = round((sum(np.trace(restrict(g @ g, V)) for g in G) / len(G)).real, 6)
    # Sym^2 of the fermionic block versus the bosonic block
    if Vf.shape[1] == 2 and Vb.shape[1] == 3:
        chi_f = np.array(out["chars_ferm"])
        chi_f2 = np.array([np.trace(restrict(g @ g, Vf)) for g in G])
        sym2 = (chi_f ** 2 + chi_f2) / 2
        out["bos_equals_sym2_ferm"] = bool(np.allclose(sym2, out["chars_bos"], atol=1e-7))
        if p == 5:
            gal = np.array([galois5(c.real) for c in chi_f])
            gal2 = np.array([galois5(c.real) for c in chi_f2])
            out["bos_equals_sym2_galois_ferm"] = bool(np.allclose((gal ** 2 + gal2) / 2, out["chars_bos"], atol=1e-7))
    # commutant on C^p
    eqs = []
    I = np.eye(n)
    for h in (S, T):
        eqs.append(np.kron(I, h.T) - np.kron(h, I))   # vec(X h - h X) for row-major vec
    A = np.vstack(eqs)
    sv = np.linalg.svd(A, compute_uv=False)
    out["commutant_dim"] = int(sum(sv < 1e-8))
    # traceless element of the span {I, Z}: Y = I - (n / tr Z) Z when tr Z != 0
    trZ = np.trace(Z).real
    out["trace_Z"] = round(trZ, 9)
    if abs(trZ) > 1e-9 and out["commutant_dim"] == 2:
        Y = I - (n / trZ) * Z
        ev = np.round(np.linalg.eigvalsh((Y + Y.conj().T) / 2), 9)
        scale = 1 / 12 if p == 5 else 1.0
        out["Y_eigen_scaled"] = sorted(set(np.round(ev * scale, 9)))
        out["Y_ratio_bos_over_ferm"] = round((1 - n / trZ) / (1 + n / trZ), 9)
    # spectrum of the Fourier quarter turn
    ev = np.linalg.eigvals(S)
    out["S_order"] = next(k for k in range(1, 9) if np.allclose(np.linalg.matrix_power(S, k), I, atol=TOL))
    out["S_spectrum"] = sorted([complex(round(e.real, 6), round(e.imag, 6)) for e in ev],
                               key=lambda c: (c.real, c.imag))
    out["S2_is_Z"] = bool(np.allclose(S @ S, Z, atol=TOL))
    return out


def main():
    print("Weil-Darstellung von SL(2,p) auf C^p\n")
    rows = []
    for p in (3, 5, 7, 11, 13):
        r = analyse(p)
        rows.append(r)
        print(f"p = {p}: |Gruppe| = {r['order']} (= p(p²-1) = {p*(p*p-1)}), W(-I) = {r['W(-I)']}")
        print(f"   bosonisch dim {r['dim_bos']}, Kern {r['kernel_bos']}, det=1: {r.get('det1_bos')}, "
              f"Frobenius-Schur {r.get('fs_bos')}")
        print(f"   fermionisch dim {r['dim_ferm']}, Kern {r['kernel_ferm']}, det=1: {r.get('det1_ferm')}, "
              f"Frobenius-Schur {r.get('fs_ferm')}")
        if "bos_equals_sym2_ferm" in r:
            print(f"   bosonischer Teil = Sym² des fermionischen: {r['bos_equals_sym2_ferm']}; "
                  f"= Sym² des galoiskonjugierten Dubletts (√5 → −√5): {r.get('bos_equals_sym2_galois_ferm')}")
        print(f"   Kommutant auf C^p: Dimension {r['commutant_dim']}, tr W(-I) = {r['trace_Z']}")
        if "Y_eigen_scaled" in r:
            print(f"   einzige spurfreie Größe im Kommutanten: Eigenwerte {r['Y_eigen_scaled']}"
                  f"{' (Skala 1/12)' if p == 5 else ''}, Verhältnis bos/ferm {r['Y_ratio_bos_over_ferm']}")
        spec = ", ".join(f"{c.real:+.0f}{c.imag:+.0f}i" for c in r["S_spectrum"])
        print(f"   Fourier-Vierteldrehung S: Ordnung {r['S_order']}, S² = W(-I): {r['S2_is_Z']}, Spektrum {spec}, "
              f"Spur {sum(r['S_spectrum']).real:+.0f}")
        print()
    print("Auswahl:")
    for r in rows:
        su2 = r["dim_ferm"] == 2 and r["kernel_ferm"] == 1 and r.get("det1_ferm")
        so3 = r["dim_bos"] == 3 and r["kernel_bos"] == 2 and r.get("det1_bos") and r.get("fs_bos") == 1.0
        print(f"   p={r['p']:2d}: fermionisch treues SU(2) Dublett: {bool(su2)}; "
              f"bosonisch treue reelle 3D Drehung von SL(2,p)/±1: {bool(so3)}")
    print("\nZweite Weil-Darstellung bei p = 5 (additiver Charakter mit Nichtrest a = 2):")
    r5b = analyse(5, a=2)
    spec = ", ".join(f"{c.real:+.0f}{c.imag:+.0f}i" for c in r5b["S_spectrum"])
    print(f"   bos dim {r5b['dim_bos']} (FS {r5b['fs_bos']}, det=1 {r5b['det1_bos']}), ferm dim {r5b['dim_ferm']} "
          f"(Kern {r5b['kernel_ferm']}, FS {r5b['fs_ferm']}, det=1 {r5b['det1_ferm']})")
    print(f"   bos = Sym² ferm: {r5b['bos_equals_sym2_ferm']}; = Sym² galoiskonjugiert: {r5b['bos_equals_sym2_galois_ferm']}")
    print(f"   S Spektrum {spec}, Spur {sum(r5b['S_spectrum']).real:+.0f}")
    r5 = next(r for r in rows if r["p"] == 5)
    same = np.allclose(sorted(np.round(r5["chars_ferm"], 6), key=lambda c: (c.real, c.imag)),
                       sorted(np.round(r5b["chars_ferm"], 6), key=lambda c: (c.real, c.imag)))
    print(f"   Charaktermenge des Dubletts gleich wie bei a = 1: {same}")
    # check the expected exact facts for p = 5
    assert r5["order"] == 120 and r5["dim_bos"] == 3 and r5["dim_ferm"] == 2
    assert r5["kernel_ferm"] == 1 and r5["det1_ferm"] and r5["kernel_bos"] == 2 and r5["det1_bos"]
    assert r5["fs_bos"] == 1.0 and r5["fs_ferm"] == -1.0
    assert r5["commutant_dim"] == 2 and r5["Y_eigen_scaled"] == [-1 / 3, 0.5] or \
        np.allclose(r5["Y_eigen_scaled"], [-1 / 3, 0.5])
    print("\nAlle Prüfungen für p = 5 bestanden.")


if __name__ == "__main__":
    main()
