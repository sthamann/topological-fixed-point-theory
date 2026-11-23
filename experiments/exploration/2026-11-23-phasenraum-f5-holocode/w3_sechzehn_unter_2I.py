"""W3: Zerlegung der 16 = Λ^gerade C^5 unter der Weil-Wirkung der binären Ikosaedergruppe.

Die Weil-Darstellung von SL(2,5) = 2I auf C^5 zerfällt in 3' (bosonisch) und 2 (fermionisch).
Auf Λ^k C^5 setzen wir die induzierte Wirkung an und zerlegen jeden Hyperladungssektor
(Y = (I - 5P)/12) in irreduzible Darstellungen über die Charakternorm <χ, χ>.
"""
import itertools
import numpy as np

exec(open(__file__.replace("w3_sechzehn_unter_2I.py", "w1_weil_darstellung.py")).read().split("def main")[0])

S, T, G = weil(5)
P = np.zeros((5, 5))
for x in range(5):
    P[(-x) % 5, x] = 1
Y = (np.eye(5) - 5 * P) / 12
# basis in which Y is diagonal: 3 even then 2 odd
w, V = np.linalg.eigh(Y)
Gd = [V.conj().T @ g @ V for g in G]
yv = np.round(w, 9)


def wedge_matrix(g, k):
    idx = list(itertools.combinations(range(5), k))
    M = np.zeros((len(idx), len(idx)), dtype=complex)
    for a, I in enumerate(idx):
        for b, J in enumerate(idx):
            M[a, b] = np.linalg.det(g[np.ix_(I, J)])
    return idx, M


print("Λ^gerade C^5 unter 2I, nach Hyperladung:")
total = 0
for k in (0, 2, 4):
    idx, _ = wedge_matrix(Gd[0], k)
    charges = [round(sum(yv[i] for i in I), 6) for I in idx]
    for q in sorted(set(charges)):
        sel = [a for a, c in enumerate(charges) if c == q]
        chars = [np.trace(wedge_matrix(g, k)[1][np.ix_(sel, sel)]) for g in Gd]
        norm = round((sum(abs(c) ** 2 for c in chars) / len(Gd)).real, 6)
        triv = round((sum(chars) / len(Gd)).real, 6)
        total += len(sel)
        print(f"   Λ^{k}, Y = {q:+.4f}: Dimension {len(sel)}, <χ,χ> = {norm} "
              f"({'irreduzibel' if abs(norm - 1) < 1e-6 else 'zerfällt'}), Anteil trivial {triv}")
print(f"   Summe der Dimensionen {total}")

# Gegenprobe: mit dem ungedrehten Spin-1 Vektor Sym^2(2) statt 3' zerfiele der Q-Block.
f = np.array([abs(y - 0.5) < 1e-9 for y in yv])
chi_f = np.array([np.trace(g[np.ix_(f, f)]) for g in Gd])
chi_f2 = np.array([np.trace((g @ g)[np.ix_(f, f)]) for g in Gd])
chi_vec = (chi_f ** 2 + chi_f2) / 2
norm_q_untwisted = round((sum(abs(c) ** 2 for c in chi_vec * chi_f) / len(Gd)).real, 6)
print(f"   Gegenprobe Sym²(2) ⊗ 2: <χ,χ> = {norm_q_untwisted} (2 = zwei irreduzible Teile, 2 ⊕ 4)")
assert abs(norm_q_untwisted - 2) < 1e-6
