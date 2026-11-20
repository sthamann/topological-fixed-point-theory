"""Auxiliary (not one of the two decisive enumerations): numbers quoted for hypothesis 4."""
import numpy as np
me, mm, mt = 0.51099895, 105.6583755, 1776.86   # MeV, remembered PDG values (tau +-0.12)
s = np.sqrt([me, mm, mt]); K = (me + mm + mt) / s.sum() ** 2
print('Koide K = %.6f  (2/3 = %.6f), K-2/3 = %.1e' % (K, 2 / 3, K - 2 / 3))
u = np.ones(3) / np.sqrt(3); c = (s @ u) / np.linalg.norm(s)
print('angle between sqrt-mass vector and (1,1,1): %.4f deg (45 deg <=> K=2/3); trivial-character weight = %.5f' % (np.degrees(np.arccos(c)), c * c))
for dt in (-0.12, 0.12):
    s2 = np.sqrt([me, mm, mt + dt]); print('  m_tau%+.2f: K = %.6f' % (dt, (me + mm + mt + dt) / s2.sum() ** 2))
print('rephasing-invariant CP phases of an n x n unitary mixing matrix, (n-1)(n-2)/2:', {n: (n - 1) * (n - 2) // 2 for n in range(1, 6)})
print('groups of order n with a non-real character: n=1,2 none; n=3 (Z3) first')
w = np.exp(2j * np.pi / 3); F = np.array([[1, 1, 1], [1, w, w * w], [1, w * w, w]]) / np.sqrt(3)
# Majorana matrix that is Z3-invariant must be symmetric: retro-circulant; real eigenvectors
Mnu = np.array([[1.0, 0.3, 0.3], [0.3, 0.3, 1.0], [0.3, 1.0, 0.3]]) * 1.0
Mnu = np.array([[2.0, 0.5, 0.5], [0.5, 0.5, 2.0], [0.5, 2.0, 0.5]])
_, V = np.linalg.eigh(Mnu)
U = F.conj().T @ V
print('|U_PMNS|^2 from circulant charged leptons x retro-circulant Majorana neutrinos:\n', np.round(abs(U) ** 2, 3))
print('observed approx |U|^2 rows: e:(0.68,0.30,0.022) -> the 0 entry (theta13) is NOT reproduced exactly')
