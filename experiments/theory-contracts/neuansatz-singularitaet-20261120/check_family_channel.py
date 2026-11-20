"""NEUSING20261120: Ist der verlangte Familienkanal (1, 1/3, eta, 2/3) eine Mischung der S4 Wirkungen selbst?"""
import json, itertools, numpy as np, sympy as sp
from scipy.optimize import linprog
# S4 auf C^3 = (Funktionen auf vier Marken minus Konstante) mal Vorzeichen (Charakter von v117: 3,-1,0,1)
def perm_mat(p):
    M = np.zeros((4, 4)); 
    for i, j in enumerate(p): M[j, i] = 1
    return M
Bs = np.array([[1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]) / 2.0      # orthonormale Basis des Komplements
reps = {}
for p in itertools.permutations(range(4)):
    sign = round(np.linalg.det(perm_mat(p)))
    reps[p] = sign * Bs @ perm_mat(p) @ Bs.T
def cyc(p):
    seen, out = set(), []
    for i in range(4):
        if i in seen: continue
        n, j = 0, i
        while j not in seen: seen.add(j); j = p[j]; n += 1
        out.append(n)
    return tuple(sorted(out))
classes = {}
for p, U in reps.items():
    classes.setdefault(cyc(p), []).append(np.kron(U, U.conj()))
S = {c: sum(v) / len(v) for c, v in classes.items()}
# Sektoren des Operatorraums: 1, 2, 3 (symmetrisch), 3 (antisymmetrisch = Drehungen)
I3 = np.eye(3).reshape(-1) / np.sqrt(3)
def sector_vals(M):
    vals = {}
    vals["1"] = I3 @ M @ I3
    d = np.diag([1., -1, 0]).reshape(-1) / np.sqrt(2); vals["2"] = d @ M @ d
    s = np.zeros((3, 3)); s[0, 1] = s[1, 0] = 1; s = s.reshape(-1) / np.sqrt(2); vals["3s"] = s @ M @ s
    a = np.zeros((3, 3)); a[0, 1] = 1; a[1, 0] = -1; a = a.reshape(-1) / np.sqrt(2); vals["3a"] = a @ M @ a
    return [round(float(np.real(vals[k])), 12) for k in ("1", "2", "3s", "3a")]
table = {str(c): sector_vals(M) for c, M in S.items()}
order = list(S)
A_eq = np.array([[1] * 5, [table[str(c)][1] for c in order], [table[str(c)][3] for c in order]])
res = linprog(np.zeros(5), A_eq=A_eq, b_eq=[1, 1/3, 2/3], bounds=[(0, None)] * 5)
# gleiche Frage fuer den Kanal des Neuansatzes: halb bleiben, halb einer der drei Schritte
V4 = [M for p, M in ((p, np.kron(U, U.conj())) for p, U in reps.items()) if cyc(p) == (2, 2)]
mine = 0.5 * np.eye(9) + sum(V4) / 6
out = {"sektorwerte je Klasse (1,2,3s,3a)": table,
       "Mischung der S4 Wirkungen mit Sektor 2 = 1/3 und Drehsektor = 2/3 existiert": bool(res.status == 0),
       "Neuansatz (halb bleiben, halb Schritt) als Quantenkanal, Sektorwerte": sector_vals(mine),
       "verlangt": [1, "1/3", "eta in [2/9,4/9]", "2/3"]}
json.dump(out, open("check_family_channel.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(out, ensure_ascii=False, indent=0))
