"""Z9: Ladungsaustausch über das Volumen (Zufluss) und die Frage nach dem Defekt im Grundzustand.

A. Zufluss. Zwei Komponenten des reellen 16-Pletts bilden eine exakte SO(2)-Ladung Q (Drehung zwischen
   Komponente 1 und 2; Massen und Wilson-Glied sind in den Komponenten diagonal). Im Sektor Q = +1
   (Teilchen von u = eta_1 + i eta_2, Löcher von v) ist eine gleichförmige Verschiebung von k_z genau ein
   SO(2)-Eichfeld entlang der Naht. Gemessen wird die Ladung der besetzten Zustände im Kern und am äußeren
   Rand über k_z.
B. Explizite Brechung. Das Wilson-Glied, das die Verdoppler entfernt (Nielsen-Ninomiya), bricht die U(1)
   der Paarungsphase. Die Energiedichte hängt dann von der relativen Phase phi0 der windenden Masse ab.
C. Folgerung für den Grundzustand.
"""
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

_src = Path(__file__).with_name("z7_ursprung_naht.py").read_text()
_ns = {"__file__": str(Path(__file__).with_name("z7_ursprung_naht.py"))}
exec(_src[:_src.index("L = 14")], _ns)
build = _ns["build"]
OK = []


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


print("A. Ladungsaustausch unter Fluss entlang der Naht (Sektor Q = +1)\n")
L, R = 12, 3.5
H0, (rad, ns) = build(L, 0.3, 1, 0, 0.25, copies=2)
nsite = L * L
qs = np.kron(np.array([[0, -1j], [1j, 0]]), np.eye(2))
q = np.kron(np.eye(nsite), np.kron(qs, np.eye(2)))
QB = np.block([[q, np.zeros_like(q)], [np.zeros_like(q), -q.T]])
print(f"   [H, Q] = {np.max(np.abs(H0 @ QB - QB @ H0)):.1e} {check('Q erhalten', np.allclose(H0 @ QB, QB @ H0))}")
wq, U = np.linalg.eigh(QB)
Up = U[:, wq > 0]
core = np.tile(np.repeat(rad <= R, ns * 2), 2).astype(float)
outer = np.tile(np.repeat(rad >= L / 2 - 2.0, ns * 2), 2).astype(float)
Pc, Po = Up.conj().T @ (core[:, None] * Up), Up.conj().T @ (outer[:, None] * Up)
ks = np.unique(np.concatenate([np.linspace(-np.pi, np.pi, 25), [-0.1, -0.05, 0.05, 0.1]]))
rows, wlow = [], {}
for kz in ks:
    H, _ = build(L, kz, 1, 0, 0.25, copies=2)
    E, V = eigh(Up.conj().T @ H @ Up, driver="evr")
    occ = E < 0
    fc = np.real(np.einsum("ij,ik,kj->j", V.conj(), Pc, V))
    fo = np.real(np.einsum("ij,ik,kj->j", V.conj(), Po, V))
    rows.append((kz, fc[occ].sum(), fo[occ].sum(), int(occ.sum())))
    lowi = np.argsort(np.abs(E))[:2]
    wlow[round(kz, 3)] = max(fc[lowi])
rows = np.array(rows)
f = dict(zip(np.round(rows[:, 0], 3), rows[:, 1]))
g = dict(zip(np.round(rows[:, 0], 3), rows[:, 2]))
left = 2 * f[-0.05] - f[-0.1]                                   # linear extrapolation to k -> 0-
right = 2 * f[0.05] - f[0.1]
jump_c = right - left
jump_o = (2 * g[0.05] - g[0.1]) - (2 * g[-0.05] - g[-0.1])
core_w = wlow[0.05]
smooth_c = rows[-1, 1] - rows[0, 1] - jump_c
print(f"   Zahl der besetzten Zustände über k_z: {sorted(set(int(x) for x in rows[:, 3]))} (Ladung erhalten)")
print(f"   Kernladung: Sprung bei k_z = 0 (beidseitig extrapoliert) {jump_c:+.4f}, stetig über die Periode {smooth_c:+.4f}")
print(f"   Außenrand:  Sprung bei k_z = 0 {jump_o:+.4f}")
print(f"   Kerngewicht der chiralen Kernmode bei k_z = 0,05: {core_w:.4f}; Verhältnis Sprung/Gewicht = {abs(jump_c) / core_w:.3f}")
check("Kern verliert genau eine Mode", abs(abs(jump_c) / core_w - 1) < 0.15)
check("Kern und Rand tauschen aus", abs(jump_c + jump_o) < 0.15)
check("stetiger Ausgleich", abs(jump_c + smooth_c) < 0.02)
print("   Unter Fluss wechselt die eine geladene chirale Mode am Kern ihre Besetzung, die Gegenmode am Rand die")
print("   umgekehrte; dazwischen liegt nur lückenhaftes Volumen, durch das die Ladung stetig umverteilt wird.")
print("   Das ist der Zufluss für die SO(2)-Ladung. Der gravitative Zufluss (c = 8) ist damit nicht gerechnet.")

print("\nB. Explizite Brechung der Paarungsphase durch das Wilson-Glied\n")
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1.0, -1.0]).astype(complex)
isy = 1j * sy


def e_density(phi0, m=0.6, rw=1.0, n=12):
    k = 2 * np.pi * (np.arange(n) + 0.5) / n
    tot = 0.0
    for kx in k:
        for ky in k:
            for kz in k:
                hs = np.sin(kx) * sx + np.sin(ky) * sy + np.sin(kz) * sz
                W = rw * (3 - np.cos(kx) - np.cos(ky) - np.cos(kz))
                D = np.kron(np.array([[m * np.exp(1j * phi0), W], [W, m]]), isy)
                Hk = np.block([[np.kron(np.eye(2), hs), D], [D.conj().T, -np.kron(np.eye(2), -hs).T]])
                E = np.linalg.eigvalsh(Hk)
                tot += -0.5 * E[E > 0].sum()
    return tot / n ** 3


phis = np.linspace(0, 2 * np.pi, 9)
V = np.array([e_density(p) for p in phis])
for p, v in zip(phis, V):
    print(f"   phi0 = {p:5.3f}: Energiedichte {v:+.6f} je Gitterzelle (relativ {v - V.min():.4f})")
barrier = V.max() - V.min()
print(f"   Periode 2 pi in phi0 = pi in der Massenphase (J-Symmetrie), Minimum bei phi0 = {phis[V.argmin()]:.3f},"
      f" Barriere {barrier:.4f}")
check("Energie hängt von der Phase ab", barrier > 0.01)
check("Periode 2 pi in phi0", abs(V[0] - V[-1]) < 1e-9 and abs(V[2] - V[6]) < 1e-9)
print("   Die U(1) der Paarungsphase ist auf dem Gitter keine Symmetrie (Nielsen-Ninomiya zwingt das Wilson-")
print("   Glied). Der verdrillte Defekt dreht phi0 einmal ganz herum und muss über die Barriere: Er berandet")
print("   eine Domänenwand. Die Wand ist der Verzweigungsschnitt. Topologisch bleibt der Schnitt unsichtbar")
print("   (Anomaliephase, c = 8, Z5), energetisch ist er eine Wand mit Spannung auf Gitterskala.")
print("   Weil phi0 und phi0 + 2 pi dasselbe Vakuum sind, hat die Wand Domänenwandzahl 1 wie bei Axionen:")
print("   Wände enden auf Nähten, und Nähte beranden Wände.")

print("\nC. Defekt im Grundzustand\n")
print("   Im unendlichen Vakuum kondensiert der Defekt nicht: Er kostet Kern- und Wandenergie. Eine geschlossene")
print("   Naht berandet eine Wand; ihre euklidische Wirkung ist S = sigma_W V + sigma_N A - log chi_E8(i), mit")
print(f"   log chi_E8(i) = log 12 = {np.log(12):.3f} aus Z4. Nähte mit Wirkung unter log 12 sind nicht unterdrückt.")
print(f"   Die Barriere {barrier:.3f} je Zelle legt sigma_W auf die Größenordnung 0,1 Gittereinheiten; damit sind nur")
print("   Nähte auf Gitterskala zugelassen, große Nähte sind durch das Wandvolumen unterdrückt. Lesart: Der")
print("   Grundzustand enthält Nähte als Fluktuationen auf Gitterskala, nicht als makroskopische Objekte.")
print("   sigma_N und die Steifigkeit der Phase sind nicht gerechnet; die Aussage ist eine Größenordnung.")

print("\nZusammenfassung")
for n, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
Path(__file__).with_name("z9_ergebnis.json").write_text(json.dumps({
    "checks": {n: c for n, c in OK}, "core_jump": float(jump_c), "core_mode_weight": float(core_w),
    "phase_potential": [float(v - V.min()) for v in V], "barrier": float(barrier)}, indent=1))
assert all(c for _, c in OK)
