"""Z7: Ursprung der Naht — Defekt des Spiegelmassenfelds.

Auf dem H6 Gitter hat jede chirale Generation psi einen Spiegelpartner chi (Nielsen-Ninomiya). Ein
Massenfeld m e^{i theta} psi chi macht die Spiegel massiv. Hypothese: Die Naht ist ein Defekt dieses
Felds, um den psi und chi mit der Vierteldrehung J: psi -> chi -> -psi vertauscht werden, während die
Massenphase um pi windet. Das ist eindeutig, weil J die Masse psi chi umpolt. J^2 = -1 auf Fermionen,
also Ordnung 2 auf allem Bosonischen: zwei Blätter, und auf Fermionen mu4.

In der Basis xi_pm = (psi -+ i chi)/sqrt 2 ist J diagonal (Phasen +-i). Mit eta_pm = e^{-+i phi/4} xi_pm
(eindeutig) bekommt eta_+ eine Majoranamasse mit Windung 1, eta_- eine ohne Windung, und beide sehen
einen Aharonov-Bohm-Fluss von einem Viertel. Erwartung nach Jackiw-Rossi: eine chirale Majorana-Weyl
Mode je Paar entlang des Defekts. Kontrolle: Gewöhnlicher Wirbel ohne J (Windung 1 für beide) gibt zwei.

Rechnung: 3D Gitter, Defektlinie entlang z, Querschnitt L x L offen, k_z erhalten. Bogoliubov-de-Gennes
mit Weyl-Kinetik, Wilson-Majoranamasse W eta_+ eta_- gegen die Verdoppler (unter J invariant).
Gezählt werden am Kern lokalisierte Zweige, die E = 0 kreuzen, mit ihrer Steigung.
"""
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

OK = []
MASS = {}


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
isy = 1j * sy


def build(L, kz, n_plus, n_minus, flux, m=0.6, rw=1.0, copies=1):
    m = MASS.get("m", m)
    """BdG Hamiltonian; species index runs over copies x (eta_+, eta_-); flux in units of 2 pi for eta_+."""
    c = (L - 1) / 2.0
    xs, ys = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")
    phi = np.arctan2(ys - c, xs - c)
    rad = np.hypot(xs - c, ys - c)
    prof = np.tanh(rad / 1.5)
    ns = 2 * copies
    dim = L * L * ns * 2
    idx = lambda x, y, a, s: ((x * L + y) * ns + a) * 2 + s
    h = np.zeros((dim, dim), complex)
    hz = np.zeros((dim, dim), complex)
    D = np.zeros((dim, dim), complex)
    charge = [flux if a % 2 == 0 else -flux for a in range(ns)]
    winding = [n_plus if a % 2 == 0 else n_minus for a in range(ns)]
    for x in range(L):
        for y in range(L):
            for a in range(ns):
                # on-site kinetic sigma_z sin kz (kept apart: the hole block needs h(-kz))
                for s in range(2):
                    for t in range(2):
                        hz[idx(x, y, a, s), idx(x, y, a, t)] += np.sin(kz) * sz[s, t]
                # on-site Majorana masses: eta_a eta_a with winding, and Wilson eta_+ eta_- on-site part
                Mdiag = m * prof[x, y] * np.exp(1j * winding[a] * phi[x, y])
                for s in range(2):
                    for t in range(2):
                        D[idx(x, y, a, s), idx(x, y, a, t)] += Mdiag * isy[s, t]
                if a % 2 == 0:
                    Wos = rw * (2 + 1 - np.cos(kz))
                    for s in range(2):
                        for t in range(2):
                            D[idx(x, y, a, s), idx(x, y, a + 1, t)] += Wos * isy[s, t]
                            D[idx(x, y, a + 1, s), idx(x, y, a, t)] += Wos * isy[s, t]
            for (dx, dy, sig) in ((1, 0, sx), (0, 1, sy)):
                x2, y2 = x + dx, y + dy
                if x2 >= L or y2 >= L:
                    continue
                dphi = np.angle(np.exp(1j * (phi[x2, y2] - phi[x, y])))   # short difference
                for a in range(ns):
                    ph = np.exp(-1j * charge[a] * dphi)
                    # kinetic (1/2i) sigma (psi_r^dag psi_{r+e} - h.c.)
                    for s in range(2):
                        for t in range(2):
                            val = sig[s, t] / (2j) * ph
                            h[idx(x, y, a, s), idx(x2, y2, a, t)] += val
                            h[idx(x2, y2, a, t), idx(x, y, a, s)] += np.conj(val)
                    if a % 2 == 0:
                        # Wilson pairing hopping -rw/2 (eta_+ at r, eta_- at r+e) and (eta_- at r, eta_+ at r+e)
                        for (aa, bb, q) in ((a, a + 1, +1), (a + 1, a, -1)):
                            pw = -rw / 2 * np.exp(1j * q * flux * (phi[x2, y2] - phi[x, y] - dphi) * 0) \
                                * np.exp(-1j * q * flux * dphi)
                            for s in range(2):
                                for t in range(2):
                                    D[idx(x, y, aa, s), idx(x2, y2, bb, t)] += pw * isy[s, t]
                                    D[idx(x2, y2, bb, t), idx(x, y, aa, s)] += pw * isy[t, s]
    H = np.block([[h + hz, D], [D.conj().T, -(h - hz).T]])   # hole block = -h(-kz)^T
    return H, (rad.ravel(), ns)


def core_branches(L, n_plus, n_minus, flux, copies=1, R=3.5, window=0.15, dk=1e-3):
    """Core-localized zero modes at k_z = 0 and their velocities along z.

    The low-energy subspace mixes core and outer-edge modes of equal energy; diagonalizing the core
    projector inside it separates them. Velocity = <dH/dk_z> in each core state.
    """
    H, (rad, ns) = build(L, 0.0, n_plus, n_minus, flux, copies=copies)
    E, V = eigh(H, driver="evr")
    S = V[:, np.abs(E) < window]
    w = np.abs(S) ** 2
    core_mask = np.tile(np.repeat(rad <= R, ns * 2), 2)
    C = S.conj().T @ (core_mask[:, None] * S)
    cval, cvec = np.linalg.eigh(C)
    Hp, _ = build(L, dk, n_plus, n_minus, flux, copies=copies)
    Hm, _ = build(L, -dk, n_plus, n_minus, flux, copies=copies)
    dH = (Hp - Hm) / (2 * dk)
    core_states = S @ cvec[:, cval > 0.7]
    edge_states = S @ cvec[:, cval < 0.3]
    vc = np.real(np.einsum("ij,ik,kj->j", core_states.conj(), dH, core_states))
    ve = np.real(np.einsum("ij,ik,kj->j", edge_states.conj(), dH, edge_states))
    return len(S[0]), cval, vc, ve


L = 14
print(f"Querschnitt {L} x {L}, Kernradius 3.5; BdG-Zweige am Kern bei k_z = 0\n")
cases = [
    ("verdrillter Defekt (J, Phase windet um pi), ein Paar", 1, 0, 0.25, 1),
    ("gewöhnlicher Wirbel (ohne J, Windung 1), ein Paar", 1, 1, 0.0, 1),
    ("nur J, keine Massenwindung", 0, 0, 0.25, 1),
    ("verdrillter Defekt, zwei Paare", 1, 0, 0.25, 2),
]
summary = {}
for name, npl, nmi, fl, cp in cases:
    nlow, cval, vc, ve = core_branches(L, npl, nmi, fl, copies=cp)
    chiral_core = int(np.sum(vc > 0.2) - np.sum(vc < -0.2))
    chiral_edge = int(np.sum(ve > 0.2) - np.sum(ve < -0.2))
    summary[name] = {"low_states": nlow, "core_velocities": [round(float(x), 3) for x in vc],
                     "edge_velocities": [round(float(x), 3) for x in ve], "net_core": chiral_core, "net_edge": chiral_edge}
    print(f"   {name}:")
    print(f"      Nullnahe Zustände {nlow}; am Kern v = {np.round(vc, 3)} (netto {chiral_core}); am Rand v = "
          f"{np.round(ve, 3)} (netto {chiral_edge})")
check("verdrillt: eine Mode je Paar", summary[cases[0][0]]["net_core"] == 1 and summary[cases[3][0]]["net_core"] == 2)
check("gewöhnlicher Wirbel: zwei", summary[cases[1][0]]["net_core"] == 2)
check("ohne Windung: keine", summary[cases[2][0]]["net_core"] == 0)
check("Kern mit Lichtgeschwindigkeit", all(abs(v - 1) < 0.02 for v in summary[cases[0][0]]["core_velocities"]))
print("\nRobustheit (verdrillter Defekt, ein Paar):")
robust = []
for L2, m2 in ((18, 0.6), (14, 0.4), (14, 0.9)):
    MASS["m"] = m2
    nlow, cval, vc, ve = core_branches(L2, 1, 0, 0.25)
    net = int(np.sum(vc > 0.2) - np.sum(vc < -0.2))
    robust.append(net)
    print(f"   L = {L2}, m = {m2}: Kern v = {np.round(vc, 3)}, netto {net}")
MASS.clear()
check("robust", robust == [1, 1, 1])
print("\nGegenprobe: Liefern die Verdoppler eines naiven Feldes die Moden? (ohne Wilson-Glied, gleiche Masse für alle)")
MASS.clear()
_build = build


def build_nowilson(L, kz, n_plus, n_minus, flux, m=0.6, rw=1.0, copies=1):
    return _build(L, kz, n_plus, n_minus, flux, m=m, rw=0.0, copies=copies)


build = build_nowilson
nlow, cval, vc, ve = core_branches(14, 1, 0, 0.25)
net_tw = int(np.sum(vc > 0.2) - np.sum(vc < -0.2))
nlow2, cval2, vc2, ve2 = core_branches(14, 1, 1, 0.0)
build = _build
print(f"   verdrillt: Kern v = {np.round(vc, 3)} (netto {net_tw}); gewöhnlich: Kern v = {np.round(vc2, 3)}")
print("   Ohne Wilson-Glied haben die acht räumlichen Verdoppler abwechselnde Chiralität; der verdrillte Defekt")
print("   trägt dann keine Mode, beim gewöhnlichen Wirbel mischen sich Moden beider Richtungen. Die 16 Moden")
print("   der Naht kommen also nicht aus den Verdopplern; sie brauchen 16 echte Paare.")
check("Verdoppler geben keine Moden (verdrillt)", net_tw == 0)
J = np.array([[0, -1], [1, 0]])
print(f"\nJ^2 = {(J @ J).tolist()} = -1 auf Fermionen, J^4 = 1: Ordnung 4 (mu4); auf Bilinearen Ordnung 2 (zwei Blätter).")
check("J^2 = -1", np.array_equal(J @ J, -np.eye(2, dtype=int)))
print("Mit 16 Paaren (eine Generation aus 16 Weyl-Fermionen und ihre Partner) trägt der Defekt 16 Majorana-Weyl")
print("Moden, also c = 8. Die Zahl der Moden ist die Zahl der Paare, deren Masse am Defekt windet.")
print("\nZusammenfassung")
for n, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
summary["checks"] = {n: c for n, c in OK}
summary["robust"] = robust
Path(__file__).with_name("z7_ergebnis.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1))
assert all(c for _, c in OK)
