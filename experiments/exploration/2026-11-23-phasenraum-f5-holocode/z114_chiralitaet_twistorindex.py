"""Z114: Tiefenwelle T5 'Chiralität als Index: Twistor-Homogenität, Familienindex 3 und radiale Leiter'.

Leitlinie 'Rekonstruktion statt Auswahl': eine Lesart zählt nur als nativ, wenn sie aus vorhandener TFPT-Struktur folgt.
Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (jede 3+1-Aussage ist Kinematik der radialen Einpaarschicht oder Lesart),
GRAV.NONCIRCULAR.01. Gelöste Bausteine (drei Familien über den mu4-Familienindex, E_- = O(-1) + 3 O mit Index 3 ohne
Spiegel, Familiencharaktere {-1, i, -i}, die 16 je Generation, E8 = D5 + A3 über Z4, Flavour) werden benutzt, nicht neu
geöffnet. Bestand: z58 (Familienindex E_-: h^0 = 3, h^1 = 0), z61 (Lesart T: Naht = Twistorgerade; h^0(O(2)) = 3 =
Lambda^2_+), z83 (radiale Standardform, CCR-Schritte S_a, S-quer_b), z97 (Z_SU = i^(2h-2), Helizität h <-> O(2h - 2);
GSO-Regel P_GSO, vektorartig je R-Klasse; S_a ist ein (1/2, 0)-Spinor, S-quer ein (0, 1/2)-Spinor; J tauscht Original und
Partner), z102 (Kleberroute geschlossen: jeder Z4-Charakter, der h = +-1/2 trennt, trennt h = 0 von h = 1; Empfehlung:
nichtzentraler Mechanismus, Index), z109 Teil 6 (der Index eines Elternoperators folgt nicht aus Daten auf Hyperflächen).

Bezeichnungen: Twistor Z = (omega^0, omega^1, pi_0, pi_1), PT' = {pi != 0}, Geraden L_x = {omega^A = x^(AA') pi_A'} = P(pi).
Konvention z97: h = (m - mquer)/2 <-> O(2h - 2) auf PT' mit Geraden P(pi) (pi = Partnerseite, wird in Teil 2 explizit
gezeigt). Die Penrose-Konvention der Aufgabe, O(-2h - 2), ist die duale Realisierung mit Geraden P(omega) der Originalseite,
also das J-Bild. Zwei Kopplungen des Familienbündels F an die Twistorgerade: K1 'Verdrillung' (Naht = Twistorgerade, F (x)
O(n_h) auf derselben Geraden) und K2 'Produkt' (Familiennullmoden der Naht mal Kohomologie der Twistorgerade, Künneth auf
Naht x L).

Teil 1  Index auf der Naht CP^1 und den Twistorgeraden L_x: Cech-Kohomologie von O(n) (Laurent-Monome, Rang des
        Korandoperators) und unabhängig der analytische d-quer-Operator auf spingewichteten Kugelfunktionen (Kern und
        Kokern per SVD); die gelösten Familienblöcke (E_-, H^0(K(D)) mit mu4-Charakteren); Tabellen h^0, h^1, Index für
        h in {+-1, +-1/2, 0} in beiden Konventionen, für F in {O, E_-, E_+ = E_-^v, 4 O, 3 O, Spiegelfamilie O(-4) + 3 O,
        Serre-Dual E_-^v (x) K}, in K1 und K2; Asymmetriesatz, affines Gesetz, Starrheit von E_- auf allen L_x.
Teil 2  Radiale su(2,2): explizite Penrose-Transformation (Standardform z83 = Cech-H^1(PT', O(m - mquer - 2)), alle 17
        Erzeuger als Differentialoperatoren auf Cech-Klassen, Z_SU = Homogenitätscharakter), Einschränkung auf L_x (Rang je
        Helizität, J-Bild), Artinhalt mit Familienbündel (K1, K2), Nullmodenzählung je Z4-Sektor (16, 16q, neutrale
        Skalare und Vektoren) für alle 256 Linienregeln (Sektor -> Original- und/oder Partnergerade), J-Kovarianz,
        Negativkontrollen (triviales Familienbündel, Gegenorientierung, Spiegelfamilie, Verdrillung statt Produkt).
Teil 3  Bilanz, Urteil für T4, Folgen für K16 und die Klausel chiral_matter der Szenario-Suche.
"""
import itertools
import json
import math
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(114)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


fs = lambda x: str(Fr(x))                                           # Bruch als Text


# ------------------------------------------------------------------ Kohomologie auf CP^1
def cech_O(n, L=None):
    """h^0, h^1 von O(n) auf CP^1 aus der Cech-Überdeckung U_0 = {pi_0 != 0}, U_1 = {pi_1 != 0}. Kozyklen auf U_01 sind
    Laurent-Monome pi_0^p0 pi_1^p1 mit p0 + p1 = n; Koketten C^0 = Monome holomorph auf U_0 (p1 >= 0) bzw. auf U_1
    (p0 >= 0); delta(f0, f1) = f0 - f1. h^0 = dim ker delta, h^1 = dim C^1 - rang delta (Fenster |p0| <= L, randfrei für
    L > |n| + 1)."""
    L = abs(n) + 3 if L is None else L
    C1 = list(range(-L, L + 1))
    ix = {p: i for i, p in enumerate(C1)}
    sp_ = [(ix[p0], 1.0) for p0 in C1 if n - p0 >= 0] + [(ix[p0], -1.0) for p0 in C1 if p0 >= 0]
    D = np.zeros((len(C1), len(sp_)))
    for j, (i, v) in enumerate(sp_):
        D[i, j] = v
    r = int(np.linalg.matrix_rank(D))
    return D.shape[1] - r, D.shape[0] - r


def dquer_wigner(n, J):
    """Analytischer d-quer-Operator auf O(n) über S^2 = SU(2)/U(1): Schnitte = Wigner-Funktionen D^j_(mu, s) zum
    Spingewicht s = n/2 (j >= |s|), d-quer = rechtsinvarianter Heber D^j_(mu, s) -> sqrt((j - s)(j + s + 1)) D^j_(mu, s+1).
    Kern und Kokern per SVD bis j <= J."""
    s = Fr(n, 2)

    def jwerte(t):
        out, j = [], abs(t)
        while j <= J:
            out.append(j)
            j += 1
        return out
    src = [(j, j - k) for j in jwerte(s) for k in range(int(2*j) + 1)]
    tgt = [(j, j - k) for j in jwerte(s + 1) for k in range(int(2*j) + 1)]
    ti = {x: i for i, x in enumerate(tgt)}
    M = np.zeros((len(tgt), len(src)))
    for c, (j, mu) in enumerate(src):
        if (j, mu) in ti:
            M[ti[(j, mu)], c] = math.sqrt(float((j - s)*(j + s + 1)))
    sv = np.linalg.svd(M, compute_uv=False) if M.size else np.zeros(0)
    r = int(np.sum(sv > 1e-10))
    return len(src) - r, len(tgt) - r


def coh(F, n):
    """(h^0, h^1) von F (x) O(n), F = Liste der Grade einer Summe von Geradenbündeln."""
    h = [cech_O(a + n) for a in F]
    return sum(x for x, _ in h), sum(y for _, y in h)


def kuenneth(F, n):
    """(H^0, H^1, H^2) von F |x| O(n) auf Naht x L (Künneth)."""
    f0, f1 = coh(F, 0)
    l0, l1 = cech_O(n)
    return f0*l0, f0*l1 + f1*l0, f1*l1


ind = lambda hh: hh[0] - hh[1]
indK = lambda H: H[0] - H[1] + H[2]
HEL = [Fr(-1), Fr(-1, 2), Fr(0), Fr(1, 2), Fr(1)]
n_z97 = lambda h: int(2*h - 2)                                      # Helizität h <-> O(2h - 2)   (Geraden P(pi))
n_pen = lambda h: int(-2*h - 2)                                     # Helizität h <-> O(-2h - 2)  (Geraden P(omega))
BUENDEL = {"O": [0], "E_-": [-1, 0, 0, 0], "E_+ = E_-^v": [1, 0, 0, 0], "4 O": [0, 0, 0, 0], "3 O": [0, 0, 0],
           "Spiegelfamilie O(-4) + 3 O": [-4, 0, 0, 0], "E_-^v (x) K": [-1, -2, -2, -2]}
EM = BUENDEL["E_-"]

# =================================================================== Teil 1
say("Teil 1: Index auf der Naht CP^1 und den Twistorgeraden L_x (Cech, analytischer d-quer-Operator, Familienbündel)")

# --- 1.1 O(n): Cech, Riemann-Roch, Serre, analytischer Operator
tab_On = {}
ok_rr = ok_serre = ok_wig = True
for n in range(-6, 5):
    h0, h1 = cech_O(n)
    w1, w2 = dquer_wigner(n, Fr(abs(n), 2) + 5), dquer_wigner(n, Fr(abs(n), 2) + 8)
    tab_On[n] = (h0, h1)
    ok_rr &= (h0 - h1 == n + 1) and h0 == max(0, n + 1) and h1 == max(0, -n - 1)
    ok_serre &= h1 == cech_O(-n - 2)[0]
    ok_wig &= w1 == (h0, h1) and w2 == (h0, h1)
azyklisch = [n for n, (a, b) in tab_On.items() if a == 0 and b == 0]
ZAHLEN["O(n) auf CP^1 (h0, h1)"] = {str(n): v for n, v in tab_On.items()}
check("[E] O(n) auf CP^1, n = -6..4: Cech-Kohomologie (Laurent-Monome, Rang von delta) gibt h^0 = max(0, n + 1), h^1 = "
      "max(0, -n - 1), Index n + 1 (Riemann-Roch) und Serre h^1(O(n)) = h^0(O(-n - 2)); der analytische d-quer-Operator auf "
      "spingewichteten Kugelfunktionen hat dieselben Kern- und Kokerndimensionen (SVD, zwei Abschneider j <= |n|/2 + 5, + 8, "
      f"stabil). Einziges azyklisches Bündel: O({azyklisch[0]})",
      ok_rr and ok_serre and ok_wig and azyklisch == [-1])

# --- 1.2 gelöste Familienblöcke (benutzt, nicht neu geöffnet)
z, w = sp.symbols("z w")
marken = [1, sp.I, -1, -sp.I]
kd = {}
for k in range(4):
    f = z**k/(z**4 - 1)
    res = [sp.nsimplify(sp.residue(f, z, m)) for m in marken]
    g = sp.simplify(-f.subs(z, 1/w)/w**2)                           # f dz bei z = 1/w
    regulaer = sp.limit(g, w, 0).is_finite
    char = sp.simplify(sp.I*f.subs(z, sp.I*z)/f)                    # Rückzug unter z -> i z
    kd[k] = {"Residuen an den Marken": [str(r) for r in res], "Residuensumme": str(sp.nsimplify(sum(res))),
             "regulär bei unendlich": bool(regulaer), "mu4-Charakter": str(char)}
schnitte = [k for k in range(4) if kd[k]["regulär bei unendlich"]]
chars = sorted(kd[k]["mu4-Charakter"] for k in schnitte)
h0E, h1E = coh(EM, 0)
ZAHLEN["Familienblöcke"] = {"E_- (h0, h1)": (h0E, h1E), "H^0(K(D)) = p(z) dz/(z^4 - 1)": kd}
check("[E, benutzt] Gelöste Familienblöcke reproduziert: E_- = O(-1) + 3 O hat h^0 = 3, h^1 = 0 (Familienindex 3 ohne "
      f"Spiegel; gerechnet {h0E}, {h1E}). H^0(K(D)) mit D = vier Marken bei mu4: z^k dz/(z^4 - 1) ist genau für k = "
      f"{schnitte} ein Schnitt (regulär bei unendlich), mit den mu4-Charakteren {chars} = {{-1, i, -i}}; der vierte Kandidat "
      f"(k = 3, Charakter {kd[3]['mu4-Charakter']}) hat bei unendlich einen Pol (Residuensumme an den Marken "
      f"{kd[3]['Residuensumme']} != 0; für k = 0, 1, 2 ist sie {[kd[k]['Residuensumme'] for k in schnitte]})",
      (h0E, h1E) == (3, 0) and schnitte == [0, 1, 2] and set(chars) == {"I", "-1", "-I"} and kd[3]["mu4-Charakter"] == "1"
      and kd[3]["Residuensumme"] == "1" and all(kd[k]["Residuensumme"] == "0" for k in schnitte))

# --- 1.3 die Twistorgerade allein
tab_L = {}
for h in [Fr(k, 2) for k in range(-6, 7)]:
    tab_L[str(h)] = {"z97 O(2h-2)": cech_O(n_z97(h)), "Penrose O(-2h-2)": cech_O(n_pen(h))}
azyk_z97 = [h for h, v in tab_L.items() if v["z97 O(2h-2)"] == (0, 0)]
azyk_pen = [h for h, v in tab_L.items() if v["Penrose O(-2h-2)"] == (0, 0)]
lin_index = {h: ind(tab_L[str(h)]["z97 O(2h-2)"]) for h in HEL}
ZAHLEN["Twistorgerade allein (h0, h1) je Helizität"] = tab_L
check("[E] Twistorgerade allein: Index von O(2h - 2) ist 2h - 1 (h = -1, -1/2, 0, 1/2, 1: "
      f"{[lin_index[h] for h in HEL]}). Über alle |h| <= 3 hat genau eine Helizität keine Kohomologie auf der Geraden: "
      f"h = {azyk_z97} (Konvention z97) bzw. h = {azyk_pen} (Penrose-Konvention O(-2h - 2), das J-Bild). Die einzige "
      "Nullstelle des affinen Index liegt an einer Weyl-Helizität; Skalar O(-2) (h^1 = 1) und beide Vektorslots (O(0): "
      "h^0 = 1; O(-4): h^1 = 3) sind nicht azyklisch",
      azyk_z97 == ["1/2"] and azyk_pen == ["-1/2"] and all(lin_index[h] == 2*h - 1 for h in HEL)
      and tab_L["0"]["z97 O(2h-2)"] == (0, 1) and tab_L["1"]["z97 O(2h-2)"] == (1, 0)
      and tab_L["-1"]["z97 O(2h-2)"] == (0, 3))

# --- 1.4 Kopplung K1 (Verdrillung) mit E_-
tab_K1 = {name: {str(h): coh(F, n_z97(h)) for h in HEL} for name, F in BUENDEL.items()}
nullmoden_Em = [str(h) for h in HEL if tab_K1["E_-"][str(h)][0] > 0]
ZAHLEN["K1 Verdrillung F (x) O(2h-2): (h0, h1)"] = tab_K1
check("[E] K1 (Naht = Twistorgerade, E_- (x) O(2h - 2)): (h^0, h^1) für h = -1, -1/2, 0, 1/2, 1 = "
      f"{[tab_K1['E_-'][str(h)] for h in HEL]}, Index 8h - 5. Nullmoden (h^0 > 0) des Familienbündels gibt es nur bei "
      f"h = {nullmoden_Em} (Homogenität 0, der |h| = 1-Slot; passend zu Familien = Lambda^2_+ aus z61); an beiden "
      "Weyl-Slots hat E_- nur Kokern (h^1 = 1 bei +1/2, 9 bei -1/2). Der Familienindex 3 ist ein Index des Vektorslots, "
      "nicht des Weyl-Slots",
      nullmoden_Em == ["1"] and tab_K1["E_-"]["1"] == (3, 0) and tab_K1["E_-"]["1/2"] == (0, 1)
      and tab_K1["E_-"]["-1/2"] == (0, 9) and all(ind(tab_K1["E_-"][str(h)]) == 8*h - 5 for h in HEL))

# --- 1.5 Asymmetriesatz und affines Gesetz
asym, gerade, affin = {}, {}, {}
ok_asym = ok_aff = True
for name, F in BUENDEL.items():
    r_, d_ = len(F), sum(F)
    I_ = {h: ind(tab_K1[name][str(h)]) for h in HEL}
    asym[name] = {str(h): I_[h] - I_[-h] for h in HEL if h > 0}
    gerade[name] = {str(h): I_[h] + I_[-h] for h in HEL if h > 0}
    affin[name] = (I_[Fr(1, 2)] - I_[Fr(-1, 2)], I_[Fr(1)] - I_[Fr(0)])
    ok_asym &= all(I_[h] - I_[-h] == 4*h*r_ and I_[h] + I_[-h] == 2*d_ - 2*r_ for h in HEL if h > 0)
    ok_aff &= affin[name][0] == affin[name][1] == 2*r_
ZAHLEN["Asymmetrie ind(h) - ind(-h)"] = asym
ZAHLEN["gerader Teil ind(h) + ind(-h)"] = gerade
check("[E] Asymmetriesatz: für jedes F gilt ind(F, h) - ind(F, -h) = 4h rang F und ind(F, h) + ind(F, -h) = 2 grad F - "
      "2 rang F. Der Familienindex (über den Grad) steht nur im geraden Teil; die Helizitätsasymmetrie ist familienblind: "
      f"E_-, E_+ = E_-^v und das triviale 4 O haben dieselbe Asymmetrie ({asym['E_-']['1/2']}, {asym['E_+ = E_-^v']['1/2']}, "
      f"{asym['4 O']['1/2']} bei h = 1/2). Affines Gesetz (Indexform des Leiterquotienten aus z102): ind(1/2) - ind(-1/2) = "
      f"ind(1) - ind(0) = 2 rang F für alle F ({[affin[k] for k in BUENDEL]})",
      ok_asym and ok_aff and asym["E_-"]["1/2"] == asym["E_+ = E_-^v"]["1/2"] == asym["4 O"]["1/2"] == 8)

# --- 1.6 E_- am Spiegelslot, Skalarslot, Starrheit
chi_k = {k: ind(coh(EM, k)) for k in range(-12, 13)}
end_E = [b - a for a in EM for b in EM]                              # End E_- = Hom(E_-, E_-)
skalar_O, skalar_E = coh([0], -2), coh(EM, -2)
check("[negativ, E] E_- am azyklischen Slot: E_- (x) O(-1) = O(-2) + 3 O(-1) hat (h^0, h^1) = "
      f"{coh(EM, -1)}, also eine Spiegelmode aus dem O(-1)-Summanden; chi(E_- (x) O(k)) = 4k + 3 ist für k = -12..12 nie 0 "
      f"(min |chi| = {min(abs(v) for v in chi_k.values())}): keine Geradenverdrillung gibt E_- den Spiegelindex 0. Am "
      f"Skalarslot verschiebt der O(-1)-Summand die Helizität um -1/2: E_- (x) O(-2) = O(-3) + 3 O(-2) gibt {skalar_E}, die "
      f"familienneutrale Skalarleiter (F = O) bleibt mit {skalar_O} unberührt. E_- ist starr: h^1(End E_-) = "
      f"{coh(end_E, 0)[1]} (End E_- = 10 O + 3 O(1) + 3 O(-1)), also auf jeder Geraden L_x vom selben Spaltungstyp; der "
      "Index ist auf allen L_x derselbe",
      coh(EM, -1) == (0, 1) and all(v == 4*k + 3 and v != 0 for k, v in chi_k.items()) and skalar_E == (0, 5)
      and skalar_O == (0, 1) and coh(end_E, 0) == (16, 0) and sorted(end_E).count(0) == 10)

# --- 1.7 Kopplung K2 (Produkt, Künneth)
tab_K2 = {name: {str(h): kuenneth(F, n_z97(h)) for h in HEL} for name, F in BUENDEL.items()}
chiF = {name: ind(coh(F, 0)) for name, F in BUENDEL.items()}
ok_prod = all(indK(tab_K2[name][str(h)]) == chiF[name]*(2*h - 1) for name in BUENDEL for h in HEL)
ZAHLEN["K2 Produkt (H0, H1, H2) auf Naht x L"] = tab_K2
ZAHLEN["chi(F)"] = chiF
check("[E] K2 (Familiennullmoden der Naht mal Twistorgerade, Künneth): Index = chi(F) (2h - 1) exakt für alle F und h. "
      f"Für E_-: (H^0, H^1, H^2) bei h = -1/2 (chiraler Slot) = {tab_K2['E_-']['-1/2']} (drei Weyl, je 2 Komponenten, "
      f"kein H^2, weil h^1(E_-) = 0), bei h = +1/2 (Spiegelslot) = {tab_K2['E_-']['1/2']} (alle Kohomologie verschwindet). "
      f"Kontrollen: die Spiegelfamilie O(-4) + 3 O (z58) gibt am chiralen Slot {tab_K2['Spiegelfamilie O(-4) + 3 O']['-1/2']}, "
      f"Index 0 (Familienspiegel entgegengesetzter Parität); das triviale 3 O gibt {tab_K2['3 O']['-1/2']} wie E_-; das "
      f"Serre-Dual E_-^v (x) K (Gegenorientierung, chi = {chiF['E_-^v (x) K']}) gibt {tab_K2['E_-^v (x) K']['-1/2']}, "
      f"Parität gekippt. Skalarslot: E_- verdreifacht ({tab_K2['E_-']['0']}), F = O nicht ({tab_K2['O']['0']})",
      ok_prod and tab_K2["E_-"]["-1/2"] == (0, 6, 0) and tab_K2["E_-"]["1/2"] == (0, 0, 0)
      and tab_K2["Spiegelfamilie O(-4) + 3 O"]["-1/2"] == (0, 6, 6) and tab_K2["3 O"]["-1/2"] == (0, 6, 0)
      and tab_K2["E_-^v (x) K"]["-1/2"] == (0, 0, 6) and chiF["E_-^v (x) K"] == -3 and tab_K2["O"]["0"] == (0, 1, 0))
befund("Teil 1: Die Twistorgerade ist für sich chiral: der Index 2h - 1 verschwindet genau an einer Weyl-Helizität (O(-1) "
       "azyklisch), Skalar und Vektoren behalten Kohomologie. Der Familienindex 3 ist helizitätsblind: in der Verdrillung "
       "(K1) liegen seine Nullmoden am Vektorslot und die Asymmetrie ist 4h rang F, unabhängig vom Grad; E_- hat dort sogar "
       "eine Spiegelmode. Erst im Produkt (K2) gilt Index = Familienindex mal Twistorindex, mit drei Weyl am chiralen und "
       "null am Spiegelslot; dafür ist h^1(E_-) = 0 nötig, aber das triviale 3 O leistet dasselbe", True)

# =================================================================== Teil 2
say("\nTeil 2: Radiale su(2,2): explizite Penrose-Transformation, Linienrestriktion, Artinhalt, Sektoren und Linienregeln")

# --- 2.1 Standardform (wie z83/z97, CCR) und Cech-Realisierung
AMP_CCR = (lambda m, k: np.sqrt(k + 1), lambda m, k: np.sqrt(m - k))
SCHRITT = (lambda k: k + 1, lambda k: k)


def standardform(M):
    """(+)_(m, mquer <= M) H_m (x) H-quer_mquer wie z97: P[a,b] = S_a (x) S-quer_b, A, B Bilineare links/rechts, E."""
    lab = [(m, mb, k, kb) for m in range(1, M + 1) for mb in range(1, M + 1) for k in range(m) for kb in range(mb)]
    ix = {s: i for i, s in enumerate(lab)}
    n = len(lab)

    def dicht(r, c, w_):
        return sps.csr_matrix((w_, (r, c)), shape=(n, n)).toarray()
    P = {}
    for a in range(2):
        for b in range(2):
            r, c, w_ = [], [], []
            for (m, mb, k, kb), i in ix.items():
                if m < M and mb < M:
                    r.append(ix[(m + 1, mb + 1, SCHRITT[a](k), SCHRITT[b](kb))])
                    c.append(i)
                    w_.append(AMP_CCR[a](m, k)*AMP_CCR[b](mb, kb))
            P[(a, b)] = dicht(r, c, w_)

    def bilinear(c_, a_, seite):
        r, cc, w_ = [], [], []
        for (m, mb, k, kb), i in ix.items():
            mm, kk = (m, k) if seite == 0 else (mb, kb)
            k1 = kk - (1 if a_ == 0 else 0)
            if mm == 1 or not 0 <= k1 <= mm - 2:
                continue
            k2 = SCHRITT[c_](k1)
            r.append(ix[(m, mb, k2, kb) if seite == 0 else (m, mb, k, k2)])
            cc.append(i)
            w_.append(AMP_CCR[a_](mm - 1, k1)*AMP_CCR[c_](mm - 1, k1))
        return dicht(r, cc, w_)
    A = {(c_, a_): bilinear(c_, a_, 0) for c_ in range(2) for a_ in range(2)}
    B = {(c_, a_): bilinear(c_, a_, 1) for c_ in range(2) for a_ in range(2)}
    return lab, ix, P, A, B, np.diag([0.5*(s[0] + s[1]) for s in lab])


def op_anwenden(terme, mon):
    """Wendet sum_t c_t (Folge elementarer Operatoren) auf ein Laurent-Monom (a0, a1, p0, p1) an, exakt (Fraction).
    Elementar: ('dom', a) = d/d omega^a, ('mo', a) = omega^a *, ('dpi', b) = d/d pi_b, ('mpi', b) = pi_b *.
    Am Ende werden Korandterme (p0 >= 0 oder p1 >= 0: holomorph auf U_0 bzw. U_1) verworfen."""
    out = {}
    for c, folge in terme:
        cur = {mon: Fr(c)}
        for art, ind_ in folge:
            neu = {}
            for mm, v in cur.items():
                mm = list(mm)
                pos = ind_ if art in ("dom", "mo") else 2 + ind_
                if art in ("dom", "dpi"):
                    e = mm[pos]
                    if e == 0:
                        continue
                    mm[pos] -= 1
                    v = v*e
                else:
                    mm[pos] += 1
                neu[tuple(mm)] = neu.get(tuple(mm), 0) + v
            cur = neu
        for mm, v in cur.items():
            out[mm] = out.get(mm, 0) + v
    return {mm: v for mm, v in out.items() if v != 0 and mm[2] <= -1 and mm[3] <= -1}


TWISTOR = {}
for a in range(2):
    for b in range(2):
        TWISTOR[f"P_{a}{b}"] = [(-1, [("dpi", b), ("mo", a)])]                 # -omega^a d/dpi_b   (Heber)
        TWISTOR[f"K_{a}{b}"] = [(1, [("dom", a), ("mpi", b)])]                 # pi_b d/domega^a    (Senker)
        TWISTOR[f"A_{a}{b}"] = [(1, [("dom", b), ("mo", a)])]                  # omega^a d/domega^b
        TWISTOR[f"B_{a}{b}"] = [(-1, [("mpi", b), ("dpi", a)])]                # -d/dpi_a o pi_b
TWISTOR["E"] = [(Fr(1, 2), [("dom", 0), ("mo", 0)]), (Fr(1, 2), [("dom", 1), ("mo", 1)]),
                (Fr(-1, 2), [("dpi", 0), ("mpi", 0)]), (Fr(-1, 2), [("dpi", 1), ("mpi", 1)])]

M = 5
LAB, IXS, PS, AS, BS, ES = standardform(M)
NR = len(LAB)
zust = lambda s: (s[2], s[0] - 1 - s[2], s[3], s[1] - 1 - s[3])    # (m, mb, k, kb) -> (alpha0, alpha1, beta0, beta1)
ZI = {zust(s): i for i, s in enumerate(LAB)}
norm_c = lambda al: math.sqrt(math.factorial(al[2])*math.factorial(al[3])/(math.factorial(al[0])*math.factorial(al[1])))
monom = lambda al: (al[0], al[1], -1 - al[2], -1 - al[3])            # Cech-Klasse omega^alpha pi^(-1-beta)


def twistor_matrix(terme):
    X = np.zeros((NR, NR))
    for al, i in ZI.items():
        for mm, v in op_anwenden(terme, monom(al)).items():
            al2 = (mm[0], mm[1], -1 - mm[2], -1 - mm[3])
            if al2 in ZI:
                X[ZI[al2], i] += float(v)*norm_c(al)/norm_c(al2)
    return X


Z83 = {}
for a in range(2):
    for b in range(2):
        Z83[f"P_{a}{b}"] = PS[(a, b)]
        Z83[f"K_{a}{b}"] = PS[(a, b)].T
        Z83[f"A_{a}{b}"] = AS[(a, b)]
        Z83[f"B_{a}{b}"] = BS[(a, b)]
Z83["E"] = ES
abw = {k: float(np.abs(twistor_matrix(TWISTOR[k]) - Z83[k]).max()) for k in Z83}
nichttrivial = sum(1 for k in Z83 if np.abs(Z83[k]).max() > 0)
homog = np.array([al[0] + al[1] - al[2] - al[3] - 2 for al in map(zust, LAB)])     # n = N_a - N_b - 2
hel = np.array([(s[0] - s[1])/2 for s in LAB])
J3R = (BS[(0, 0)] - BS[(1, 1)])/2
Z_SU = np.exp(1j*np.pi*np.diag(ES + 2*J3R))
zsu_fehler = float(np.abs(Z_SU - 1j**(homog % 4)).max())
EV = np.diag(ES)
weyl = {fs(hh): [int(np.sum((hel == hh) & np.isclose(EV, N + 0.5))) for N in range(1, 5)] for hh in (0.5, -0.5)}
ZAHLEN["Penrose-Transformation"] = {"Zustände (M = 5)": NR, "max |Twistor - z83| je Erzeuger": abw,
                                    "Z_SU - i^n": zsu_fehler, "Weyl-Zählung bei E = N + 1/2, N = 1..4": weyl}
check("[E] Explizite Penrose-Transformation: |m, mquer, k, kquer> -> omega_0^k omega_1^(m-1-k) pi_0^(-1-kquer) "
      "pi_1^(-1-(mquer-1-kquer)) sqrt(kquer! (mquer-1-kquer)!/(k! (m-1-k)!)) ist eine Cech-Klasse in H^1(PT', O(n)) mit "
      f"n = m - mquer - 2 = 2h - 2. Alle 17 Erzeuger der Standardform (z83/z97, CCR, M = 5, {NR} Zustände, {nichttrivial} "
      "nicht null) sind exakt die Twistor-Differentialoperatoren P = -omega^a d/dpi_b, K = pi_b d/domega^a, A = omega^a "
      "d/domega^b, B = -d/dpi_a o pi_b, E = (omega.d_omega - pi.d_pi)/2 auf Klassen modulo Korändern (max. Abweichung "
      f"{max(abw.values()):.0e}). Das Zentrum Z_SU = exp(i pi (E + 2 J3_R)) aus z97 ist genau der Homogenitätscharakter "
      f"i^n (Fehler {zsu_fehler:.0e}). Die Konvention z97 ist damit die Realisierung mit Geraden P(pi), pi = Partnerseite. "
      f"Beide Weyl-Leitern existieren: Weyl-Zählung bei E = N + 1/2 für h = +1/2 {weyl['1/2']}, für h = -1/2 {weyl['-1/2']}",
      max(abw.values()) < 1e-12 and nichttrivial == 17 and zsu_fehler < 1e-12
      and weyl["1/2"] == weyl["-1/2"] == [2, 6, 12, 20])

# --- 2.2 Einschränkung auf die Geraden L_x
def restriktion(n, x, Kmax=5):
    """Rang der Einschränkung H^1(PT', O(n)) -> H^1(L_x, O(n)), omega^A = x^(A0') pi_0 + x^(A1') pi_1, auf den Cech-Klassen
    mit omega-Grad <= Kmax. Gibt (Rang, Zahl der Klassen, dim H^1(L, O(n)))."""
    zeilen = list(range(n + 1, 0))                                  # p0 in [n+1, -1], p1 = n - p0
    zi = {p: i for i, p in enumerate(zeilen)}
    spalten = []
    for k in range(Kmax + 1):
        s_ = n - k
        for a0 in range(k + 1):
            for p0 in range(s_ + 1, 0):
                spalten.append((a0, k - a0, p0, s_ - p0))
    R = np.zeros((len(zeilen), len(spalten)), complex)
    for c, (a0, a1, p0, p1) in enumerate(spalten):
        for i in range(a0 + 1):
            for j in range(a1 + 1):
                e0, e1 = p0 + i + j, p1 + (a0 - i) + (a1 - j)
                if e0 <= -1 and e1 <= -1:
                    R[zi[e0], c] += (math.comb(a0, i)*x[0, 0]**i*x[0, 1]**(a0 - i)
                                     * math.comb(a1, j)*x[1, 0]**j*x[1, 1]**(a1 - j))
    if R.size == 0:
        return 0, len(spalten), len(zeilen)
    sv = np.linalg.svd(R, compute_uv=False)
    return int(np.sum(sv > 1e-9*max(1.0, sv[0]))), len(spalten), len(zeilen)


xs = [np.zeros((2, 2), complex)] + [rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)) for _ in range(3)]
HEL3 = [Fr(k, 2) for k in range(-3, 4)]
rang_p, rang_o, klassen = {}, {}, {}
for h in HEL3:
    rr = [restriktion(n_z97(h), x) for x in xs]
    rang_p[h] = {r[0] for r in rr}
    klassen[h] = rr[0][1]
    rang_o[h] = {r[0] for r in (restriktion(n_z97(-h), x) for x in xs)}       # J-Bild: Originalgerade P(omega)
sichtbar_p = {h: max(rang_p[h]) > 0 for h in HEL3}
sichtbar_o = {h: max(rang_o[h]) > 0 for h in HEL3}
ZAHLEN["Linienrestriktion"] = {"Rang auf P(pi) (Partner)": {fs(h): sorted(v) for h, v in rang_p.items()},
                               "Rang auf P(omega) (Original, J-Bild)": {fs(h): sorted(v) for h, v in rang_o.items()},
                               "Klassen auf PT' (omega-Grad <= 5)": {fs(h): v for h, v in klassen.items()}}
check("[E] Linienrestriktion H^1(PT', O(2h - 2)) -> H^1(L_x, O(2h - 2)) an x = 0 (die Naht) und drei zufälligen komplexen "
      f"x: Rang je h = -3/2..3/2 {[sorted(rang_p[h]) for h in HEL3]}, unabhängig von x und gleich h^1(L, O(2h - 2)), "
      f"surjektiv; die Klassen selbst sind für jedes h vorhanden ({[klassen[h] for h in HEL3]} bis omega-Grad 5). Feld"
      "daten auf der Geraden haben genau die Helizitäten h <= 0 (ein Halbraum): der Spiegel h = +1/2 und auch der "
      "Vektor h = +1 haben auf keiner Geraden Daten (ihre Felder kommen erst aus Normalenableitungen der Ordnung 2h). Das "
      f"J-Bild (Geraden P(omega) der Originalseite) sieht genau h >= 0 ({[sorted(rang_o[h]) for h in HEL3]})",
      all(rang_p[h] == {cech_O(n_z97(h))[1]} for h in HEL3) and all(sichtbar_p[h] == (h <= 0) for h in HEL3)
      and all(sichtbar_o[h] == (h >= 0) for h in HEL3) and all(v > 0 for v in klassen.values()))

# --- 2.3 Artinhalt mit Familienbündel
def leiter_dim(h, E):
    """Dimension der Leiter der Helizität h auf der Energie E (K-Typ (j1, j2), j1 - j2 = h, j1 + j2 + 1 = E)."""
    j1, j2 = (E - 1 + h)/2, (E - 1 - h)/2
    return int((2*j1 + 1)*(2*j2 + 1)) if j1 >= 0 and j2 >= 0 and (2*j1) == int(2*j1) else 0


ENERGIEN = [Fr(k, 2) for k in range(2, 13)]
voll = {E: {h: leiter_dim(h, E) for h in HEL3} for E in ENERGIEN}
cnt_cech = {E: {h: int(np.sum(np.isclose(EV, float(E)) & np.isclose(hel, float(h)))) for h in HEL3}
            for E in ENERGIEN if E <= 3}
leiter_ok = all(cnt_cech[E][h] == voll[E][h] for E in cnt_cech for h in HEL3)


def inhalt(F, basis):
    """Helizitätsinhalt von F (x) (Leitern der Basis-Helizitäten): Summand O(a) verschiebt h um a/2 (K1)."""
    out = {}
    for h0, mult in basis.items():
        for a in F:
            out[h0 + Fr(a, 2)] = out.get(h0 + Fr(a, 2), 0) + mult
    return out


basis_weyl = {Fr(1, 2): 1, Fr(-1, 2): 1}
basis_voll = {h: 1 for h in [Fr(k, 2) for k in range(-12, 13)]}
inh_weyl = {name: inhalt(F, basis_weyl) for name, F in (("E_-", EM), ("4 O", [0]*4), ("E_+ = E_-^v", [1, 0, 0, 0]))}
inh_voll = {name: inhalt(F, basis_voll) for name, F in (("E_-", EM), ("4 O", [0]*4), ("E_+ = E_-^v", [1, 0, 0, 0]))}
netto = {name: [sum(inh[h]*leiter_dim(h, E) for h in (Fr(1, 2),)) - inh[Fr(-1, 2)]*leiter_dim(Fr(-1, 2), E)
                for E in ENERGIEN] for name, inh in inh_voll.items()}
mult_mitte = {name: sorted({inh[h] for h in HEL3}) for name, inh in inh_voll.items()}
ZAHLEN["Artinhalt"] = {"Weyl-Basis, K1": {k: {fs(h): v for h, v in sorted(d.items())} for k, d in inh_weyl.items()},
                       "volle Standardform, K1: Vielfachheiten |h| <= 3/2": mult_mitte,
                       "netto Weyl-Index je Energie (volle Standardform)": netto}
check("[negativ, E] Artinhalt (Leitern auf PT', K-Typen aus der Cech-Zählung bis E = 3 bestätigt): mit E_- verdrillt wird "
      f"ein Weyl-Paar zu {{h: Vielfachheit}} = {{{', '.join(f'{fs(h)}: {v}' for h, v in sorted(inh_weyl['E_-'].items()))}}}"
      ": drei Familien je Helizität +-1/2 und eine um -1/2 verschobene Kopie mit ganzer Helizität (der O(-1)-Summand, "
      "unter P_GSO im Sektor R = +-i verworfen). Auf der vollen Standardform ist die Verdrillung eine Umbenennung: jede "
      f"Helizität kommt rang F mal vor ({mult_mitte}), der Netto-Weyl-Index ist auf jeder Energie 0 ({netto['E_-'][:6]}...). "
      "Im Produkt K2 ist der Inhalt chi(F) mal die Standardform, ebenso vektorartig. Familienbündel und triviale Kontrollen "
      "unterscheiden sich nur in der Vielfachheit",
      leiter_ok and inh_weyl["E_-"] == {Fr(0): 1, Fr(1, 2): 3, Fr(-1): 1, Fr(-1, 2): 3}
      and all(m_ == [len(F)] for m_, F in ((mult_mitte["E_-"], EM), (mult_mitte["4 O"], [0]*4)))
      and all(all(v == 0 for v in netto[k]) for k in netto))

# --- 2.4 Sektoren und Linienregeln
# Z4-Sektoren (z97, z102): R = 1 neutral ((45,1) + (1,15); Skalarleiter und Vektoren, Familienbündel O), R = i (16, E_-),
# R = -1 ((10,6), Skalare), R = -i (16q, konjugierter Familienraum H^0(E_-)^* = H^1(E_-^v (x) K)). Vielfachheit im Produkt K2.
GSO = lambda R, h: (h.denominator == 1) == (R in (1, -1))                 # P_GSO: ganz mit R = +-1, halb mit R = +-i
MULT = {1: 1, 1j: chiF["E_-"], -1: coh([-1, -1, -1, 0, 0, 0], 0)[0], -1j: -chiF["E_-^v (x) K"]}
SEITEN = ("o", "p")                                                       # Originalgerade P(omega), Partnergerade P(pi)
TEILM = [frozenset(s) for s in ((), ("o",), ("p",), ("o", "p"))]
Jseite = lambda S: frozenset({"o": "p", "p": "o"}[x] for x in S)
RS = (1, 1j, -1, -1j)
RNAME = {1: "R = 1", 1j: "R = i", -1: "R = -1", -1j: "R = -i"}
Rinv = {1: 1, 1j: -1j, -1: -1, -1j: 1j}
sicht_rest = {"o": sichtbar_o, "p": sichtbar_p}
sicht_index = {"o": {h: ind(cech_O(n_pen(h))) != 0 for h in HEL3}, "p": {h: ind(cech_O(n_z97(h))) != 0 for h in HEL3}}


def bewerte(regel, sicht, mult16=None):
    m16 = MULT[1j] if mult16 is None else mult16
    wz = lambda R, h, m_: m_ if GSO(R, h) and any(sicht[s][h] for s in regel[R]) else 0
    N16 = wz(1j, Fr(1, 2), m16) - wz(1j, Fr(-1, 2), m16)
    N16q = wz(-1j, Fr(-1, 2), m16) - wz(-1j, Fr(1, 2), m16)
    spiegel = min(wz(1j, Fr(1, 2), m16), wz(1j, Fr(-1, 2), m16))
    return {"N16": N16, "N16q": N16q, "chiral": N16 != 0 and N16 == N16q and spiegel == 0,
            "skalar": wz(1, Fr(0), 1) > 0, "vektor": wz(1, Fr(1), 1) > 0 and wz(1, Fr(-1), 1) > 0,
            "alle": all(len(regel[R]) > 0 for R in RS),
            "J": all(regel[Rinv[R]] == Jseite(regel[R]) for R in RS),
            "universell": len({regel[R] for R in RS}) == 1}


REGELN = [dict(zip(RS, t)) for t in itertools.product(TEILM, repeat=4)]
erg = {}
for krit, sicht in (("Restriktion", sicht_rest), ("Index", sicht_index)):
    ev_ = [(r, bewerte(r, sicht)) for r in REGELN]
    volle = [(r, e) for r, e in ev_ if e["chiral"] and e["skalar"] and e["vektor"] and e["alle"]]
    erg[krit] = {"Regeln": len(REGELN), "chiral": sum(e["chiral"] for _, e in ev_),
                 "chiral, Skalar, beide Vektoren, alle koppeln": len(volle),
                 "davon J-kovariant": sum(e["J"] for _, e in volle),
                 "universelle Regeln chiral": sum(e["chiral"] for _, e in ev_ if e["universell"] and e["alle"]),
                 "J-kovariante Lösungen": [{RNAME[R]: sorted(r[R]) for R in RS} | {"N16": e["N16"]}
                                           for r, e in volle if e["J"]]}
uni_p = bewerte({R: frozenset({"p"}) for R in RS}, sicht_rest)
uni_p_idx = bewerte({R: frozenset({"p"}) for R in RS}, sicht_index)
ZAHLEN["Linienregeln"] = erg
ZAHLEN["universelle Partnergerade"] = {"Restriktion": uni_p, "Index": uni_p_idx}
check("[negativ, E] Universelle Nahtgerade (alle Sektoren über dieselbe Gerade P(pi), P_GSO, Familienvielfachheit aus "
      f"K2): N16 = {uni_p['N16']}, N16q = {uni_p['N16q']}: 16 und 16q liegen bei derselben Helizität, also Eichchiralität "
      "0 (CPT verlangt N16q = N16); unter dem Restriktionskriterium sieht die Naht den R-neutralen Vektor nur mit h = -1 "
      f"(beide Vektoren: {uni_p['vektor']}), unter dem Indexkriterium (h^0 - h^1 != 0) beide ({uni_p_idx['vektor']}), aber "
      f"ebenfalls ohne Eichchiralität ({uni_p_idx['chiral']}). Über alle 4 universellen Regeln: chiral "
      f"{erg['Restriktion']['universelle Regeln chiral']} bzw. {erg['Index']['universelle Regeln chiral']}. Die "
      "Linienchiralität ist eine reine Helizitätsauswahl, gleich für 16 und 16q; die Halbraumauswahl, die h = +-1/2 trennt, "
      "trennt auch h = +-1 (Indexform des Leiterquotienten), der Skalar liegt auf dem Rand und bleibt",
      uni_p["N16"] == -3 and uni_p["N16q"] == 3 and not uni_p["chiral"] and uni_p["skalar"] and not uni_p["vektor"]
      and uni_p_idx["vektor"] and not uni_p_idx["chiral"] and erg["Restriktion"]["universelle Regeln chiral"] == 0
      and erg["Index"]["universelle Regeln chiral"] == 0)

prj = [{R: frozenset(x) for R, x in zip(RS, t)} for t in ((("o", "p"), ("o",), ("o", "p"), ("p",)),
                                                         (("o", "p"), ("p",), ("o", "p"), ("o",)))]
e_prj = [bewerte(r, sicht_rest) for r in prj]
lsg = erg["Restriktion"]["J-kovariante Lösungen"]
check("[E, bedingt] Alle 256 Linienregeln (jeder Z4-Sektor koppelt über die Originalgerade P(omega), die Partnergerade "
      "P(pi), beide oder keine), Feldsichtbarkeit aus den gerechneten Restriktionsrängen, P_GSO, Vielfachheit aus K2: "
      f"chiral (16 bei einer Helizität, Spiegel 0, 16q CPT-gegengleich), mit R-neutraler Skalarleiter und beiden "
      f"Vektorhelizitäten, alle Sektoren koppelnd: {erg['Restriktion']['chiral, Skalar, beide Vektoren, alle koppeln']} "
      f"Regeln (Indexkriterium: {erg['Index']['chiral, Skalar, beide Vektoren, alle koppeln']}); J-kovariant (J: R -> R^-1, "
      f"Original <-> Partner) davon genau {erg['Restriktion']['davon J-kovariant']} (Index: "
      f"{erg['Index']['davon J-kovariant']}): R = +-1 über beide Geraden, R = i und R = -i über entgegengesetzte Geraden "
      f"(P_RJ und seine Gegenorientierung), N16 = {[e['N16'] for e in e_prj]}. J-Kovarianz erzwingt für die J-selbstgepaarten "
      "Sektoren R = +-1 beide Geraden; einseitig dürfen nur die Uhrsektoren R = +-i sein, und genau sie sind die Fermionen "
      "(R^2 = (-1)^F)",
      erg["Restriktion"]["davon J-kovariant"] == 2 and erg["Index"]["davon J-kovariant"] == 2
      and erg["Restriktion"]["chiral, Skalar, beide Vektoren, alle koppeln"] == 6
      and sorted(e["N16"] for e in e_prj) == [-3, 3] and all(e["chiral"] and e["skalar"] and e["vektor"] and e["J"]
                                                              for e in e_prj)
      and sorted(x["N16"] for x in lsg) == [-3, 3])

# --- 2.5 Nullmodenzählung unter P_RJ und Negativkontrollen
P_RJ = prj[0]
o_line = {h: kuenneth(EM, n_pen(h)) for h in (Fr(1, 2), Fr(-1, 2))}       # 16 auf der Originalgeraden, Produkt K2
kontr = {"triviales 4 O": bewerte(P_RJ, sicht_rest, mult16=chiF["4 O"])["N16"],
         "Spiegelfamilie O(-4) + 3 O": bewerte(P_RJ, sicht_rest, mult16=chiF["Spiegelfamilie O(-4) + 3 O"])["N16"],
         "Gegenorientierung (P_RJ gespiegelt)": e_prj[1]["N16"],
         "Serre-Dual E_-^v (x) K": bewerte(P_RJ, sicht_rest, mult16=chiF["E_-^v (x) K"])["N16"]}
k1_spiegel = coh(EM, n_pen(Fr(-1, 2)))                                    # Verdrillung K1 am Spiegelslot (Originalgerade)
k1_chiral = coh(EM, n_pen(Fr(1, 2)))
ZAHLEN["P_RJ"] = {"Regel": {RNAME[R]: sorted(v) for R, v in P_RJ.items()}, "Bewertung": e_prj[0],
                  "Familienvielfachheit je Sektor (K2)": {RNAME[R]: v for R, v in MULT.items()},
                  "16 auf der Originalgeraden, K2 (H0, H1, H2)": {fs(h): v for h, v in o_line.items()},
                  "Kontrollen N16": kontr, "K1 statt K2 (h0, h1) Spiegel/chiral": (k1_spiegel, k1_chiral)}
check("[E unter P_RJ] Nullmodenzählung je 16 unter P_RJ (R = i über die Originalgerade, Penrose-Konvention O(-2h - 2), "
      f"Produkt K2): chiraler Slot h = +1/2 (H^0, H^1, H^2) = {o_line[Fr(1, 2)]} (drei Weyl, Netto-Index "
      f"{e_prj[0]['N16']} = Familienindex 3 mal Twistorindex), Spiegelslot h = -1/2 {o_line[Fr(-1, 2)]} (Index 0, keine "
      "Mode); die 16q CPT-gegengleich über die Partnergerade; R-neutrale Skalarleiter und beide Vektorhelizitäten über "
      f"beide Geraden. Kontrollen: triviales 4 O gibt {kontr['triviales 4 O']} statt 3, die Spiegelfamilie "
      f"{kontr['Spiegelfamilie O(-4) + 3 O']}, die Gegenorientierung {kontr['Gegenorientierung (P_RJ gespiegelt)']} "
      f"(Index kippt mit der Orientierung), das Serre-Dual {kontr['Serre-Dual E_-^v (x) K']}. In der Verdrillung K1 statt "
      f"des Produkts hat E_- am Spiegelslot (h^0, h^1) = {k1_spiegel} (eine Mode aus O(-1)) und am chiralen {k1_chiral}: "
      "Spiegelindex 0 gibt nur K2",
      o_line[Fr(1, 2)] == (0, 6, 0) and o_line[Fr(-1, 2)] == (0, 0, 0) and e_prj[0]["N16"] == 3
      and kontr == {"triviales 4 O": 4, "Spiegelfamilie O(-4) + 3 O": 0, "Gegenorientierung (P_RJ gespiegelt)": -3,
                    "Serre-Dual E_-^v (x) K": -3}
      and k1_spiegel == (0, 1) and k1_chiral == (0, 9))
befund("Teil 2: Die radiale Standardform ist exakt die Cech-Kohomologie von PT' (Penrose-Transformation explizit, Z_SU = "
       "Homogenitätscharakter); ihr Artinhalt bleibt mit jedem Familienbündel vektorartig. Chiral ist nur, was die Naht "
       "über ihre Gerade sieht: ein Halbraum der Helizitäten, mit universeller Gerade gleich für 16 und 16q "
       "(Eichchiralität 0). Chiralität mit Skalarleiter und beiden Vektorhelizitäten verlangt eine sektorabhängige "
       "Geradenwahl; J-kovariant ist genau P_RJ (Uhrsektoren R = +-i auf entgegengesetzten Tomita-Geraden, R = +-1 auf "
       "beiden), dann Netto-Index 3 je 16 mit azyklischem Spiegelslot", True)

# =================================================================== Teil 3
say("\nTeil 3: Bilanz und Urteil für T4")
STATUS = {
    "Index auf der Twistorgeraden": ("E", "2h - 1, einzige Nullstelle h = +1/2 (O(-1)), Skalar und Vektoren nicht azyklisch"),
    "Penrose-Transformation der Leiter": ("E", "Standardform = Cech-H^1(PT', O(m - mquer - 2)), 17 Erzeuger exakt, "
                                               "Z_SU = i^n; Konvention z97 = Geraden P(pi) der Partnerseite"),
    "Familienindex in der Verdrillung K1": ("E (negativ)", "Nullmoden nur am Vektorslot, Asymmetrie 4h rang F "
                                                          "familienblind, Spiegelmode aus O(-1)"),
    "Familienindex im Produkt K2": ("E", "Index = chi(E_-) (2h - 1), drei Weyl am chiralen, null am Spiegelslot; "
                                        "h^1(E_-) = 0 nötig, 3 O leistet dasselbe"),
    "Artinhalt (Spektrum)": ("E (negativ)", "vektorartig mit jedem Familienbündel; Spiegelleiter existiert"),
    "universelle Nahtgerade": ("E (negativ)", "Halbraum h <= 0: 16 und 16q gleich (Eichchiralität 0), "
                                              "Vektor nur selbstdual gesehen"),
    "Linienregel P_RJ": ("C (neue Prämisse)", "J-kovariant, einzig (bis Orientierung) mit Chiralität, Skalarleiter, "
                                              "beiden Vektoren; Netto-Index 3 je 16"),
    "Kontrollen": ("E", "triviales 4 O: 4; Spiegelfamilie: 0; Gegenorientierung: -3; K1 statt K2: Spiegelmode"),
    "Chiralität als Index gesamt": ("bedingt", "nativ: Twistor-Azyklizität und Familienindex; Prämissen T, P_GSO, "
                                               "Produktkopplung K2, P_RJ; als Spektralaussage nein"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
ZAHLEN["T5_urteil_kurz"] = {
    "Urteil": "bedingt (Prämissen T, P_GSO, K2, P_RJ); nicht nativ; als Spektralaussage gescoptes Negativ",
    "T4": "Indexroute trennt die Weyl-Helizitäten ohne die Skalarleiter zu verlieren; Eichchiralität nur mit "
          "sektorabhängiger Geradenwahl P_RJ, das Spiegelpaar bleibt als nahtsterile Leiter; offen: wählt die Uhr R die "
          "Tomita-Seite nativ?",
    "K16": "Linienindex holomorph und exakt, kein Spiegelproblem auf der Geraden; K16 bleibt L",
    "chiral_matter": "neuer Wert 'Twistorlinienindex mit P_RJ' (bedingt) neben P_A3; universelle Gerade: nicht chiral"}
urteil = ("T5 ist bedingt, nicht nativ, und als Aussage über das Spektrum ein gescoptes Negativ. (1) Die Twistorgerade "
          "ist für sich chiral: Der Dolbeault-Index von O(2h - 2) ist 2h - 1 und verschwindet nur bei h = +1/2, denn O(-1) "
          "ist das einzige azyklische Bündel auf CP^1 (Cech und analytischer d-quer-Operator gleich). Die Nullstelle liegt "
          "an einer Weyl-Helizität; Skalar (h^1 = 1) und beide Vektorslots (h^0 = 1, h^1 = 3) behalten Kohomologie. In der "
          "Penrose-Konvention O(-2h - 2) liegt sie bei h = -1/2; die Konventionen sind J-Bilder (Geraden der Partner- bzw. "
          "Originalseite). (2) Der Familienindex 3 ist helizitätsblind. Verdrillt (K1: Naht = Twistorgerade, E_- (x) "
          "O(2h - 2)) liegen seine Nullmoden nur am Vektorslot h = +1 (Homogenität 0, passend zu Familien = Lambda^2_+), "
          "die Asymmetrie ind(h) - ind(-h) = 4h rang F hängt nicht vom Grad ab (E_-, E_+ und 4 O gleich), und E_- hat am "
          "Spiegelslot eine Mode aus O(-1) (chi = 4k + 3 nie 0). Erst im Produkt (K2: Familiennullmoden mal Twistorgerade) "
          "gilt Index = chi(E_-) (2h - 1): drei Weyl am chiralen Slot (h^1 = 6, h^2 = 0), am Spiegelslot verschwindet alle "
          "Kohomologie. Dafür ist die Spiegelfreiheit h^1(E_-) = 0 nötig (O(-4) + 3 O gibt Index 0), das triviale 3 O "
          "leistet aber dasselbe. (3) Explizite Penrose-Transformation: Die radiale Standardform ist exakt Cech-H^1(PT', "
          "O(m - mquer - 2)); alle 17 Erzeuger sind die Twistor-Differentialoperatoren auf Klassen modulo Korändern, und "
          "Z_SU aus z97 ist der Homogenitätscharakter i^n. Beide Weyl-Leitern existieren (n(n+1) beidseitig), mit jedem "
          "Familienbündel bleibt der Artinhalt vektorartig (Netto-Weyl-Index 0 auf jeder Energie). Chiral ist nur, was eine "
          "Gerade sieht: Die Restriktion hat Rang h^1(L, O(2h - 2)), also Felddaten genau für h <= 0, ein Halbraum. Er "
          "trennt mit h = +-1/2 auch h = +-1 (Indexform des Leiterquotienten aus z102); der Skalar liegt auf dem Rand und "
          "bleibt, von den Vektoren fehlt nur eine Helizität (beim Zentrumskleber aus z102 fiel die Skalarleiter weg, mit "
          "dem Viertelwinkel die R-neutralen Vektoren ganz). (4) Mit einer universellen Nahtgeraden liegen 16 und 16q bei "
          "derselben Helizität: Eichchiralität 0, und der R-neutrale Vektor wird nur selbstdual gesehen. Von allen 256 "
          "Linienregeln (Sektor -> Original-, Partnergerade, beide, keine) geben Chiralität, R-neutrale Skalarleiter und "
          "beide Vektorhelizitäten zusammen 6, wenn jeder Sektor an die Naht koppelt; J-kovariant sind genau 2, die Regel P_RJ und ihre Gegenorientierung: die "
          "Uhrsektoren R = +-i (die Fermionen, R^2 = (-1)^F) koppeln über entgegengesetzte Tomita-Geraden, die "
          "J-selbstgepaarten Sektoren R = +-1 über beide. Unter P_RJ, P_GSO und K2 hat jede 16 Netto-Index 3 (Familienindex "
          "mal Twistorindex), der Spiegelslot ist azyklisch, Skalarleiter und beide Vektorhelizitäten bleiben. Kontrollen: "
          "triviales 4 O gibt 4, die Spiegelfamilie 0, die Gegenorientierung -3, die Verdrillung K1 statt K2 eine "
          "Spiegelmode. Folge für T4: Die Indexroute löst den Konflikt aus z102 auf der Geraden (sie trennt die "
          "Weyl-Helizitäten, ohne Skalare zu verlieren), gibt Eichchiralität aber nur mit der neuen, J-kovarianten Prämisse "
          "P_RJ; das Spiegelpaar (16 bei h = -1/2, 16q bei h = +1/2) bleibt als Leiter im Spektrum und ist nur nahtsteril "
          "(am Spiegelslot verschwindet alle Kohomologie der Sektorgeraden). Offen ist, ob die Uhr R die Tomita-Seite "
          "nativ wählt; das ist die reduzierte Form von T4. Folge für K16: Auf der Geraden gibt es kein Spiegelproblem, der "
          "Index ist holomorph und exakt; die Frage verlagert sich in die Kopplung des Spiegels im Spektrum, K16 bleibt L. "
          "Folge für chiral_matter: Die Szenario-Suche kann neben P_A3 den Wert 'Twistorlinienindex mit P_RJ' führen, bedingt; "
          "mit universeller Geraden ist er nicht chiral. Nichts ist befördert; 3+1-Aussagen bleiben Kinematik der "
          "Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z114_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
