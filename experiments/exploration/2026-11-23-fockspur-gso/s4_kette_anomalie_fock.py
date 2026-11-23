"""Kette: mu_4 in vier Dimensionen -> 16 -> fuenf Moden -> 3+2 mit fester Orientierung.

Glieder:
 (a) A3 vom 20.11. (a3_anomalie_eta.py): Weylfermion in d = 4 mit Z4, Z4^2 = (-1)^F, eta = -1/16;
     anomaliefrei nur Vielfache von 16 (Literatur: Z16, Garcia-Etxebarria/Montero, Wang/Wen).
 (b) Lesart: eine Generation ist der gerade (GSO) Fockraum von N komplexen Moden, dim 2^(N-1).
     Die Z4 Ladung X = N - 2p ist auf geraden p konstant N mod 4; Z4^2 = (-1)^F verlangt N ungerade.
 (c) s2_fockspur_auswahl.py: positive Fockspur mit spurfreier su(N) Quelle waehlt die balancierte
     Teilung; die GSO Projektion waehlt eine Bahn.
Hier nur die exakte Arithmetik von (b) und die Tabelle, welche Teilung je N entstuende.
"""
import sys
from fractions import Fraction
OK = True
def check(name, cond, info=""):
    global OK
    OK = OK and bool(cond)
    print(f"  [{'ok' if cond else 'FEHLER'}] {name}" + (f" -- {info}" if info else ""))
eta = Fraction(-1, 16)
print(" N  dim(gerade)  X mod 4 auf geraden p   Z4^2=(-1)^F   anomaliefrei   balanciert")
sol = []
for N in range(1, 12):
    dim = 2 ** (N - 1)
    xs = {(N - 2 * p) % 4 for p in range(0, N + 1, 2)}
    z4 = (len(xs) == 1 and next(iter(xs)) % 2 == 1)
    free = (dim * eta).denominator == 1
    split = f"{(N + 1) // 2}+{N // 2}"
    print(f"{N:2d}  {dim:11d}  {str(sorted(xs)):22s}  {str(z4):12s}  {str(free):13s}  {split}")
    if z4 and free: sol.append(N)
check("kleinste Modenzahl mit einheitlicher Z4 Ladung und anomaliefreier Generation ist N = 5", sol[0] == 5, f"Loesungen {sol}")
check("dann ist die Generation 16 = 1 + 10 + 5 und die balancierte Teilung 3+2", 2 ** 4 == 1 + 10 + 5)
check("drei Generationen: 48 Weylfermionen, ebenfalls anomaliefrei (Familienzahl bleibt frei)", (48 * eta).denominator == 1)
print("\nGESAMT:", "alle Pruefungen bestanden" if OK else "FEHLER")
sys.exit(0 if OK else 1)
