Gebiete: Quantenstatistik und Randfermionen (GSO, Spinstrukturen) · Zufallsmatrizen und Hubbard–Stratonovich · Kobordismus und Anomalien · allgemeine Relativitätstheorie in Spinorform (Petrov) · Kristallographie (Oktaedergruppe) · algebraische Quantenfeldtheorie (modulare Inklusionen)

# Exploration 23. November 2026: Fockspur statt Felddeterminante

Anlass: zwei Berichte vom 23. November (F1 Quellenphase und Ladungszerlegung des Ursprungskanals)
und die neuen Atlasbefunde vom selben Tag (origin/main, Commits 38362a47e, 459ae74d1, a9d827064:
fünf neue Problemfelder, 27 neue Voraussetzungen). Der Nutzer bat, beides zu prüfen und daraus
die fundamentale Lösung zu suchen. Das ist Explorationsmodus mit ausdrücklich vorgegebenem
Material; PLAN, RESULTS und progress wurden nicht gelesen.

Vielfaltsprüfung: Auf main liegt nur `2026-11-20-kleinstes-verbundsystem` (Operatoralgebren,
Quanteninformation, Codes, Kobordismus, 𝔽₂, Extremalprinzip, Automaten, Kausalmengen). Neu hier:
Randfermionen/GSO, Hubbard–Stratonovich, Petrov-Spinoren, Kristallographie, modulare Inklusionen.
Der perfekte Fünf-Qubit-Code (16 = 1 + 3·5) wurde zuerst ebenfalls erwogen; er ist dort als A2
bereits gerechnet und wird hier nicht wiederholt.

## Ansätze

| Nr | Gebiet | Grundprinzip | Soll liefern | Freie Parameter | Erste freie Wahl | Gerechnet |
|---|---|---|---|---|---|---|
| F1 | Quantenstatistik, Randfermionen | Die positive Quelle ist eine Zustandssumme über den Fockraum von N Randmoden, nicht eine Determinante über Feldern. Gerader Teil = Summe beider Spinstrukturen des Kreises. | Welche Darstellung (Frage aus Bericht 1); 3+2; Auswahl einer der beiden Bahnen ±(2,2,2,−3,−3) | ρ (Kopplung), T | Spurfreie su(N) Hubbard–Stratonovich Quelle | `s2_fockspur_auswahl.py` |
| F2 | Zufallsmatrizen | Vorstufe mit reellen Moden: so(2N) Quelle auf dem Spinorfockraum | komplexe Struktur J (SO(2N) → U(N)) | ρ | Quelle in so(2N) | `s3_majorana_stufen.py` |
| F3 | Kobordismus | ℤ₁₆ Anomalie mit μ₄ in d = 4 (A3 vom 20.11.) und Lesart „Generation = gerader Fockraum von N Moden“ | N = 5, dann über F1 3+2 | keine | Minimalität der Generation; d = 4 | `s4_kette_anomalie_fock.py` |
| F4 | Allgemeine Relativität, Spinoren | Träger = Sym⁴ℂ² = Weylspinor Ψ_ABCD (Newman–Penrose Ψ₀…Ψ₄); 3+2 = Parität des Spingewichts unter der π-Drehung; vier Marken = Hauptnullrichtungen; j = 1728 | P2 und 3+1 aus einem Objekt | keine | Lesart Träger = Krümmung | nein |
| F5 | Kristallographie | Sechs Ereignislabels des Ursprungskanals = ±e_i, Matching = Raumspiegelung, Zentralisator B₃ ≅ O_h (Ordnung 48); 3+2 = T₁ᵤ ⊕ E_g; 5 = 2d − 1 | Verbindung Träger ↔ drei Raumachsen | keine | Lesart Label = orientierte Achse | nein |
| F6 | Algebraische QFT | Halbseitige modulare Inklusionen der markierten Quellalgebra erzeugen Translationen (Vorschlag aus Bericht 2) | Zeit, Raumrichtung | keine | welche Inklusion | nein |

## Begründung der Auswahl zum Rechnen

F1 beantwortet eine in Bericht 1 ausdrücklich offene Frage („Welche Darstellung?“) mit einer
kleinsten entscheidenden Rechnung. F2 und F3 prüfen, wie weit dieselbe Quelle nach vorn und
hinten trägt. F4 und F5 sind Lesarten mit bekannten Einwänden (F4: Spin–Statistik, die 16 trüge
Lorentzspins 0 bis 3, und der Atlasbefund `pauli-lorentz-kinematics-and-spin-two-hypercharge-exclusion`;
F5: dieselbe Doppellesart Farbe/Raum, die A1 am 20.11. als stärksten Einwand hatte). F6 braucht
erst eine konkrete Inklusion aus der Quelle.

## Kleinste entscheidende Tests

- F1: Ist die vollständige antiperiodische Paarquelle aus Bericht 1 eine Fockspur? Wählt die Spur
  über die 16 dieselbe Form wie die über 32, und trennt sie die beiden Bahnen? Gegenprobe:
  Felddeterminante über 16 Gewichten. Abbruch, wenn die Spur über 16 nicht 3+2 wählt.
- F2: Liefert die so(10) Fockquelle allein schon 3+2? Abbruch, wenn ja (dann wäre F1 überflüssig).
- F3: Ist N = 5 die kleinste Modenzahl mit einheitlicher ℤ₄ Ladung und anomaliefreier Generation?
- Basislinie für F1: welche Teilung wählt dieselbe Spur bei N = 2 … 9 Moden?
