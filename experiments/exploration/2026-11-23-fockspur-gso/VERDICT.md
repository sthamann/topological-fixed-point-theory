# Urteil

Gerechnet am 23. November 2026 mit `.venv/bin/python` (sympy 1.14, numpy 2.4, scipy 1.17).
Vier Skripte, alle mit Rückgabewert 0; Ausgaben in `ausgabe_s*.txt`. PLAN, progress, Ledger und
Atlas sind unberührt.

## Übersicht

| Nr | Urteil | Stufe | Erste freie Wahl |
|---|---|---|---|
| Nachprüfung der Berichte | bestätigt, eine Präzisierung | – | – |
| F1 Fockspur und GSO | **offen**, exaktes Teilergebnis | 1 gegeben N = 5; an der Basislinie | spurfreie su(5) Quelle, Lesart „16 = Fockraum“ |
| F2 Majoranavorstufe | **offen** | 1 (komplexe Struktur J) | Quelle in so(2N) |
| F3 Kette μ₄ → 16 → 5 → 3+2 | **Kandidat** als Reduktion | 2 für P2 unter benannter Lesart | d = 4 mit μ₄; Minimalität der Generation |
| F4, F5, F6 | **offen**, nicht gerechnet | – | siehe IDEAS.md |

Keine vollständige Lösung. Was sich ändert: P2 hängt nach dieser Runde an weniger Eingaben als
vorher, und die in Bericht 1 offene Darstellungsfrage hat eine exakte Antwort.

## 1. Nachprüfung der beiden Berichte (`s1_pruefe_bericht.py`)

Unabhängig nachgerechnet und bestätigt:

- **Ursprungskanal.** Aus den 15 Transpositionen mit Gewichten (32, 5, −4)/648 folgt das Spektrum
  1(1), 3/4(5), 11/36(9), 1/4(10). Ladungssektoren von ad(6Y): Blocknormen 317/108 (0 → 0),
  53/54 (±5 → ±5), 41/216 (±5 → ∓5), alle Mischungen mit dem Nullsektor exakt null.
  Singulärwerte des Übergangs (1/4)×3 und (1/36)×3. ‖[Φ, ad_Y]‖² = 1025/972,
  ‖[Φ, Ad_R]‖² = 41/27, Verhältnis 36/25, beide aus demselben Block. Bericht 2, Abschnitt 1 stimmt.
- **Darstellungsspuren.** Tr₁₀, Tr₁₆, Tr₂₄, Tr₂₄₈ wie angegeben, einschließlich Tr₂₄₈Σ⁴ = 36p₂²
  und Tr₂₄₈Σ⁶ = 30p₂³. Quartikauswahl: 5 und 24 bevorzugen 3+2, 10 und 16 bevorzugen 4+1.
- **Determinantenquelle.** ρ tr Σ² − tr log(1+Σ²) hat für ρ ∈ {0,4; 0,5; 0,7; 0,9; 0,99} globale
  Minima genau auf 3+2 mit der angegebenen v²-Formel; bei ρ > 1 der Ursprung.
- **E₈ Fenster.** ρ₈(10) = 3075/247 und ρ₈(15) = 52660/5643 exakt. Die Cartan-Hessian bei q = 10
  hat genau die Eigenwerte 0,408990 (A₁), 2,659621 (A₄, vierfach), 3,181899 (A₂, zweifach),
  34,081245 (radial). Positiv bei q = 10; 12,5; 15; 20, negativ bei 25. 28 Wurzeln stehen senkrecht
  auf h, der Zentralisator hat Dimension 36.
- **Freie Stromdarstellung.** ‖rΩ‖² = 3 und b₃b₁rΩ = Ω sind elementar richtig.

**Präzisierung.** Die 3+2 Auswahl der Determinantenquelle gilt numerisch auch unterhalb von
ρ = 1/3 (geprüft bei 0,2; 0,3; 0,34 mit Übereinstimmung zur 3+2-Formel auf acht Stellen).
Die Schranke 1/3 ist eine Grenze des Beweises, keine der Aussage.

**Nicht nachgerechnet** (fehlende Archivdaten oder Aufwand): die vollständigen Wardidentitäten
von W für SU(3), SU(2), U(5), Spin(10), SU(4) (480 Kopplungen; der Atlasbefund
`rr-w-fixed-hypercharge-ten-doublets-20261123` deckt nur Y ab), die globale E₈ Suche, der
Starrheitssatz für D, die F_a Nichteindeutigkeit und die 2π-Gegenprobe zu R066.

## 2. F1: Fockspur statt Felddeterminante (`s2_fockspur_auswahl.py`)

**Exakt gezeigt.**

1. Für einen Einteilchenoperator Σ auf ℂ⁵ mit Eigenwerten 2yᵢ/T, Σyᵢ = 0:
   Tr_{Λℂ⁵} e^{−TΣ̂} = 2⁵ Π cosh yᵢ, Tr_{Λℂ⁵}(−1)^F e^{−TΣ̂} = 2⁵ Π sinh yᵢ,
   also **Tr₁₆ = 16 Π cosh yᵢ (1 + Π tanh yᵢ) > 0** auf Λ^gerade ℂ⁵.
2. Das vollständige antiperiodische Paarprodukt der Quelle aus Bericht 1,
   Π_n det(1 + g²Σ²/λ_n²) mit λ_n = (2π/T)(n + ½), ist det cosh(gTΣ/2). Die Quelle
   det(I + Σ²) aus Bericht 1 ist die Einmoden-Abschneidung dieser Fockspur. Die Spiegelungs-
   kopplung τ ⊗ Σ und der antiperiodische Kreis aus P1 erzeugen also genau die Zustandssumme
   der fünf Moden.
3. log(1 + Π tanh yᵢ) = y₁⋯y₅ + O(y⁷). Bis zur vierten Ordnung sind die Spuren über 32 und über
   16 gleich; die Projektion auf die 16 setzt erst mit **det Σ in fünfter Ordnung** ein.

**Numerisch (globale Minimierung, 300 Starts je Punkt, ρ = 0,12 … 0,49).**

| Quelle | Auswahl |
|---|---|
| Tr₃₂ (antiperiodisch, ganzer Fockraum) | 3+2, beide Bahnen gleichwertig |
| **Tr₁₆ (gerader Fockraum = TFPT 16)** | **3+2, immer dieselbe Bahn: Σ ∝ +(2,2,2,−3,−3) = −6Y** |
| Tr₁₆′ (ungerader Teil) | 3+2, die Gegenbahn |
| Felddeterminante über die 16 Gewichte (cosh oder log(1+w²)) | 4+1 |

**Antwort auf „Welche Darstellung?“** Wird die 16 so gelesen, wie TFPT sie definiert, als
Fockraum der fünf Moden, dann ist die positive Quelle eine Spur und hängt nur von det₅ ab: 3+2.
Werden die 16 als unabhängige Felder integriert, entsteht det₁₆ und 4+1. Die Umkehr aus der
Tabelle in Bericht 1 ist also kein Widerspruch, sondern unterscheidet zwei physische Lesarten
der 16. Zusätzlich beseitigt die Projektion auf die 16 die Restentartung der beiden
Grassmannbahnen, die Bericht 1 ausdrücklich offen ließ: Die Chiralität des Focksektors legt die
Orientierung von Σ relativ zu Y fest (Teilchen-Loch-Tausch vertauscht 16 ↔ 16′ und Σ ↔ −Σ).

**Basislinie.** Dieselbe Spur wählt für jede Modenzahl N die balancierteste Zweierteilung
(1+1, 2+1, 2+2, 3+2, 3+3, 4+3, 4+4, 5+4). Die 3+2 enthält also keine Information über die Zahl
fünf hinaus. Die Orientierungsauswahl ist für jedes ungerade N generisch (det Σ ≠ 0).

**Stärkster Einwand.** In einem endlichen invarianten Maß bleibt ⟨Σ⟩ = 0; die Auswahl betrifft
den Sattel, nicht einen Zustand. Gebraucht wird ein Konzentrationsgrenzwert (viele Kopien,
Volumen). Dazu ist die spurfreie su(5) Quelle selbst gewählt: Eine Hubbard–Stratonovich
Zerlegung einer SU(5) Strom-Strom-Kopplung liefert sie, aber diese Kopplung folgt nicht aus P1.

## 3. F2: reelle Moden (`s3_majorana_stufen.py`)

Für Σ ∈ so(2N) auf dem Spinorfockraum von 2N Majoranamoden (N = 5 und 8, ρ = 0,2 … 0,45) sind
im globalen Minimum alle |μ_k| gleich: **Σ ∝ J mit J² = −1, SO(2N) → U(N)**. Die chirale Spur
wählt eine Orientierungsklasse (gerade Anzahl negativer μ_k). Am J-Sattel sind alle su(5)
Richtungen stabil (Krümmung 2ρ − sech²μ* > 0): Eine so(10) Quelle allein bricht nicht weiter
nach 3+2. Die komplexe Struktur J („Multiplikation mit i“ der Gitterkarte) und die Teilung 3+2
sind zwei Stufen mit zwei Quellen, nicht eine.

## 4. F3: Kette μ₄ → 16 → 5 → 3+2 (`s4_kette_anomalie_fock.py` mit A3 vom 20.11.)

- A3 (20.11.): Ein Weylfermion in d = 4 mit ℤ₄, ℤ₄² = (−1)^F, hat η = −1/16; anomaliefrei sind
  Vielfache von 16 (Literatur ℤ₁₆). A3 schloss: „erzeugt weder 5 noch 3+2“.
- Neu: Liest man die Generation als geraden Fockraum von N Moden (dim 2^(N−1)), so ist die ℤ₄
  Ladung X = N − 2p auf allen Zuständen gleich N mod 4, also ℤ₄² = (−1)^F genau für ungerades N.
  Die kleinste anomaliefreie Lösung ist **N = 5** (weitere: 7, 9, 11, …).
- Mit F1 folgt aus N = 5 die Teilung 3+2 mit festgelegter Orientierung, also die Hyperladung bis
  auf die Skala.

**Urteil: Kandidat als Reduktion, Stufe 2 für P2 unter benannter Lesart.** P2 (fünf Plätze,
3+2) wird zurückgeführt auf: vier Dimensionen mit μ₄ (Spin-ℤ₄) als Struktur, die Lesart
„Generation = gerader Fockraum von Randmoden“ mit Minimalität und eine positive Fockspur-Quelle
mit spurfreiem Feld. Keine angepasste Konstante.

**Basislinie.** N = 5 ist das erste Glied einer unendlichen Reihe; erst die Minimalität wählt
aus. Die balancierte Teilung ist generisch. Der Treffer ist strukturell (5, 16, 3+2, Orientierung
aus einer Kette), nicht numerisch.

**Erste freie Wahl.** Die Minimalität der Generation und d = 4. Damit ist 3+1 hier Eingabe,
nicht Ausgabe.

**Stärkster Einwand.** Jedes Glied für sich ist bekannt oder generisch (ℤ₁₆, 2^(N−1), balancierte
Teilung). Neu ist nur die Verbindung über die Fockspur, die in F1 die Darstellungsfrage
entscheidet. Die Familienzahl, D₅ ⊕ A₃ (also die Teilung 16 = 10 + 6 reeller Moden), c₃, φ₀, 41,
Raumzeit und Gravitation liefert die Kette nicht.

## 5. Lesart des Ladungsblocks (zu Bericht 2)

Die beiden Sektoren ±5 von End(ℂ⁵) sind genau die Richtungen (3̄, 2)_{+5/6} und (3, 2)_{−5/6},
also die X,Y Richtungen der SU(5). Der Ursprungskanal mischt X- und X̄-artige Operatoren mit
Stärke 41/216; die geforderte Referenzkohärenz ΔJ = 10 entspricht einem X X Paar. Das ist eine
Benennung, kein Ergebnis. Sie stützt die Lesart aus Bericht 2, dass der Kanal ein Auslese- oder
Compilerkanal ist und keine physische Zeitentwicklung geladener Felder.

## 6. Was jetzt noch Eingabe ist

Nach dieser Runde stehen für P2 und die Fermionstruktur noch:

1. **N = 5** über die Minimalität (oder gleichwertig ein Prinzip, das genau einen kleinsten
   anomaliefreien Focksektor verlangt);
2. **10 + 6** innerhalb der 16 reellen Moden, also D₅ ⊕ A₃ beziehungsweise das Element der
   Ordnung 4 (Kac-Marke 4); die Fockspur selbst wählt balanciert (8 = 4 + 4) und gibt das nicht;
3. **die spurfreie Quelle** und ein Konzentrationsgrenzwert für einen ausgewählten Zustand;
4. **d = 4 mit μ₄**, das in dieser Kette Eingabe ist.

## Kleinste entscheidende Tests für eine Folgesitzung

1. **Zwei Quellen aus einer.** Gibt es eine einzige positive Fockspur über 16 reelle Moden
   (E₈ Niveau eins, 8 komplexe Moden), deren globales Minimum 10 + 6 statt 8 + 8 wählt, etwa durch
   die Summe über Spinstrukturen getrennt für zwei Modengruppen? Abbruch, wenn jede solche Spur
   balanciert bleibt.
2. **Konzentration.** Für K Kopien (Volumen) die exakte Integration ∫dΣ e^{−K Γ(Σ)} rechnen und
   prüfen, ab welchem K die Bahn −6Y einen Zustand auswählt. Erwartet: Konzentration ∝ 1/√K.
3. **Normierung.** Fixiert die Kombination T, ρ aus P1 (antiperiodischer Kreis, c₃ = 1/(8π)) den
   Punkt im Fenster? Abbruch, wenn ρ frei bleibt.
