# TFPT: Was funktioniert, was fehlt, was eine vollständige Lösung wäre

Kompakter Stand vom 20. November 2026, ergänzt am 23. und 24. November (Befunde 20.–24. November, Kandidat H6 bis
H6h47, Leitlinie „Rekonstruktion statt Auswahl“ und Abschnitt 5), für Explorations- und Lösungssitzungen. Seit 24. November
liegen die Explorationszweige vom 20. November in main (Quaternionrahmen, Fibonacci-Uhr, Vertauschung und Auslese, fremde
und operative Ursprünge, Projekt Genesis, die Vertragsserie vom 14./15. Oktober); ihre Befunde stehen in ihren Verzeichnissen
unter `experiments/`. Nicht zusammengeführt ist `codex/c3-p1-p2-abgleich-20261105` (121 Konflikte mit dem heutigen Stand). Neuere Befunde stehen in den Atlasproblemen mit Datum
ab 20261120 und unter `experiments/exploration/`; dieser Stand wird nicht automatisch
fortgeschrieben. Er ersetzt
keine Originalquelle. Statuskürzel: **E** exakt unter genannten Voraussetzungen, **C** bedingte
physikalische Identifikation, **N** numerischer Vergleich, **O** offene Herleitung. Maßgeblich
bleiben `verification/status_ledger.csv`, die Papiere unter `documents/` und der Atlas.

## 1. Das Ziel

TFPT erklärt bereits sehr viel aus zwei gesetzten Eingaben. Gesucht ist die **vollständige
Lösung**: ein einziges Prinzip oder ein einziger Prozess, aus dem

1. die beiden Eingaben P1 und P2 selbst folgen, ohne zusätzliche freie Wahl,
2. die vorhandenen Ergebnisse (Abschnitt 3) als Folgerungen wieder herauskommen,
3. die offenen Teile (Abschnitt 5) mitgeliefert werden: physischer Zustand und Zeit,
   3+1 dimensionale Raumzeit, Kontinuum und Dynamik, Gravitation.

Teillösungen zählen, wenn sie mindestens eine Eingabe ohne neue freie Wahl erzeugen und
die Ergebnisse unverändert tragen. Atlasziel: `map:toe-goal`, `map:ur-to-tfpt`.

**Leitlinie seit 24. November: Rekonstruktion statt Auswahl.** Bedingungen an Kinematik, Symmetrie und Topologie legen
eine Theorie nicht fest. Der bisherige Weg fand zu ihnen einen passenden Operator, stellte dessen Mehrdeutigkeit fest und
ergänzte das nächste Auswahlprinzip. Das wird nicht fortgesetzt. Stattdessen gilt:
- **Dynamik:** Gravitative und materielle Dynamik werden aus den vollständigen Antworten des originalen Prozesses
  rekonstruiert (GNS, Osterwalder–Schrader, Tomita–Takesaki), nicht aus frei entworfenen positiven Propagatoren.
- **Chiralität:** Chirale Materie darf aus Indizes kommen. Sie muss nicht aus gegappten Spiegeln eines voll diskretisierten
  Fermionoperators entstehen.

Erster Schritt (Z58):
- Aus der Calderón/Szegő-Kovarianz allein (quasifrei wie im Original) folgen Modulgenerator, CCR, gaußsche Ströme und die Leiter statt des Bands.
- Der Spurkanal-Propagator ist die Fisher-Information des Originals.
- Die Chiralität folgt aus Indizes, ohne Spiegellücke.
- In 2D ist auch der Tensorkanal rekonstruiert (Z59, c = 1 je Quarter).

Offen sind der Tensorkanal und die Rekonstruktion von 3+1 aus 4D-Antworten. Z60 schärft das: Die Zahl der Richtungen
ist die Spektraldimension (Weyl-Gesetz) einer Antwort. Die Quelle hat d = 1, die E₈-Naht kompakt und starr d = 8 = 4 + 4.
Gebraucht wird eine originale Antwort mit d = 4, nichtkompakt und reflexionspositiv.

**Gesamtstand und Strategien (Stand H6h47):** [`GESAMTSTAND-UND-STRATEGIEN.md`](GESAMTSTAND-UND-STRATEGIEN.md), mit
Tor-Karte T1–T8, Kriterien K12–K20 und sieben Strategien (S1–S7).

**Gesamtbild (Vogelperspektive, Stand H6h35):** [`VOGELPERSPEKTIVE.md`](VOGELPERSPEKTIVE.md).
- H6h35, die Zusammenführung:
  - Die 3+1-Kinematik hängt an einer Aussage I_pos, gleichwertig zu T′.
  - Der Feldkern ν = 1 ist durch zwei Auslesungen erzwungen (Z74).
  - Die kontinuierliche Zeit ist die parabolische Wiesbrock-Translation im Vakuumsektor; im thermischen Original gibt es sie
    nicht (Z72).
  - Der 4D-Inhalt ist gezählt, die Gruppe ist Lesart (Z73).
  - Offen: I_pos, die Oszillatornormierung (Toeplitz → CCR), Vakuum gegen thermisch, die nichtlineare Gravitation, die
    Mehrpunktantworten und der Restflavor.
- Gelöste Bausteine (Abschnitt 3: Teilchen, Flavour, Massen, α, E₈) werden übernommen, nicht neu aufgemacht.
- Holomorph ist die Naht kanonisch eine Twistorgerade: π_*(η ⊗ π*O(2)) = O(1) ⊕ O(1), mit η = (A,A) − (P,P) als einziger
  vierteldrehungsinvarianter 2-Torsion (Z64, Z70).
- Lorentzisch wird sie nur über das konjugierte Partnerdublett. Das ist die Prämisse T′: Naht = P(V) ⊂ P(V ⊕ V̄*).
- Der Dublettzustand ist ρ_m, eindeutig je Windungszahl m; m selbst ist offen (Z63).
- Linear gibt es beide Gravitonhelizitäten; nichtlinear ist die Lage googly-dual (Z65).
- Der globale Modulfluss ist keine physikalische Zeittranslation (Z70). Die lokale Modulzeit auf Nahtbögen ist kontinuierlich
  (Z71). Das eigene Spektrum legt keine Skala fest (Z68).
- Vorbestand und Firewalls stehen in `experiments/exploration/2026-11-23-phasenraum-f5-holocode/PRUEFUNG-BESTAND-H6h34.md`.
- Offen: T′, m, die physikalische Zeit, die nichtlineare volle Gravitation, die Flächennormierung und der Restflavor.

## 2. Die gesetzten Eingaben

- **P1, Randkern (Naht):** orientierter, reflexionspositiver Randkern mit Einheitswindung.
  Seine operative Saat enthält Zustand und Zeit. Normierung c₃ = 1/(8π). [gesetzt; c₃ unter
  Normierungsverträgen hergeleitet, C]
- **P2, Träger:** fünf Plätze mit markierter Aufteilung 3 + 2 als physische Schnittstelle.
  [gesetzt; algebraische Reduktionen vorhanden, C]

## 3. Was funktioniert

**Teilchen.** Die gerade äußere Potenz von C⁵ hat 1 + 10 + 5 = 16 Komponenten. Mit den
Platzhyperladungen (−1/3, −1/3, −1/3, 1/2, 1/2) ist das genau eine Standardmodellgeneration
einschließlich rechtshändigem Neutrino. Erhaltungsgruppe (SU(3) × SU(2) × U(1)_Y)/Z₆. [E]

**Vier Marken.** Derselbe Divisor vom Grad 4 auf der Sphäre liefert H¹ vom Rang 3
(drei Familienrichtungen) und h⁰(O(D)) = 5 (die fünf Trägerplätze). Die Nahtnormierung
benutzt dieselbe Viererstruktur. Vier Marken → Doppelüberlagerung → Torus, μ₄ Deckstruktur,
Genus 3 Überlagerung mit gehobener Clock. [E für die Beziehung, C für die physische Lesart]

**E₈.** Trägergitter D₅ und Familiengitter A₃ kleben über Z₄ Glue zum geraden unimodularen
E₈: q(D₅) + q(A₃) = 2, Wurzeln 52 + 60 + 64 + 64 = 240, Dimension 248. Unter positiver
Standardnorm, einseitiger Integralität und primitiver Projektion erzwingt eine primitive
integrale Ladung n = 4 und ganz E₈ (R102). Hammingcode [8,4,4] → Konstruktion A → E₈;
Clifford-, Code- und Zwei-Qubit-Beschreibungen derselben Struktur. [E]

**Clocks.** Coxeterzahl 30, μ₄ als Automorphismus, Netz C₅/C₆/C₃₀, explizite
Invariantenabbildung E₈ ↔ G₃₁. Clocks greifen über konkrete Wirkungen ineinander. [E]

**Seed und Massen.** φ₀ = 1/(6π) + 3/(256π⁴) = 1/(6π) + 48 c₃⁴, λ_Y = √(φ₀(1−φ₀)),
Massenformel m̂ = (v/√2) λ_Y^L Λ = (vπ/√2) c φ₀^K. Transportlängen L und Potenzen K für alle
neun geladenen Fermionen aus Residuen-, Gradierungs- und Blattoperatoren; Leptonkoeffizienten
(16/7, 4/3, 7/6); Quellverhältnis m̂_u/m̂_d = 55/117. [E innerhalb der Annahmen; Vergleich mit
Pol- oder laufenden Massen braucht die Transfervorschrift, C]

**Alpha.** Teilcheninhalt und Transportbudgets geben den Materiekoeffizienten 41. Die
Wurzel der Randgleichung ist α⁻¹ = 137,03599921684, global eindeutig bewiesen. Gemessen
137,035999177(21). [E für Wurzel und Eindeutigkeit, C für die Lesart als Elektromagnetismus]

**Mischung und Kosmologie.** sin²θ₁₂ = 1/3 − φ₀/2 = 0,30675 (NuFIT 0,307), Koeffizient
ε = q(A₃) φ₀ = (3/4) φ₀ aus derselben Gluenorm; CKM mit s₁₂ = λ_C, δ = π/3 + 3λ_C²;
Ω_b = (1 − 1/(4π)) φ₀ = 0,04894 (Planck 0,04897). Koide: Quelle 0,664, Vermutung
Quelle → Pol mit Faktor 53/54 gibt 0,66666610. [N, C; Koide Vermutung]

**Quanten- und Prozessstruktur.** Vollständige Antworten bestimmen den minimalen
Prozessquotienten; ein positiver Wortkern rekonstruiert eine zyklische Darstellung (GNS) [E].
Im ursprünglichen L/J Baustein reichen genau 68 innere Richtungen: 5 sichtbar, 63 Gedächtnis
(R022/R023) [E]. Endliches E₈/Code/Clifford Gerüst; exakte stochastische Interpolation des
Transferspektrums auf endlicher Stufe [E].

**Jüngste Befunde (18.–20. November).** τ = i ist der symmetrieerzwungene kritische Punkt
jeder Spektralfunktion der antiperiodischen Struktur (KERNSYM). Bei τ = i 32 chirale Nullmoden
in der Verklebung (NEUSING). Der Schritt 1 + i ist der Frobenius der Kurve y² + y = x³ + x + 1
über F₂; E₈ als Gitter der 120 Einheitsquaternionen von SL₂(F₅). Eine blinde Regelsuche über
Vorzeichen liefert 5 und 16 mit eindeutigem Partner A₃ (GENESIS). Alles Hypothese oder
bedingt, nichts befördert.

### Befunde 20.–23. November, konsolidiert

Quellen: Atlasprobleme mit Datum 20261120–20261123, `experiments/exploration/2026-11-20-*`
und `2026-11-23-*`. Nichts davon ist im Ledger befördert.

- **P2 ist weitgehend auf „4“ und „minimal“ zurückgeführt.**
  - Kette μ₄ → 16 → 5 → 3+2 (Exploration 23.11., F3): Ein Weylfermion in d = 4 mit ℤ₄,
    ℤ₄² = (−1)^F, ist anomaliefrei in Vielfachen von 16. Der gerade Fockraum von N Moden trägt
    eine konsistente ℤ₄ Ladung genau für ungerades N; das kleinste ist N = 5. [E für die
    Glieder, Kandidat Stufe 2 für P2 unter der Lesart „Generation = gerader Fockraum“]
  - 3+2 als Minimum: Die positive Fockspur über Λ^gerade ℂ⁵ wählt 3+2 und die Orientierung
    Σ ∝ −6Y; Γ = ρp₂ − tr log(1+Σ²) wählt global 3+2 für 1/3 < ρ < 1. Basislinie: dieselbe Spur
    wählt für jedes N die balancierteste Teilung, also steckt in 3+2 nur die Information „5“.
    Felddeterminante über die 16 Gewichte wählt dagegen 4+1. [E für die Minima]
  - Klebegleichung m = n+1, n² = 2^(m−1), m/4 + (n−1)/n = 2 hat nur (n, m) = (4, 5); die
    Korpuskarte ist die Standardform von M₂(ℂ), σ dreht 3 Produktbasen und fixiert 2
    verschränkte; alle 18 σ-invarianten Lagen von D₅ + A₃ in E₈ drehen genau 3 der 5 D₅
    Achsen. [E; Lesart als P2 bedingt, Ganzzahligkeit offen]
  - Vier Marken aus quadratischer Vollständigkeit: Genau eine quadratische Relation zwischen
    den logarithmischen Formen erzwingt n = 4 Marken; ein quadratisches Differential ohne
    Nullstellen mit nur einfachen Polen hat g = 0, n = 4 und Windung +1 am Pol. Das
    Produkt liefert fünf Trägermoden mit Viertelcharakteren; R² teilt 3+2, Y = (I−5R²)/12.
    [E für die Geometrie, bedingt als Ursprung]
- **c₃ und Seed haben eine Strukturlesart.** Viertelinvarianz und fehlende Nullmoden wählen
  auf dem quadratischen Torus die Spinstruktur (−,−); ein Viertel der Dirac Wärmespur ergibt
  s·H → 1/(8π) = c₃ (die Viertelzählung ist gewählt). Für X ∼ N(0, c₃) ist κ₄(X²) = 48c₃⁴,
  genau der quartische Seedterm; die freie fermionische Kreisquelle trägt ihn, und die
  Einheitsladung W₁ liefert relativ Δκₙ = n!·2ⁿ⁻³c₃ⁿ. Der Basisterm ist (1 − tr Y²)/π =
  1/(6π). [bedingt]
- **Lokale Eichstruktur.** Das Original W erfüllt alle Wardidentitäten (SU(3), SU(2), Y,
  U(5), Spin(10), SU(4)) und ist lokal parallel; die Krümmung von Gr₃(ℂ⁵) spannt
  s(u(3)⊕u(2)); Berryfelder geben keinen quadratischen Eichterm; vollständige Polarisation
  fixiert den lokalen Dirac im festen Rahmen. Das volle E₈ Modell ist im Fenster 10 ≤ q ≤ 15
  lokal stabil zertifiziert, global offen. [E im Geltungsbereich]
- **Dynamik ist nicht autonom endlich.** Der ursprüngliche Kanal relaxiert exakt mit (3/4)ⁿ
  (Spektrum 1, 3/4, 11/36, 1/4); sein Hauptlogarithmus ist kein GKSL Generator; jede exakte
  autonome Erweiterung braucht eine unendlichdimensionale Umgebung; additive Y Erhaltung
  scheitert. [E im Geltungsbereich]
- **Skala und Abstand.** Umläufe und Ortsabstände bestimmen denselben Dämpfungsfaktor; nur
  y = 1/3 trägt Einheitswindung; Nullenergiedaten bestimmen die Kragenlänge nicht; gleiche
  frühe Energieantwort lässt verschiedene Innenräume zu. [E/bedingt]
- **Methodik.** R113 ist eine eingeschränkte Route, kein Türsteher (Entscheidung 23.11.).

**Muster dieser Runde.** Fast jede verbleibende Eingabe läuft auf zwei Fragen hinaus: woher
die **Vier** kommt (μ₄, vier Marken, d = 4, Viertelzählung in c₃) und welches Prinzip
**Minimalität** verlangt (kleinstes anomaliefreies N, kleinste Klebelösung). Ein neuer
Ursprung sollte beides zugleich liefern.

### Kandidat H6: Ursprung von P1 aus einer hintergrundfreien Naht (Exploration 23.11.)

Quelle: `experiments/exploration/2026-11-23-phasenraum-f5-holocode/VERDICT.md` (Stand H6h47),
Skripte `n1`, `v1`, `y1`, `z1`–`z99`, Lean-Vertrag `experiments/theory-contracts/h6-exakte-teile-20261123/`; Atlasproblem `problem:c3-p1-p2-shared-origin-reconciliation`,
Befunde mit Präfix `h6`. Nichts davon ist im Ledger befördert. Prüfstand mit 20 Kriterien (Stand H6h47): 13
hergeleitet, 7 unter benannter Lesart, keines gescheitert, offen oder gesetzt. K12 (3+1) zählt seit der Leitlinie
„Rekonstruktion statt Auswahl“ nur als Lesart: d = 4 folgt kinematisch, nicht aus den Antworten (Z60).

- **Prämissen:**
  - **A\*:** Die Naht hat keine Hintergrundgeometrie; Basiswechsel sind Symmetrien, die Form
    erfüllt ⟨T⟩ = 0.
  - **C″:** Die Familien sind der chirale Index der vier Marken.
  - **D:** Gravitation als Zustandsgleichung (Jacobson).
  - **E:** Die Teilung 3+2 ist lokal.
  - **G:** Die Naht ist ein verdrillter Defekt eines Paarungsmassenfelds.
- **Was folgt:**
  - Basiswechsel als Symmetrie geben die GSO-Summe. Nur 16 Majorana-Weyl haben Twistgewicht 1,
    also E₈ (c = 8) mit flachem (P,P) Sektor. [E]
  - ⟨T⟩ = 0 heißt E₆(τ) = 0, also τ = i als einzige stationäre Form. Stabil ist sie, weil der
    Sattel Verdrillung (Impuls) und Temperatur (Energie, positive Wärmekapazität) trennt.
    Zeitneutralität und KMS-Zustand bei β = L folgen, ebenso Kragenbreite L/2 und μ₄ als
    Stabilisator. [E]
  - Die gelifte Vierteldrehung gibt c₃ = 1/(8π) je Spinorsektor exakt, ohne Fixpunktkorrektur. [E]
  - Die vier Marken bei μ₄ geben den Fünferträger H⁰(O(D)) und den chiralen Familienindex
    H⁰(K(D)) = 3 (ker D₋ = 0) mit Charakteren {−1, i, −i}. [E] Der physische 3+2-Schnitt ist der
    V₄-Schnitt der Marken (GESAMT2, Z34): schwach = V₄-invariant {1 − z⁴, z²}, Farbe = Familien als
    V₄-Darstellung und = Λ²F als D₄-Darstellung. Halbdrehung (Z10) und T-Kern (Z16) geben dieselben
    Dimensionen, aber nicht den äquivarianten Teilraum (4 von 5 Hyperladungen verschieden). [E]
  - Λ²₊ = Λ²(Familien), also d = 4; Lipschitzgitter ℤ⁴ mit gleichwertigen Achsen. [E/C]
  - Lorentz bei Dimension ≤ 4 ohne Feinabstimmung; Reflexionspositivität über T². [E, perturbativ]
  - Yang–Mills aus Fermionschleifen über die lokale 3+2-Verbindung. [E im Modell] Die Lesart, die
    Kopplungen würden an der Gitterskala induziert und seien dort gleich (sin²θ_W = 3/8), fällt
    gegen die Daten (Z25, zwei Schleifen): Das SM spreizt an 1/a um Δα⁻¹ = 15–19, E₈-Zwischeninhalt
    mit vollständigen Multipletts vereinigt höchstens bei 1,1·10¹⁶ GeV, und die induzierte Kopplung
    ist an 1/a stark (α⁻¹ ≤ 23 gegen 60–74 aus der RG). Überleben kann sie nur mit SO(10) als
    nackter Eichgruppe (45 Eichbosonen, nicht aus dem 3+2-Split); dann liegt p → e⁺π⁰ unter
    3,5·10³⁴ a, also bei Hyper-K. [N/X]
  - Unter D liegt der Gitterabstand oberhalb der Planckskala: robust a/ℓ_P ≥ 1,65 (bosonischer
    Sockel ohne Randterm); die genaue Zahl hängt am Regulator. [N]
  - Die Naht ist ein verdrillter Defekt mit einer chiralen Mode je Fermion mit Spin-ℤ₄-Ladung 1.
    16 schwere reelle Fermionen (10,1)+(1,6) geben c = 8; der Schnitt ist genau bei ν ∈ 16ℤ
    topologisch unsichtbar. Die Wand am Schnitt ist unter G kein Gitterartefakt: Die Phasenbarriere
    ist der relevante Operator Λ² Re(Φ²) (ΔV = 0,254 m², auch für die lokale GJ-Form). Satz: Eine
    exakte ortslokale U(1) mit geladenem Φ lässt einem gegappten Wirbel keinen Index
    (Nielsen–Ninomiya). Geschützte Naht, lokale positive Transfermatrix und flache Phase ohne
    Abstimmung schließen sich daher aus. Der lokale Weg sind zwei abgestimmte Gegenterme
    (Re Φ², Re Φ⁴) wie die kritische Masse bei Wilson-QCD; U(1)_J ist flach, schützt den Index
    aber nicht (Z22, Z24). Der nichtlokale Weg ist eingeschränkt positiv (Z26): Lüschers projizierte
    Weylfelder sind auf der Ankeralgebra reflexionspositiv, bei gedrehter Phase unter der Reflexion
    auf χ = (1 − D/2)ψ, auf der vollen lokalen Feldalgebra aber nicht (wie Lu–Shao,
    arXiv:2611.20084); die Phase
    ist flach, der Index gleich der Windung, und H6s J ist in diesen Variablen exakt. Preis: Die
    Wirkung in χ ist potenzartig nichtlokal und der Transfer hat ein Kontinuum, unter der strengen
    Gitterlesart (lokale hermitesche H) also unzulässig. [N/E, frei]
  - Gerechnet ist der Zufluss der SO(2)-Ladung und die Gravitationsanomalie: c = 0,4994 und
    h_R = 1/16 je Paar, also c = 8 und h_R = 1 für 16 Paare. Die Volumenwand trägt die Chernzahl 1
    je Windung, also c_Wand = c_Naht. Der Fluss in Krümmung ist nicht gerechnet. [N/E]
  - TFPTs Einheitswindung [u_Σ] = 1 ist die relative Chernzahl der Massenpfaffianlinie auf der
    Normalscheibe; im Sektor 1 erzwingt sie einen Rangverlust, wählt aber den Sektor nicht. Unter G
    kommen verdrillte Nähte durch geschlossene Flächen paarweise. Ein Ordnungsfeld mit Ladung 1
    unter einer Hopf-U(1) auf S³ erzwingt genau eine Naht; offen ist die Ladungswahl. [E/N]
  - Geometrischer Ursprung der Marken: Eine verzweigte Doppelüberlagerung w² = s(z) der Sphäre mit
    Torus hat genau 4 Verzweigungspunkte; mit Vierteldrehung ist s = a z⁴ + e, also j = 1728,
    τ = i. Die μ₄-äquivariante Abbildung T: H⁰(O(D)) → H⁰(K(D)) macht die Familiencharaktere zum
    nichttrivialen Teil des Fünferträgers [E]; ihr Kern {1, z⁴} ist aber nicht das V₄-Dublett (Z34).
  - Naht-16 und Generations-16 lassen sich nicht vereinigen (Z23): Die Naht-16 liegt in keinem
    E₈-Multiplett, und ein reelles Multiplett mit Paarstruktur J gibt Familien = 3 × konj(Naht),
    also vektorartig. Konsistent ist nur die direkte Summe: Familien 3 × 16 (leicht), Eltern
    (10,1)+(1,6) (schwer), anomaliefrei mit gemeinsamem Spin-ℤ₄ J (ν ≡ 0 mod 16). Ein Singulett-Φ
    kann die chiralen Familien nicht paaren (16 ⊗ 16 ohne Invariante); damit entfällt die
    Familienmassen-Obstruktion des gemeinsamen Operators aus Z18. [E]
  - Spinorordnung (Z28): Der reine Spinor s im Fockraum Λ^gerade ℂ⁵ hat den Stabilisator SU(5),
    mit 3+2-Teilung G_SM (33 von 45 Richtungen gebrochen). P2 ist damit als Ordnungsparameter
    lesbar. Die kleinste zentrumerhaltende Potenz s^⊗4 hat den Stabilisator G_SM × ℤ₄, also
    π₁ = ℤ₄: Die interne Viererholonomie folgt aus dem Symmetriebruch, falls s^⊗4 kondensiert. Ein
    J-ungerades ⟨Φ⟩ lässt nur ℤ₂ (Wand); ein Ladung-4-Kondensat ist ohne Wand verträglich. Die
    Rangschranke der zusammengesetzten Verbindung (47 < 72) fällt mit einer unabhängigen
    Elternverbindung. [E algebraisch, Dynamik offen]
  - Gravitationsstruktur (Z29): Die Familiengruppe S₄ lässt auf Λ²F nur δ_ij zu, also die
    isotrope Plebański-Simplizität; die Urbantke-Metrik ist lorentzsch nur mit Re(a b̄) = 0. Der
    positive Kegel von Herm(2) legt die Minkowskiform fest. Die Wardidentität erzwingt die
    Spin-2-Struktur, aber nicht G (ξRφ² ist kovariant); auf dem festen Gitter verletzt schon eine
    reine Koordinatenänderung die Wardidentität (C_long = 0,2513 bei m = 0). [E/N]
  - Toeplitz-Ursprung (Z30): Eine irreduzible nichtunitäre Isometrie ist der einseitige Shift mit einem
    Defekt vom Rang 1 (Index −1). Ihre zweiseitige Vervollständigung trägt eine Spiegelung, die die
    Entwicklungsrichtung umkehrt, und J = RC mit J² = −1. Die vier Marken ergänzen sich als
    Qubit-Zustände über ihren Kommutator zu Herm(2); das Oktaeder trägt S₄ und die Tetraeder-
    Nullvektoren die Signatur 1 + 3. [E unter der Prämisse Irreduzibilität] Einschränkung (Z32): Ein
    festes Markenquadrat hat nur D₄; geometrisch isotrop wird es erst über die drei Rahmen.
  - Anomalienbilanz (Z31): c = 8 auf der Naht verlangt die axiale, gravitativ anomale Ladung; eine
    direkte J-Fortsetzung gibt c ∈ 16ℤ und ein gerades Niveau. Ein Quartett-Kondensat ohne anomale
    Antwort trägt keine Nahtmoden. Der Bruch Re Φ² hat axiale Ladung 4 und trägt auch im
    Ladung-4-Weg die Wand. Entscheidung 1 (eine Naht) und 2 (ein Block) sind beide ein trivialer
    Kommutant, also eine einzige Prämisse (Irreduzibilität). Ein Halbgruppenregler macht
    S_ent = A/(4G_ind) für minimal gekoppelte Felder exakt, aber nur als Eigenzeitregularisierung;
    als Propagatorschnitt verletzt er die Reflexionspositivität, und ξ bleibt frei (Z33). [E/N]
  - Lokale E₈-Naht ohne Gitter (Z32): Der Hardy-Kozyklus Tr(P U^m (1−P) U^{−m} P) = m ist der
    Schwinger-Term; er gibt Stromniveau 1 und c = 1 je komplexes Fermion, also c = 8 für 16 reelle.
    D₈ + s ist gerade und unimodular, und es gibt das lokale konforme E₈-Netz. Die Kette
    Index → Niveau → Anomalie → E₈-Netz ist damit exakt, unter den H6-Eingaben 16 und GSO. Der
    globale Shift und ein diskreter Gibbs-Zustand tragen keine Raumzeittranslationen (keine echte
    halbseitige modulare Inklusion). [E]
  - Originalbestand (Z34): Die (ℤ₂)³-Überlagerung trägt r_* r^* O(1) = O(1) ⊕ 6O ⊕ O(−1), mit den
    Deckhälften 5 + 3 (schwach, Farbe | Familien), E₋ ≅ E₊^∨ und O(−1) als Gegenpart ohne Nullmode;
    relative Selbstdualität zwingt L = r^*O(1). Im Randkern ist ε_car = ι_C|car = +1, wenn der
    Träger in Ran(P_prim) liegt (Typfehler; Träger braucht eine eigene Auslesung, etwa 2P_{V₄} − 1).
    Die Stresssättigung (12.11.) macht „ein E₈-Block“ zur Schranke c_W < 17/2. (c₁(L₂), c₁(L₃)) =
    (1, 0) gibt w₂ = 1, also keinen Spin(10)-Lift; E₈ enthält (Spin(10) × SU(4))/ℤ₄. [E/N]
  - Randoperator und Nahtinhalt (Z35): Die acht Schnitte von L = r^*O(1) (Deckhälften 5 + 3) sind
    reell genau die 16 Nahteltern (10,1) + (1,6); ihre Ströme so(10)₁ + so(6)₁ = D₅ + A₃ mit
    c = 5 + 3 und Kleber 5/8 + 3/8 = 1 erfüllen die Voraussetzungen der Stresssättigung. Der originale
    Randkern (ungerade, invertierbar, flach) verlangt den Spin-Dirac auf K^{1/2} ⊗ E mit einem
    Rang-1-Kanal zwischen O(1) und O(−1). Das einfache λ₁⁺ = 0,0852 des nativen Operators ist nur
    Konsistenz: Die Fixture hat diagnostische Zufalls-Yukawas (rng252), die zweite Vielfachheit weicht
    ab, und der Kanal bricht die Decksymmetrie (Z37). Calderón- und Hardy-Projektor haben denselben
    Toeplitz-Index −1. [E/N]
  - V₄-Umstellung (Z36): Ladungen, sin²θ_W = 3/8 und G_SM bleiben; die V₄-äquivariante Abbildung
    Träger → Familien hat genau das V₄-Dublett als Kern; V₄ wird eine diskrete Farbdrehung in SU(3),
    die Vierteldrehung wirkt auf dem Dublett mit det −1. [E]
  - Auswertungsgeometrie (Z37): Der Rang-1-Kanal ist die kanonische Erweiterungsklasse der
    Euler-Folge 0 → O(−1) → ℂ² ⊗ O → O(1) → 0; das volle Dublett ist trivial, die Windung sitzt in der
    Polarisation (Spin(10)-Lift möglich). Yang–Mills plus Fermiondeterminante (log det = 4 log t − 2t²)
    wählen stets t = 1: Die gewöhnliche Wirkung wählt den Kanal nicht. [E/N]
  - Flussnaht (Z38): Vierteldrehung und Twistgewicht 1 erzwingen V₄, acht Schnitte und 16 Majorana;
    5 + 3 lebt im periodischen Spinsektor, J² = −1 im (A,A)-Thetamodell. Eine Flussnaht auf der Zigarre
    (Ricci-Soliton) trägt mit den Flüssen 2 + 6 × 1 + 0 genau acht chirale Moden (flach nur 1), mit
    Niveau 8 und c₋ = 8; auf dem Gitter hält ihr Index ohne exakte chirale Symmetrie (Wilson
    Chiralitätssumme N, Overlap N exakte Nullmoden). Damit entfällt für die Naht die Wand aus
    Entscheidung 4; offen sind Zufluss, Gitterzigarre und Rückwirkung. [E/N]
  - Zufluss und Gitterzigarre (Z39): Die acht Moden tragen Normaldrehimpulse b = (3/2, 1/2⁷), die
    Anomaliematrix ist [[8, 5], [5, 4]] (I₄ = 4x² + 5xy + 2y² − p₁/3). Ganzzahlige b-Ladung auf allen
    240 Wurzeln hat nur die ungerade Spinorverklebung (dort ist q primitiv); die Auswahl hängt am
    Spinlift. Ein E₈-Chern–Simons-Gitter mit ganzzahligen l_q, l_b reproduziert die Anomaliematrix
    (Zufluss ohne Paarungsfeld). Die Rückwirkung f = (1 − y)(R² + Cy) löst die gekoppelten
    Dilaton-Maxwell-Gleichungen exakt (Modellklasse, nicht aus TFPT). Eigene Gitterzigarre: Overlap-Index
    N in 36 Fällen mit exakten Nullmoden der Chiralität +1, Wilson-Massen linear gegen null, auch
    rückgewirkt. APS: Eine 4D-anomaliefreie U(1) kann bei gebrochenem Fluss Nahtindex 1 haben. q-neutral
    sind nur 20 der 40 D₅-Wurzeln: Der Flusshintergrund lässt von Spin(10) nur U(5). [E/N]
  - Gemeinsamer Ursprung (Z40): Komplexe chirale Spin(10)-Multipletts haben gerade Stromniveaus
    (2, 16, 56, 4, 4), also kein Spin(10)₁ aus komplexem Fluss (dafür Majorana-Eltern oder eine
    Nahtsymmetrie). Zwei Schnitte von O(1) tragen keine lokale SU(2), erst O² der Euler-Folge. Das
    Triplett von ∇_t ist nur bei t = 1 masselos; ein fermionischer Polarisationsterm (Γ = 2t²)
    stabilisiert t* = 1/√2 (a = 1), dort ist es massiv: Rand- und schwache Verbindung sind verschieden.
    5k₅ + 3k₃ = 8 hat genau (1, 1). Der Hardy-Kern ist der kleinste positive Wallach-Punkt ν = 1 von
    det(I − ZW†)^{−ν}, genau dort gilt die masselose Wellengleichung (Cayley); vier tetraedrische
    Translationen brauchen die gemeinsame Bedingung (Σp)² ≥ 3Σp². E₋ = O(−1) + 3O gibt Familienindex 3
    ohne Spiegel, eine Generation mit V₄-Hyperladung ist anomaliefrei. Spin 2: Positivität verlangt
    Q²C ≥ Z; der ganze Formfaktor e^{−p²/M²}/p² der GRAV-/QGAMB-Ledgerzeilen verletzt das. Die sieben
    Zeilen tragen seit 24.11. einen Positivitätsvermerk (v1037). [E/N]
  - Altbestand und Graviton der Spektralwirkung (Z41): Der Altbestand hat eine reelle 16-Majorana-Naht mit
    Zufluss (SEAM.S3, p+ip × 16, c = 8); sie trägt Spin(10)₁ ohne U(1) und ersetzt die komplexe Flussnaht.
    Sein Masterfunktional, Tr e^{−D²/Λ²} mit Naht-KMS-Abschneider, gibt für Spin 2 den Kinetikfaktor
    a(x) = 1 + 12x g(x) = 1 − 3x/10 + x²/28 − … (Barvinsky–Vilkovisky), nicht e^x: tachyonischer Geist bei
    p_E² ≈ 6,77 Λ² (Residuum −2,59) und komplexe Polpaare; jeder getestete vollständig monotone Abschneider
    zeigt dasselbe. Positiv ist der Wallach-Kern (Wallach-Menge = 4D-Unitaritätsschranke Δ ≥ 1, Helizität 2 bei
    Δ = 3). M(P) ist die Yukawa-Kopplung an die Higgsrichtung, Γ_pol ihr Fubini–Study-Term. [E/N]
  - Richtungsentscheid und Prozesssuche (Z42, Z43): Die Dirac-Spektralwirkung erfüllt für jeden Abschneider die
    Summenregel ∫H_f/x^{3/2} = 0 und scheidet als fundamentale Gravitationsdynamik aus (auch mit κx; schon
    b = 1/28 < a² verletzt die Positivität); positiv ist ein Schur-Gedächtnis C = 1/x + 1/(x + 10/3). Die zweite
    Blatthälfte muss orientierungsumkehrend sein (8 − 8 = 0); zwei Bosonenpaare geben den Matrixkern exakt, der aber
    eine weitere Regel braucht; die Hopf-Windung der Higgsrichtung braucht ein verdrilltes Bündel; der SM-Inhalt
    hebt die quadratische G-Divergenz auf (48 − 48). Die neue OriginLab-Lane `process_rules` sucht über 53 760 ganze
    Prozesse: Acht lasttragende Prinzipien wählen Kern ν = 1, Spin-2-Gedächtnis, reelle 16er-Naht mit umgekehrtem
    Blatt und Familienpaar Z₃ × Z₄; kein exaktes Restsymmetriemuster passt zu den Neutrinodaten. [E/N]
  - Übergangsfamilie und ausgebaute Prozesssuche (Z44, Z45): a_q(z) = z⁻¹(1 − q z⁴) gibt für jedes q das Bündel
    O(1) + 6O + O(−1), die Hardy-Kompression S†(I − qS⁴) mit Index +1 und den positiven Randkern C_q (Pol 1 − q²
    plus Marchenko-Pastur-Band). Ihr Spin-2-Gedächtnis ist exakt H/x = 1 − a x + b x² mit a = q²/(1 − q²)²,
    b/a² = (1 + q²)/q²; a = 3/10 verlangt q² = (8 − √55)/3, dann b ≈ 0,553 (Schur-Kern 9/50, Spektralwirkung 1/28).
    Kein Rang-eins-Auslesen ist masselos und verbunden oder universell gekoppelt; das Graviton braucht eine
    Ward-Struktur. Die Paarkanäle schließen zu sp(4, ℝ), erst mit einseitiger Blattaktion zu su(2, 2); 5 + 3 ist
    eine Spin^c-Zählung (F = E ⊗ O(1), c₁ = 8); die Z-Masse sitzt in der Hopf-Phase, also bleibt M(H). Die Lane
    `process_rules` hat Übergangs- und Ausleseuniversen: Struktur legt d = 1, k = 4 fest, q nie; neun lasttragende
    Prinzipien (mit „eine Übertragung“ und Orientierungsumkehr, ohne Spiegelblatt) wählen über 67 200 Prozesse
    Kern ν = 1, Transfergedächtnis, reelle 16er-Naht und Z₃ × Z₄. [E/N]
  - q-Regel und Ward-Struktur (Z46, Z47): a = 3/10 ist eine Skalenkonvention (a = (3/10) f(0)/∫f), skalenfrei ist nur
    b/a² (25/63, 25/42, 0; für alle drei Abschneider < 1), also legt die Dirac-Zahl q nicht fest (Korrektur zu Z45). Die originale
    Quarter-Gaußquelle (Calderón/Szegő, b = 4π) hat ω(J₋₁J₁) = 1/(e^{4π} − 1); als Leiterbesetzung des Nullvektors
    gelesen folgt q = e^{−2π}, der Hardy-Zustand gibt q = 0. Spin 2 ist dann ein Einsteinpol mit Gewicht 1 − e^{−4π}
    plus eine schmale massive Linie mit Gewicht 3,5·10⁻⁶. Ohne die Skalenkonvention wählen acht Prinzipien genau
    diesen Prozess (Korrektur H6h30: im zertifizierten Modus zusammen mit dem Spiegelblatt bei q = 1/√2 zwei Prozesse). Ward-Struktur: Der masselose Pol ist durch den Fredholmindex +1 geschützt (G nicht); Weinbergs
    Softbedingung verlangt eine einzige Randquelle; Ward plus Positivität lassen das Gedächtnis nur im Spin-2-Kanal
    (W F(□) W, linearisiert ½k⁴hP²h), der Spurkanal bleibt Einstein. Gegenprüfung eines externen Laufs: 32er-Naht nur
    über die Größenordnung ausgeschlossen; Z₃ × Z₄ exakt (2 − √3)/6 = TM2 bei 15°; die Daten passen zu TM1 bei 15°
    (sin²θ₁₃ = (2 − √3)/12), Herkunft offen. [E/N]
  - Zuordnung, Quelle, Bimetrik, Familiengruppen (Z48, Z49): Die Quarterzwiste α = 1/4, 3/4 sind die Vierteldrehungsgrade
    der Hardy-Marken k = 4(r + α); die originale Kovarianz ist 1/(1 + e^{πk}), und J₋₁ wirkt auf den Einteilchen-Marken als
    Übertragungsschritt S⁴ (Korrektur Z54: auf der Besetzungsleiter ist der Schritt die Polarphase J₋₁(N+1)^{−1/2}).
    Damit ist die Zuordnung hinter q = e^{−2π} hergeleitet (Prämisse bleibt die Gaußquantisierung; seit Z58 folgt die
    gaußsche Struktur der Ströme aus der quasifreien Originalkovarianz, die Quasifreiheit selbst bleibt Prämisse). Ward verlangt
    Minimalkopplung über die Noetherladung |0⟩; Quarter-Harmonische gäben eine EP-Verletzung e^{−2π}|Δh|. Nichtlinear
    passt eine Bimetrik mit (M_f/M_g)² = 1/(e^{4π} − 1); kubische und quartische Kopplungen bleiben frei. Δ(6n²) mit
    Majorana-Klein-Resten gibt bis n = 16 nur TM2 bei diskreten Winkeln, keinen im θ₁₃-Fenster; Δ(96) macht die
    Familienuhr Majorana-verträglich; TM1 bei 15° mit μ-τ-Reflexion passt (δ = ±π/2), Mechanismus offen. [E/N]
  - Metrischer Stern und quartischer Rest (Z50, Z51): Der Vierteltransfer ist exakt der Massenoperator eines endlichen
    metrischen Sterns (Zentrum plus vier Arme, M_(j,n) = (M_c/2) q^(n+1), Potential (−4 + e₁ + e₃ − 4e₄)/2; B U = V A_rect,
    ein Nullvektor); minimale Knotenquellen koppeln universell (η_v/M_v = 1/M_eff). Korrigiert: Die zentrierte
    Feldkovarianz ist verbunden und masselos; „keine Nullstelle“ ist kein Positivitätszertifikat (Schranke xC ≤ 1/κ
    schließt die Dirac-Profilklasse auch mit Einstein-Zusatz aus); mit freier Skala trifft jedes q die Zahl 3/10.
    Ungerichtete Kanten legen die kubische Ordnung fest; frei bleibt eine quartische Zahl t, ein Viererkontakt auf den
    Transferkanten, der in der Tiefe wie q^{−2n} wächst; schwach gekoppelt bis zu endlicher Tiefe. [E/N]
  - Gravitation bis zur Nahtskala, Spiegelbilanz, Gaußzustand (Z52, Z53): Jedes flache, ungerichtete
    Fierz-Pauli-Potential ist ε(D∧D∧Q(e_p, e_k)); t misst die Form von Q (V 0, Mittel 1, e_p e_k 4, Quadrate −2), nur
    e_p e_k ist kubisch in Y = I − X. Mit Nahtskala ~ Planck (GRAV.GAPCORR.01) koppeln die Arme stark; darunter bleibt
    Einstein mit G aus GRAV.NONLINEAR.01 plus Scalaron (3,06·10¹³ GeV), der ein gesunder Spinnull-Pol und damit
    Ward-verträglich ist. Der KMS-Zustand einer freien Majorana-Naht ist exakt gaußsch: Die Gaußprämisse ist die
    Freie-Naht-Prämisse. Eine Familie mit ν_R hat 16 Weyl mit Spin-Z₄ (X ≡ 1 mod 4), alle Anomalien null, 16 ≡ 0 mod 16
    (ohne ν_R 15): notwendige Bedingung für Spiegel ohne Brechung erfüllt; Spielzeug 0+1 D: 4 Majoranas entartet,
    8 gegappt. Δ(6n²) gibt auch mit Dirac-Neutrinos kein TM1. [E/N]
  - Prüfung des externen Zweigs zu H6h27 (Z54, Z55; Zweig `research/h6h27-source-dynamics-20261124`, inzwischen in main):
    Volle KMS legt den gaußschen Stromzustand fest (gleiches Mittel genügt nicht); J₋₁ ist auf der Besetzungsleiter
    keine Isometrie, der Schritt ist die Polarphase (Korrektur zu Z48). Der Thermofeldoperator ist im Paarsektor
    A⁺(N+1)A: diskrete Leiter Δn statt Band; der gewichtete Stern hat endlich viele leichte Moden. Die t-abhängige
    Vierpunktantwort ist eine erzeugende Funktion Σ t_n (rP)^{n+1} (flach rP/(1 − rP), gewichtet Li₃(r) bei p → 0); die
    Vierteldrehung setzt nur die Arme gleich, je Tiefe bleibt ein t_n. G skaliert ohne Normregel frei. TM1 bei 15° mit
    Z₃, Majorana-Klein und μ-τ-Reflexion ist in jeder endlichen Gruppe ausgeschlossen (Spuren −(2 ∓ √3)/4 nicht
    ganz); TM2 verletzt die Solarsummenregel (> 3,5σ gegen NuFIT 6.1). [E/N]
  - Zusammenführung, Spiegellücke, Flächenregel, Clock-Familien (Z56, Z57; der externe Zweig ist in main, Merge
    a7b1ef933 mit zertifiziertem Modus run(certified=True); dort wählt die skalenfreie Suche zwei Prozesse, originaler
    Zustand und Spiegelblatt bei q = 1/√2, Korrektur zu Z46): Der lokale Spiegeloperator H = −(J/4)Σ(B_A + B_A⁺)² auf 16 Spin(10)-Fermionen × 2 ist
    J(C₂(Spin(10)) − C₂(Spin(11))); für Spin(7) und Spin(9) voll diagonalisiert, für Spin(11) über eine exakte
    Polynomrealisierung: eindeutiger Singulett, E₀ = −136J, Lücke 9J, Fermionlücke 29J/2, ⟨X_A X_B⟩ = 68/5, nichtgaußsche
    Vierpunktantwort (verbunden −2822/25), E(h) bis h⁴; die starke Projektion lässt sieben naive Knoten. Die
    Plebański-Flächenregel (Kante quadratisch in B_p − B_k) gibt t = 1 (Epsilon-Kontraktion = −64 × HR-Potential). Die
    originalen Clocks F_C, F_J tragen die TM1-Spalte (2/3, 1/6, 1/6), der Winkel hängt an einem freien Gewicht; die
    eingefrorenen Werte brechen TM1 und liegen näher an NuFIT 6.1; G = v_geo⁻². [E/N]
  - Rekonstruktion statt Auswahl (Z58, Leitlinie aus Abschnitt 1): Aus der originalen quasifreien Kovarianz
    c(x) = 1/(1 + e^{4πx}), x = r + α, folgt ohne weitere Eingabe der modulare Einteilchengenerator 4πx. Die Stromantworten m/(e^{4πm} − 1),
    Kommutator m und KMS-Verhältnis e^{4πm} geben CCR und die lineare Dispersion 4πm. Sugawara, Wick und das
    Jacobi-Tripelprodukt zeigen: Der Zustand ist Thetaverteilung der Ladung mal gaußsche Oszillatoren, der Modulgenerator ist
    quadratisch in den Strömen plus eine Ladungsfunktion, sein Spektrum liegt in πℕ. Damit ist Band gegen Leiter
    entschieden: Das Band A⁺A ist nicht spurklasse und keine Einschränkung eines quadratischen Zweimodenoperators. Unter
    Modulfluss und modularer Konjugation bleibt genau die Leiter, 4π je Quant (OriginLab `reconstructed_dynamics`). Der
    Spurkanal-Propagator ist die Fisher-Information κ₂ = 0,798 der relativen Entropie, die Kopplungen sind κ₃, κ₄; der
    Thetasektor dominiert (Faktor 724) und ist stark nichtgaußsch. Die E₈-Naht hat bei τ = i die Modulantworten
    0, E₄/6, E₄/(4π), 3E₄/(8π²) − E₄²/27. Chiralität ohne Spiegellücke: ein chiraler Zweig, Spektralfluss 1,
    Toeplitz-Index +1 mit Restmischung q^Tiefe, Familienindex 3 ohne Spiegel, Z₁₆ erfüllt. In 2D ist auch der
    Tensorkanal rekonstruiert (Z59): Die Virasoro-Antworten ω([L_m, L₋ₘ]) der Kovarianz geben c = 1 je Quarter und das
    Twistgewicht 1/32, also die ganze induzierte Gravitation der Naht (Polyakov-Wirkung mit c). Die Quelle trägt aber nur
    eine lichtartige Translation; die vier Marken geben die Nullrichtungen von R^{1,3}, nicht die Dynamik. Offen:
    4D-Tensorkanal, 3+1-Rekonstruktion, Familiengewichte. [E/N]
  - Spektraldimension (Z60): Die Zahl der Richtungen, die Antworten rekonstruieren, ist ihre Weyl-Dimension. Die
    Quarterquelle hat d = 1; ihre Viererstruktur ist eine Graduierung derselben Kette (direkte Summe), kein gemeinsames
    Spektrum vier vertauschender Generatoren (das gäbe d = 4). Das E₈-Ladungsgitter hat d = 8 und enthält D₄ ⊕ D₄
    (Hurwitz-Quaternionen, Index 4), ist aber kompakt und starr (248 Ströme vom Gewicht (1,0), keine Moduli). Kein
    vorhandenes Original rekonstruiert ein nichtkompaktes 3+1; K12 zählt deshalb seit H6h32 als Lesart. [E/N]
  - Gemeinsamer Elternoperator (Paket `2026-11-24-h6-common-parent-synthesis`, andere Sitzung, 40/40, 187/187, 201/201):
    - Der Wilson-Halbraum A_q(k) = (S⁺ − w(k))(I − qS⁴) über drei transversale Impulse hat genau eine masselose
      Weyl-Ecke. Er ist der nächste Kandidat für die d = 4-Antwort, aber seine k-Abhängigkeit ist eingesetzt, nicht
      rekonstruiert.
    - Der gewichtete Stern hat einen expliziten Verflechter B U = V A_q.
    - Der naive Familienlinien-Index ist 8 mod 16 (keine volle 16).
    - Die γ₅-Spur gibt t = 1, ist aber indefinit.
    - Zur endlichen Spiegellücke gibt es keinen Satz bei endlichem J. [E/N]
  - Vogelperspektive (Z61, Z62, `docs/research/VOGELPERSPEKTIVE.md`; korrigiert in H6h34: Vorbestand CELEST.SEAM.01 und
    WOIT.OS.TWISTOR.01; das schwache Dublett ist nicht H⁰(O(1)); Familien ≅ Λ²₊ ⊗ χ_q):
    - h⁰(O(4)) = 5 (Träger = Weylspinor), h⁰(O(2)) = 3 (Familien = Λ²₊), h⁰(O(1)) = 2 (Spinor), h⁰(O(1) ⊗ ℂ²) = 4
      (Raumzeit).
    - Die vier Marken sind die Nullrichtungen eines Spin-2-Zustands im Träger. Die Weierstrass-Invarianten der
      Doppelüberlagerung sind g₂ = I und g₃ = J; τ = i heißt also J = 0, und die H6-Marken sind (|2⟩ − |−2⟩)/√2.
    - Die Standardform eines Dublettzustands hat Herm(2) = R^{1,3} mit dem Zukunftskegel als natürlichem Kegel. Nur n = 2
      gibt einen Lorentzkegel.
    - Unter T (Naht = Twistorgerade) sind die Orte die Lagen der Naht. Holomorphe Twistorfunktionen geben exakt □φ = 0, die
      Weyl-Gleichung und das selbstduale Graviton.
    - Lücken: Die verzweigten Blätter geben O(1) ⊕ O(−1) statt O(1) ⊕ O(1), der Dublettzustand ist nicht aus dem Original
      hergeleitet, und der nicht-selbstduale Teil fehlt. [E/L]
  - Parallele Agenten (Z63–Z70, H6h34):
    - Der Dublettzustand ist ρ_m = diag(1, e^{−πm})/(1 + e^{−πm}), eindeutig, treu, KMS und basispunktfrei. Neutral ist nur
      m = 4, getrennte Marken verlangen ungerades m.
    - Satz: π_*M ist balanciert, außer für Pullbacks; z61s O(1) ⊕ O(−1) ist der Sprungfall, der das schwache Dublett trägt.
    - Das kanonische η = (A,A) − (P,P) gibt O(1) ⊕ O(1); die Spinbrücke der externen Analyse ist bestätigt.
    - Die Signatur (1,3) kommt nur über V̄*, die Naht selbst gibt (2,2) oder (4,0). Damit ist T äquivalent zu T′.
    - Linear tragen beide Blätter beide Helizitäten; nichtlinear gilt die googly-Dualität. Penroses Deformationen tragen Ψ.
    - Der periodische Modulfluss ist keine Zeittranslation. Aus dem eigenen Spektrum folgt keine Skala; am KMS-Punkt gilt
      ⟨T⟩ = 0, S = ln 2/4 je Quarter.
    - Der H6-Spin-2-Zustand gibt keine Flavourgewichte; Flavour bleibt Bestand.
    - OriginLab `reconstruction_search` (14 Tests).
    - Zweite externe Analyse (Z71):
      - Die lokale Modulzeit auf Nahtbögen ist dicht und nicht periodisch.
      - Die 4D-Fortsetzung stimmt auf der Linie, ist aber nicht erzwungen (der Nahtwinkel gibt die 2D-Identität).
      - Zwei maximal inkompatible binäre Beobachtungen geben τ = i und das kanonische η. [E/C]
  - Zusammenführung (Z72–Z74, H6h35):
    - L₀ ist elliptisch (die Uhr), P₀ parabolisch mit Spektrum [0, ∞). Im Vakuum der Uhr bilden Nahtbögen halbseitige
      modulare Inklusionen, und die Wiesbrock-Differenz ist P₀; im Original (β = 4π) nicht (Leckage 9·10⁻⁴).
    - D_m = D₁·sin(mγ)/sin γ ist eine Teleskopsumme. Σm²qᵐ ist nativ, eine su(2) auf den Paarräumen nicht (Toeplitz statt
      CCR).
    - Die Standardform liefert J = †, Herm(2), V ⊕ V̄* und die Form (2,2); T′ ⇔ I_pos.
    - R ist Ladungsdrehung mal Parität, also schließen Spinbrücke und Neutralität sich aus; mit Uhrenquotient m = 1.
    - Der Kern det(I − ZW†)⁻¹ ist durch die Hardy- und Stromauslesung erzwungen. [E/C]
  - Zweipfadkanal (Z76, H6h36, vierte externe Analyse):
    - Die Erhaltung des Originalzustands erzwingt zwischen den Paarschalen die Schwinger-CCR-Gewichte p_k = (k+1)/(m+1).
    - Kanal mit Petz-Rückkanal gibt exakt den Sphären-Casimir ℓ(ℓ+1).
    - K/π ≡ Q_a − Q_b mod 4.
    - Damit ist die Oszillatornormierung bedingt gelöst, als zustandsbalancierte Auslesung. [E/C]
  - Thermischer Kanal (Z77, H6h37, fünfte externe Analyse; korrigiert in Z80, H6h40):
    - [J₋₁, B_{r,m}] = B_{r,m+1} − B_{r−1,m+1} gilt exakt.
    - Z77 hatte die Übergänge vertauscht. Richtig zugeordnet (die Formel der Analyse ist das Gewicht „Loch runter“) ist der
      thermische Zweipfadkanal bei jeder Temperatur exakt zustandserhaltend, eindeutig unter den beschränkten Lösungen. [E]
    - Kantenaufgelöst (Mischung verschobener Vakua der Einteilchendaten) gibt Petz ℓ(ℓ+1) bei jeder Temperatur; der
      quasifreie Zustand selbst ist keine solche Mischung. [E]
  - Araki-Woods und Orte (Z78, H6h38, sechste externe Analyse):
    - Die Oszillatoren d, e aus Strom und Tomita-Partner sind kanonisch, und ihr Vakuum ist die Purifikation des
      thermischen Originals.
    - c = κ(a + b⁺) gibt vertauschende P_μ mit P² = 0.
    - Die Lagen X sind exakt aus der gemischten Korrelation auslesbar: (M_X − M_X⁺)/i = X für jedes q.
    - Damit ist „Vakuum gegen thermisch“ strukturell beantwortet, und I_pos hat einen Kandidaten. [E/C]
  - Mehrpunkt und Spin 2 (Z75, H6h39):
    - Die nichtabelschen Nahtströme sind nicht gaußsch (f-Dreipunkt, Nullvektoren).
    - Eine skalare 4D-Lesart ist ausgeschlossen, möglich sind nur Vektorkomposita.
    - T(z) bildet auf Δ = 2 und (j,j) ab, nicht auf Spin 2.
    - Helizität ±2 gibt es nur als Vierquantenkompositum im Araki-Woods-Dubleton, ohne kovariante T^μν-Kopplung
      (Weinberg–Witten), mit einem Speziesproblem.
    - Gescoptes Negativ für „Graviton aus der Virasoro-Antwort“. [E]
  - Quarterindex als Dublett (Z79, H6h40):
    - Kinematisch ist der Quarterindex ein vollwertiges Spinordublett. Beide Quarter haben dieselbe Stromantwort, der
      Oszillatorzustand ist U(2)-invariant.
    - Mit dem Tomita-Partner entsteht der Twistorraum mit CCR-Form (2,2). Der Kommutant von τ = (Ad J)⁺ in SU(2,2) ist die
      Lorentzgruppe mit Dilatation; das Zentrum ist der Modulfluss.
    - Die Orte sind Twistortranslationen zum Röhrenpunkt κ²(I − iXᵀ), mit dem ν = 1-Kern. [E]
    - Das Nahtdublett ist er nicht: Die Vierteldrehung der Marken hat auf V₁ Ordnung 4, auf dem fermionischen Quarterpaar 2
      und auf den Oszillatoren 1. Auf Fermionebene ist der Quarterindex genau das η-Dublett der Spinbrücke. [E]
    - Auf H⁰(N) = V* ⊗ η wirkt die Markensymmetrie euklidisch, nie lorentzsch; T′ (η = V̄*) ist nicht
      markenverträglich. I_pos für die markierte Naht bleibt offen. [E/O]
  - Analysen 7 und 8 (Z80, H6h40): Überlapp der Relativzustände 1/det(I − (i/2)(Y − X)), daraus s² = det(Y − X). Die Wurzel
    des Kerns ist nicht positiv (Wallach), zusammengesetzter Spin 2 hat keinen Pol, und Cl₂ folgt aus der
    Vorzeichensymmetrie (bedingt). [E/C]
  - Vorprüfung (Z81, H6h41):
    - R = i^(Q_a − Q_b), R² = (−1)^F: Die Uhr ist ein Spin-ℤ₄-Erzeuger. [E]
    - Die Einteilchendaten sind die Saatfamilie {ρ_m, m ungerade} auf Fermipaaren. [E]
    - Die Marken liegen über einem Punkt der physischen Windung. [E]
  - Strategien parallel gerechnet (Z82, Z83, Z86, Z87, Z88, H6h42):
    - S2 (Z83): su(2,2) schließt mit den Tomita-Wurzeln des nativen paarerhaltenden Hebers L^nat = Σ(k+1)ψ*_{k+1}ψ_k.
      Die Moden sind Energiestufen der ν = 1-Leiter, und die Minkowskienergie hat das Spektrum [0, ∞). Das Speziesproblem
      verschwindet auf Einteilchenebene. [C, eine Identifikation]
    - S1 (Z82): R ist auf dem Träger die ℤ₄ der 16 (Zentrum von Spin(10), bis auf Eichung), Spin und Statistik gelten innen
      exakt. Das Original erreicht nur die „2“, und die Quarter-Oszillatoren geben Fermionen kein 2π-Vorzeichen.
      [E, teilweise]
    - S5 (Z86): keine Flächennormierung aus der Verschränkung Original|Partner (Typ I, kein Horizont); eine Uhr für alle
      Sektoren. [E negativ]
    - S6 (Z87): Saat ρ₁ + Spin-ℤ₄ + CAR; der Ursprung ist auf die Prämissen P_sym und P_clock eingegrenzt, die
      Gaußprämisse folgt aus P_clock in allen gerechneten Abschneidungen. [E/C]
    - Analyse 9 (Z88): Der volle 4D-Fockraum passt energietreu nicht in eine Naht (ab E = 36 mehr Zustände). Die lineare
      Leiter folgt aus der Hardy-Invarianz des Modulflusses. [E]
    - S7 (OriginLab `scenarios`): Von 960 Szenarien bleiben 8 mit Schrittart-Spinor, inneren Marken, parabolischer Zeit,
      Energiestufen und Saat; die Gravitation ist offen.
  - Gravitation als Deformation (Z85, H6h43): Die Homogenität der Twistorfunktionen ist der Modulgrad, und die Tomita-
    Konjugation ist die Twistorkonjugation. Aus den nativen Quellen (gaußsch, Vertexoperatoren, Ströme, Ladung) folgt keine
    universelle Grad-2-Klasse; die zweite Ordnung ist Heisenberg, nicht Einstein. Die Marken sind keine Weylkrümmung
    (Quadrupoltest, Bianchi IX), sondern Hitchins vier Punkte, ein Symmetriedatum. [E negativ, E/L]
  - Osterwalder-Schrader (Z84, H6h44): Native reflexionspositive Spiegelungen gibt es; Θ_n ist die Bahn von τ aus Z79.
    Für das thermische Original gibt die Rekonstruktion genau den Modulfluss zurück. Positive kontinuierliche Zeit entsteht
    nur im Vakuumsektor (wie Z72). T′ gilt nur punktweise; R hat auf U(2) keinen Fixpunkt. [E/C]
  - Bilanz der parallelen Runde (H6h44): Raumzeit = radiale Leiter mit Schrittart-Spinor, Zeit parabolisch im Vakuumsektor,
    Marken und Uhr innen, Ursprung Saat ρ₁ mit P_sym und P_clock; Gravitation offen (radialer Kanal oder räumlicher Schnitt).
  - Folgerunde (Z89–Z93, Z95, H6h45):
    - Der Infinity Twistor ist bis auf seine Phase nativ und gibt die Einstein-Klammer exakt. Der radiale Kanal gibt nur
      AdS, und G bleibt frei. [C/E]
    - Ein räumlicher Schnitt fixiert G nicht; das radiale Feld induziert kein G. [E negativ]
    - Die Identifikation hängt nur noch an P_koh. Der freie 4D-Fockraum folgt aus bosonischer Komposition (P_ident), aber
      keine native Fusion gibt eine 4D-Wechselwirkung. [E/C]
    - Die „3“ ist nur algebraisch. Es gibt keine native Verklebung von Spin und Uhr; die Statistik steckt im Q²/2.
      [E negativ]
    - Die Kugel D(V, η) hat keinen Kern aus dem Ladungssektor und gehört nach innen. [E negativ]
    - Die elfte Analyse hält in den Rechnungen, auch die Schwelle N = 6 für identische Saaten (Korrektur H6h46). [E]
  - Starrheit und Ursprung (Z94, H6h46):
    - Die Starrheit des Originals trägt die Nichtentartung, nicht die Ladungsauswahl. S(2) ist bewiesen, für drei und vier
      Paare erschöpfend in Fenstern. [E/O]
    - P_clock folgt nicht aus P_sym. Es ist gleichwertig zur Zeitneutralität des Gesetzes P_tn (Modulzeit-Torus
      selbstdual). Zwei unabhängige Ursprungsprämissen bleiben: P_sym und P_tn. [E]
  - Folgeschritte G1–G3 und zwölfte Analyse (Z96–Z99, H6h47):
    - Zwölfte Analyse (MacDowell–Mansouri, Plebański, Clifford, Λ–G) hält in den Rechnungen. Die Clifford-Projektionen sind
      Instanton und Anti-Instanton der Twistorfaserung; die Λ–G-Relation ist COSMO.LAM.01 in neuen Variablen. [E]
    - Gravitation: Die Kinematik ist nativ – MacDowell–Mansouri-Algebra (AdS, Radius halbe Vakuumgröße), τ als
      Punktspiegelung, Σ± aus Original und Partner, Tomita-Familie als flaches Vakuum (F = 0). Die Phase von I ist eine
      Eichrichtung, Positivität fixiert g nicht, der Plebański-Sektor ist nur bis Z₄ fest, nativ ist AdS statt dS. [E/C]
    - Spin-Statistik: Die Regel (−1)^{2h} = R² schließt unter su(2,2), Fusion und gekleideten Vertizes (Schrittart ×
      Ladungsschritt), ist aber nicht nativ (Prämisse P_GSO) und vektorartig. Chiralität nur über den E₈-Kleber mit
      P_A3, im Konflikt mit der Skalarleiter und dem Familiengitter A3. [E/C]
    - Wechselwirkung: Vertexoperatoren sind auf Niveau 1 die su(2)₁-Ströme; Δ = 1 ist starr; eine verbundene 1→2-Amplitude
      braucht einen gewählten Zweikörperterm (Doppeltransfer, h = 4). [E negativ/O]
    - Szenario-Suche mit acht Achsen (13 824): 12 Szenarien erfüllen alles außer Graviton und Chiralität. [E]
  - Instantonen (Z33): Unter Spin(10) hat ein Einheitsinstanton 4 Eltern- und 12
    Familiennullmoden; mit gesättigten Familiennullmoden entsteht Re[e^{iθ}Φ² C_fam], die Wand ist
    also auch im Kontinuum möglich, bei nackter SO(10) aber ~10⁻¹⁶⁴ unterdrückt. [N]
  - TOE-Lückenkarte (Z21): TFPT.TOE.COMPLETE.01 hat acht Tore, alle [O]; H6 trägt fast nur zu T1
    und T2 bei. Am meisten blockieren die lokale wechselwirkende 3+1D-Theorie (T3/T5) und chirale
    Materie ohne Spiegel (T4).
  - Satzzahl: Gleichförmig drei Sätze widersprechen α unter der Lesart EM1. Offen bleibt der Zweig
    mit zwei verborgenen Sätzen; entscheidend ist c₋/k. [N/C]
  - a/ℓ_P mit vollem Feldinhalt: im H6-Schema etwa 10,8; über Regulatoren hinweg nicht festgelegt. [N]
    Sakharov gegen Jacobson im selben Regler entscheidet nicht (Z27): Der universelle m² ln m-Teil
    stimmt auf beiden Seiten, der Flächenkoeffizient hängt aber stärker von der Metrikkopplung des
    Gitters ab als vom Feldinhalt (R = κ/σ für H6 1,18–2,29 kanonisch, −4,6 bis 2,3 über alle
    Regeln). Robust: Eichfelder brauchen unter D den Randterm, und das Gitter erzeugt einen
    gradientenfreien TT-Term (Gravitonmasse auf Gitterskala). D ist nicht widerlegt. [N]
  - Der KMS-Zustand bei β = L ist eindeutig der Zustand größter Entropie bei mittlerer Energie null
    (Gibbs). [E]
  - Lean (87 Theoreme, ohne sorry): Viertelsektoren, Familiencharaktere, 3+2, ℤ₁₆-Arithmetik,
    X-Ladungen, E₆(i) = 0. Offen ist die Eindeutigkeit, weil die Valenzformel fehlt. [E]
- **Zurückgenommen oder korrigiert:**
  - Ein „Gitter in Planckeinheiten“ gibt es nicht.
  - Die Gabel c = 8 gegen 24 ist nur gegen die Operatorlesart entschieden.
  - Die Familiencharaktere legen die Satzzahl nicht fest.
  - J folgt nicht aus der Modenzahl, sondern aus der Blattstruktur.
  - Die Einheitswindung erzwingt nicht genau eine Naht.
  - Ein Elternoperator mit einem Multiplett liefert kein E₈.
  - Die Wand am Schnitt hängt nicht am Wilson-Regulator allein: Unter G ist sie für jede lokale
    Form ein relevanter Operator (Z24).
  - Ein gemeinsames Multiplett für Naht und Familien (Z18) liefert keine chiralen Familien (Z23).
  - Craigs eichinvariantes Maß (arXiv:2610.09093) deckt den SU(5)-Zweig von TFPT ab, Spin(10) mit
    der 16 aber nicht; eine Übertragung über die Reduktion auf SU(2) benutzt Craig nicht (Z28).
  - Die Wardidentität allein macht a/ℓ_P nicht zur Vorhersage (Z29).
  - Z35s nativer Vergleich ist nur Konsistenz mit einer Diagnosefixture, und sein Rang-1-Kanal bricht
    die Decksymmetrie (Z37).
  - Der 3+2-Schnitt der Halbdrehung (Z10) und der T-Kern (Z16) sind nicht der äquivariante
    V₄-Schnitt des Originals (Z34).
  - Die geometrische S₄-Isotropie aus Z29/Z30 gilt nicht für ein festes Markenquadrat (nur D₄), und der
    Halbgruppenregler ist kein physischer Propagatorschnitt (Z32, Z33).
  - Ein Ladung-4-Kondensat schützt den Nahtindex nicht automatisch; es braucht die anomale axiale
    Antwort, und deren Gitterbruch trägt auch dort die Wand (Z31).
  - Wandfrei ist nicht nur „nichtlokal oder positivitätsverletzend“ (Z19): Lüschers Weylfelder
    sind wandfrei und auf einer eingeschränkten Algebra positiv, wenn man auf χ = (1 − D/2)ψ
    reflektiert (Z26); lokale Reflexionspositivität auf der vollen Algebra haben sie nicht (Lu–Shao).
  - Die ℤ₁₆-Bedingung erzwingt die 16 nur in der partnerfreien Fassung.
  - Die Untergrenze a/ℓ_P ≥ 3,3 gilt nur für die gerechneten Regulatoren.
  - sin²θ_W = 3/8 ist bei H6 nicht induziert mit den Daten verträglich, sondern nur von einer nackten
    SO(10) geerbt (Z25). Nebenbefund: v266 zählt (15,1,1) komplex; mit der einen reellen 45 aus E₈
    und zwei Schleifen sind beide v266-Inhalte von Super-K ausgeschlossen. Unabhängig bestätigt
    durch v1036 (`verification/v1036_ps_proton_real45.py`, Bericht in
    `experiments/theory-contracts/ps-proton-real45-20261124/`); PS.THRESHOLD.02 und PS.PROTON.02
    sind im Ledger korrigiert (24.11.), die abhängigen Zeilen noch nicht.
  - Die Gitterzigarren-Tabellen der achten externen Analyse waren unbelegt; ersetzt durch Z39. Die
    GSO-Auswahl durch die Normaldrehung ist bedingt (Spinlift). „Ein E₈-Block im freien Flussmodell“
    gilt als c = 8, nicht als Spin(10)₁ auf dem Elternoperator (Z39, Z40). Dass eine anomaliefreie U(1)
    die reduzierte Anomalie löscht (H6h20), gilt nur bei ganzzahligem Fluss (Z39, APS).
  - GRAV.KMS.HESSIAN.01 identifiziert den resummierten Spin-2-Hessian mit e^u [C]; die tatsächliche
    Wärmeleitungs-Resummation von Tr e^{−D²/Λ²} gibt a = 1 − 3x/10 + … mit euklidischer Nullstelle (Z41).
    Ledger unverändert; ein Vermerk ist eine eigene Entscheidung.
- **Offen:**
  - Was die primitive Ladung bzw. den Sektor 1 und genau einen windenden Block wählt: Keine
    Grundzustandsregel gefunden; Lesarten sind ein Ordnungsfeld ohne Mittelwert oder N_fam = 3k
    mit LEP. TFPTs Hopf-Lift (v957) mit Deck im Raum verbietet eine einzelne Naht; im Deck-internen
    Zweig erzwingt eine physische Deck-Holonomie auf S³ einen Flussschlauch, um den Φ verschwinden und
    um π winden muss (G folgt), minimal genau einen auf einer Polfaser (Z20). Einen Block
    wählt nur die Forderung „E₈ auf der Naht“ (sonst D₂₄⁺).
  - Warum es neben den drei Familien die 16 schweren Eltern gibt (die direkte Summe ist erzwungen,
    ihr Grund nicht).
  - c₋/k der Naht (Satzzahl).
  - Der gravitative Fluss in Krümmung.
  - Eine Gitter-Wardidentität für die Metrikkopplung bis O(a²∂²): Sie entscheidet Sakharov gegen
    Jacobson und damit a/ℓ_P (Z27); dieselbe Struktur braucht ein Spin-2-Pol (T7).
  - Die Valenzformel in Lean.
  - Nichtperturbative Universalität und eine chirale Gitterformulierung mit Eichfeldintegration;
    Entscheidung 4 ist dieses Problem selbst (Z24). Frei ist ein positiver chiraler Regler gezeigt
    (Z26); offen sind der wechselwirkende Fall und ob H6 eine nichtultralokale Wirkung zulässt.
  - θ₁₃ und δ aus TFPTs Kanälen (sin²θ₁₃ = 0,0231; 1,1σ bis 2,2σ je nach Datensatz laut externer
    Prüfung; δ = 240° ist ein zugewiesener Texturwert).
  - Ein gemeinsamer Elternoperator.
  - Ob H6s Dynamik das Kondensat s^⊗4 und eine unabhängige Spin(10)-Elternverbindung erzeugt, und ob
    ein Ladung-4-Kondensat den Nahtindex c = 8 erhält (das hängt mit Entscheidung 2 zusammen).
  - Für die Flussnaht G′ (Zufluss, Rückwirkung und Gitterzigarre sind gerechnet, Z39): woher Spin(10)₁
    bzw. reelle Nahteltern kommen, wenn der Fluss komplex ist; die Herleitung von E₈-Zufluss und
    Dilatonwirkung aus der zweiten TFPT-Blatthälfte; die Wahl des Spinlifts; ein TFPT-Sektor mit
    M(P) = m(2P − I), der den Rang-1-Kanal stabilisiert, neben einer eigenen schwachen Verbindung (Z40).
  - Ein Operator ∂̄ + ∂̄† auf O(1) ⊕ 6O ⊕ O(−1), geprüft am originalen Randoperator (Calderón-
    Hauptsymbol, relativer Index, Nullmoden mit Ladungen, c_W < 17/2, relative Determinante), und die
    Entscheidung, wie der Spin(10)-Lift (w₂ = 1) aufgelöst wird.
  - Vier positive, kommutierende Translationen aus halbseitigen modularen Inklusionen lokaler
    Typ-III-Algebren (das E₈-Netz liefert eine lichtartige); sie müssen einen gemeinsamen Impuls in
    Herm(2)₊ auslesen ((Σp)² ≥ 3Σp², Z40), freier Kandidat ist der Matrixkern bei ν = 1. Dazu der
    't-Hooft-Vertex mit den
    TFPT-Yukawas (Sättigung der Familiennullmoden, Bezug zu θ).
  - Die Prozessalgebra (Toeplitz/Hardy) als Regler statt eines festen Gitters: Das feste Gitter bricht
    sowohl die axiale Symmetrie (Wand) als auch die Diffeomorphismen (Z24, Z29, Z31). Dazu gehören
    die Identifikation des Hardy-Kerns mit dem P1-Kern und die Unimodularität der Indexabbildung.
  - Eine Regel für den Gesamtkoeffizienten G, Gitterkanten, die die Metrik selbst definiert, und das
    Soldering von Herm(2) bzw. Λ²F an die Raumzeitverschiebungen.
  - Ein positiv-spektraler Spin-2-Sektor (Q²C ≥ Z) statt eines bloß polfreien Formfaktors; die
    GRAV-/QGAMB-Zeilen mit dem ganzen Formfaktor tragen dazu seit 24.11. einen Vermerk (v1037, Z40). Das
    Spektralwirkungs-Graviton mit vollständig monotonem Abschneider ist nicht positiv (Z41); gesucht ist ein
    Spin-2-Sektor in Wallach-Form oder ein anderer Abschneider.
  - Die 4D-Einbettung der reellen Majorana-Chern-Naht (Blatt mit Kante) und ihre μ₄-Kondensation.
  - Die Herkunft der Prinzipien, die die Prozesssuche verlangt (Z43, Z45): konforme Kovarianz des Kerns,
    Orientierungsumkehr, eine gemeinsame Übertragung und die Familienorientierung (TM1 statt TM2 bei 15°); dazu die
    Gaußquantisierung der Quelle als letzte Prämisse der Zustandsregel (seit Z58 für die Ströme aus der quasifreien Originalkovarianz
    gefolgert, die Quasifreiheit selbst bleibt Prämisse; J₋₁ = S⁴ gilt auf den Einteilchen-Marken, auf der Besetzungsleiter ist der Schritt die Polarphase, Z48, Z54),
    die Skala μ, M_c (also G), die quartische Zahl t des metrischen Sterns (die kubische Ordnung legen ungerichtete
    Kanten fest, Z51; t ist die Form der Endpunkt-Quadratform Q, Z52; je Tiefe ein t_n, über die Vierpunktantwort
    bestimmbar, Z55), die Wahl zwischen Band- und Leiteroperator (Z54; seit Z58 durch Rekonstruktion entschieden: Leiter) und die
    Dynamik an der Nahtskala; die Dynamik der
    Spiegellücke in 4D (die mod-16-Bedingung ist erfüllt, Z53); der Mechanismus für TM1 bei 15° (Restsymmetrie allein
    reicht weder mit Majorana- noch mit Dirac-Neutrinos, Z49, Z53).
  - Die eine Prozessregel Z[J, A, g, Φ], aus der Naht, Eichwirkung und metrische Variation gemeinsam
    folgen.

## 4. Wiederkehrende Muster und Verbindungen

- **Dieselbe Vier überall:** vier Marken, μ₄, Z₄ Glue, Nahtnormierung, 4 = 3 + 1.
- **3 und 5 aus einem Objekt:** Familien (3) und Träger (5) aus demselben Divisor.
- **2/3 dreimal:** Recovery, Koide, Nariai (Atlas `pattern:two-thirds-recovery-koide-nariai`).
- **Ein Seed, korrelierte Ausgaben:** φ₀ speist Flavor, Mischung, Ω_b; Flavorleiter und
  Alphabudget treffen sich am Anker (1, 1, 2). Korreliert heißt nicht unabhängig bestätigt.
- **Code ↔ Gitter ↔ Clock:** Hamming → E₈ → G₃₁ → Clocknetz als konkrete Abbildungen.
- **Abschluss ist kein Ursprung:** Closure verdichtet Struktur, beweist aber keine Herkunft.
- **Lücke endlich → physisch:** endliche Algebra braucht eine physische Auswahlbrücke.
- **Auslesung verliert Information:** Momentaufnahmen brauchen Gedächtnis (68 = 5 + 63).

## 5. Was fehlt

- **Herkunft von P1 und P2:** warum ein Randkern mit Einheitswindung, warum fünf Plätze
  mit 3 + 2, warum c₃ = 1/(8π). Das ist der Kern der Suche (Aufgabe F1/F2, `map:universalraum`).
- **Die gemeinsame physische Quelle:** volles wechselwirkendes Feldintegral, sein Spektrum,
  Zustands- und Zeitwahl (`map:state-selection`).
- **Raumzeit:** Metrik, 3+1, lokale Kausalstruktur, Kontinuumslimes (`map:geometry-selection`,
  `map:continuum`, `map:local-net`).
- **Gravitation:** Flächennormierung und lokaler Boost/Horizontfluss (`map:area-normalization`,
  `map:local-boost`).
- **Restflavor:** absolute Quarkamplitude, Neutrinozweig, Phasen.

**P2 folgt aus P1.** Aus der markierten Viermarkennaht von P1 folgen die drei logarithmischen
Richtungen, der Fünferträger H⁰(O(D)) und seine Teilung 3 + 2 exakt; der Anker a = (1,1,2) ist
die einzige Zerlegung von 4 in drei positive Teile und gibt g_car = 5, |ℤ₂| = 2, N_fam = 3
(P2.PARTITION.01, AX.P2.01, Atlas `pattern:fourmark-three-five-normalization`). Die Zuordnung
zur physischen Eichwirkung ist eine eigene Frage innerhalb dieser Herleitung. **Der Kern der
Suche ist damit allein der Ursprung von P1.**

**Was an P1 noch Eingabe ist.** Die flache Sphäre mit zwei Blättern erzwingt vier Marken
(Gauß–Bonnet, Defizit π je Marke) und c₃ = 1/(|ℤ₂|·2π·χ(S²)) = 1/(8π). Offen sind: die
quadratische Lage der Marken (μ₄, τ = i; im Original ausdrücklich zusätzlich genannt; KERNSYM:
einziger symmetrieerzwungener kritischer Punkt, aber Sattel, kein Minimum), die Viertelzählung
in der Wärmespurlesart von c₃, Zustand und Zeit (die operative Saat, bisher bedingt) und warum
die Naht überhaupt eine flache zweiblättrige Fläche ist. Kandidat H6 (Abschnitt 3) leitet diese
Punkte unter den Prämissen A\*, C″, D, E, G her; offen bleibt dort vor allem, was den Defekt im
Grundzustand erzwingt.

**Was nicht fehlt:** Die vier R113 Endpunktwerte (Z0, ZU, ZY, ZUY) der ursprünglichen
Quellverträge sind ein Test, ob *diese* Originalquelle ihre markierten Antworten liefert.
Diese Route ist blockiert; das ist ein eingeschränkter Befund über sie, keine Lücke der
Theorie und keine Voraussetzung für einen anderen Ursprung (Entscheidung vom
23. November 2026, `AGENTS.md`).

## 6. Was ausgeschlossen ist (immer mit Geltungsbereich)

Der Atlas führt rund 345 eingeschränkte Ausschlüsse, der Prüfkatalog das Blatt
„Gegenprüfungen“. Beispiele: zustandserhaltende gaußsche oder Bogoliubov-Quellenabbildung;
klassische konstante 4D Wärmeproduktklasse (scheitert am Vakuum/Einstein-Verhältnis);
projektive GHZ Nahtkonsistenz; freie Standardmodellmatrizen als Suchgröße für Yukawa. Ein
Ausschluss gilt nur unter seinen genauen Voraussetzungen; eine neue Voraussetzung öffnet
den Weg wieder.

## 7. Wo mehr steht

- **Atlas** (MCP `tfpt-knowledge`): `theory_overview` (Karte, 21 Muster), `research_state`,
  `search_knowledge`, `list_attempts`, `get_dossier`; Web http://127.0.0.1:8731.
- **Prüfkatalog** `docs/catalog/TFPT_Pruefkatalog_2026-11-20.xlsx`: Verifikationen (707),
  Lean Beweise, Lean Kernelexport, Messungen und Suchziele, Atlasbefunde, Versuchsprojekte,
  Versuchsskripte, Einzelprüfungen, Wolfram, Gegenprüfungen, frühere TFPT Versuche.
- **Papiere:** `documents/tfpt_1_architecture_e8` (Eingaben, E₈, Alpha),
  `tfpt_2_standard_model` (Teilchen, Massen, Mischung), `tfpt_4_frontier`, `tfpt_5_redteam`,
  `origin_theory`. Ledger: `verification/status_ledger.csv`; Vorhersagen:
  `verification/predictions_frozen.json`.
- **Gesamtbild der Prozessidee:** `docs/research/fundamental-plan/GESAMTBILD.md`.


### Integration of the certified research branch and further exploration (24 November)

The H6h27 research branch (252bac7a) is integrated with versioned certified replay, preserving the existing H6h29 main lane and its historical outputs. See `experiments/exploration/2026-11-24-h6h26-source-dynamics/` and OriginLab README. New independent exploratory checks in `experiments/exploration/2026-11-24-pair-and-refinement/` reran with 43/43 and 87/87 grouped checks. Local quartic mirror parents, clock compressions, and diagonal geodesic refinement are candidates, not physical T1–T8 closure. Ledger and frozen predictions unchanged; generated registry and Atlas not refreshed in the integration session.
