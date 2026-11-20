Gebiete: Operatoralgebren und Selbstbezug · Quanteninformation und endliche Geometrie · Quantenfehlerkorrektur · topologische Phasen und Kobordismus · arithmetische Geometrie über endlichen Körpern · Extremalprinzipien · Codes als Fixpunkte der Dualität · Quantenzellautomaten · Kausalmengen

# Exploration 20. November 2026: das kleinste zusammengesetzte Quantensystem

Erste Sitzung im Explorationsordner. Kein Ansatz startet bei Randkern, vier Marken, D₅ + A₃,
3 + 2, Clocknetz oder φ₀. Der Atlas wurde erst nach der Ideenphase befragt (Abschnitt
„Abgleich“ in `VERDICT.md`).

## Ansätze

| Nr | Gebiet | Grundprinzip | Soll liefern | Messgröße | Freie Parameter | Erste freie Wahl | Gerechnet |
|---|---|---|---|---|---|---|---|
| A1 | Operatoralgebren, Quanteninformation, endliche Geometrie | Die kleinste nichtkommutative Algebra M₂(ℂ) beschreibt sich selbst. Ihre Standardform ist ℂ² ⊗ ℂ²: System und Spiegel, mit Zustand, Spiegelung J und modularer Zeit. Gefragt wird, was die vollständige unverzerrte Beschreibung dieses Systems erzwingt. | 5 Plätze mit markierter Teilung 3 + 2, Λ^gerade ℂ⁵ = 16 mit Hyperladungen, μ₄, E₈, Cl(1,3) mit 3 + 1, Zustand und Zeit | 16 Fermionen je Familie, Eichgruppe, 3 + 1 | Gewichte des treuen Zustands (ein Verhältnis) | Lesart „Platz = antikommutierender Erzeuger“ | `a1_zwei_qubits.py` |
| A2 | Quantenfehlerkorrektur | Perfektion: ein Quantenfreiheitsgrad wird gegen jeden lokalen Einzelfehler geschützt und kein Syndrom bleibt ungenutzt. | 5 Plätze, 16 = 1 + 3·5, 4 Prüfungen, Schwelle 3 + 2, Träger über 𝔽₄ | 16 Fermionen je Familie | keine | Markierung einer der zehn gleichwertigen 3 + 2 Teilungen | `a2_perfekter_code.py` |
| A3 | Topologische Phasen, Kobordismus | Die Theorie ist Rand einer invertierbaren reflexionspositiven Theorie; mit einer ℤ₄, deren Quadrat die Fermionparität ist, muss die η Anomalie verschwinden. | 16 Weylfermionen genau in vier Dimensionen, ν_R, μ₄ als Zentrum von Spin(10) | 16 mit ν_R, 3 + 1 | keine | Dimension 4 oder die Zahl 16 (eine von beiden ist Eingabe) | `a3_anomalie_eta.py` |
| A4 | Arithmetische Geometrie über 𝔽₂ | Alles ist über dem kleinsten Körper definiert; Plätze sind Punkte, Zählungen sind Frobenius Fixpunkte. | 5 = 3 + 2 (ℙ¹(𝔽₄) unter Frobenius; fünf elliptische Kurven über 𝔽₂), 5 und 41 als Punktzahlen, 68 | Materiekoeffizient 41 | keine | Auswahl der Kurve und des Erweiterungsgrades | `a4_endliche_koerper.py` |
| A5 | Extremalprinzip | Der Grundzustand minimiert jede vollständig monotone Paarenergie zugleich (universelle Optimalität, Cohn, Kumar, Miller, Radchenko, Viazovska). | E₈ in Dimension 8, Leech in 24 = 3 · 8 | keine direkt | keine | Dimension 8 statt 24 | nein |
| A6 | Selbstbezug, Codes | Eine binäre Beschreibung, die ihre eigene Prüfung ist (C = C^⊥) und konsistente Phasen hat (Gewichte in 4ℤ), in kleinster Länge. | Länge 8, eindeutiger Code [8,4,4], E₈ | keine direkt | keine | keine bis E₈; danach die Zerlegung in D₅ ⊕ A₃ | `a6_selbstdual.py` |
| A7 | Quantenzellautomaten | Lokalität, Homogenität, Isotropie und Unitarität bei kleinster innerer Dimension (D'Ariano, Perinotti): der Weyl Automat, zwei gekoppelt geben Dirac. | Dynamik und Kontinuumslimes der vier Komponenten aus A1 | 3 + 1 Dispersion | keine nach Wahl des Gitters | Raumdimension 3 des Gitters | nein |
| A8 | Kausalmengen | Raumzeit ist eine lokal endliche Halbordnung; Dimension und Volumen folgen aus Zählungen. | 3 + 1, Kausalstruktur, Flächengesetz | 3 + 1 | Streudichte | Dimension der Einbettung | nein |

## Der Zusammenhang der gerechneten Ansätze

A1, A2, A3 und A6 beschreiben dasselbe Objekt von vier Seiten. Das ist das Ergebnis der
Rechnungen und war nicht die Ausgangsannahme.

1. **Selbstbeschreibung.** M₂(ℂ) wirkt auf sich selbst; der Darstellungsraum ist ℂ⁴ = ℂ² ⊗ ℂ².
   Die Zerlegung in zwei Faktoren ist dabei gegeben (Algebra und Kommutant), ebenso ein treuer
   Zustand, die Spiegelung J und der modulare Fluss. Die Nebendiagonalen tragen modulare Ladung
   ±1 (A1, Teil 5).
2. **Vollständige unverzerrte Beschreibung.** Die 15 Paulioperatoren zerfallen in fünf
   kommutierende Tripel (fünf unverzerrte Basen) oder dual in Sätze von fünf antikommutierenden
   Operatoren. Beide Fünfen teilen sich zwangsläufig in 3 + 2: drei Produktbasen und zwei
   verschränkte, beziehungsweise drei nichtlokale und zwei lokale Erzeuger (A1, Teil 1).
3. **Perfekter Code.** Die Gleichung D^(n−k) = 1 + (D² − 1)n hat mit einem logischen
   Freiheitsgrad nur die Lösung D = 2, n = 5. Die Syndrome der fünf Plätze bilden genau einen
   solchen Satz von fünf Basen eines Zwei-Qubit-Registers (A2, Teil 2).
4. **Operatoralgebra.** Die geraden Monome der fünf Erzeuger sind eine Basis von M₄(ℂ):
   1 + 10 + 5 = 16 = Λ^gerade ℂ⁵. Die spurfreie Klassenfunktion der 3 + 2 Teilung gibt die
   Hyperladungen einer Standardmodellgeneration (A1, Teil 2).
5. **Dirac Algebra.** Dieselben fünf Erzeuger sind γ⁰, γ⁵ (lokal) und γ¹, γ², γ³ (nichtlokal).
   Die 3 + 2 Teilung der Plätze ist die Teilung in drei Raumrichtungen und (Zeit, Chiralität)
   (A1, Teil 3).
6. **Gitter.** Die Cliffordgruppe über ℚ(i) hat Ordnung 46080 und Zentrum μ₄; ihr invariantes
   Gitter hat 240 kürzeste Vektoren, das sind die 60 Stabilisatorzustände mal μ₄. Die fünf
   Erzeuger zerlegen sie in 5 × 48, die Hopfabbildung schickt sie auf die zehn Vektorgewichte
   von D₅ mit Faser 24-Zell (A1, Teil 4). Derselbe Gittertyp entsteht aus dem eindeutigen
   kleinsten selbstdualen doppelt geraden Code (A6).
7. **Anomalie.** Die Phase i wirkt als Skalar auf Λ^gerade ℂ⁵. In vier Dimensionen ist die
   η Invariante eines solchen Weylfermions −1/16; sechzehn heben sich auf, fünfzehn nicht (A3).

## Fortsetzung am selben Tag: die drei entscheidenden Tests

Nach der ersten Runde blieben drei Tests offen (Wörterbuch, Glue, Familien). Sie sind gerechnet.
Der Abgleich mit dem Korpus wurde dabei vertieft und korrigiert die erste Runde: TFPT führt die
Zwei-Qubit-Struktur des Trägers bereits selbst (v752, v774, v775, v783, v798, v799, v845),
einschließlich eines gefeuerten Negativbefunds. Neu ist deshalb nur, was unten steht.

| Nr | Frage | Ergebnis | Skript |
|---|---|---|---|
| B1 | Warum M₂(ℂ) und warum D₅ ⊕ A₃? | Drei Forderungen an ein n-Niveau-System mit m Plätzen (m = n + 1; n² = 2^(m−1); Klebenorm m/4 + (n−1)/n = 2) haben genau die gemeinsame Lösung n = 4, m = 5. Je zwei davon, die (iii) enthalten, reichen. n = k² gibt k = 2. Die Gitterbedingung allein lässt D₅, D₁₃, D₂₁ mit A₃ zu. | `b1_klebegleichung.py` |
| B2 | Ist die TFPT Karte die Standardform von M₂(ℂ)? Wo liegt D₅ ⊕ A₃? | Ja, Koordinate für Koordinate: Familienpaare = X, Y, Z Komponenten des Operators, Ankerpaar = Einskomponente (Zustandsgerade ℂΩ), J = Multiplikation mit i, σ = Drehung X → Y → Z als u ⊗ ū, äußere Spiegelung = a ↦ a*. Von 7560 Lagen von A₃ sind 18 σ-invariant, in genau zwei Bahnen (12 + 6). In allen 18 dreht σ drei der fünf D₅ Achsen und lässt zwei fest. Die 6 enthalten eine komplexe Gerade und sind die σ-stabilen Teilungen des Korpus. | `b2_woerterbuch.py` |
| B3 | Feuern die Gegenproben? | Ja, alle acht: falsche Ankerachse (112 von 240 Wurzeln), Rechenbasis statt Paulikomponenten (0 von 240), einseitige Drehung (240 von 240 weichen ab, Ω nicht fest), vier falsche Paare (m, n). | `b3_gegenproben.py` |
| B4 | Sind die fünf Plätze fünf unverzerrte Basen, und ist 3 + 2 die Teilung Produkt / verschränkt? | Ja in dieser Karte: σ ist lokal, lässt genau drei Kontexte fest (alle verschränkt, darunter die Bellbasis = Ankerklasse) und dreht in jedem der drei σ-invarianten Sätze die 3 Produktbasen, während die 2 verschränkten fest bleiben. Klasse in L/(1+i)L = Kontext; ħ = 0 genau bei gemeinsamem Pauli. | `b4_fuenf_basen.py` |

### Die Kette nach der Fortsetzung

1. Die kleinste nichtkommutative Algebra M₂(ℂ) wirkt auf sich selbst: ℂ⁴ = ℂ² ⊗ ℂ² mit Zustand
   Ω, Spiegelung a ↦ a* und der Phase i.
2. In Paulikomponenten ist das genau die Karte, die TFPT für E₈ benutzt (B2). Deck σ,
   Orientierung J, Ankerpaar und die Länge 8 sind damit Teile eines Objekts und keine getrennten
   Eingaben mehr.
3. Die vollständige unverzerrte Beschreibung hat n + 1 = 5 Basen; sie zerfallen erzwungen in
   3 Produktbasen und 2 verschränkte. Das ist die Wirkung von σ auf den fünf Trägerplätzen (B4).
4. Operatoralgebra (16) mal Zustandsraum (4) ist die Klebeklasse von D₅ ⊕ A₃; (4, 5) ist die
   einzige Lösung der drei Forderungen (B1). Unter σ ist 4 = 1 + 3 (B2, Teil 5).
5. Der Anker (1, 1, 2) aus AX.P1.01 sind die Normen von 1, i, 1 + i auf der Zustandsgeraden.

## Was die Kette nicht liefert

Die Ganzzahligkeit (warum das Gitter der Stabilisatorzustände), den Faktor 2π in
c₃ = 1/(2 · 4 · π), eine nichttriviale modulare Zeit (für den Spurzustand Ω ist der Fluss
trivial; ein anderes Gewicht wäre ein freier Parameter), die kontinuierliche Eichgruppe, φ₀, 41,
α, Massen, Metrik, Kontinuum und Gravitation. Die physische Lesart „Klasse = Materiezustand“
bleibt durch ARF.ROOTCLASS.01 (v775) ausgeschlossen; die Hyperladungstabelle aus A1 ist die
Identität aus v774 und kein neuer Treffer.
