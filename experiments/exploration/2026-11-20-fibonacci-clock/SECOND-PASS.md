# TFPT: Ursprung aus minimalen binären Geschichten — unabhängige Folgerunde

**Stand: 20.11.2026, keine vollständige TFPT, F1 offen.**  Alle exakten Programme in diesem Paket geprüft, ohne Parameterfits. Diese Sitzung beginnt nicht bei P1/P2 oder den TFPT Clocks, sondern bei minimalen zulässigen Übergangsregeln und deren unitärer Kategorifizierung.

## 1. Minimalregel und ehrliche erste Wahl

Hypothese: Das ursprüngliche Etwas ist eine binäre Unterscheidung, die bei iterativer Wiederholung eine **primitive** (zusammenhängende und aperiodische), positive-Entropie-Geschichte bildet. Wähle unter solchen Regeln diejenige mit minimalen erlaubten Übergängen. Vollständige Enumeration aller 16 binären 2×2 Matrizen: drei primitive gelabelte Matrizen, zwei Isomorphieklassen. Minimal (drei Kanten) ist `A_F=[[0,1],[1,1]]` bis auf Zustandsumbenennung; die einzige Hinzunahme der verbotenen Kante gibt `A_2=[[1,1],[1,1]]`.

Die Wahl „zwei Zustände, primitiv, minimal drei Kanten“ ist **kein** von der Natur abgeleiteter Satz. Sie ist die erste freie *Regel*, aber keine frei angepasste Zahl und keine TFPT-Eingabe. Die physische Realisierung auch der einen möglichen Komplettierung wird nicht schon von der Minimalregel erzwungen. Unter dieser Regel erscheinen die Eigenwerte φ und 2 ohne Fit.

## 2. Unendliche Geschichte und modularer Zustand

Für die jeweils zusätzlich gewählte Cuntz–Krieger-Algebra O_A mit kanonischer Wortlängen-Gaugeaktion gilt nach Enomoto–Fujii–Watatani (primitive Graphen): einzigartiger KMS-Zustand bei β=ln ρ(A). Die positive normierte Eigenvektorlösung des Fibonacci Graphen ist `p=(φ⁻²,φ⁻¹)`, mit exakter Konsistenz aller Zylindermaßen. Der vollständige Binärgraph hat β=ln2.

**Wichtig:** die Produktzustände sind *nicht* KMS bei einer gemeinsamen β für die unskaliert-diagonale Gaugeaktion, weil lnφ ≠ ln2. Für eine KMS_1-Dynamik muss man `α_t=γ^Fib_{(lnφ)t} ⊗ γ^2_{(ln2)t}` wählen. Diese relative Skalierung stammt aus den Spektralradien, ist aber eine zusätzliche Normierungsentscheidung. Der Produktmodularfluss ist dann automatisch nichtperiodisch. Die GNS-Faktoren sind Typ III_{1/φ} (primitive Periode 1) und III_{1/2}; ihr Tensorprodukt ist Typ III_1, da lnφ/ln2 irrational. Exakter Beweis: φ^n=F_n φ+F_{n-1} ist für n>0 irrational, sodass φ^n ≠ 2^m für ganze positive m,n. Damit ist `{φ^m2^n:m,n∈Z}` dicht in R_{>0}.

**Grenze:** Ein Typ-III_1-Faktor allein ist keine lokalisierte relativistische Feldtheorie. Seine Einbettung in ein lokal-kausales Netz, positiver Energiefluss, nichtperiodischer physischer Lorentzboost, echte Raumzeit und Gravitation sind weiterhin offen. Die einzelne Gaugewirkung ist periodisch und ihr Liouville-Spektrum ist zweiseitig, keine bereits physische positive Zeittranslation.

## 3. Dieselbe Fibonacci Fusion und fünf Quantenereignisse

Zusätzliche freie Wahl: unitäre braidgefaltete Kategorifizierung `τ⊗τ=1⊕τ` sowie fünf τ-Objekte mit Gesamtladung τ. Das Fusionsraum-Hilbertsystem hat exakt Dimension 5. Wird ein Ereignis als Umgebungspunkt behandelt und nur die anderen vier verflochten, zerfällt das System in Gesamtladungsräume 2+3. Die exakten Jones/Temperley–Lieb Projektoren mit `a=φ^-1`, `a²+a−1=0` erfüllen die gesamten Idempotenz-, Fernkommutativitäts- und TL-Relationen.

Der Rang der durch die ersten drei Vertauschungen erzeugten komplexen Algebra ist 13, entsprechend M2⊕M3; mit der vierten Vertauschung ist er 25, entsprechend M5. Die Ränge sind **exakte Charakteristik-Null-Zertifikate** über zulässige Spezialisierungen `a mod 11=3`, `a mod19=4`, `a mod29=5`. Überall Rang (13,25). Jeder Nichtnull-Minor nach Reduktion ist ein Nichtnullpolynom im algebraischen Feld und belegt eine Charakteristik-Null-Untergrenze; 2²+3² bzw. 5² geben die matching Obergrenzen. Diese Rechnung nutzt keinerlei E8, P1 oder P2 als Input.

Die algebraisch erlaubte determinant-eins Gruppe der Blockzerlegung ist `S(U3×U2) ≅ (SU3×SU2×U1)/Z6`, mit spurfreier relativer Ladung `diag(-1/3×3,+1/2×2)`. **Das ist der Stabilisator eines 3+2-Projektors, nicht bereits eine physikalisch lokale Eichtheorie.** Die U1 ergibt sich nicht als dichter Braidgenerator. Fünf τ-Fusionsgeschichten sind keine fünf lokalen CAR-Fermion-Einteilchenmoden. Die 16 Spinor-Komponenten bleiben eine nur unter zusätzlichem CAR-Funktor gültige Folgerung.

Auch **drei Familien dürfen nicht mit dem Farbtripel identifiziert werden**: auf einer irreduziblen SU3-Fundamentaldarstellung sind alle unabhängigen mit SU3 kommutierenden Operationen nach Schur skalar, während der braidgewonnene A3-Zyklus nicht skalar ist. Eine Familie und ein Farbtripel brauchen grundsätzlich unterschiedliche Freiheitsgrade.

## 4. E8 ohne D5+A3 als Start

Mit einem **zusätzlich begründungsbedürftigen** zweiten positiven A4-Ladungsgitter kann man A4⊕A4 über den Z5-Diskriminanten gluen. Exakt zwei der vier Isomorphismen `k→mk` mit m=1..4 erlauben gerade unimodulare Verklebung: m=2 und m=3. Für jede ist die Wurzelzählung `240=20+20+50+50+50+50`. Explizite Gram-Matrix von acht ganzzahligen Generatoren hat alle positiven Hauptminoren `(2,3,4,8,12,16,20,1)`, Determinante 1 und gerade Diagonale; damit ist das Gitter E8.

Das ist eine **bekannte A4+A4-Realisation von E8**, neu hier nur als alternativer Anschluss an die Fibonacci Fünferalgebra. Es ist nicht gezeigt, dass Reflexion zwei positive A4-Gitter produziert, die tatsächliche Quelle primitiv integral gluet oder die resultierenden Gitter zu einer chiralen E8_1-Quellenalgebra quantisiert werden. Aus beliebigen Anyondaten kann die absolute chiral-gravitative Ladung c=8 nicht folgen, da invertibles E8-Stapeln die Anyonkategorie unverändert lässt und c um 8 verschiebt.

## 5. No-gos: Zeitwahl und positive Antwort

**Schur-Obstruktion:** Bei vollständiger B5-Invarianz muss jeder treue stationäre Zustand auf M5 mit allen Generatoren kommutieren, also ρ=I5/5. Der Tomita–Takesaki-Modularfluss auf M5 ist dann identisch. Bei bloßer B4-Invarianz bleibt eine freie kontinuierliche Sektorengewichtung: `ρ=(p2/2) I2 ⊕ ((1-p2)/3) I3`. Auch das bestimmt keine physische Clock.

**Unabhängiger Readout:** Gleichwahrscheinliche lokale Braidbewegung und Projektion auf denselben Fusionbaum liefern einen verbundenen symmetrischen 5-Knoten Markovkanal `T=I−(1/(4φ))L`. Sein Laplacepolynom lautet `x(x−2)(x³−8x²+18x−10)`. Sein einziger stationärer Zustand ist I5/5, mit Sektorwahrscheinlichkeiten 2/5 und 3/5. Diese sind **nicht** die unendlichen KMS-Gewichte φ⁻² und φ⁻¹. Eine physische Auswahl und Konsistenzbrücke ist notwendig.

**Reine Finitalgebra kann c3 nicht geben:** Alle aus den algebraischen F/R-Daten mit endlich vielen algebraischen Operationen gebildeten Skalarwerte bleiben algebraisch. `1/(8π)` und `φ0` sind transzendent; eine echte analytische Kontinuums-/Winkel-/Integralkonstruktion ist notwendig. Der bekannte WZW-Level-1-Koeffizient `k/(8π)` trifft bei k=1 dieselbe Zahl, ist aber unter anderer Funktionalrolle und Trace-Normierung kein Beweis für den originalen TFPT Randkern. Der originale TFPT c3-Text führt zudem andere konditionale Normierungsargumente an (`documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex`), die nicht durch das Fibonacci Modell ersetzt werden.

**Kein automatischer 3+1D Lift:** Ortsfeste Fibonacci Anyonen sind natürliche 2+1D Objekte. Walker–Wang oder verwandte 3+1D Bulk-Modelle erlauben topologische Erzeugung bei **vorausgesetzter** Dimension, aber keine bereits emergente dynamische Metrik und keinen Graviton.

## 6. Nullmodelle und Freiheitsbilanz

| Zahl/Struktur | Basislinie im vorab deklarierten Universum | Einordnung |
|---|---|---|
| Min. binärer primitiver Graph | 3/16 primitive Matrizen, 2 Isomorphieklassen; genau 1 minimale | Präzise endliche Klassifikation, Minimalitätswahl extra |
| Unabhängige Clock-Paarung | 2/4 geordnete Paare der zwei Graphklassen sind gemischt und inkommensurabel | Physische Produktbildung extra |
| Fünferalgebra | Rang M2⊕M3=13 und M5=25 exakt aus F-Projektoren; keine freie Ganzzahl | Braiding der 5 τ gesondert gewählt |
| 3+2-Eichgruppe | jeder Rang-3-Projektor auf C5 hat abstrakt denselben Stabilisator | Keine unabhängige Physikbestätigung |
| E8 A4-Doppel | 2/4 diskriminante automorphe Kleber liefern gerades Gitter | Doppelung, Positivität, Integralkleber extra |
| c3 1/(8π) | standardmäßig auch WZW k=1/GHY Planck-Konventionen | ohne gleiche Quellidentifikation Zahlentreffer schwach |

## 7. Urteil und nächster Härtetest

**Stufe 1: KANDIDAT** (neue unabhängige algebraische Ursprungskette, mit exakten Nullmodellen). **Stufe 2 bis 4: offen/nicht erreicht.** Gegenüber den bereits katalogisierten Versuchen ist kein universeller Überlegenheitsbeweis möglich. Wir haben zwei zusätzliche echte Zwischensätze: exakte Fünfer-Sektorverbindung und aus *vollständiger* binärer 2×2-Primitivität folgende Paarung zweier inkommensurabler KMS-Clocks. Das ist keine vollständige physische TFPT.

Entscheidendes neues Experiment: **eine einzige Zustands-/Lokalitätsbrücke** muss (i) die hier differierenden endlichen Readout- und unendlichen KMS-Zustände als Auslesungen *derselben* positiven Historiendynamik realisieren, (ii) eine nichttriviale halbseitige modulare Inklusion und einen Generator positiver Energie liefern, (iii) den relevanten lokalen CAR-Index der 3+1D-Materie inklusive drei unabhängiger Familien erzeugen, (iv) die zweite A4-Faser und integrale E8-Quelle bestimmen und (v) die absolute ursprüngliche P1-Response samt c3, Raumzeit und Gravitation wirklich rekonstruieren. Scheitert eine dieser verlangten Identifikationen, ist die Behauptung eines einzigen Urprozesses zu verwerfen.

## Quellen

* https://www.aidansims.com/papers/aHLRS2012.pdf (Gaugedynamik und KMS)
* https://arxiv.org/abs/1409.6481 (Faktortypen höherer Graphalgebren)
* https://academic.oup.com/nsr/article/3/1/68/2460241 (E8-Stapelung, zentrale Ladung mod 8)
* https://arxiv.org/abs/math-ph/0411058 (Type III1 in relativistischer Quantenfeldtheorie)
* https://link.springer.com/article/10.1007/JHEP08%282024%29084 (Fibonacci Zöpfe)
* https://www.math.auckland.ac.nz/~hekmati/Students/MarcLescano.pdf (Fusions- und Pentagon-Daten)
* Repo: `research/originlab/P1-FROM-PHYSICAL-NECESSITY-20261116.md`, `research/originlab/P1-EQUIVARIANT-RESIDUE-ORIGIN-20261115.md`, `research/originlab/F1-E8-ANOMALY-INFLOW-ORIGIN-20261116.md`, `documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex`.