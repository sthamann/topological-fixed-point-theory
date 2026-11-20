# Neuheitsprüfung der Ursprungsideen

Stand: 20. November 2026. Geprüfter Hauptzweig: `a1273187b82488d34d2a474470152a41d42d4055`, Commitzeit `2026-11-20T14:02:49Z`.

## Zugriff und Reichweite

Die aktuellen Dateien `docs/research/TFPT-STAND.md`, `docs/research/EXPLORATION.md`, `AGENTS.md`, `README.md` und `docs/registry/INDEX.md` wurden direkt aus dem angegebenen Repository gelesen. Der Hauptzweig wurde über die GitHub Branch API verifiziert. Das örtliche ältere Checkout `/workspace/scratch/da205a085916/tfpt` steht dagegen auf `ef937aeb916bd7f58b85c7fc2575d9adfe265f42` vom 16. November; es ist nur als Suchhilfe für historische Originalquellen verwendet worden.

In der verfügbaren Werkzeugliste gibt es keine Werkzeuge `tfpt-knowledge`, `theory_overview`, `search_knowledge`, `list_attempts` oder `get_dossier`. Der aktuelle lokale Atlas wurde daher nicht live konsultiert. Seine im Prüfkatalog exportierten Befunde sind verfügbar. Die Aussage „kein Treffer“ bedeutet hier „kein passender Treffer in den geprüften Dateien und Suchbegriffen“, niemals „dies wurde noch nie versucht“.

Der aktuelle Ordner `experiments/exploration/` enthält nur `README.md`. Dort liegen noch keine älteren Explorationsordner; relevante frühere Arbeit liegt vor allem unter `experiments/theory-contracts/`, `research/originlab/` und `verification/`.

Der aktuelle Prüfkatalog wurde vollständig heruntergeladen, als XLSX mit XML gelesen und auf allen zwölf Blättern durchsucht. Dateigröße: 3.665.827 Byte. Die folgenden Zeilenzahlen schließen Überschriften und Einleitungszeilen ein und sind deshalb keine Zahlen wissenschaftlicher Resultate.

| Blatt | XML Zeilen |
|---|---:|
| Überblick | 23 |
| Verifikationen | 710 |
| Messungen und Suchziele | 130 |
| Versuchsprojekte | 485 |
| Lean Beweise | 869 |
| Lean Kernelexport | 797 |
| Atlasbefunde | 2020 |
| Versuchsskripte | 2230 |
| Einzelprüfungen | 4864 |
| Wolfram | 7 |
| Gegenprüfungen | 28 |
| Frühere TFPT Versuche | 260 |

Die vollständigen lexikalischen Treffer stehen in `catalogue_hits.json` und `catalogue_extra_hits.json`. Sie enthalten auch falsche Treffer, beispielsweise Jordan Wigner statt Jordan Algebra, Frobenius Norm statt endlicher Körper und M2/M3 als Massen oder Momente. Die Beurteilungen unten trennen diese ausdrücklich.

## Ergebnis nach Ansatz

| Ansatz | Tatsächlicher Bestand | Neuheitsurteil und genaue Grenze |
|---|---|---|
| Kleinster perfekter Quantencode mit Korrektur eines beliebigen lokalen Fehlers | Die Begriffe `[[5,1,3]]`, perfect code, Singleton und entsprechende deutsche Varianten liefern keinen passenden Eintrag. Zwei Fünfqubit Treffer betreffen die Kraus Dilatation des schon vorgegebenen Erholungskerns, Projekt `qc-recovery-kernel`. | Die vorgeschlagene Minimalität aus Fehlerkorrektur ist in der geprüften Erfassung nicht als Ursprungstest dokumentiert. Sie ist von der vorhandenen Fünfqubit Schaltung zu unterscheiden. Keine weltweite Neuheit des bekannten Codes behaupten. |
| Operationaler Kegel und SO(d) äquivariante antisymmetrische Dynamik | Lorentzkegel und Choi Abbildungen sind Bestand: `compiler-cone-object-audit` und `compiler-process-reconstruction-20260918`. Allgemeine Komposition, Kausalität und lokale Tomographie wählen die Algebra nach einem expliziten Gegenmodell nicht aus. | Die zusätzliche Auswahl durch eine nichtverschwindende, voll SO(d) äquivariante alternierende Abbildung von zwei Observablen zu einem Generator ist in den geprüften Treffern nicht dokumentiert. Diese stärkere Voraussetzung muss genannt werden. Der Lorentzkegel allein ist kein neuer Befund und keine Raumzeit. |
| Kleinste bosonische invertierbare chirale Phase | Exakt derselbe Ansatz ist in `F1-E8-ANOMALY-INFLOW-ORIGIN-20261116.md`, Abschnitt 3, formuliert: positive gerade unimodulare Ladungsgitter, kleinster nichtnull Rang, E8 als eindeutiges Rang acht Gitter. Der Text nennt ausdrücklich die Unabhängigkeit von der vorherigen Wahl sechzehn Majoranas. | Direkte Wiederholung eines bekannten bedingten Ansatzes. Ein neuer Cartanbaum Zensus verbessert seine Basislinie oder Reproduktion, ist aber kein neuer Ursprungsmechanismus. Die Quelle warnt ausdrücklich vor der zusätzlichen 2+1 dimensionalen Phasenannahme und dem fehlenden Übergang zu 3+1. |
| Kleinste Fixpunktalgebra mit zwei nichtisomorphen nichtkommutativen einfachen Sektoren | Bedingte Erwartungen und Fixpunktalgebren kommen vor, unter anderem bei μ4 Mittelung, E6+A2 und Familienrückwegen. Kein Treffer beschreibt die Minimierung zu M2(C) ⊕ M3(C). Breite M2/M3 Suchen liefern überwiegend Massen, Momente und Meilensteinbezeichnungen. | Als Selektionsvorschlag in der geprüften Erfassung neu. Die zentrale freie Wahl ist jedoch gerade die Forderung zweier nichtisomorpher nichtkommutativer Summanden. Zwei lediglich verschiedene Sektoren erlauben bereits M2 ⊕ M2. Die minimale Dimension fünf ist deshalb nur bedingt erzeugt. |
| Kausal invariante Automaten und Umschreibungen | `v298_e8_network_attractor` hat den affinen E8 Graphen bereits als Eingabe. `v299_hypergraph_0d_autonomous` erzeugt E6→E7→E8 unter dem ausdrücklich gesetzten Dreiarmsaattyp (2,3,·) und Wachstum des längsten Arms. `automaton_multiway_tower_probe.py` ist registriert. `foreign_universes.py` durchsucht inzwischen kleine Graphen, elementare Automaten und reguläre Polytope. | Der Bereich ist bekannt. Neu wäre eine quellenfreie Regel mit einer belegten Auswahl von Dimension, Materie und Dynamik. Die früheren negativen A_n und D_n Kontrollen schließen nur den bezeichneten Wachstumsselektor aus, keine beliebigen kausalen Quantensysteme. |
| Selbstbeschreibung oder minimale Beschreibungslänge | Keine passende wissenschaftliche Ursprungsrechnung in den geprüften Treffern. Der einzige Selbstreferenztreffer betrifft die technische Integrität des Ledgers. | Kein dokumentierter passender Treffer. Ohne feste universelle Beschreibungssprache, Kodierung, Maß und Auswahlregel besteht aber kein eindeutiges Extremalprinzip. Ein fehlender Treffer erzeugt noch keinen Kandidaten. |
| Zustandsummen mit Pachner Konsistenz | Keine passenden Treffer zu Pachner, Turaev Viro, Crane Yetter oder Triangulationsunabhängigkeit. Es existieren zahlreiche andere Frobenius und Klebekonstruktionen. Ein lexikalischer `state_sum` Treffer ist nur eine Erwartungswertsumme. | Der konkret formulierte Ansatz ist in der geprüften Erfassung nicht gefunden. Die schon gesetzte Dimension der Mannigfaltigkeiten, die Kategorie und die Gewichte müssten als freie Eingaben offengelegt werden. Allgemeine Klebeverträglichkeit liefert diese Auswahl noch nicht. |
| Endliche Körper und elliptische Frobeniusdynamik | Der aktuelle TFPT Stand nennt bereits den F2 Befund y²+y=x³+x+1 und 1+i. `check_h10_frobenius.py` vom 19. November prüfte außerdem Spuren einer CM Kurve gegen Zufallskurven und fand keinen Vorsprung. | Die spezielle F2 Kurve erneut zu verwenden wäre eine Fortsetzung des Bestands. Der frühere negative CM Spurtest schließt nur die dort gewählten Spuren als Ereignisprotokoll aus. Ein anderes endliches Körpermodell ist dadurch nicht generell verworfen, benötigt aber eine unabhängige Kurvenauswahl. |

## Vier entscheidende Originalquellen

### 1. E8 aus einer invertierbaren Phase ist bereits vorhanden

Die Originalnotiz vom 16. November sagt unter den expliziten mikroskopischen Axiomen einer bosonischen abelschen Chern Simons Phase mit positivem geradem unimodularem Gitter und minimalem nichtnull Rang genau E8 voraus. Sie hält zwei noch zusätzliche Entscheidungen fest: die physische Quelle wählt diese Phase und sie wählt zur Modularitätsdiagnose eine minimale doppelte Überlagerung der Sphäre. Vier Marken werden dadurch bedingt gezählt, ihre Lage und die physische Zeit werden nicht bestimmt. Drei gleich orientierte E8 Faktoren besitzen weiterhin eine lokale Gravitationsanomalie, auch wenn ein skalarer T Multiplikator trivial wird.

Quelle: [F1-E8-ANOMALY-INFLOW-ORIGIN-20261116.md](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/research/originlab/F1-E8-ANOMALY-INFLOW-ORIGIN-20261116.md).

### 2. Komposition und lokale Tomographie genügen nicht

Im expliziten Gegenmodell wirken eine kommutative Algebra C^4 und die linke reguläre Darstellung von M2(C) auf demselben Hilbert Schmidt Raum. Beide besitzen die genannten zyklischen und trennenden Eigenschaften und erfüllen die allgemeinen Kompositionsbedingungen. Der Unterschied ihrer Zentren und ihrer Kommutativität bleibt bestehen. Das Ergebnis schließt diese allgemeinen Anforderungen als hinreichenden Algebra Selektor aus. Zusätzliche Postulate wie Purifikation oder eine besondere Observable Dynamik Korrespondenz liegen außerhalb dieses Ausschlusses.

Quelle: [PROCESS_RECONSTRUCTION.txt](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/experiments/theory-contracts/seam-joint-source-origin-audit-20261016/PROCESS_RECONSTRUCTION.txt), insbesondere Abschnitt 3. Katalog: Atlasbefunde Zeile 31.

### 3. Ein Lorentzkegel ist bereits mathematisch vorhanden

Die bestehende Prüfung `compiler-cone-object-audit` beschreibt eine exakte algebraische Verbindung des positiven Kegels eines integralen Ankerkommutanten zu dem verwendeten E8 Operationsgitter. Sie lässt Zustand, erlaubte Prozesse, lokale Felder und gemeinsame Dynamik ausdrücklich offen. Die spätere Prozessprüfung enthält eine komplexe Kompositionsvorschrift sowie Choi, Pauli und Lorentz Identitäten, ohne diese physische Lücke zu schließen.

Quellen: [compiler-cone-object-audit/contract_index.json](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/experiments/theory-contracts/compiler-cone-object-audit/contract_index.json), [compiler-process-reconstruction-20260918/contract_index.json](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/experiments/theory-contracts/compiler-process-reconstruction-20260918/contract_index.json). Katalog: Versuchsprojekte Zeilen 107 und 127.

### 4. Der GHZ Ausschluss hat einen bestimmten Reduktionsbegriff

Die reine GHZ Familie ist unter gewöhnlicher partieller Spur keine projektiv konsistente Familie: Beim Entfernen eines Faktors verschwinden die globalen Kohärenzen. Die Originalnotiz behandelt zusätzlich drei klar genannte Anforderungen an die vergrößerte Quelle. Das Ergebnis darf nicht auf andere grobe Abbildungen, bedingte Quantencodes oder nichtfaktorisierende Algebren übertragen werden. Genau diese Alternativen werden in der Quelle als außerhalb des Ausschlusses benannt.

Quelle: [F1-P1-PROJECTIVE-CONSISTENCY-AND-COHERENT-GLUING-20261116.md](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/research/originlab/F1-P1-PROJECTIVE-CONSISTENCY-AND-COHERENT-GLUING-20261116.md).

## Technische Warnung zur neuen Kanalidee

Der Bildraum eines beliebigen unitalen vollständig positiven idempotenten Operators ist nicht automatisch eine Unteralgebra für die ursprüngliche Matrixmultiplikation. Die Choi Effros Konstruktion kann ein anderes Produkt benötigen. Für eine direkte Klassifikation als Summe von Matrixalgebren muss man daher einen bistochastischen Kanal, einen geeigneten treuen invarianten Zustand mit entsprechendem Satz oder ausdrücklich eine bedingte Erwartung voraussetzen. Die Auswahl der Teilräume darf nicht erst in frei gesetzten Krausoperatoren versteckt werden.

Aus M2 ⊕ M3 folgt eine zentrale spurfreie Richtung mit Größenverhältnis 3:(−2). Eine physikalische Hyperladungsnormierung folgt daraus allein noch nicht. Die nichtisomorphen Größen unterscheiden die beiden Summanden kanonisch; bei M2 ⊕ M2 fehlt diese Unterscheidung.

## Fazit für die neue Sitzung

Am deutlichsten unterscheidbar vom geprüften Bestand sind die **Minimalität des perfekten Fehlercodes** und die **zusätzliche äquivariante Observable Dynamik Korrespondenz**. Die **Algebra mit zwei ungleichen nichtkommutativen Sektoren** ist ein klarer bedingter Test mit einem besonders sichtbaren Risiko, die gewünschte Zerlegung schon in der Voraussetzung festzulegen. Der **invertierbare E8 Phasenansatz ist eine direkte Wiederholung** und sollte als Kontrolle des Neuheitsfilters erscheinen.

Weder die Repository Recherche noch der Katalogvergleich belegt einen Abstand zu sämtlichen früheren Versuchen. Eine gerechtfertigte Aussage kann lauten: Hier wurden bestimmte neu formulierte Voraussetzungen exakt ausgewertet und von vorhandenen Konstruktionen getrennt. Ein Beweis, dass diese Vorschläge physikalisch besser sind als alle bisherigen, liegt damit noch nicht vor.
