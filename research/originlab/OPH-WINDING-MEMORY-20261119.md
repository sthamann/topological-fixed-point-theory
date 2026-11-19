# Vollständig positive Viermarkenmischung · 19. November 2026

OPHCP20261119. Die Suche vergleicht eine autonome vierdimensionale Gedächtnisalgebra mit einer Fortsetzung, die nur die diskrete Deckwirkung respektiert. Die physische P1 und R113 Quellenidentifikation bleibt offen.

## Klasse und Schranke

Die Stationen 0 bis 3 haben den Windungsoperator N = diag(0,1,2,3)/4. Die kontinuierliche, vollständig positive, spurerhaltende Halbgruppe Φ respektiert die gesamte Konjugationsgruppe von N. Ihre Populationszeit hat Generator m/2 je gerichteter nächster Kante des Viererrings, mit Gesamtverlustrate m > 0. Die Annahme betrifft den autonomen Kanal auf M₄. Im gemeinsamen affinen Tensorraum können andere Dynamiken und Kompensationen vorliegen.

Der Übertragsoperator E₀₃ hat die eindeutige Frequenz −3/4. Damit ist Φₜ(E₀₃) = λₜ E₀₃ und |λₜ| = exp(−γt). Der 2×2 Hauptminor der Choi Matrix auf den Vektoren |00〉, |33〉 erzwingt |λₜ| ≤ p₀₀(t). Die Populationskurve p₀₀(t) = (1+2exp(−mt)+exp(−2mt))/4 hat Anfangssteigung −m. Daher gilt γ ≥ m. Populationsmischung mit unveränderter Übertragskohärenz ist in der gesamten Klasse ausgeschlossen. Der spezielle diagonalenreine Mischgenerator aus dem Vergleich hat bei Einheitsrate einen expliziten negativen bedingten Choi Zeugen mit Wert −1.

Mit normierter Spurmessung C_vac(t) = (2 exp(−mt)+exp(−2mt))/16 gilt damit auf demselben t und derselben Referenz
```
|tr_4(E_30 Φ_t(E_03))| ≤ (sqrt(1+16 C_vac(t))-1)/4 .
```
Das ist eine Schranke für die gemessene Übertragsantwort. Sie setzt eine tatsächliche Bindung der Standortprojektion und des Übertragsoperators an dieselbe Quelle voraus. Die R100 Energieratio e^(−3t/4) ist eine andere Antwort, solange die Operatoridentifikation fehlt.

## Scharfer positiver Kandidat und weitere gebundene Antwort

Setze J = E₁₀+E₂₁+E₃₂ und B = E₀₃. Die vier Lindblad Richtungen J, J†, B, B† mit gleichen Gewichten m/2 sind zeitkovariant und liefern dieselbe Populationsmischung. Der Generator ist bzgl. der normierten Spur selbstadjungiert. Sein Kern ist eindimensional, die bedingte Choi Matrix ist ein expliziter positiver Gram. Der Übertrag erreicht γ = m.

Auf den drei gleichfrequenten Kanten gilt der exakte Generatorblock −m I + (m/2) A(P₃). Für jeden Generator der erklärten Klasse sind die Beträge der benachbarten Kopplungen dieses Blocks höchstens m/2. Das folgt aus Cauchy Schwarz für die positive Lindblad Koeffizientenmatrix bei feststehenden diagonalen Sprungraten. Andere Beiträge auf der Diagonale haben nichtpositive reelle Teile über den unvermeidbaren Verlust −m hinaus. Die komponentenweise Betragsdynamik ist daher durch den angegebenen reellen Block dominiert. Für S = J+B ergibt die Summation der beiden Frequenzblöcke die scharfe Hüllkurve
```
|tr_4(S† Φ_t(S))| ≤ exp(−mt)(1+3 cosh(mt/√2)+2√2 sinh(mt/√2))/4 .
```
Der konstruierte Generator erreicht sie. Die notwendige Anfangsdämpfung ist mindestens m/2 für diese normierte Spurkorrelation. Der Operator S trägt zwei Frequenzen, weshalb seine gesamte Antwort nicht die reine Übertragsrate misst.

## Kontrolle und Originalstatus

Fordert man nur den diskreten Deckschritt, liefert die Mischung der zyklischen Konjugationen durch S und S† einen exakten Gegenfall. Sie hat dieselben Populationsraten und lässt S fest. Sie ist unter exp(i t N) nicht kovariant, und ihr stationärer Raum auf vollem M₄ hat Dimension vier. Diese Kontrolle begrenzt den Ausschluss auf die erklärte Modellklasse.

`python3 -B research/originlab/f1_oph_winding_memory.py` prüft die Matrizen, die verglichenen Kovarianzen und die Choi Zeugen. Der reguläre OriginLab Aufruf `--campaign f1_oph_coherent_histories` enthält das Ergebnis. Der kombinierte Lauf mit `--campaign r113 --repo <ursprünglicher Checkout>` meldet R113 als BLOCKED. Fehlend sind `original_blocked_fiber_packet`, die vollständige gewichtete physische Randform und die gebundene absolute Leerreferenz. Das endliche Modell erzeugt daher noch keinen originalen Wert oder eine vollständige fundamentale Lösung.

Quellen zur Ausgangsdynamik: [Kohärente OPH Randwege](OPH-COHERENT-HISTORY-20261119.md) und [OPH Nahtoperator](https://github.com/FloatingPragma/observer-patch-holography/blob/main/model/polynomial_sv/SPEC.md). Allgemeiner Hintergrund: [Chruściński, vollständige Positivität und Relaxationsraten](https://refubium.fu-berlin.de/handle/fub188/50247).
