# Schwache Eichwirkung, normaler Higgsmodus und gewichteter Übergang

`TFPTNORMALWEAKTIME20261110` · 10. November 2026

Die direkte Identifikation der zwei holomorphen Moden der Determinantenlinie mit dem physikalischen Higgsdublett scheitert an einem einzigen Eichschritt. Der bekannte Ausweg über die Familie schwacher Richtungen rekonstruiert das bereits gegebene Dublett. Auf dem bereits hergeleiteten minimalen normalen Bündel lässt sich ein geladener Übergang samt relativem Gewicht explizit auswerten. Die dafür benutzte Produktzeit ist eine Voraussetzung der Rechnung.

## Der erste fehlende Schritt im Paper

Im [Originalabschnitt zur Bosonauswahl](../../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/001-subsection-from-bridge-to-geometric-completion.tex#L247) folgen auf `L₂ = det E₂ = O(1)` die zwei holomorphen Moden und ihre Deutung als ein schwaches Dublett. Die Dimension und die Indexfrage wurden bereits im [Clutch und Pushforward Audit](../original-admissible-source-closure-20261025/NORMAL-CLUTCH-TO-SOURCE-PUSHFORWARD.txt) getrennt.

Der direkte Vergleich der Wirkungen entscheidet die verbleibende Identifikation. Die globale vertikale Bündeleichung `z = −I₂` wirkt auf `det E₂` mit `det(z)=+1`. Auf jeder Sektion dieser Linie wirkt sie daher mit `+1`. Auf dem schwachen Fundamentalfeld und seiner konjugierten Yukawaspalte wirkt sie mit `−1`. Für jede eichäquivariante lineare Abbildung T der ersten Darstellung in die zweite folgt `T = −T = 0`. Selbst eine nichtlineare äquivariante Abbildung in das Dublett müsste jedes Bild im Fixraum von `−I₂` haben, also im Nullraum. Die zwei holomorphen Moden können unter einer separaten Rotation der Normalkugel ein geometrisches Dublett bilden. Das vertikale Eichzentrum lässt die Normalkugel punktweise fest. Der Gegenfall betrifft daher genau das natürliche Determinantenbündel mit seiner gegebenen Eichwirkung.

Der unveränderte 96er Operator v252 bestätigt den Unterschied. Mit `U_z = π(z)Jπ(z)J⁻¹` ist der Majoranablock gerade und der ganze Yukawablock ungerade. Die schwache positive Wurzel hat Rang 24 und erfüllt `[Q_w,T₊]=T₊`. Das gilt für die gesamte v252 Trägerdarstellung, bevor ein physischer Zustand ausgewählt ist.

Die [bereits vorhandene Familie](../tfpt-unified-source-derivation-20261029/aps_source_review/CLAUDE17-FAMILY-INDEX-REVIEW.txt#L23) löst ein anderes Problem. Ihre Indexlinie liegt über dem Raum W der Richtungen eines schon gegebenen schwachen Moduls V. Die Dualisierung `H⁰(W,L*)* ≅ V` erhält die volle Eichwirkung einschließlich Hyperladung. Damit kann die Familie V transportieren. Ihre Definition setzt V ein. Auch die Bewegung einer gemeinsamen Eichrichtung des Bündels und aller Felder erzeugt für sich keine Zeitdynamik. Auf einem Zeitintervall hat jede solche Kurve einen Lift U(t). Für eine statische Ausgangskonfiguration liefert die gemeinsame Eichtransformation den zeitlichen Anschluss `A_t = −(∂_t U)U⁻¹`. Dann gilt `F_tN=0` und `D_t H=0`. Derselbe geometrische Orbit lässt sich so mit verschwindender zeitlicher Energie durchlaufen. Ein physischer Übergang braucht eine Änderung eichinvarianter relativer Daten oder gebundene Endpunktdaten samt wirklichem Anschluss.

## Geladene Normalmode und klassisches Minimum

Die [normale Krümmungsrechnung](../chat-triangulation-20261108/TFPT-NORMALKRUEMMUNG-UND-SCHWACHER-ORBITTRANSPORT-20261109.txt#L58) liefert am vertikalen Yang–Mills Minimum die Spaltung `E₂=L⊕1` mit `deg L=1` und `F_L−F_1=2πi dA/A_N`. Sie lässt die Lage des Higgs relativ zur Krümmung zunächst frei. Der obere schwache Wurzelanteil ist ein Abschnitt von `Hom(1,L)=L`, der untere liegt auf der dualen Linie. Für einen tangentialen transversal polarisierten Eichanteil auf einem lokalen Produkt ist der vertikale kinetische Operator die kovariante skalare Laplace auf L. In der positiven Flusskonvention ergibt die Kähleridentität
```text
∇*∇ = 2 ∂̄*∂̄ + 2π/A_N .
```
Die untere Schranke wird von den zwei holomorphen Sektionen erreicht. Die neutrale Richtung erlaubt im gleichen Produkt eine konstante Normalmode mit Eigenwert null. Diese Aussagen beziehen sich auf das genannte tangentiale Quadrat und seine feste Hintergrundverbindung. Für die Normalanteile einer Eichform oder eine dynamisch integrierte Verbindung wäre ein eigener Operator nötig. Die runde Spektralform stimmt mit dem Ladungseinsfall der [Monopolharmonischen von Wu und Yang](https://www.sciencedirect.com/science/article/pii/0550321376901437) überein.

Die [lokalen Terme im Original](../../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/002-subsection-closed-texorpdfstring-41-41-chain-after-famil.tex#L140) und die volle v252 Wärmeauswertung R031 motivieren die nächste präzise bedingte Rechnung. Sei die normale Metrik rund, die Determinantenkrümmung gleichmäßig oder mit positiver zentraler Yang–Mills Energie frei minimierbar. Sei die lokale Energie die Summe der positiven Eichenergie, `Z|∇H|²` und `q|H|²+λ|H|⁴` mit konstantem q und `Z,λ>0`. Die separate untere Schranke aller Summanden ist erreichbar. Man setzt die bekannte minimale Eichverbindung ein, nimmt H parallel auf der flachen Grad null Linie und `|H|²=max(0,−q/(2λ))`. Jeder Minimierer muss alle Schranken sättigen. Bei nichtverschwindendem H liegt es deshalb auf der flachen Linie, falls diese im Konfigurationsraum zugelassen ist. Am R031 Zweig mit wörtlichem c=−6, d=1 und Raumkrümmung null ist der dortige quadratische Koeffizient durch die positive Majoranaantwort ebenfalls nichtnegativ. Die Rechnung setzt den Vorzeichenwechsel des physischen Potentials nicht voraus.

Falls der Selector nur die Grad eins Sektionen zulässt, ändert sich dieselbe Minimierung. Für eine normierte holomorphe Sektion s gilt in normiertem Flächenmaß `⟨|s|²⟩=1`, `⟨|s|⁴⟩=4/3` und `⟨|∇s|²⟩=2π/A_N`. Die Aussage gilt für jede normierte Sektion dieser zweidimensionalen geometrischen Mode durch Rotationssymmetrie. Mit `u=|h|²` und `H=hs` wird die Energie `(q+2πZ/A_N)u+4λu²/3`. Das Minimum liegt bei `max(0,−3(q+2πZ/A_N)/(8λ))`. Bei gleichzeitig erlaubtem konstantem Feld und Grad eins Mode wird die integrierte Quartik `λ(u₀²+2u₀u₁+4u₁²/3)`. Bei fester Gesamtnorm kostet jede Beimischung der Grad eins Mode `2πZ u₁/A_N + λu₁²/3`. Die beiden Bereiche sind also physikalisch unterscheidbar. Der Projektor auf die Grad eins Linie ist der Spektralprojektor der vertikalen spurfreien Krümmung, weshalb die relative Orientierung zur Krümmung eichinvariant gelesen werden kann.

Im direkten lokalen Produkt mit flacher Higgsrichtung und den üblichen kanonischen elektroschwachen Termen folgt im klassischen Quadrat `m_W²−cos²(θ_W)m_Z²=2π/A_N`. Der normale Anteil addiert sich zum geladenen Higgsbeitrag, während er beim neutralen konstanten Modus null ist. Die [aktive Standardmodellprojektion im Paper](../../../documents/tfpt_2_standard_model/sections/006-section-neutrinos.tex#L332) hat diese Verschiebung auf der Baumebene nicht und führt anschließend separate Matchings ein. Die Differenz ist eine Prüfbedingung für die physische Reduktion. Ohne deren Ableitung liefert die normale Lücke weder einen W Pol noch eine physische Kalibrierung von A_N.

## Übergang mit allen Normalmoden

Auf der runden Sphäre mit Radius R setze `x=(1+cos θ)/2` und normiere die Fläche durch `dx dφ/(2π)`. Die beiden normierten niedrigsten Ladungsmoden sind lokal `s₊=√(2x)` und `s₋=√(2(1−x)) exp(±iφ)` in einer passenden Bündelkarte. Für den Grad eins Higgszweig `H=hs₊` liefert das lokale Eichquadrat die Multiplikation `b·2x` mit b proportional zur Higgsnorm und b≥0. Der benutzte tangentiale Anteil ist räumlich und transversal. Bei tangential räumlich homogenem H und verschwindendem räumlichen Tangentialimpuls koppelt er im betrachteten Quadrat weder an den zeitlichen noch an die normalen Eichanteile. Mit einer zusätzlichen lokalen Produktzeit liefert das gewöhnliche tangentiale kinetische Quadrat `−∂²_t`.

Die axialen Teilräume tragen die vollständigen, normierten Jacobi Basen `s₊P⁽⁰˒¹⁾ₙ(2x−1)` beziehungsweise `s₋P⁽¹˒⁰⁾ₙ(2x−1)`. Nach Abtrennen von s wirkt die Laplace mit Radius eins als `−x(1−x)∂²_x−(β+1−3x)∂_x+1/2`. Jacobi Differentialgleichung und Dreitermrekursion geben
```text
λₙ = (n(n+2)+1/2)/R²
(J±)ₙₙ = 1 ± dₙ,  dₙ = 1/((2n+1)(2n+3))
(J±)ₙ,ₙ₊₁ = aₙ = √((n+1)(n+2))/(2n+3)
K± = diag λ + bJ± .
```
Der erste Übergang nach außen hat auf beiden Seiten den Faktor `b√2/3`. Jeder weitere Nachbarübergang ist für b>0 ebenfalls ungleich null. Damit kann kein fester endlicher Normalmodenausschnitt, der die Grundmode enthält, diese Operatoren reduzieren. Schon die lokale Lie Klammer zweier konjugiert geladener Grundmoden erzeugt einen neutralen Anteil proportional zu `|s|²`. Seine nichtkonstante Komponente hat Normquadrat 1/3 bei s₊ oder s₋, der gekreuzte Anteil hat Normquadrat 2/3. Wiederholtes Klammern erzeugt `|s|^{2j}s`. Da `|s|²` ein Intervall von Werten annimmt, sind diese Funktionen linear unabhängig. Der Grundmodenraum allein ist folglich auch unter der gegebenen Eichklammer nicht abgeschlossen.

Die Rückkehr zur niedrigsten Mode lässt sich trotzdem vollständig rechnen. Für z>0 und den hinteren Anteil ab n gilt
```text
mₙ(z) = 1/(z+λₙ+b(1±dₙ)−b²aₙ² mₙ₊₁(z)).
1/(z+λ_N+2b) ≤ m_N(z) ≤ 1/(z+λ_N).
```
Dies ist die Schur Rekursion mit einer garantierten Startschranke durch `0≤J≤2`. Der konstante b Fall hat positiven tangentialen euklidischen Zeitkern
`C±(t)=⟨e₀, K±⁻¹ᐟ² exp(−|t|√K±) e₀⟩/2`.
Der positive Spektralanteil ist die konkrete Laplace Darstellung. Bei den Testeinheiten R=b=z=t=1 sind die beiden projizierten Resolventen ungefähr 0,3580060594 und 0,4704568663. Die Kovarianzen sind ungefähr 0,0983973510 und 0,1643444586. Die Verwendung dieser Testeinheiten fixiert keine physische Skala.

## Relatives Gewicht ganzer Normverläufe

Für periodische Testzeit der Länge T und beschränktes `b(t)≥0` gilt `ℒ±=−∂²_t+diag λ+b(t)J±` und `Δ=ℒ₊−ℒ₋=2b(t)diag dₙ≥0`. Der relative Operator `ℒ₋⁻¹ᐟ²Δℒ₋⁻¹ᐟ²` ist spurklassig. Vergleich mit dem freien Operator begrenzt seine Spur durch `2||b||∞ Σ dₙ T coth(T√λₙ/2)/(2√λₙ)`. Das normierte Verhältnis der beiden Gaußgewichte eines komplexen Feldes ist damit
```text
Z₊[b]/Z₋[b] = det_F(I+ℒ₋⁻¹ᐟ²Δℒ₋⁻¹ᐟ²)⁻¹ .
```
Es ist für nichtverschwindendes b strikt kleiner als eins. Die inverse Potenz zählt zwei reelle Komponenten. Die Spurform der Logarithmusdifferenz ist `∫₀¹ Tr(Δ(ℒ₋+sΔ)⁻¹)ds`. Beim Test R=b=T=1 ergibt der endliche Teil der Logarithmusdifferenz bei Tiefe 96 ungefähr 0,58385734. Die analytische Obergrenze multipliziert ihn höchstens mit `1/(1−2ε)`, `ε=1/(4·96²)`, und addiert höchstens `1/(8·95²)+1/(6·95³)`. Hierfür ist der hintere Schur Beitrag höchstens ε und `ℒ₋+sΔ≥1/2`; die Spur der übrigen Diagonale wird mit dem freien Zeitpropagator beschränkt. Zusammen entspricht das einem Verhältnis von ungefähr 0,5577. Die ausgewiesenen endlichen Eigenwerte wurden numerisch berechnet. Die Tail Schranken sind analytisch.

Auf der unendlichen tangentialen Testzeit mit kompakt getragenem b ist die zweite funktionale Ableitung von `log(Z₊/Z₋)` bei b=0 für zwei Zeitschnitte exakt
```text
Σₙ dₙ/λₙ exp(−2√λₙ|t−s|).
```
Die gleichen Nebendiagonalen heben sich in der Differenz weg, und die Differenz der quadrierten Diagonalen ist `4dₙ`. Der freie Zeitpropagator trägt `1/(2√λₙ)`. Am Kontakt und bei R=1 ergibt partielle Bruchzerlegung `−1−(π/√2)cot(π/√2)=0,6910120633…`. Eine Variation in der Norm ist dabei von einer linearen Higgsvariation um H=0 zu trennen. Da b quadratisch in H ist, erzeugt dieser relative Zweig um H=0 keinen gemischten Term linear in jeder von zwei Higgsvariationen. In der Normantwort ist das Spektralmaß diskret. Dies liefert im zuvor ausgewerteten APS Vergleich keinen kontinuierlichen Ausgleich auf einem Frequenzband.

## Prüfbarer Anschluss

[Gezielter reproduzierbarer Check](check_normal_weak_higgs_transfer_20261110.py) · [R031 und Originale im Ergebnisledger](../../../docs/research/fundamental-plan/RESULTS.md#r031-vorhandener-ganzer-higgsoperator-verbindet-bewegung-quartik-und-neutrinomasse)

Der erste physische Anschluss verlangt eine wirkliche äquivariante Realisierung des Higgsfelds gemeinsam mit der schwachen Verbindung. Die natürliche Determinantenlinie allein ist durch das Eichzentrum als Träger ausgeschlossen. Die bekannte Familie über schwachen Richtungen bietet die passende Wirkung aus einem schon gegebenen V. Ob der Originalselector diese Familie und die zugehörige relative Higgsrichtung auswählt, muss am tatsächlichen ursprünglichen Operator entschieden werden. Danach muss derselbe Anschluss die Endpunktdaten und das absolute gemeinsame Gewicht einschließlich Leerreferenz liefern. In der obigen Rechnung wurden die feste minimale Normalverbindung, die runde Normalmetrik und die tangentiale Produktzeit ausdrücklich benutzt. Das relative Gaußgewicht betrifft zwei axiale Teilräume bei vorgegebener Higgsnorm und liefert keine vollständige gemeinsame Feldintegration. `Z₀`, `F(0)`, `C_UY` und `H_tr,pos` bleiben physisch unberechnet.

## Nullmodenprüfung ohne Wahl der minimalen Verbindung

Die Lücke für den geladenen Zweig ist auch außerhalb des runden Minimalhintergrunds strukturell sichtbar. Für jede glatte unitäre Verbindung auf E₂ über S² mit ungeradem Determinantengrad gilt
```text
dim_R ker(∇*∇ auf su(E₂)) ≤ 1 .
```
Das setzt für den tangentialen Vektoranteil dieselbe lokale Produktreduktion voraus. Die Aussage über den vertikalen Operator selbst braucht weder eine runde Metrik noch das Yang–Mills Minimum. Der Kern besteht wegen `⟨a,∇*∇a⟩=||∇a||²` aus parallelen Abschnitten. Gäbe es zwei reell unabhängige parallele traceless antihermitesche Abschnitte, wären ihre Werte in jeder Faser unabhängig. Mit ihrer Klammer ergäben sie ein paralleles Dreibein. Daher würde F mit ganz su(2) kommutieren und wäre skalar. Schon ein nichtverschwindender paralleler Abschnitt liefert zwei parallele komplexe Eigenlinien. Bei skalarer Krümmung haben beide Linien dieselbe Krümmung und deshalb denselben ganzzahligen Grad. Ihre Gradsumme wäre gerade. Das widerspricht dem ursprünglichen Grad eins. Auf der kompakten Sphäre ist das Spektrum diskret. Somit fehlen mindestens zwei der drei schwachen vertikalen Nullmoden auf jedem solchen Hintergrund. Das bekannte Minimum realisiert genau eine neutrale Nullmode.

Ein zusätzlicher nichtnegativer lokaler Higgs Hesseterm kann die fehlenden Nullmoden in diesem Quadrat nicht herstellen. Der Kern einer Summe nichtnegativer quadratischer Formen liegt im Kern jedes Summanden. Damit scheidet auch die Wahl einer anderen glatten Verbindung oder Metrik auf demselben normalen Bündel als Reparatur der direkten lokalen Produktprojektion aus. Die physische Reduktion müsste erklären, weshalb ihre schwachen Felder einen anderen Operator tragen oder wie der ursprüngliche gemeinsame Beitrag die reine vertikale Laplace Lesart verändert. Die Aussage ist keine absolute Lückenschranke über alle Verbindungen und keine Behauptung über die vollständige eichfixierte Quantenquelle.
