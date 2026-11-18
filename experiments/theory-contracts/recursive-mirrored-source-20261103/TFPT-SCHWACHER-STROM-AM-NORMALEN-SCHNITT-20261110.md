# Der schwache Strom am normalen Schnitt

`TFPTWEAKEQUATOR20261110` · 10. November 2026

Die fehlenden schwachen Richtungen lassen sich im ganzen E₈ bei festgehaltener Hyperladung nicht durch andere Wurzeln ersetzen. Auf dem Äquator gibt es dagegen einen vollständigen schwachen Strom. Die Normalkappen und der gegebene Viertelkreis unterscheiden, welche Fortsetzung und welcher Zustand diesen Strom tragen. Das ergibt eine unnormierte, direkt ausrechenbare Rückkehrprobe für den ursprünglichen Selektor.

## Der ganze E₈ Träger und die wirklichen Ladungen

Der [ursprüngliche Trägerabschnitt](../../../documents/archive--tfpt-45--source_extracts--04_qft_source/sections/003-source.tex#L44) gibt E=E₃⊕E₂ und Y=−P₃/3+P₂/2. Der [normale Bosonabschnitt](../../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/001-subsection-from-bridge-to-geometric-completion.tex#L247) gibt c₁(det E₂)=1. Wir verwenden diese Deklarationen gemeinsam mit dem vorhandenen [vollen E₈ Bündellift](../original-admissible-source-closure-20261025/JOINT-LIFT-NATIVE-CLOCK-CONNECTION.txt#L53). Der Vergleich verwendet tatsächliche Ladungen. Die Wahl der physisch zulässigen Zustände wird dabei nicht vorweggenommen. Zur Unterscheidung von lokalen Eichalgebren und globalen Flussdaten siehe [Tong, Line Operators in the Standard Model](https://arxiv.org/abs/1705.01853).

Im ursprünglichen achtkoordinatigen Wurzelbild sind
```text
6Y  = (−2,−2,−2,3,3,0,0,0),
2T₃ = ( 0, 0, 0,1,−1,0,0,0).
```
Unter den 112 ganzzahligen Wurzeln mit Y=0 liegen die sechs Farbwurzeln eᵢ−eⱼ, die zwölf Familienwurzeln ±eᵢ±eⱼ mit i,j≥6 und ±(e₄−e₅). Bei den Halbwertwurzeln verlangt Y=0 die Gleichung 2(s₁+s₂+s₃)=3(s₄+s₅). Die Farbsumme ist ungerade. Daher sind die ersten fünf Vorzeichen alle gleich. Die acht zulässigen Halbwertwurzeln haben T₃ Gewicht null. Insgesamt sind 28 Wurzeln neutral unter Y. Die einzige neutrale Wurzel mit T₃ Gewicht +1 ist e₄−e₅; beim negativen Gewicht ist es ihr Negatives.

Die beiden Wurzeln und ihr Cartan bilden ein direktes schwaches A₁ Ideal in der Lie Algebra des Y Zentralisators. Alle anderen neutralen Wurzeln sind dazu orthogonal. Weder ihre Summe noch Differenz mit der schwachen Wurzel hat Normquadrat zwei. Sie kommutieren daher mit diesem Ideal. Bei jedem gegebenen toralen Lift mit schwachen Graden a und 1−a trägt die eindeutige schwache Wurzel den Normalgrad 2a−1. Der Zusatz der Familienkoordinaten ändert diese Paarung nicht.

Das Ergebnis gilt auch für nichttorale lokale Trivialisierungen desselben fest markierten Bündels. Sein schwaches reelles adjungiertes Bündel besitzt w₂=c₁(E₂) mod 2. Das folgt im toralen Bild direkt aus der ungeraden Windung des reellen Rotationspaars, und beide Seiten sind Bündelinvarianten. Wäre ein glattes punktweises su(2) im ganzen Y Zentralisator mit der ursprünglichen schwachen Wirkung vorhanden, hätte seine Projektion auf das schwache Ideal nichttriviale Wirkung. In jeder Faser ist eine solche Projektion eine Abbildung einfacher reeller dreidimensionaler Lie Algebren, also entweder null oder ein Isomorphismus. Bei glatter Abhängigkeit auf zusammenhängender Basis ist die Wahl konstant, weil der Normbetrag eines Lie Algebra Isomorphismus in den Killing Normen fest ist. Ein nichtverschwindender Zweig gäbe eine globale schwache Dreierbasis und würde w₂=0 erzwingen. Eine Beigabe aus Familie oder anderen Y neutralen E₈ Richtungen umgeht das nicht. Eine Kontraktion der ganzen E₈ Verklebung, die Y samt seinem schwachen Ideal festhält, scheitert aus demselben Grund.

## Normale Fermionmoden und der wirkliche Schnitt der zwei Kappen

Auf dem gegebenen runden Normalprodukt mit Monopolverbindung sind die beiden Komponenten eines links schwachen Dubletts in benachbarten Liniengraden d+1 und d. Der bereits geprüfte [global verklebte v252 Träger](../../../semantic_explorer/research/program.json#L38056) verwendet bei trivialer Farbe den linken Twist m_L=−1, also die Grade 0 und −1. Der allgemeine Fall bleibt hier sichtbar. Dieser v252 Bündellift und der zuvor geprüfte E₈ Lift sind zwei bereits benannte Realisierungen. Die determinantabhängige v252 Korrektur ist auf den alten E₈ Gewichten keine einzelne torale Umwertung. Für die Rechnungen wird keine neue Identität zwischen den beiden Liften angenommen.

Die beiden Spinorkomponenten des mit O(d) gekoppelten Operators haben auf der Sphäre die homogenen Gewichte d−1 und d+1. Der gekoppelte Diracoperator ist auf SU(2) durch die passende Leiterwirkung zwischen beiden Gewichten gegeben. Auf einem Spin j Block liefert das Produkt der Leiterfaktoren
`(j+1/2)²−d²/4`. Mit `j=(|d|−1)/2+n` ergibt das für positive Eigenwerte des Quadrats
```text
R² D_d² = n(n+|d|), n = 1,2,…
```
Bei |d|>0 hat der unterste einseitige Block n=0 Dimension |d| und normale Chiralität sign(d). Bei d=0 gibt es keinen solchen Nullblock. Das entspricht der homogenen Spinrealisierung bei [Jante und Schroers, §2](https://arxiv.org/html/1312.4879#S2). Für aufeinanderfolgende Absolutgrade a,a+1 gilt streng
```text
n(n+a) < n(n+a+1) < (n+1)(n+a+1), n≥1 .
```
Die positiven Spektren sind disjunkt. Also kann auch ein nichtlokaler beschränkter Operator, der beide vollen freien Quadrate miteinander verflicht, außerhalb der Nullräume keinen geladenen Wechsel tragen. Für d=−1, d+1=0 fehlt zudem auf einer Seite der Nullraum. Es gibt dann überhaupt keinen solchen Wechsel. Der Test setzt die freie normale Kinetik an beiden Enden und einen gemeinsamen Radius voraus.

Die Identifikation mit geometrischen Drehungen hat einen weiteren, einfachen Gegenfall. Das skalare H⁰(O(1)) trägt in der runden Geometrie Spin 1/2. Die normalen Fermionnullmoden des tatsächlichen linken Splitzweigs liegen dagegen nur auf O(−1). Sie sind unter der geometrischen Drehung Singuletts. Allgemein ist der normale Spinindex für den linken Leptonensektor aus drei Dubletts ungerade. Denn `c₁(E₂⊗F₃)=3c₁(E₂)+2c₁(F₃)`, und der gewöhnliche normale Spinindex ist dieser Grad. Die Hyperladung −1/2 isoliert diesen Sektor im unveränderten v252. Falls beide ganzen Chiralitätskerne nur physische Dubletts tragen sollen, müssten ihre Dimensionen gerade sein. Auch ohne Spinquantenzahlen kann daher die Auswahl des ganzen normalen Nullraums kein vollständiges schwaches Leptonendublett liefern. Der Beweis betrifft die physische Lesart dieses deklarierten Normaloperators, ohne daraus eine Aussage über andere Arten ursprünglicher Reduktion abzuleiten.

Der Äquator der normalen Sphäre hat nun zwei echte Cauchyräume für die Gleichung D⁺=0. Im gemeinsamen Nordrahmen der Linie O(d−1) sind sie
```text
N = span zⁿ, n≥0,
S_d = span zⁿ, n≤d−1 .
```
Die Nordform folgt aus holomorpher Fortsetzung in die Nordscheibe. Eine in der Südkoordinate w=1/z holomorphe Funktion wird in den Nordrahmen durch z^{d−1}f(1/z) geschrieben. Der normierte Äquatorfaktor der runden Monopollinie ist konstant und ändert die geschlossenen Fourierunterräume nicht. Die Fredholmpaarung hat Schnittdimension max(d,0) und Kodimension der Summe max(−d,0). Ihr Index ist d. Dies ist der konkrete normale Schnittfall der [Calderón und Bojarski Verklebung](https://arxiv.org/abs/0803.4160).

Der geladene Randoperator z bildet S_d unitär auf S_{d+1} ab. Auf N bildet er dagegen nur auf zN ab. Mit beiden Randoperatoren auf N komprimiert ist der geladene Block der einseitige Shift S. Sein Kommutator mit der adjungierten schwachen Richtung ist
```text
[W₊,W₋]|_(N⊕N) = diag(I−|0⟩⟨0|,−I).
```
Die fehlende Mode ist exakt rang eins. Auch ein anderer unitärer Operator, der beide Kappen gleichzeitig erhält und die Grade vertauscht, kann den Unterschied der Fredholmindizes d und d+1 nicht beseitigen. Dasselbe ist mit vertauschter Chiralität im rechten Randraum sichtbar. So lokalisiert sich das Problem der durchgehenden schwachen Wirkung auf eine präzise Randlinie, bevor ein vollständiges Maß festgelegt wird.

Auf derselben Nordhälfte ist p₀=I−SS* zugleich der Defektprojektor für die gerichtete relative Fredholm Determinante des Shiftpaars. Für 0≤q<1 ist `(I−q S*S)(I−q SS*)⁻¹=I−q p₀`, also ist die Determinante genau `1−q`. Es ist im gegebenen Normalrahmen derselbe Projektor wie in der schwachen Kommutatorabweichung. Die [bereits bewiesene Rang eins Formel für die ursprüngliche Naht](../original-admissible-source-closure-20261025/PAARINDEX-EINFACH.txt#L18) hat genau diese Operatorform. Eine physische Zuordnung der beiden Randlinien oder von q folgt daraus noch nicht.

## Ungekürzter schwacher Strom und zwei überprüfbare Zustände

Auf dem ganzen Randkreis vor der Kappenkompression ist die Lösung der Lie Algebra einfach. Wenn U(z) die gegebene Grad eins Verklebung ist, dann sind `W₊=zT₊`, `W₋=z⁻¹T₋` und `T₃` genau die durch U transportierten konstanten Generatoren. Ihr Einteilchenkommutator ist exakt `[W₊,W₋]=2T₃`. Das verwendet dieselbe Ladung und den vorhandenen v252 Operator, keinen neuen E₈ Generator.

Für die folgenden Rückkehrwerte wird der Äquatoroperator zusätzlich mit der bereits vorliegenden freien Viertelkreiskinetik `h=γ_F(n+α)/R` getestet, α=1/4 oder 3/4. Dies ist ein konkret benannter Anschluss als Test. Der geometrische Äquator und der tangentiale ursprüngliche Viertelkreis sind damit nicht als physisch identisch gesetzt. Bei H=0 wirkt der Majoranablock nur auf schwachen Singuletts und verändert diesen freien schwachen Unterblock nicht.

Im v252 Matrixaufbau hat T₊ exakt 24 Einträge der Betragsnorm eins, mit je zwölf in der Teilchen und in der konjugierten Hälfte. Die Ränge nach (γ_F,6Y) sind (1,−3):3, (1,1):9, (−1,3):3 und (−1,−1):9. Bei festem freiem See P=1_(h<0) erzeugt W₊ auf der positiven γ_F Hälfte in jeder der zwölf Richtungen genau einen orthogonalen Teilchen Loch Zustand: n=−1 in der unteren Komponente wird nach n=0 in der oberen verschoben. Die umgekehrte Ordnung entsteht auf der anderen γ_F Hälfte. Der Energieabstand ist jedes Mal 1/R, unabhängig vom Viertel. Auf dem signierten vollständigen komplexen CAR Test liefert das
```text
⟨Ω|J(W₋) e^(−t(H−E₀)) J(W₊)|Ω⟩ = 12 e^(−t/R),
⟨Ω|J(W₊) e^(−t(H−E₀)) J(W₋)|Ω⟩ = 12 e^(−t/R).
```
Diese Formeln betreffen den unveränderten See im genannten vollen Testträger. Eine zusätzliche reale physische Reduktion, Eichpaarung oder Auswahl ist nicht eingerechnet. Pro Dublett hat der normalgeordnete Kommutator beim See n<0 den zentralen Zusatz −1, weil die Besetzung nur beim Weg von n=−1 nach n=0 gekreuzt wird. Bei umgekehrter Besetzung ist der Zusatz +1. Die Beiträge sind damit −12 und +12 in den beiden γ_F Hälften. Ihre Summe verschwindet, während beide positiven Rückkehrnormen 12 bleiben.

Transportiert man dagegen denselben See und seine freie Kinetik mit U, dann sind `P'=UPU*` und `h'=UhU*`. In beiden schwachen Richtungen werden nun gleiche verschobene Fourierenergien verbunden und die Ströme kommutieren mit P' und h'. Die normalgeordneten Ströme erzeugen dort keine Teilchen Loch Zustände, die beiden Rückkehrwerte sind null. Das ist ein exakter Vergleich zweier benannter Referenzen auf demselben Randstrom. Für seine Übertragung muss die Quellenabbildung entscheiden, ob und wie der ursprüngliche See, die beiden Kappen und die Viertelzeit zusammengehören.

Bei endlicher inverser Temperatur β lässt sich derselbe Vergleich ohne Abschneiden des Fouriersees rechnen. Schreibe ω=1/R und f(E)=(1+exp(βE))⁻¹. Die Wick Antwort eines Übergangs mit Energieänderung δ verwendet f(E)(1−f(E+δ)). Für δ≠0 gilt exakt
```text
f(E)(1−f(E+δ)) = (f(E)−f(E+δ))/(1−exp(−βδ)).
```
Die Differenz teleskopiert über den ganzen Viertelkreis. Die beiden signierten Hälften liefern zusammen für jede Ordnung und 0≤t≤β
```text
K_β(t)=12[exp(−ω t)+exp(−ω(β−t))]/[1−exp(−βω)].
```
Das gilt für beide ursprünglichen α. Bei gemeinsam transportierter Kinetik und thermischem Zustand ist die Energieänderung null und die Antwort stattdessen zeitkonstant, nämlich `24 Σₙ f((n+α)/R)(1−f((n+α)/R))`. Der Grenzwert bei β→∞ ist null. Der Vergleich unterscheidet beide Referenzen auch bei endlicher Temperatur über die wirkliche Zeitabhängigkeit. Er benutzt weiterhin den freien signierten v252 Anschluss.

## Rücktest am ursprünglichen Selektor

Der [aktive Paperabschnitt zur admissiblen Zeit](../../../documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex#L48) verweist auf den endlichen positiven Transfer. Die [unveränderte v330 Verifikation](../../../verification/v330_qgamb_admissible_measure.py#L42) hat genau zwei verschiedene positive angeregte Energien. Soll ein Transfer die zwei gleich schnell abklingenden geladenen Randantworten oben samt ihrer Ladung tragen, müssen die durch W₊ und W₋ erzeugten Vektoren am selben positiven Eigenwert liegen. Die Eindeutigkeit der positiven Laplacetransformation liefert dies aus der einzelnen Exponentialrate. Die entgegengesetzten T₃ Ladungen machen die beiden Vektoren orthogonal. Ein solcher Eigenwert braucht also Multiplizität mindestens zwei. Der tatsächliche skalare v330 Transfer hat keinen solchen Eigenwert. Keine Wahl einer gemeinsamen Zeiteinheit ändert die Multiplizität. Dies schließt seine direkte Lesart als Träger dieses geladenen Tests aus und lässt den v330 Beweis seines endlichen Maßes bestehen.

Der [ursprüngliche vierteilige Selektor](../../../documents/archive--tfpt-45--source_extracts--04_qft_source/sections/002-section-closure-theorems-from-the-admissibility-complex.tex#L78) fixiert Nahtvorzeichen, Farbe, Determinanteninvolution und die wirkliche Spektrallücke von D_adm. Die Rang eins Normalverschiebung gibt weder die Wirkung dieser vier Faktoren auf die Zwischenlinie noch die ursprüngliche Zeitkarte vor. Der [installierte Quellenaufruf](../original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/original_source_adapter.py#L309) verlangt für die originale Faser eine tatsächliche Selektorform, Randpaarung, Zeit, Referenz und Maßbindung. Die erneute direkte Prüfung des ursprünglichen QΘ Wortaufrufs meldet in der Faser `original_blocked_fiber_packet` und im ganzen Endpunkt `evaluate_original_full_blocked_endpoint_source`. Die vier registrierten Z Werte und die Maßmasse bleiben ohne Wert. Das Randergebnis bestimmt dort weder die ganze Auswahl noch die Integrationsmasse. Die kleinste bindende Fortsetzung ist konkret: denselben markierten schwachen Hin und Rückweg mit den beiden tatsächlichen Randpaarungen lesen. Auf dem alten See die Rang eins Linie pro Dublett und die beiden positiven Normen wiederfinden oder die genau zugehörige Änderung von Referenz und Zeit nachweisen. Danach den originalen unnormierten gemeinsamen Beitrag mit der vorhandenen leeren Referenz auswerten. Ein normiertes Einzelbeispiel kann den fehlenden Quellenaufruf nicht ersetzen.

[Reproduzierbarer Prüfcode](check_native_weak_equator_20261110.py). `v1` 17/17, `v252` 11/11 und `v330` 6/6 lokal erneut bestanden. Die Algebra, die Kappenrechnung und der gegebenenfalls mitgeführte See wurden getrennt überprüft. Das physische (Z_0), (F(0)), (C_{UY}) und (H_{mathrm{tr,pos}}) bleiben unberechnet.
