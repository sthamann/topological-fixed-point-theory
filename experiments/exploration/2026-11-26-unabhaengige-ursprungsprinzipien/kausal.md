# Kausale Rekonstruktion vor der TFPT Naht

## Ergebnis

Der unabhängige Ansatz liefert die lokale kausale Bühne vor den TFPT Objekten: Vier Dimensionen, ein Lorentzkegel und die Maxwell Dynamik folgen gemeinsam aus einer prämetrischen Antwortstruktur. Der exakte Test besteht samt Gegenfall.

Er ist trotzdem noch keine TFPT Lösung. Die Antwortstruktur selbst ist gewählt. Sie erzeugt weder den Zustand von P1 noch den Fünferträger von P2. Der Befund ist deshalb ein tragfähiger Kandidat für die offene Raumzeitbrücke, aber kein Kandidat der TFPT Stufe 1 bis 4.

## Die neue Voraussetzung

Der Ausgangspunkt ist eine differenzierbare Mannigfaltigkeit unbekannter Dimension. Eine Metrik, ein Lichtkegel und eine Verbindung werden nicht vorausgesetzt.

Es gelten vier physische Annahmen:

1. Elektrische Ladung ist erhalten. Der Strom ist eine verdrehte Form $J$ vom Grad $n-1$, lokal also $dH=J$ mit einer Anregung $H$ vom Grad $n-2$.
2. Magnetischer Fluss ist erhalten. Die Feldstärke ist eine Zweiform $F$ mit $dF=0$.
3. Die Vakuumantwort ist lokal, linear, verlustfrei und reziprok. Sie vertauscht $F$ und $H$ punktweise, ohne vorher eine Metrik oder einen anderen Gradwechsler einzusetzen.
4. Scharfe Wellenfronten besitzen einen einzigen hyperbolischen Kegel. Das ist die Forderung nach fehlender Doppelbrechung im Vakuum.

Die dritte Annahme ist die erste zusätzlich gewählte Struktur. Gerade sie macht den Ansatz prüfbar. Weil $F$ Grad 2 und $H$ Grad $n-2$ hat, kann eine solche interne Vertauschung nur gelten, wenn

\[
2=n-2,\qquad n=4.
\]

Die Literatur zur prämetrischen Elektrodynamik setzt vier Dimensionen gewöhnlich bereits voraus. Die Gleichung oben ist deshalb keine Behauptung aus diesen Arbeiten. Sie ist die elementare Folgerung aus der hier ausdrücklich stärkeren Forderung, dass Reziprozität ohne zusätzliches gradänderndes Objekt wirkt.

In vier Dimensionen wirkt die konstitutive Antwort $κ$ auf dem sechs dimensionalen Raum der Zweiformen. Reziprozität verlangt bis auf eine positive Skala

\[
κ^2=-λ^2 I.
\]

Das ist die Closure Relation. Sie hat genau die algebraische Form einer Hodge Operation auf Zweiformen. Die allgemeine Wellenfläche einer lokalen linearen Antwort ist die quartische Tamm Rubilar Fläche

\[
\mathcal G^{ijkl}(χ)q_iq_jq_kq_l=0.
\]

Fehlt Doppelbrechung, wird dieses Polynom zum Quadrat einer quadratischen Form. Hyperbolizität wählt deren Signatur $1+3$. Damit ist die konforme Lorentzmetrik rekonstruiert. Eine Skala, ein möglicher Axionanteil und die Dynamik von $κ$ sind dadurch noch nicht bestimmt.

## Kleinstes unterscheidendes Experiment

Das Skript [kausaltest.py](kausaltest.py) prüft eine exakt lösbare Klasse. Für eine uniaxiale Antwort

\[
ε=\operatorname{diag}(a,a,b)
\]

ist nach Entfernung des longitudinalen Maxwell Faktors die physische charakteristische Fläche

\[
P=
(-a\omega^2+k_x^2+k_y^2+k_z^2)
(-ab\omega^2+a k_x^2+a k_y^2+b k_z^2).
\]

Im Testfall $b=a$ erfüllt $κ$ exakt die Closure Relation und beide Faktoren fallen zusammen:

\[
P=(-a\omega^2+\lvert\mathbf k\rvert^2)^2.
\]

Die Hesse Matrix der quadratischen Form besitzt für $a=1$ die Eigenwertsignatur $(-1,+1,+1,+1)$.

Die Kontrolle setzt $a=1$ und $b=2$. Dann ist der Closure Defekt ungleich null und

\[
P=(-\omega^2+k_x^2+k_y^2+k_z^2)
(-2\omega^2+k_x^2+k_y^2+2k_z^2).
\]

Es entstehen zwei verschiedene Kegel. Das Experiment unterscheidet daher genau die tragende Annahme. Es zählt keine Dimensionen nach und passt keine Zahl an.

Erfolgskriterium: Dimension vier, verschwindender Closure Defekt und ein perfektes Quadrat im Testfall, dazu nichtverschwindender Defekt und zwei verschiedene Faktoren in der Kontrolle.

Abbruchkriterium: Wenn die Kontrolle ebenfalls ein perfektes Quadrat liefert oder Closure die beiden Fälle nicht trennt, trägt die vorgeschlagene Identifikation nicht.

## Weg bis Quelle, Zustand, Materie und Gravitation

### Quelle und erste Dynamik

Die Quelle ist zunächst konkret: $J$ ist der erhaltene elektrische Strom. Die Gleichungen $dF=0$ und $dH=J$ mit $H=κ(F)$ bilden ein geschlossenes hyperbolisches Feldsystem, sobald $κ$ gegeben ist. Raumzeit und Dynamik entstehen hier zusammen, weil derselbe Hauptteil sowohl die Entwicklung als auch den Lichtkegel bestimmt.

### Freier Fall und volle Metrik

Der Lichtkegel liefert nur die konforme Klasse. Beobachtete freie Fallbahnen ergänzen nach Ehlers, Pirani und Schild eine projektive Struktur. Ihre Kompatibilität mit den Lichtstrahlen liefert eine Weyl Verbindung. Die experimentelle Abwesenheit des zweiten Uhreneffekts reduziert diese auf die Levi Civita Verbindung einer Lorentzmetrik, eindeutig bis auf einen konstanten Maßstab.

### Zustand

Hier liegt die erste echte Lücke nach der Kausalrekonstruktion. Der Hauptteil legt keine globale Quantenvakuumfunktional fest. Auf einer allgemeinen gekrümmten Raumzeit gibt es aus dem Kegel allein keinen ausgezeichneten Grundzustand. Benötigt wird ein mikroskopischer Prozess, der zugleich $κ$, einen positiven lokalen Zustand und dessen Entwicklung erzeugt. Genau das darf nicht durch die nachträgliche Wahl eines Hadamard Zustands ersetzt werden.

### Materie und TFPT

Die prämetrische Algebra liefert sechs reelle Zweiformrichtungen, nach Wahl des Lorentzkegels drei komplexe selbstduale Richtungen. Sie liefert keinen natürlichen Fünferträger und keine Teilung $3+2$. Der vorhandene TFPT Weg von P1 über vier Marken zu P2 kann erst angeschlossen werden, wenn P1 als Zustand oder Defekt des kausalen Feldsystems hergeleitet ist.

Das ist ein harter Typentest: Eine Zerlegung $3+3$ der elektromagnetischen Zweiformen darf nicht wegen ähnlicher Zahlen als $3+2$ gelesen werden. Ohne einen Operator, der den P1 Randkern samt Einheitswindung und Normierung aus dem kausalen Zustand erzeugt, bleibt Materie offen.

### Gravitation

Die Geometrie allein liefert noch keine Einstein Dynamik. Zwei bekannte Anschlusswege zeigen genau, welche Zusatzannahme fehlt:

* Ehlers, Pirani und Schild rekonstruieren die Geometrie aus Licht und freiem Fall, nicht ihre Feldgleichung.
* Lovelocks Satz reduziert in vier Dimensionen eine lokale, natürliche, symmetrische, divergenzfreie Metrikgleichung zweiter Ordnung auf Einstein Tensor plus kosmologischen Term. Diese Eigenschaften sind neue Voraussetzungen.
* Jacobsons lokale Horizontthermodynamik liefert die Einstein Gleichung aus $δQ=T dS$, wenn Entropie proportional zur Fläche ist. Der Flächenkoeffizient und der zugrunde liegende Zustand bleiben Eingaben.

Für TFPT ist der sinnvolle Anschluss deshalb klar: Der gesuchte gemeinsame Zustand muss den konstitutiven Kegel und die lokale Horizontentropie aus derselben Antwort erzeugen. Erst dann wären Zustand, $G$ und Gravitation nicht bloß an die fertige Raumzeit angehängt.

## Vergleich mit der bisherigen Exploration

Der bisherige H6 Zweig beginnt bei der Naht und benötigt für die Lorentzlesart das konjugierte Dublett $T'$ beziehungsweise $I_{pos}$. Die prämetrische Route beginnt vorher bei Erhaltung und Antwort. Sie ersetzt $T'$ durch zwei andere, offen benannte Voraussetzungen: reziproke lokale Vakuumantwort und genau einen hyperbolischen Wellenkegel.

Der Atlas führt die Auswahl von Metrik, $3+1$ und lokaler Kausalstruktur weiterhin als offen. Die Suche nach „Doppelbrechung“, „Fresnel“, „Ehlers Pirani Schild“ und „konstitutiver Tensor“ lieferte in den jeweils ersten zwanzig Treffern keinen direkten geprüften Vorläufer. Das ist ein Neuheitshinweis, kein Beweis, dass der gesamte Korpus den Gedanken nie berührt hat.

Der nächste verwandte Atlasbefund ist `finding:tfpt-joint-time-20261111-shared-light-cone-test`. Er prüft einen gemeinsamen Lichtkegel erst nach Wahl der tatsächlichen OS Zeit und einer kanonischen lokalen Zweitableitungskinetik. Der neue Test setzt diese Zeit nicht ein. Er rekonstruiert den Kegel aus dem konstitutiven Hauptsymbol und zeigt mit der anisotropen Kontrolle, wann genau die Rekonstruktion scheitert.

## Urteil und nächster entscheidender Test

**Urteil: offen als kausaler Gesamtweg, exakt positiv für die Raumzeitstufe, kein TFPT Kandidat.**

Der Ansatz beseitigt einen realen Engpass, weil er vier Dimensionen und Lorentzsignatur aus einer einzigen physikalischen Antwortregel gewinnt. Er löst den Ursprung nicht, solange die Regel und ihr Zustand gewählt sind. Die erste offene Brücke lautet:

> Erzeugt ein gemeinsamer mikroskopischer Zustand eine reziproke, nicht doppelbrechende Antwort $κ$ und zugleich den P1 Randkern als seinen topologischen Defekt?

Der kleinste nächste Test ist daher nicht ein weiteres Medium. Aus dem vollständigen TFPT Quellhauptsymbol muss, ohne eingesetzte Raumzeitmetrik, der effektive konstitutive Tensor extrahiert werden. Dann sind drei Dinge auf demselben Zustand zu prüfen: Closure, perfektes Quadrat der Tamm Rubilar Fläche und der Fredholmindex des behaupteten P1 Defekts. Scheitert einer dieser drei Tests, verbindet diese Route Raumzeit und TFPT nicht.

## Primärquellen

* G. F. Rubilar, Y. N. Obukhov und F. W. Hehl, [Generally covariant Fresnel equation and the emergence of the light cone structure in linear premetric electrodynamics](https://arxiv.org/abs/gr-qc/0109012), insbesondere die quartische Fresnel Gleichung und die Closure Relation.
* C. Lämmerzahl und F. W. Hehl, [Riemannian light cone from vanishing birefringence in premetric vacuum electrodynamics](https://arxiv.org/abs/gr-qc/0409072), zur Reduktion auf einen hyperbolischen Lorentzkegel.
* A. Favaro und L. Bergamin, [The non birefringent limit of all linear, skewonless media and its unique light cone structure](https://arxiv.org/abs/1008.2343), zur allgemeinen skewonfreien Klasse und der Rolle der Hyperbolizität.
* J. Ehlers, F. A. E. Pirani und A. Schild, [The geometry of free fall and light propagation](https://bpfoundations.org/plans/ehlers2012.pdf), Originalarbeit von 1972 in der Wiederveröffentlichung von 2012.
* D. Lovelock, [The Einstein Tensor and Its Generalizations](https://doi.org/10.1063/1.1665613), Journal of Mathematical Physics 12, 498 bis 501, 1971.
* T. Jacobson, [Thermodynamics of Spacetime: The Einstein Equation of State](https://arxiv.org/abs/gr-qc/9504004), zur bedingten Gravitationsbrücke über lokale Horizonte.

Keine Atlasdaten, Ledgerwerte oder bestehenden TFPT Dateien wurden durch diese Exploration verändert.
