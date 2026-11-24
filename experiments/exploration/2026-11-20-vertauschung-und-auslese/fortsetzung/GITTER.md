# Vom Fano Code zu einem markierten Teilgitter von E₈

**Ergebnis.** Sobald das ursprüngliche Prozessmodell eine dreielementige Marke mit ihrer Permutationswirkung und den sieben Fano Linien liefert, lässt sich eine konkrete Einbettung von A₃ und D₅ berechnen. Die drei Familienrichtungen entstehen als Differenzen von Komplementpaaren. Die übrigen fünf Achsen tragen D₅. Zwei Orientierungen eines aus denselben Marken gebildeten Sechsecks geben dem Fünferträger eine echte Zweierbahn. Die physische Bedeutung der gewählten Auslesung und die Identifikation der Prüfkoordinate mit dem Orientierungsgedächtnis bleiben Zusatzannahmen.

## Der Code und seine Vergleichsklasse

Indiziere die sieben Punkte mit den von null verschiedenen Vektoren aus F₂³. Für verschiedene Punkte a,b ist die zugehörige Linie {a,b,a+b}. Die sieben Indikatorwörter spannen einen binären Code mit Gewichtsverteilung 1+7z³+7z⁴+z⁷ auf. Es gibt 16 lineare Regeln für ein zusätzliches Prüfbit auf diesem vierdimensionalen Code. Genau eine macht alle Gewichte gerade. Sie ergänzt die Gesamtparität und führt sogar auf Gewichte 0,4,8 mit Verteilung 1+14z⁴+z⁸. Die neue Koordinate trägt den Index 000.

Die Kontrolle betrachtet alle 6 724 520 Familien aus sieben verschiedenen Dreiermengen auf derselben beschrifteten Siebenpunktmenge. Nach Paritätserweiterung haben die sieben Erzeuger Gewicht vier. Ihr gesamter Spann ist genau dann doppelt gerade, wenn je zwei Erzeuger gerades Skalarprodukt besitzen. Zwei verschiedene ursprüngliche Dreiermengen müssen sich daher jeweils in genau einem Punkt schneiden. Dann erscheinen unter den sieben Linien alle 21 Punktpaare genau einmal. Die erschöpfende Zählung liefert **30 von 6 724 520**. Die 30 sind Beschriftungen derselben Fano Ebene. Die Quote beschreibt diese spezifische Kombinationsklasse.

Im reellen Punktkoordinatenraum ist das Gitter
$$
L=\left\{z/\sqrt{2}\ \middle|\ z\in\mathbb{Z}^{8},\ z\bmod 2\in C\right\}.
$$
Der Urbildindex 16 und der Skalierungsfaktor 1/16 im Achtvolumen geben Kovolumen eins. Die doppelte Geradheit des Codes liefert ein gerades Gitter. Bei Normquadrat zwei bleiben 16 Achsenwurzeln und 14·16=224 Wurzeln auf Gewichtsviererträgern. Die direkt berechnete einfache Wurzelbasis hat die drei E₈ Astlängen 1,2,4.

## Das Teilgitter aus der Komplementauslesung

Fasse einen Bitvektor als Teilmenge der drei Marken auf. Das Komplement ist Translation um u=111. Das Paar {000,111} ist die Referenz. Die sechs anderen Punkte bilden drei Paare aus einem Einzelpunkt x und seiner komplementären Zweiermenge y. Wähle auf jedem Paar die orthonormalen Richtungen
$$
f_C^+=\frac{e_x+e_y}{\sqrt{2}},\qquad f_C^-=\frac{e_x-e_y}{\sqrt{2}}.
$$
Ein Vertauschen der beiden Punkte ändert allein die Orientierung der Differenzlinie. Die Auslesung der drei **ungeordneten** aktiven Paare summiert die Komponenten auf jedem Paar. Ihr Kern hat als Basis genau die drei aktiven Differenzen. Diese nennen wir F. Die orthogonale Ergänzung V umfasst die drei aktiven Mittelrichtungen und die beiden Referenzrichtungen.

Die Koordinaten einer Gitterwurzel im Paarenrahmen lauten t/2, wobei die beiden Einträge des jeweiligen Paares t=(z_x+z_y,z_x−z_y) sind. Die exakte Rechnung ergibt alle 112 ganzzahligen D₈ Wurzeln ±f_i±f_j sowie 128 halbzahlige Gewichte einer Spinorchiralität. Ein Gittervektor mit Träger ganz in F oder ganz in V muss ganzzahlig sein, weil seine übrigen Koordinaten verschwinden. Im ganzzahligen D₈ Anteil ist allein die gerade Gesamtsumme vorgeschrieben. Daher sind die beiden Schnitte exakt die primitiven Teilgitter D₃≅A₃ und D₅. Ihre Wurzelsysteme haben 12 und 40 Elemente.

Der Homomorphismus χ(t/2)=Σ_{j∈F}t_j mod 4 hat Kern D₃⊕D₅. Der Kern besitzt Index vier, denn die halbzahligen Gewichte liefern beide ungeraden Werte. Die Klassen 0,1,2,3 enthalten **52,64,60,64** Wurzeln. Die Diskriminantengeneratoren der zwei Koordinatenblöcke haben Normquadrate 3/4 und 5/4. Die Summe zwei erklärt, weshalb die Verklebung ein gerades Gitter ergibt.

**Basislinie.** Jede der 56 Auswahlen von drei Paarenrahmenachsen liefert abstrakt denselben Typ D₃⊕D₅. Unter den reinen Permutationen der drei Marken bleiben genau zwei der 56 Dreierauswahlen invariant, die drei Mittelrichtungen und die drei Differenzen. Erst die Forderung, dass Familienrichtungen bei der Auslesung ungeordneter Paare unsichtbar sind, zeichnet die Differenzen aus. Ohne den Komplementrahmen bewahren die 168 Fano Kollineationen keines der 7 560 eingebetteten A₃ Wurzelsysteme. Diese Aussagen betreffen unterschiedliche Symmetriegruppen.

## Die beiden Orientierungen

Die sechs nichtleeren echten Teilmengen bilden unter Inklusion einen Sechserkreis. Seine ganzzahlige nullte Homologie ist eine triviale Linie. Die erste Homologie ist ebenfalls eine Linie. Gerade Markenpermutationen erhalten ihren primitiven Kreiszyklus, ungerade kehren ihn um. Damit sind die beiden möglichen primitiven Orientierungen eine kanonische Zweiermenge. Die drei Komplementpaare und diese beiden Orientierungen bilden einen markierten Fünferträger. Das Paar aus i und {j,k} entspricht der Transposition (jk). Eine zyklische Orientierung der drei Marken entspricht dem zugehörigen Dreizyklus. Die Zuordnung respektiert jede Umbenennung der Marken.

Auf den ursprünglichen **positiven Punktvektoren** permutieren die drei Marken die drei aktiven Mittelrichtungen. Sie fixieren beide Referenzlinien. Wer zusätzlich alle sieben ursprünglichen Fano Punktvektoren positiv erhält, hat genau zwei reelle diagonale Hebungen der Wirkung auf die zusätzliche Prüfkoordinate. Bleibt diese trivial, gibt es auf den Referenzrichtungen zwei Einpunktbahnen. Trägt e₀₀₀ das Signum, so vertauscht jede ungerade Permutation die zwei Vektoren (e₁₁₁±e₀₀₀)/√2 mit positivem Koeffizienten. Sie entsprechen den beiden Orientierungszuständen. Die Prüfkoordinate des Fano Codes ist bei dieser Konstruktion eine eigene Linie. Ihre Identität mit einem Vakuum in einer Hecke Darstellung ist unbewiesen.

Lässt man ortsabhängige Vorzeichen der alten Punktvektoren zu, existieren 256 diagonale Hebungen, beziehungsweise 16 Klassen unter diagonalen Vorzeichenwechseln. 32 der Hebungen erhalten den Raum der aktiven Differenzen. Von diesen vertauschen 16 bei ungeraden Permutationen die beiden Referenzlinien. Die Festlegung der sieben alten Punktvektoren ist somit für die genannte Eindeutigkeit wesentlich.

Das Signum auf der Prüfkoordinate erhält das binäre Prüfbit, denn +1 und −1 stimmen modulo zwei überein. Es ist eine echte Gitterisometrie. Ein zusätzlich geforderter äquivarianter ganzzahliger Prüfwert scheitert bereits am Siebenpunktwort aus lauter Einsen. Dieses Wort bleibt unter allen Permutationen fest und hat Parität eins. In einer ganzzahligen Signumlinie ist der einzige unter einer ungeraden Permutation feste Wert null. Ein orientiertes Gedächtnis kann daher eine Vorzeichenlinie und eine Paritätsklasse bereitstellen, aber keinen bevorzugten ungeraden ganzzahligen Repräsentanten.

## Reproduktion und physische Reichweite

`python3 paired_frame.py` reproduziert Wurzeln, Verklebung, alle 56 Rahmenkontrollen, die beiden invarianten Dreierauswahlen, sämtliche Vorzeichenhebungen und die Wirkung auf den primitiven Sechserzyklus. Das Skript benötigt nur die Python Standardbibliothek. Die Komplementbildung wird durch die vorliegende Dreimarkenstruktur ermöglicht. Ein nativer Messprozess, der gerade über diese Paare summiert und ihre beiden Kreisorientierungen bewahrt, würde die physische Auswahl stützen. Die Gitterrechnung allein bestimmt weder die Randamplitude c₃ noch Zustand, Zeitentwicklung oder Gravitation.
