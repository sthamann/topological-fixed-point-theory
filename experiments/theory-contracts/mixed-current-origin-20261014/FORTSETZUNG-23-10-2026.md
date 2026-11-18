# TFPT: Was nichtunitäre Endpunktgrenzen leisten können

Stand: 23. Oktober 2026. Vertrag `UR.SEAM.MIXED_CURRENT_ORIGIN.01`.
**Ergebnis: PARTIAL. Kein mikroskopischer E₈-Spinorsektor hergeleitet.**

Die bisher offene Möglichkeit war konkreter als ein weiterer Ramond-Ansatz:
Kann das tatsächlich vorhandene neutrale Achtpräbild-Paritätspaar durch
nichtunitäre, verschmierte und renormierte Operatorfolgen einen einzelnen
Spinorsektor hervorbringen? Der frühere Energieausschluss für glatte
**unitäre** Phasen beantwortete diese Frage nicht.

Jetzt ist eine größere Klasse entschieden. Auf dem festgehaltenen
ursprünglichen Quellträger bleibt jeder vollständige Grenzzustand in der
alten Darstellung, wenn seine Feldnorm endlich und nicht null bleibt,
seine ursprüngliche Zeitantwort stetig ist und die vollständige Wirkung
der ursprünglichen Observablen erhalten bleibt. Eine beschränkte mittlere
Energie braucht der Beweis ausdrücklich nicht.

## Warum die Zeitantwort entscheidet

Für einen beliebigen präparierten Vektor \(v_N=T_N\Omega_N\) sei
\(n_N=\|v_N\|^2\). Die ursprüngliche positive Quellenergie ist \(H_N\).
Dann gilt für jedes \(E,t>0\):

\[
\|\mathbf1_{(E,\infty)}(H_N)v_N\|^2
\leq\frac{n_N-\langle v_N,e^{-tH_N}v_N\rangle}{1-e^{-tE}}.
\]

Anschaulich: Zustandsgewicht, das zu immer höheren Energien entweicht,
fehlt in jeder Antwort bei positiver imaginärer Zeit. Bleibt dagegen die
Zeitantwort bei null stetig und stimmt dort mit der Feldnorm überein,
kann kein endliches Gewicht auf diese Weise verloren gehen.

Die bereits bewiesene vollständige Spektralausschöpfung des ursprünglichen
Streifens liefert in jedem festen Energiefenster nur endlich viele
Zustände. Zusammen ergeben beide Aussagen Normkompaktheit entlang einer
Teilfolge. Mit den tatsächlichen Observablen-Matrixelementen folgt: Der
Grenzzustand ist in der alten globalen Darstellung normal.

Dasselbe folgt direkt aus der **reellen ursprünglichen Zeit**: Das Mittel
des Korrelationskerns über ein kleines Zeitintervall liefert die Schranke
mit \(1-\sin(aH_N)/(aH_N)\). Ein neuer euklidischer Zeitprozess wird für
den Ausschluss also nicht vorausgesetzt.

Im festgehaltenen D₈/even-CAR-Wörterbuch enthält die alte NS-Darstellung
die Klassen \(\mathbf1\oplus v\), aber weder \(s\) noch \(c\). Damit
erzeugt die untersuchte Klasse keinen neuen Spinorsektor. Dieser Schluss
setzt das benannte Wörterbuch voraus; er leitet weder acht unabhängige
Materiekanäle noch dessen physische P1-Auswahl her.

## Zwei konkrete Anschlussprüfungen

**Register und gemeinsame Energie.** Die lokale gerade Registerobservable
\(Z\) zwingt bei einem entfernten lokalen Feld sämtliche registerändernden
Anteile seiner Vakuumwirkung in Norm gegen null. Das gilt auch für Summen
verschiedener Registergrade. Zusätzlich sind die tatsächlichen skalierten
Vakuumenergieunterschiede aller vier Flussblöcke gleichmäßig beschränkt.
Der analytische Beweis benutzt die einzige Betragskante der vollständigen
Bloch-Grundenergiefunktion und deren Fourierkoeffizienten. Damit benötigt
der Schluss keine vier getrennt gewählten Energiesubtraktionen. Eine neue
Gaußregel, die \(Z\) global macht, verändert diese Voraussetzung.

**Ein tatsächlicher nichtunitärer Operator.** Berechnet wurde
\(e^{-sH_N}P_E\Omega_N\) am vollständigen ursprünglichen QWZ-Streifen,
einschließlich massiver Moden und beider Ränder. Seine Slater-Kovarianz
bestimmt alle endlichen CAR-Wörter; der Zeitkern folgt aus einem exakten
Determinantenquotienten. Eine zweite Berechnung durch direkte
Slater-Entwicklung und QR-Zerlegung bestätigt Kovarianz und Zeitkern.

Die endliche Kontrolle umfasst beide Viertelsektoren bei \(N_x=32,N_y=8\),
alle 512 Einteilchenmoden und alle 256 besetzten Moden. Beide Ränder
reagieren weiterhin; die acht **gemeinsamen Bloch-Restklassenladungen**
ändern sich nicht. Diese Ladungen werden nicht mit acht einzelnen
links-chiralen Cartanladungen gleichgesetzt. Das ist kein analytischer
Kontinuumsnachweis für diesen speziellen Filter.

In einem ausdrücklich getrennten chiralen Weyl-Vergleich mit der bereits
gegebenen Faltungszeit \(K=H_{\mathrm{cover}}/8\) kann man den Zeitkern
geschlossen summieren. Dieser feste Zeitmaßstab teilt keinen primitiven
Spektralfluss auf. Für gegenüberliegende Endpunkte und
\(r=e^{-2s}\) ergibt sich

\[
F_s(t)=\left[\frac{(1+r e^{-t})(1-r)}{(1-r e^{-t})(1+r)}\right]^4.
\]

Bei jedem \(s>0\) bleibt das ein neutrales Paar im alten Fockraum.
Beim Entfernen des Filters gilt \(F_s(0)=1\), aber für jedes feste
\(t>0\) geht \(F_s(t)\) gegen null. Genau dort geht die globale
Zeitstetigkeit verloren. Diese Vergleichsformel ist **keine** Identifikation
des scharfen ursprünglichen Gitter-Paritätsoperators mit einem unregulierten
bosonischen Vertex.

## Was als Herkunftsaufgabe bleibt

Das tatsächlich konstruierte neutrale Defektpaar bleibt erhalten. Sein
korrektes Vorzeichen auf allen 16 gefalteten Komponenten ist ein positiver
Operatorbefund. Auch die bekannten lokal normalen Spinordarstellungen und
die bedingte E₈-Erweiterung bleiben gültig.

Die neue Grenze lautet enger: Ein zusätzlicher Spinorsektor kommt nicht
durch raffiniertere Operatorfolgen auf dem **gleichen globalen Träger mit
dem gleichen vollständigen Zeit- und Observablenvertrag** hinzu. Jetzt
benötigt man eine ursprünglich begründete geladene Sektorabbildung, die
auch die Zeitantwort und lokale Observablenwirkung nachweist, oder eine aus
P1 abgeleitete andere lokale Netz-/Gaußstruktur. Die bereits bekannte
bedingte Formel für die geladene Energie ist kein neuer frei zu wählender
Hamiltonoperator; ihre mikroskopische Implementierung fehlt.

Das ist kein Ausschluss aller P1-Grenzprozesse. Insbesondere folgt aus
einem bloßen Darstellungswechsel nicht automatisch, dass die
Observablen-Matrixelemente auf alten niedrigen Energien geändert werden.
Bei unveränderten Operatoren kann stattdessen Zustandsgewicht ins hohe
Spektrum entweichen; dann scheitert die hier verlangte Zeitstetigkeit.
Umgekehrt beweist eine divergierende mittlere Energie allein noch kein
solches Entweichen.

Claude Opus 4.8 hat die tragenden Formeln unabhängig geprüft. Seine erste
Interpretation der Observablenvoraussetzung wurde korrigiert; die Review
ist mit diesen Grenzen dokumentiert. FE-GEN, ALG-EXH, chirales Maß,
P1-Quellenauswahl und die vollständige TFPT sind weiterhin nicht bewiesen.

Ausführlicher Beweis: `CONTINUOUS-TIME-NONUNITARY-SECTOR-TEST.txt`.
Ausführbare Kontrollen: `check_continuous_time_sector_test.py` und
`continuous-time-sector-results.json`. Review: `CLAUDE-CONTINUOUS-TIME-REVIEW.txt`.
