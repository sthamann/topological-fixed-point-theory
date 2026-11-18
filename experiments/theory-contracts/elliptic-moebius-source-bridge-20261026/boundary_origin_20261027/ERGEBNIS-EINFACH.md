# Was der neue Anschluss erklärt

Die bisher fehlende einmalige Randkorrektur lässt sich jetzt aus einer konkreten mathematischen Nahtoperation gewinnen. Dazu verwenden wir die ursprüngliche Nahtfunktion auf der bereits vorhandenen dreiblättrigen Kurve. Die Identifikation dieser Operation mit dem vollständigen physischen Operatorpaar bleibt gesondert zu beweisen.

## Das Bild dazu

Stell dir ein Kabel mit drei Adern vor. Der ursprüngliche Nahtumlauf wird einmal ausgeführt und auf alle drei Adern übertragen. Auf jeder Ader entsteht genau eine kleine Verschiebung. Wir vergleichen denselben Umlauf mit Spintransport und mit dem zugehörigen skalaren Referenztransport. Die gemeinsamen Beiträge fallen heraus. Übrig bleibt genau der Faktor, den der bereits hergeleitete vierdimensionale Symmetrieraum braucht.

Die drei Adern sind keine neu gewählten Bauteile. Sie sind die drei Blätter der bereits aus der ursprünglichen Gleichung gewonnenen Kurve.

## Die vollständige mathematische Verbindung

Die Originalquelle liefert

\[
w^3=\frac{z^2-1}{z^2+1},\qquad \eta=z(1-w^3).
\]

Damit ist die ursprüngliche Nahtfunktion

\[
z=\frac{\eta}{1-w^3},\qquad \operatorname{div}(z)=D_q-D_p.
\]

Ihre drei Nullstellen sind genau die drei bereits bestimmten q Punkte. Der alte einmalige Kreisumlauf wird auf der Kurve zur Multiplikation mit dieser Funktion. Im zugehörigen holomorphen Randraum ist der Quotient durch diese Multiplikation genau der Raum der drei Punktwerte.

Der relative Defekt von Spintransport und gleich verkabeltem Referenztransport ist deshalb

\[
\delta(S_0\otimes L,z)\otimes\delta(L,z)^{-1}
=\operatorname{Nm}_{D_q}(S_0)^{-1}.
\]

Das gilt für jede gemeinsam mitgeführte Grad null Koeffizientenlinie L. Auch die ursprünglichen Vorzeichenlinien müssen dafür nicht weggelassen werden. Ihre Beiträge und die Vertauschungszeichen der drei Blätter heben sich im relativen Vergleich auf.

Für n Umläufe ergibt dieselbe Rechnung die n te inverse Potenz. Der ursprüngliche ganzzahlige Spektralfluss eins wählt damit genau eine Kopie innerhalb dieses Randmodells. Die endlichen Clockrelationen erkennen dagegen nur Restklassen der Potenz. Daher hatten sie diese Frage allein nicht entscheiden können.

Zusammen mit dem bereits bewiesenen Theta Anschluss erhalten wir die Kette

\[
\text{ursprünglicher einmaliger Umlauf}
\longrightarrow
\text{relative Randlinie auf drei Blättern}
\longrightarrow
N_q(S_0)^{-1}
\longrightarrow
H^0(\mathcal O(2\Theta))\otimes N_q(S_0)^{-1}.
\]

Der letzte Raum besitzt bereits die vollständig geprüfte komplexe N64 und Clockwirkung. Diese alte Rechnung wurde wiederverwendet.

## Was die zusätzliche Prüfung korrigiert

Ein einzelner Fourierdefekt und der Defekt aller drei Blätter sind verschiedene Operatoren. Ihr Index ist minus eins beziehungsweise minus drei. Die richtige Viertelphase entsteht beim vollständigen relativen Dreiervergleich.

Ebenso darf man die Punktwerte nicht beliebig als neue Nullbedingungen in einen gewöhnlichen Diracraum schreiben. Auf dessen üblicher H1 Domäne sind solche Punktauswertungen nicht stetig. Im hier benutzten holomorphen Hardyraum sind sie stetig und die Multiplikation mit z liefert den Quotienten tatsächlich. Diese Unterscheidung ist für den Beweis entscheidend.

## Wo die fundamentale Herleitung jetzt steht

Die mathematische Herkunft des Faktors ist im beschriebenen Hardy Vergleich bestimmt, einschließlich seiner Potenz und sämtlicher benötigter endlicher Phasen. Offen ist die nächste physische Identifikation: Erzeugen die ursprünglichen Randoperatoren tatsächlich diesen gemeinsamen Spin und Referenzvergleich, mit derselben Polarisation, Reflexion, Zeit und demselben Zustand?

Die bisher gelesenen Originaltexte legen den Nahtumlauf und die dreiblättrige Karte fest. Ihre Definition von D_rel und D_ref als ungeradem und geradem Operatoranteil beweist den benötigten Vergleich noch nicht. Deshalb lautet die nächste Aufgabe, dieses konkrete Operatorpaar an seiner Quelle zu identifizieren. Eine neue Suche nach passenden Matrizen oder einer passenden Potenz ist dafür nicht nötig.

Auch danach müssen die integrale E8 Gitterherkunft und der ursprüngliche lokale Feldtransport mit seinen Produkten und seinem Zustand erhalten bleiben. Der vorliegende Anschluss ersetzt diese Schritte nicht durch einen komplexen Basiswechsel.

## Belege

Der vollständige Beweis steht in `source/ORIGINAL-BOUNDARY-NORM-ORIGIN.txt`, der unabhängige Review in `uniqueness/THREE-CABLED-HARDY-SEAM-REVIEW.txt`. `originals/INVERSE-Q-POINT-NORM-ORIGIN.txt` dokumentiert die gelesenen Originalstellen und die analytische Grenze der Punktrestriktion. `uniqueness/NORM-POWER-UNIQUENESS.txt` enthält den vollständigen Test der scheinbar verschiedenen Potenzen. `ETA-CONVENTION-AND-BOUNDARY-FACTOR.txt` klärt, weshalb die vorhandene Eta Phase den fehlenden Operatorvergleich nicht allein beweist.

Status: neuer mathematischer Anschluss mit unabhängiger Prüfung. Keine Behauptung einer vollständig hergeleiteten physischen TFPT Theorie. Frühere ausgeschlossene Quellenidentifikationen bleiben ausgeschlossen.
