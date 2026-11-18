# Der Anschluss lässt sich weiter vereinfachen

Die ursprüngliche Nahtfunktion erzeugt auf der bekannten dreiblättrigen Kurve den benötigten Randfaktor. Für die bereits klassifizierte komplexe Symmetriestruktur genügt sogar der Vergleich des verschobenen mit dem unverschobenen Zustand desselben mathematischen Spinraums. Ein zusätzlicher physischer skalarer Referenzträger ist dafür keine notwendige Voraussetzung.

Diese Ergänzung präzisiert die nächste Aufgabe aus `ERGEBNIS-EINFACH.md`. Der dort bewiesene Vergleich zweier Kabel bleibt richtig. Er ist eine mögliche exakte Normierung, aber nicht die einzige Route zur bereits bestimmten Modulklasse.

## Einfach vorgestellt

Der ursprüngliche Umlauf läuft durch ein Kabel mit drei Adern. Verglichen werden das einmal verschobene Kabel und dasselbe Kabel vor der Verschiebung. Aus den drei zurückbleibenden Randwerten entsteht der gesuchte Faktor. Die gesamte Rechnung benutzt den alten Umlauf und die alte dreiblättrige Kurve.

Die endliche Clock zeigt nur die Stellung eines Zeigers. Der ursprüngliche ganzzahlige Umlauf zählt zusätzlich, wie oft er gedreht wurde. Dadurch wird aus einer zuvor nur bis auf Restklassen bestimmten Potenz genau die einmalige Einfügung.

## Was jetzt mathematisch feststeht

Auf der Originalkurve gilt

\[
z=\frac{\eta}{1-w^3},\qquad \operatorname{div}(z)=D_q-D_p.
\]

Der Nahtshift durch z hat auf den drei holomorphen Blättern genau drei fehlende Randwerte. Seine relative Determinantenlinie gegenüber der Identität ist für die bereits bestimmte Spinfamilie E

\[
\delta(E,M_z)\otimes\delta(E,I)^{-1}
=N_q(S_0)^{-1}\otimes B_{\mathcal O}^{-1}.
\]

Der letzte Faktor trägt genau das bereits zugelassene Vorzeichen des Quarterlifts. Deshalb entsteht dieselbe vollständige komplexe Theta, N64 und Clockmodulklasse. Für n Umläufe entsteht die n te Potenz; der ursprüngliche Spektralfluss eins wählt einen Umlauf.

Der Vergleich mit einem getrennten skalaren Kabel entfernt zusätzlich die Ordnungsdeterminante und liefert den wörtlichen Repräsentanten Nq hoch minus eins. Diese stärkere Normierung muss nicht als Voraussetzung für die schon bewiesene komplexe Modulklasse eingeführt werden.

## Welche physische Information erhalten bleiben muss

Die Ordnungsdeterminante ist als fermionische Linie ungerade. Die drei fehlenden Modi tragen diese Parität ebenso wie der ursprüngliche einzelne Modus: drei und eins sind beide ungerade. Beim Vergleich zweier verschobener Kabel würden sich zwei solche Paritäten aufheben.

Für gewöhnliche komplexe Symmetriematrizen bleibt davon das bekannte Vorzeichen. Bei einer fermionischen Feldwirkung oder einer Pfaffianlinie ist die Parität zusätzliche Information. Sie darf nicht stillschweigend entfernt werden.

## Der nächste entscheidende Schritt

Am tatsächlichen ursprünglichen Randproblem ist jetzt der Shift gegen die Identität auf demselben Spinraum zu prüfen. Dabei müssen der wirkliche Cauchydatenraum, die Polarisation, der Zustand, die Reflexion, die Zeit und die fermionische Orientierung zusammenpassen. Der allgemeine Split D_ref/D_rel liefert diesen Nachweis allein noch nicht.

Der direkte Wechsel vom Spinträger zu einem skalaren Referenzträger ist damit keine vorgeschaltete Pflicht mehr. Die Prüfung verfolgt die einfachere vorhandene Nahtoperation und behält ihre ungerade Defektlinie bei.

Die vollständige fundamentale TFPT Herleitung ist weiterhin nicht abgeschlossen. Neu bestimmt ist die konkrete mathematische Herkunft des Randfaktors und seiner einmaligen Potenz, zusammen mit der einfacheren Referenzwahl im bisherigen komplexen Modulumfang. Die ursprüngliche physische Feldabbildung und die integrale Gitterherkunft bleiben eigene Beweisschritte.

## Nachweise

`source/ORIGINAL-BOUNDARY-NORM-ORIGIN.txt` enthält den vollständigen Naht und Jetbeweis. `uniqueness/THREE-CABLED-HARDY-SEAM-REVIEW.txt` prüft ihn unabhängig. `uniqueness/SAME-SPIN-REFERENCE-SIMPLIFICATION.txt` enthält die einfachere Referenz auf demselben Spinraum und die Grenze der Gradierung. `originals/ACTUAL-REFERENCE-SPIN-TEST.txt` vergleicht die tatsächlichen Originaloperatoren mit diesen Kandidaten.
