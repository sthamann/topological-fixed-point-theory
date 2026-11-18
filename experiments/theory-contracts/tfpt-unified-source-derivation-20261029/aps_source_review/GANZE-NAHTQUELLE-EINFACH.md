# Was die ganze gekoppelte Quelle mit den Nahtmoden macht

Stand vom 29. Oktober 2026. Der vollständige fundamentale Abschluss ist noch nicht bewiesen.

Die Rechnung aus dem einzelnen skalaren Sektor ist jetzt auf die ganze vorhandene konstante Quelle übertragen. Das Ergebnis lässt sich einfach beschreiben: Eine Nahtmode bleibt genau dort erhalten, wo die tatsächliche Kopplungsmatrix eine Nullrichtung besitzt. Die Kopplung kann die vorher getrennten Richtungen miteinander verbinden und ihre Nullmoden aufheben.

Dafür gibt es eine vollständige Formel und einen unabhängigen Beweis auf denselben Randbedingungen. Der Nachweis berücksichtigt auch den Unterschied zwischen einem Operator und seinem Adjoint. Eine Reihe abgeschnittener Computersimulationen ist nicht nötig.

Für die schon vorhandenen exakt geprüften Eingabefälle folgt unmittelbar:

| Vorhandener Eingabefall | Kernmoden | Kokernmoden |
| --- | ---: | ---: |
| Invertierbare Kopplung | 0 | 0 |
| Erste rangdefiziente Kontrolle | 2 | 2 |
| Zweite rangdefiziente Kontrolle | 1 | 1 |

Diese Zahlen betreffen die angegebenen mathematischen Kontrollen. Die Eingabematrizen wurden nicht als physische Massen aus TFPT hergeleitet. Der volle Operator hat in dieser Klasse gleich viele positive und negative interne Richtungen. Sein gewöhnlicher Index ist null. Das entfernt weder die alte geometrische Nahtrechnung noch alle relativen Feldantworten.

Die Klasse ist damit entschieden: Mehr Zufallsversuche oder ein genetischer Suchlauf würden an dieser Stelle dieselben Kernbedingungen wiederfinden. Ein anderer ursprünglicher Mechanismus muss seine tatsächlichen anderen Voraussetzungen zeigen, zum Beispiel ein aus der Quelle hergeleitetes variables Feld oder eine andere Nahtwirkung. Wir wählen keine Matrix so, dass drei gewünschte Lösungen herauskommen.

Zwei konkrete Zuordnungen aus Claudes letzter Antwort wurden außerdem korrigiert. Erstens wirkt die Konjugation C0 auf die ursprünglichen komplexen Zahlen anders, als dort angenommen wurde. Sie kehrt die Nahtphase nicht automatisch um. Der tatsächliche Produktreflektor transportiert die Nahtbedingung in ihre inverse Bedingung; die bereits bewiesene Reflexionspositivität der offenen Quelle bleibt gültig.

Zweitens muss die Phase einer Kreisbewegung von der Phase eines geometrischen Paralleltransports unterschieden werden. Die vorgegebene Viertelphase der Kreiswirkung bleibt richtig. Der natürliche Paralleltransport der ausdrücklich berechneten Gaußmode hat im untersuchten Fall einen anderen Wert. Eine zusätzlich festgelegte Verbindung darf diese beiden Größen nicht stillschweigend gleichsetzen.

Der nächste ursprüngliche Anschluss ist konkret. Die alte Kurvenrealstruktur E0=Q tau0 kehrt auch die Basisrichtung um: z wird zu i mal dem komplex konjugierten z. Auf dem Einheitskreis entspricht das theta→pi/2−theta. Wir prüfen diese wirkliche geometrische Operation zusammen mit den vorhandenen Spinfaserlifts, beiden Randseiten und dem vollständigen Feldoperator. Ein reiner Vergleich der Nahtphasen genügt nicht; auch die Verbindung und der Definitionsbereich müssen mitgehen.

Der alte Beweis der einmaligen Dreipunktnorm bleibt abgeschlossen. Die gemeinsame Quelle erhält ihre bisherigen Zustands und Zeitdaten. Noch nicht bewiesen ist, dass der jetzt geprüfte gesamte Konstruktionstransport genau die ursprüngliche physische TFPT Erzeugungsregel ist.

Die geprüften Beweise heißen FULL-COUPLED-SEAM-INDEX.txt, WHOLE-CLUTCH-ENERGY-INDEPENDENT.txt, CLAUDE25-CLUTCH-REAL-TYPE-REVIEW.txt und DEFECT-TIME-PHASE-INDEPENDENT.txt. Sie sind im gemeinsamen Atlas verzeichnet. Alle zehn zugehörigen Quellenintervalle wurden nach dem Speichern zurückgelesen; ATLAS-WHOLE-CLUTCH-READBACK.json enthält die Verweise.
