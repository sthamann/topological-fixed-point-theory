# TFPT: Windung und gemeinsame Clock als Auswahloperation

Stand 30. September 2026. Forschungs-ID UR.SOURCE.RULE_SEARCH.01, PARTIAL. Die vollständige fundamentale Lösung ist nicht gefunden. Dieses Ergebnis ist ein bedingter Herkunftskandidat für die bereits bekannte neutrale Antwort R, kein neues Naturgesetz und kein neuer Beweis ihrer physischen Realisierung.

Die Methode wurde gewechselt: keine Anpassung eines Filters an Zielzahlen, sondern Ableitung einer gemeinsamen Operation aus Windung, Clock-Wirkung und den tatsächlichen geladenen Quellenfeldern. Originalquelle: lokales Repository tfpt-theoryv4, HEAD d44f1040d733338e1c0596433b5b2578d7c00295. Vorbestehende Arbeitsänderungen bleiben separat.

## Was direkt aus der Quelle entsteht

Die adjoint-even Paaroperation auf den ursprünglichen 48 Spinorwurzelfeldern liefert B48=12T+R. Ihre Vakuumnorm ist 624=576+48. Ein unabhängiger Durchlauf der ursprünglichen Chevalley-/Ward-Rekursion liefert ebenfalls 624. Der Abzug der unverbundenen Identitätsantwort entfernt den Stressanteil 12T nicht.

Das 64er chirale Spinorpaket liefert B64=16T. Es ist ein bestimmtes Quellenpaket, nicht die gesamte Menge der 240 E8-Wurzeln. Deshalb ist

    B48 − (3/4) B64 = R

eine tatsächliche relative Operation auf denselben geladenen Paaren. Die bekannten Momentmatrizen und Normen werden dadurch nicht zu neuen Entdeckungen umetikettiert.

## Woher der Faktor 3/4 kommen kann

Auf dem nativen vierdimensionalen SU4-Familienraum wirkt die konkrete Originalclock als eine feste Linie plus ein irreduzibler dreidimensionaler Block. Verwendet werden die ursprünglichen Matrizen U=diag(1,i,−i) und M; beide haben Determinante 1. Der dreidimensionale Block ist die Rotationsdarstellung des S4, nicht die gewöhnliche Standarddarstellung mit Determinante gleich dem Permutationsvorzeichen.

Jeder mit dieser gemeinsamen Clock kommutierende hermitesche Generator hat die Form diag(a,b,b,b). Für die ausdrücklich eingeschränkte Klasse geschlossener Einparameter-Schleifen exp(2πitH) gilt: a,b sind ganzzahlig. Windung 1 erzwingt a+3b=1. Ein Rang-eins-Generator erzwingt zusätzlich b=0, a=1. Damit ist H=P_fest=diag(1,0,0,0) eindeutig.

Sein spurloser Anteil ist

    P_fest − I4/4 = diag(3,−1,−1,−1)/4 = −F.

Die relative Gegenseite F hat auf den 16 Feldern der ausgezeichneten Linie Gewicht −3/4 und auf den 48 übrigen Feldern Gewicht +1/4. Ihre tatsächliche Paarantwort ist genau R. Das Vorzeichen hängt von der gewählten relativen Orientierung ab und darf nicht nachträglich zur Übereinstimmung umgestellt werden.

F erfüllt exp(2πiF)=iI4 und Tr(F²)=3/4. Unter der zusätzlichen Vorgabe dieses zentralen Endpunkts ist es auch der minimale spurlose hermitesche Logarithmus. Der Beweis gilt für alle ganzzahligen Lifts: Eigenwerte 1/4+m_j mit Summe m_j=−1 haben Quadratsumme Summe m_j²−1/4≥3/4. Im Clock-Kommutanten ist der Minimierer eindeutig. Der zentrale Endpunkt iI4 wird hierbei nicht mit der nichtzentralen Familienclock oder dem Gaußschen Vierteldreher gleichgesetzt.

## Die erste noch nicht hergeleitete Identifikation

Die Rechnung setzt voraus, dass die ursprüngliche Nahtschleife gerade in diesem vierdimensionalen Familienraum durch eine solche Rang-eins-Einparameter-Schleife realisiert wird. Dies ist eine konkrete, prüfbare zusätzliche Identifikation. Der volle Quellenraum ist größer; Rang eins und Windung allein bestimmen dort keinen solchen Familienport. Die Originale geben den Naht-K1-Generator und den abstrakten Rang-eins-Transfer separat vor, liefern aber in den geprüften Stellen diese Feldabbildung nicht.

Auch innerhalb des Kandidaten folgen weder die physische Kopplung g=c3 noch die Gleichsetzung der konformen Sewingzeit τ mit α. Die Quellenantwort ist 48g⁴exp(−2τ); die vollständige TFPT-Zielantwort entsteht erst mit diesen zusätzlichen Zuordnungen. Die 3+1D-Realisierung, ein Primzahlgesetz und eine überlegene AI-Architektur werden nicht behauptet.

Die passende nächste Herkunftsfrage ist damit konkret: Realisiert die ursprüngliche relative Naht ihren primitiven Modus als diese clockverträgliche Familientransformation, und bestimmt dieselbe Nahtkopplung auch g und τ? Den Kandidaten weiter an die Zielzahlen anzupassen würde diese Frage umgehen.

Die verwendete Gittervertex-Expansion und Heisenberg-Konvention sind in der Primärliteratur angegeben: [Higher Airy structures, W algebras and topological recursion, §3.2](https://archive.mpim-bonn.mpg.de/2143/1/preprint_2018_69.pdf). Die TFPT-Daten stammen aus den vorhandenen Modulen neutral_source_response, source_charge_response und current_block_geometry; die physische Nahtdefinition wurde an den archivierten Originalabschnitten 01_boundary_kernel_source und 03_em_flavor_source abgeglichen.

Das beigefügte Rechenpaket enthält die Ableitung, den unabhängigen Ward-Vergleich und die Clock-/Windungsprüfung samt tatsächlichen Ausgaben. Es benötigt das vorhandene lokale TFPT-Repository und SymPy; es enthält keinen neuen Unterbau und keine neue Physiksimulation.
