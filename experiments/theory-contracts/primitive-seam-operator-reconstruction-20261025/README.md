# TFPT: Was die Prüfung der Originalpapers jetzt klärt

Stand: 25. Oktober 2026. Forschungsstand: **Teilergebnis mit vollständigen Beweisen im angegebenen Geltungsbereich. Die vollständige fundamentale TFPT Lösung ist damit noch nicht hergeleitet.**

Der neue Befund betrifft die tatsächliche Beweiskette. Im ursprünglichen elektromagnetischen Abschluss wird aus einer Windung ein einzelner Nahtmodus. Im Quantenfeldteil werden verbundene Korrelationen als vollständige Korrelationen eingesetzt. Beide Übergänge brauchen eine Korrektur. Für den ersten liegt jetzt ein präziser Rekonstruktionssatz vor. Für den zweiten ist die richtige Rekonstruktion angegeben und der fehlerhafte Schritt mit einer exakten Rechnung widerlegt.

## Das Bild: ein Rundweg und eine Tonleiter

Eine Windung sagt: Du gehst einmal um den Rundweg. Daraus folgt noch keine Tonleiter, kein bestimmter Anfangston und keine Regel, die jedem Schritt den nächsten Ton zuordnet.

TFPT besitzt bereits zwei passende Objekte: den Randoperator `B` und die Nahtoperation `u`. Wenn dieselbe ursprüngliche Konstruktion zeigt, dass `u` die von `B` geordneten Stufen immer um genau eins verschiebt, und das zugehörige tatsächliche Indexpaar den Wert eins trägt, wird die Tonleiter zwingend. Der Anfangsraum besteht dann aus genau einer Richtung. Alle weiteren Richtungen entstehen durch wiederholtes Anwenden von `u`, vorwärts und rückwärts.

**Dieser bedingte mathematische Satz ist bewiesen. Die Originalpassagen liefern die erforderliche Beziehung zwischen ihren beiden Operatoren bisher nicht.** Gleiche Windung, gleiche Viertelphase und gleicher Index erzwingen diese Beziehung nicht; dafür steht ein explizites Gegenbeispiel im Beweis.

## Der positive Reparatursatz

Die genaue Voraussetzung ist

\[
e^{itB}u e^{-itB}=e^{it}u,
\qquad 0\notin\sigma(B),
\qquad \operatorname{ind}(P_BuP_B)=-1,
\quad P_B=1_{(0,\infty)}(B).
\]

Dann wird der vorher frei benannte Nahtmodus kanonisch:

\[
\boxed{\Pi_\Sigma=P_B-uP_Bu^*=1_{(0,1]}(B),
\qquad\operatorname{rank}\Pi_\Sigma=1.}
\]

Aus einem normierten Vektor dieser eindeutig bestimmten Linie entsteht die vollständige Basis `uⁿψ`. Damit folgt auf dem gesamten ursprünglichen Hilbertraum

\[
B\simeq-i\partial_\theta+\epsilon,\qquad
u\simeq e^{i\theta},\qquad0<\epsilon<1.
\]

Eine zusätzliche tatsächliche Umlaufphase `exp(2πiB)=iI` bestimmt `ε=1/4`. Der Index allein bestimmt diese Phase nicht. Ein zusätzlicher freier Hardy Raum, eine neue Quelle oder eine Gleichsetzung des Calderón Projektors mit dem Quantenzustand werden für den Satz nicht gebraucht.

Im Original stehen `BΣ` und die Windung von `uΣ` nebeneinander. Für die Anwendung des Satzes muss der Spektralfluss ausdrücklich zu **diesem selben Operatorpaar** gehören. Ein anderer Operator mit derselben topologischen Klasse reicht dafür nicht. Der Satz beschreibt einen ausreichenden direkten Reparaturweg; er schließt andere, aus den Originalprinzipien hergeleitete Quantisierungen nicht aus.

## Der zweite Fehler: Es fehlen vollständige Korrelationen

Das Quantenfeldkapitel definiert verbundene Korrelationen als Ableitungen von `log Z`. Anschließend behandelt es diese Größen direkt als vollständige Schwinger Funktionen. Bereits eine gewöhnliche Gaußquelle zeigt den Fehler: Ihre verbundene Vierpunktfunktion ist null, ihre vollständige Vierpunktfunktion enthält die drei Paarungen. Lässt man diese Paarungen weg, kann die rekonstruierte Norm negativ werden.

Der richtige Weg lautet

\[
\text{verbundene Korrelationen}\;\longrightarrow\;
Z/Z[0]=\exp W\;\longrightarrow\;
\text{vollständige Korrelationen}\;\longrightarrow\;
\text{positive Hilbertraumrekonstruktion}.
\]

Für Fermionen gehören die richtigen Vorzeichen und sämtliche Wick Paarungen dazu. Die Positivität muss für das vollständige Funktional gelten. Die Korrektur konstruiert das in den Originalen noch nicht aus den primitiven Randdaten definierte `Z_adm` nicht automatisch. Sie repariert die logische Form des Schritts.

Eine frühere zu starke Forderung ist ebenfalls korrigiert: Geladene ungerade Fermionenfelder dürfen in einer größeren graduierten Feldalgebra leben. Erst deren tatsächlich zulässige beziehungsweise eichinvariante Operatoren müssen als physische Observablen wirken. Alle ungeraden Felder unmittelbar in die Observablealgebra einzubauen ist keine notwendige Bedingung.

## Wie das mit dem bereits Erreichten zusammenhängt

| Verbindung | Quellengeprüfter Stand |
|---|---|
| Ursprüngliche P1 und P2 Daten → Träger, Compiler und E₈ Struktur | Die vorhandenen algebraischen Resultate und ihre Voraussetzungen bleiben erhalten. E₆ und E₇ beschreiben gemeinsame Schnitte und Transportfenster. |
| Konkrete ursprüngliche mikroskopische Quelle → ein CAR Kanal → neutraler Achteroperator | Im vorigen Vertrag hergeleitet, einschließlich kontrolliertem iteriertem Grenzwert und nichtverschwindender Antwort. Kein neuer Beweis dieses bereits geschlossenen Schritts erforderlich. |
| Ursprüngliches Paar `(BΣ,uΣ)` → derselbe Kanal | Neuer vollständiger Satz unter den oben angegebenen Operatoridentitäten. Ihre Herleitung aus dem Originalkern bleibt offen. |
| Vollständiges ursprüngliches Quellfunktional → positive Feldrekonstruktion | Mathematische Form korrigiert. Das volle Funktional und die kompatible lokale Abbildung bleiben nachzuweisen. |
| Gemeinsame Randstruktur → vollständige wechselwirkende Raumzeitphysik | Durch diese Ergebnisse nicht abgeschlossen. Die Originalverträge verlangen dafür weitere konkrete Dynamik und kontrollierte Grenzwerte. |

Die aktive Fassung von Paper 2 nennt selbst Voraussetzungen wie Lokalität, Positivität nach Integration, einen gleichmäßigen Gap, Tightness und den zulässigen Eichquotienten. Die dort angegebene endliche `C₆` Spektralrechnung ersetzt den Nachweis dieser Voraussetzungen für eine Raumzeitquelle nicht. Außerdem wird dort der Dirac Operator aus einer bereits gegebenen Randkovarianz konstruiert. Diese Richtung darf nicht umgekehrt als Herleitung derselben Kovarianz gelesen werden.

Der ursprüngliche Vertrag `TFPT4D.LATTICE.ACTION.01` lässt eine konsistente lokale Gittertheorie mit thermodynamischem Grenzwert und kontrolliertem Verhalten bei niedriger Energie zu. Ein exakter kontinuierlicher OS Grenzwert ist für diesen Weg eine zusätzliche Verstärkung. Wir müssen die vollständige Lösung deshalb nicht unnötig an eine stärkere Bedingung binden, als das Original verlangt.

## Der jetzt entscheidende Anschluss

Der nächste Schritt ist die tatsächliche gemeinsame Wirkung von `BΣ` und `uΣ` auf einem ursprünglichen dichten Operatorbereich samt der zugehörigen Zeit und des Indexpaares. Trägt diese Identifikation, liefert der neue Satz die ganze Modenstruktur und den Rang eins ohne freie Auswahl. Trägt sie nicht, muss die ursprüngliche Quelle eine andere konkrete Quantisierung begründen.

Daneben muss das vollständige graduierte Quellfunktional denselben lokalen Feldern, derselben Reflexion und derselben Zeit zugeordnet werden. Bei einer vollen CAR Abbildung mit derselben positiven Vakuumzeit folgt der freie Quellzustand bereits aus dem vorigen Reduktionssatz. Bei einer nur geraden Abbildung muss die Zustandsverträglichkeit eigens gezeigt werden.

## Belege und Prüfung

- [Positiver Rekonstruktionssatz und expliziter Unterschied zwischen Windung und Zeitwirkung](SPECTRAL-SEAM-RECONSTRUCTION.txt)
- [Alternative direkt mit dem ursprünglichen Calderón Projektor](TOEPLITZ-SEAM-RECONSTRUCTION.txt)
- [Korrektur der vollständigen Korrelationen und der Feldabbildung](OS-SOURCE-REPAIR.txt)
- [Unabhängige Prüfung des Rekonstruktionssatzes](REVIEW.txt)
- [Exakte Kontrollrechnung](check_first_identifications.py) und [Ergebnis](identification-results.json)
- [Voriger Vertrag: mikroskopische Quelle und Achterantwort](../paper-source-reconstruction-20261025/README.md)

Die Kontrollrechnung besteht normal und mit `python -OO`. Die Beweise gelten auf unendlichen Hilberträumen; die exakten Kontrollbeispiele sind ihre zusätzliche Vorzeichen und Gegenfallprüfung. Sie prüfen keine bisher fehlende ursprüngliche Quelle. Der Forschungsstatus und die offenen Voraussetzungen werden im gemeinsamen Wissensatlas dokumentiert. Es erfolgt keine Hochstufung des ursprünglichen Ledgers.
