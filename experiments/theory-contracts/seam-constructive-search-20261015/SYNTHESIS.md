# Was die konstruktive Suche jetzt entscheidet

Stand: 15. Oktober 2026. Forschungsziel: dieselbe ursprüngliche Quelle soll
Markierung, Zustand und Zeit liefern. Das ist weiterhin nicht bewiesen.
Die folgende Fortsetzung schließt konkrete algebraische Fragen und verhindert,
dass ihre Lösungen als Herleitung einer unveränderten physischen Quelle zählen.

## 1. Der tatsächliche Suchraum ist abgeschlossen

Der Lauf `run-20261015T121548Z` hat auf beiden H200 zusammen alle 262.144
kanonischen Phasenbelegungen genau einmal klassifiziert. Vier Belegungen erfüllen
die festen Constraints; alle vier wurden automatisch vom lokalen Lean-Kernel
geprüft. Achtkanalquelle, Halbsektor und quadratisches Zeitgesetz waren Eingaben.
Die vier Klassen waren analytisch erwartet; neu ist hier die vollständige
ausführbare Klassifikation mit Zertifizierung, nicht eine physische Quellenauswahl.
Der unveränderte abgeschlossene Raum wird beim nächsten Start wiederverwendet.

Die Fortsetzung arbeitet mit den Originaloperatoren. Auf der globalen Erweiterung
identifiziert die bereits vorhandene Quellenparität die beiden Zeichen durch
`X_plus = Pi X_minus`, bei unverändertem Q, K und Viertelmarker. Der kanonische
Zustand wird ebenfalls erhalten. Damit bleiben zwei Klassen mit fest markiertem
P-Lift. Eine zusätzlich physisch festgelegte Operation X darf dabei nicht
umdefiniert werden; ein lokaler Feldintertwiner wurde damit nicht bewiesen.

## 2. Der Einzelwindungsschritt und sein tatsächlicher Auswahlgehalt

Der volle ursprüngliche Schritt `T1(p)=Cp+e0` erfüllt für alle Ladungen
`K_alpha(T1 p)-K_alpha(p)=(Q(p)+alpha)/8`. Die vierte Potenz liefert genau
die bereits bekannte W4-Wirkung. Universelle Cocycle- und Wardidentitäten sind
zusätzlich in Lean geprüft.

Die diskrete Erweiterung allein lässt beide P-Lifts zu. Eine vermeintliche
Auswahl aus `C X C^-1 = nu X` wäre ohne Begründung dieser skalaren Kovarianz falsch:
die explizite zulässige Gegenkonstruktion hat
`C X C^-1=(i X+X D)/sqrt(2)`. Ihr Gradoperator ist aus demselben Wort rekonstruierbar:
`D=sqrt(2) X^-1 C X C^-1-i I`. Er ist kein frei hinzugefügter Speicherzustand.
Der Torus wird dabei auf die Halbsektoren gehoben; seine alte Integerperiodizität
darf nicht unbemerkt auf die gesamte Erweiterung übertragen werden.

Innerhalb der bereits stipulierten positiven E8-Stromdarstellung gibt es hingegen
eine konkrete Auswahl: X kommutiert mit den SU(8)-Wurzelströmen, weil ihre
Ladungen zu s orthogonal sind und Q=0 besitzen. Der Zyklus wird durch diese
Stromgruppe und einen Torusweg realisiert. Daraus folgt `nu=+i` oder `-i`, also
P-Liftzeichen +1. Die Verbindung gilt auf dem gegebenen E8-Modul. Die zusätzliche
Einbettung sämtlicher ungerader CAR-Felder und die Herkunft dieses Moduls aus
P1/P2 sind damit nicht abgeleitet. Die ältere Lorentz-Gitterquelle der Signatur
(9,1) wird nicht mit diesem positiven Achtgitter gleichgesetzt.

## 3. Positive Zustandsfortsetzung ist geschlossen, Zustandsauswahl nicht

Für jedes algebraische endliche Wort `y=a0+a1 X` gilt
`E0(y* y)=a0* a0+phi^-1(a1* a1)`.
Ein positiver ursprünglicher Funktional omega0 liefert deshalb die positive
Fortsetzung `omega0 o E0`. Unter dualer Z2-Invarianz ist diese Fortsetzung auf
der Verschiebungsalgebra eindeutig. Lean prüft die allgemeinen Algebra-,
Stern-, Positivitäts- und Eindeutigkeitssätze; die konkrete Übertragung auf
endliche komplexe Gittersummen ist in `SOURCE_BRIDGE.md` analytisch bewiesen.
Eine vollständige konkrete Finsupp-Instanziierung im Lean-Kernel ist nicht behauptet.

Diese Eindeutigkeit verliert entscheidende Information: Das festgehaltene
K besitzt die beiden Grundvektoren `0` und `r0=(1/2)^4(-1/2)^4`.
Das war bereits am 14. Oktober bewiesen. Alle diagonalen Mischungen dieser
beiden Zustände liefern auf der bloßen Verschiebungsalgebra denselben Funktional.
Die Fortsetzung klärt jetzt die unterschiedlichen Auswahlannahmen:

| Zusatzbedingung im gegebenen Grundraum | Ergebnis |
|---|---|
| Energie null und duale Z2-Invarianz | ganze Familie diag(a,1-a) |
| zusätzlich Invarianz unter dem vorhandenen Austausch Rc | I2/2 |
| originale Cartan-Vakuumwerte | eindeutig der ursprüngliche Vektor 0 |
| originale Cartan-Vakuumwerte und Rc-Invarianz zugleich | unvereinbar |

Die Symmetrie eines Hamiltonoperators erzwingt keine Invarianz seines Zustands.
Auch der zweite K-Grundvektor ist kein zweites VOA-Vakuum: Er hat konformes
Gewicht eins. Cartan-Nullmoden, lokale Felder und Verschiebungswörter dürfen
deshalb nicht als dieselbe Beobachtungsalgebra behandelt werden.

## 4. Was daraus für das ursprüngliche Ziel folgt

Der vorhandene Quellenvertrag `quotient-locality-e8-20261014` beweist bereits:
Die unveränderte geprüfte CAR/QWZ-Quelle hat einen einzigen Grundzustand und
eine positive Lücke; die stipulierte E8-Erweiterung mit derselben Zeit hat zwei.
Außerdem gibt es auf der alten Hilbertraumdarstellung keinen nichttrivialen
Operator mit den erforderlichen halbintegralen Cartanladungen. Alle vier
Phasentreffer teilen diese Eigenschaften. Keiner ist daher eine äquivalente
Realisierung der unveränderten Quelle mit erhaltenem Hilbertraum, Cartan und K.
Eine positive Zustandserweiterung ist möglich; sie fügt Sektoren hinzu und
beweist deren physische Herkunft nicht. Eine Wahl des Startzustands entfernt
den zusätzlichen Nullvektor nicht.

Der nächste konstruktive Kandidat muss zuerst eine quellenbegründete
Sektor-/Hilbertraumwirkung mit dem neutralen Partner liefern, oder ein anderes
aus der Quelle hergeleitetes gemeinsames Zeit-/Operatorwörterbuch. Erst diese
geänderte, explizit begründete Voraussetzung öffnet einen neuen Suchraum.
Ein eingesetzter Sektorenergieoffset, ein weiterer Phasenscan oder das Weglassen
der Cartanobservablen erfüllt diese Bedingung nicht. Dieser Ausschluss betrifft
die benannte Quellrealisierung, nicht TFPT oder alle P1/P2-Realisierungen.

## Belege und Wiederverwendung

- Produktionsdaten: [results.json](results.json), unveränderte Laufdaten unter `runs/`.
- Äquivalenz und expliziter Gradlift: [branch_equivalence/PROOF.md](continuation/branch_equivalence/PROOF.md).
- Volle Windung: [T1/PROOF.md](continuation/T1/PROOF.md), [Lean-Beleg](continuation/t1_certificate/receipt.json).
- Stromgruppenbrücke: [current_group/PROOF.md](continuation/current_group/PROOF.md).
- Positive Fortsetzung: [SOURCE_BRIDGE.md](continuation/state_expectation/SOURCE_BRIDGE.md), [Lean-Beleg](continuation/state_expectation/receipt.json).
- Zustandsannahmen und Originalausschluss: [ground_state/PROOF.md](continuation/ground_state/PROOF.md).
- Primäre Vertexformel und Stromautomorphismen: Dong–Nagatomo, [Abschnitte 2.2–2.4](https://arxiv.org/pdf/math/9808088).

Die kleinen Fortsetzungsfragen wurden exakt auf der CPU und mit Lean geprüft.
Ein erneuter H200-Lauf derselben Phasenklasse wäre keine neue Suche.
