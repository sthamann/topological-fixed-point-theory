# Prüfung der eingefügten TFPT-Zusammenstellung

14. Oktober 2026. Geprüfte Eingabe: `Eingefügter Text.txt`, vom Nutzer in dieser Sitzung bereitgestellt.

**Urteil: Eine Mischung aus richtigen endlichen Resultaten, bedingten Theoremen, konkreten Formelfehlern und unbelegten Verknüpfungen. Die Überschrift „sämtliche geschlossenen Lösungen“ und die behauptete Gesamtvereinigung sind nicht durch den Inhalt gedeckt.**

Die Prüfung unterscheidet direkte Nachrechnung, Originalquellenvergleich und nicht nachgewiesene Behauptungen. Ein im Text nicht mit Quelle versehener großer Beweis wird nicht allein deshalb als falsch bezeichnet. Ebenso wird ein numerisch bestätigter Wert nicht zum Herkunftsbeweis seiner Eingaben.

## Direkt korrigierte Formeln

### Elektromagnetischer Fixpunkt: alpha und sein Kehrwert verwechselt

Die eingefügte Gleichung setzt q(alpha)=delta_top exp(-2/alpha). Die Originalfunktion `verification/v3_em_alpha.py`, Zeilen 17–28, verwendet dagegen

\[
q(\alpha)=48c_3^4 e^{-2\alpha}.
\]

Beide Funktionen wurden mit denselben übrigen Angaben und 60 Dezimalstellen Rechengenauigkeit gelöst:

| Funktion | Erhaltenes alpha^(-1) |
|---|---:|
| Eingefügter Text: exp(-2/alpha) | 137,001462266837407… |
| Original v3: exp(-2 alpha) | 137,035999216840712… |

Der im Text genannte Wert gehört zur Originalfunktion und **nicht** zur dort abgedruckten Formel. exp(-2/alpha) kommt an anderer Stelle im kosmologischen Readout vor; diese Funktionen dürfen nicht verwechselt werden. Die separate analytische Eindeutigkeitsaussage des vorhandenen Alpha-Contracts betrifft ihre festgelegte Originalfunktion.

### Helle Zweierblöcke von v1024: Besetzungszahl fehlt

Der Originaloperator wurde unverändert geladen; seine 26 Prüfungen im Teil `parent_internal_intertwiner` bestanden. Für k doppelt und s einfach besetzte Familien ergibt die erhaltene Ladungszerlegung den hellen Zweierblock

\[
\begin{pmatrix}kB_0+s/4&-c\sqrt{k}\\-c\sqrt{k}&kB_0+s/4+\delta\end{pmatrix},
\]

also

\[
E_{k,s,\pm}=kB_0+\frac{s}{4}+\frac{\delta\pm\sqrt{\delta^2+4kc^2}}2,
\quad k\ge1,\quad k+s\le3.
\]

Das pauschale +1/4 im eingefügten Text gilt nur für s=1. Im tatsächlich vorhandenen Vakuum-/Paarblock k=1,s=0 sind die Energien 0,5434012371… und 2,3432654295… . Die abgedruckte Formel liefert stattdessen 0,7934012371… und 2,5932654295… . Sie ist als allgemeine Spektralformel falsch. Die angegebene Übergangsamplitude 1600/72889 ist dagegen mit dem ursprünglichen k=1,s=0-Block verträglich; der fehlerhafte gemeinsame Energieoffset fällt aus diesem Übergang heraus.

### MaxEnt-Identität: gleicher reduzierter Zustand als Voraussetzung

Aus tau=E_A(rho) allein folgt nicht S(tau)-S(sigma)=D(sigma||tau). Dafür muss auch E_A(sigma)=tau gelten (mit den üblichen Trägerbedingungen). Ein direkter Gegenfall zu der unqualifizierten Form ist E_A=id, rho=tau=diag(3/4,1/4), sigma=diag(1/4,3/4). Dann ist die linke Seite null, die rechte (log 3)/2>0. Die korrekte MaxEnt-Aussage gilt auf der gemeinsamen Zwangsmenge. Ebenso darf die Takesaki/Petz-Zerlegung nicht auf eine beliebige Kompression übertragen werden; die zustandserhaltende bedingte Erwartung und ihre passende Zustandsabbildung gehören zum Satz.

### Der behauptete Fünfer-Intertwiner scheitert an der Operatorwirkung

Dies ist ein direkter algebraischer Widerspruch im Text, nicht lediglich ein fehlender Quellenbeleg. Auf dem angegebenen Quartikcode hat der genannte Vorzeichenoperator die Eigenwertaufteilung 4+1; auf den RR-Monomen hat die angegebene Involution 3+2. Eine Isometrie T zwischen beiden fünfdimensionalen Räumen ist unitär. Aus T C_v=P_v T würde deshalb dieselbe Spur folgen. Tatsächlich sind die Spuren 3 und 1. Die charakteristischen Polynome (x-1)^4(x+1) und (x-1)^3(x+1)^2 sind verschieden. Genau dieser behauptete gemeinsame Operator kann daher nicht durch T identifiziert werden. Auch seine Wiederverwendung als T1-Schritt in Abschnitt 8 ist unzulässig. Andere Aktionen oder andere Abbildungen werden dadurch nicht ausgeschlossen.

### Der allgemeine Lindblad-Ausschluss trägt — mit einem neuen Beweis

Die im Text verwendeten Choi-Ränge 6 für Phi und 16 für Phi^2 wurden an den tatsächlichen Cliffordwörtern bestätigt. Daraus folgt mehr als der ältere orbit-skalare Contract: Kein endlichdimensionaler zeitlich homogener GKLS-Generator auf demselben M4 kann Phi als positiven Zeitschritt erzeugen.

Begründung: Bei einer endlichen homogenen Quantendynamik hat die Entwicklung jedes reinen Zustands für alle t>0 dieselbe Trägerdimension. Das folgt aus der analytischen Trajektoriendarstellung: Der Träger ist e^(tG) S(psi), wobei S(psi) durch Sprungoperatoren und ihre iterierten Kommutatoren mit G aufgespannt wird, unabhängig von t; e^(tG) ist invertierbar. Wendet man dies auf den maximal verschränkten Zustand mit einer untätigen Ancilla an, ist dessen Entwicklung der Choioperator. Daher bleibt dessen Rang für t>0 konstant. Phi und Phi^2 verletzen diese notwendige Bedingung. Die zugrunde liegende Trägeraussage steht in [Fagnola–Girotti, Lemma C.2](https://link.springer.com/article/10.1007/s11785-026-01976-8#app-3).

Dies ist ein hier nachvollzogener stärkerer Anschluss, kein rückwirkend dem alten Contract zugeschriebener Beweis. Er verbietet weder einen größeren gemeinsamen Elternprozess mit Umwelt noch zeitabhängige Generatoren oder diskrete Schritte. Die zusätzlichen Rangzahlen für die Fünfer- und Zwanzigerkanäle sind in diesem Audit nicht originalbestätigt.

### Die Lüders-Folge braucht den Kanal zwischen den Messungen

Der exakte Operatorcheck ergibt

\[
\sqrt{E_i}\,\Phi(E_j)\sqrt{E_i}=M_{ij}E_i,\qquad
M=\frac1{36}\begin{pmatrix}19&7&10\\7&19&10\\10&10&16\end{pmatrix}.
\]

Damit entsteht unter dem spezifizierten Instrument tatsächlich eine autonome Markovfolge. Die Formulierung muss jedoch „Lüders-Messung, Phi-Schritt, nächste Lüders-Messung“ lauten. In M=(B+J3)/2 ist J3 die Matrix mit allen Einträgen 1/3. Ohne den dazwischenliegenden Phi-Schritt ergibt die gleiche Messung stattdessen Diagonaleinträge 2/3 und Nebendiagonaleinträge 1/6. Die ursprüngliche Quelle wählt das Instrument dadurch noch nicht aus.

### Thermische Präparation: den Sektor ausdrücklich angeben

Für den im Text angegebenen vollen Paar-Hamiltonoperator sind die Energien 0, 8g, 4g auf den Rängen 1,5,10. Bei beta g=log(7)/8 gilt daher

\[
\rho_\beta=\frac{P_1+P_5/7+P_{10}/\sqrt7}{1+5/7+10/\sqrt7}.
\]

Er hat Rang 16. Der Choi-Zustand (7P1+P5)/12 entsteht exakt erst im eingeschränkten Raum W6=P1\oplus P5. Der Text erwähnt W6; korrekt ist daher die ausdrücklich auf diesen Raum beschränkte Gibbs-Aussage. Aus dem vollen Hamiltonoperator allein folgt keine Auslöschung des Zehnersektors. Auch der ursprüngliche Paar-/Pulscontract arbeitet ausdrücklich unter dieser Sektorauswahl.

## Was an den acht Abschnitten trägt — und was daraus nicht folgt

| Abschnitt | Belastbarer Inhalt | Korrektur des behaupteten Umfangs |
|---|---|---|
| 1: Naht, vier Marken, RR und Anker | Die Aussagen über den vorgegebenen elliptischen Torus, die vier Halbperioden, H1 der vierpunktierten Sphäre, h0(O(4))=5 und die Ankerarithmetik sind konkrete mathematische Identitäten. Der Index-4-Glue D5+A3 nach E8 und der Halbspinor sind vorhandene Konstruktionen. | P1/P2 sind hier Prämissen. H1-Rang drei ist nicht ohne Feldabbildung die physische Familienzahl. RR-Dimension fünf wählt nicht automatisch fünf physische Träger. Die Quadratklasse tau=i und die Spinstruktur benötigen ihre Auswahlbedingungen. Der SM-Darstellungsinhalt ist nicht schon eine wechselwirkende chirale 4D-Theorie. |
| 2: Fünferräume und Quartikcode | Der explizite Code und sein Ein-Register-Löschschutz sind sinnvolle exakte endliche Aussagen. Clockordnungen, CRT-Zerlegung und die Theta-/Hecke-Identitäten tragen unter den angegebenen Definitionen. | Der behauptete Intertwiner ist unter den angegebenen Operatoren unmöglich: 4+1 und 3+2 haben verschiedene Spuren. Andere Verbindungen benötigen eine andere explizite Wirkung. Der Wortlaut zur goldenen Projektion muss sich auf die 240 Wurzeln und ihre projizierten Schalen beziehen, nicht auf das gesamte unendliche Gitter. Arithmetische Identitäten liefern keinen physikalischen Zeitgenerator. |
| 3: Flavor und Alpha | phi0=0,0531719521768455… und lambda=0,2243762368847218… wurden direkt nachgerechnet; die Residuenmatrix hat Determinante acht. Die Original-Alpha-Funktion liefert den genannten Root. | Die Alpha-Formel im Text ist vertauscht. Die weiteren Zahlen gehören zu verschiedenen Readout-/Transferverträgen; insbesondere N_* ist eine deklarierte Eingabe. Die Angabe „81+108 Fälle“ ersetzt den vollständigen Stabilitätsbeweis einschließlich aller Unterbündel, Grade und Fahnen nicht. Die konkrete große Bündelprüfung ist durch den eingefügten Text allein nicht reproduzierbar belegt. |
| 4: MSC und geladene Naht | D=log((1-C)C^-1) und C=(1+exp D)^-1 sind für 0<C<I exakt; sie parametrisieren denselben quasifreien Zustand. Die v1024-Paarkorrelation wurde am Originalmodell reproduziert: kappa(1)=-0,0028834107044733… und kappa'(0)=-1/300. Der angegebene Spinorvektor hat Normquadrat vier und die behauptete ganzzahlige Glueverschiebung. | Eine Kovarianz-Inversion wählt C und die Quelle nicht aus. Quasifreie Rekonstruktion und nichtgaußsches Parent sind verschiedene Verträge. Der positive geladene Randanschluss setzt im Original u.a. acht CAR-Kopien und das gerade Randnetz voraus; Breite acht ist nicht gleich acht abgeleitete chirale Arten. Zusätzliche Linkqubits sind ausdrücklich zusätzliche Struktur. Neutralität eines Gewichtes beweist keine Kondensation. |
| 5: Quell- und Prozessmodelle | Der markierte CPTP-Kanal, die explizite Lüders–Phi-Folge und bedingte Puls-/Parent-Modelle tragen. Der Choi-Rangsprung schließt sogar allgemeine homogene GKLS-Generatoren auf demselben M4 aus. Die Paarkorrelation und der konkrete Originalblock wurden neu nachgerechnet. | „Vollständig gelöst“ kann jeweils das spezifizierte endliche Modell meinen. Es ist kein Beleg, dass dessen Kopplungen, Präparation, Pulszeit, Reset oder Operatordictionary aus P1/P2 stammen. Die helle Spektralformel benötigt s. Ein endlicher Gibbszustand des vollen 16D-Paars hat vollen Rang; der rank-6-Choi-Zustand setzt die Beschränkung auf den entsprechenden 6D-Sektor voraus. |
| 6: Zukunftsverdichtung und Kopplungen | Ker(R)-Invarianz ist die richtige lineare Bedingung für eine wohldefinierte reduzierte Dynamik auf im(R). Die Dimensionen 136 und 325 sind mit den angegebenen Jordan-Zerlegungen konsistent. Eine vollständige adaptive Recordkanal-Fehlerschranke kann teleskopisch wachsen. | Die erlaubten Operationen und Effekte müssen vorher festliegen. Erlaubt man weitere Kopplungen, kann ein zuvor ausreichender Zustandsraum unzureichend werden. Die Aussage ist keine allgemeine Garantie für Hylæan-/KI-Lernen. Die Entropieformel braucht die oben fehlende Zwangsmenge. |
| 7: Gravitation und Spiegelgap | Die homogene Skalarformel löst das gewählte masselose FRW-Modell. Die Hopping-Ungleichung ist unter kontraktiven Links eine exakte Quadratsummenidentität; 3-6*(3/10)=6/5 ist richtig. | Die Eingaben M_P, räumliche Geometrie, Delta, t und Gauss-Sektor-Restschranke werden dadurch nicht hergeleitet. Freier Kubo-Anschluss, klassische ADM-Algebra und eine positive modellinterne Spiegelenergie sind keine gemeinsame nichtlineare chirale 4D-Quantentheorie. |
| 8: Gemeinsames Gesamtmodell | Die Forderung „es existiert ein M, das alle T_i gleichzeitig erfüllt“ ist genau die richtige logische Form. Halbseitige modulare Inklusionen sind eine ernsthafte bedingte Rekonstruktionsroute. | Der Text führt das erforderliche M nicht vor. Er setzt getrennte Zustände, Hilberträume, Generatoren und Zusatzgrade durch Formulierungen wie „derselbe Elternprozess“ gleich. Das ist die fehlende Konstruktion selbst, kein Beweis. Weder Polmassen noch Spinorkondensation noch die angegebene Defektkompensation folgen aus den danebenstehenden Teilrechnungen. |

Diese Tabelle bestätigt keine Quellen, die hier nicht geprüft wurden. Insbesondere die Zahlen zu minimalen Zukunftsfiltern, zur detaillierten Bündelzählung und zu allen Echo-/Doppelantwortverträgen werden ohne ihren vollständigen jeweiligen Operations- und Quellvertrag nicht pauschal als universelle TFPT-Aussagen übernommen.

## Warum Abschnitt 8 kein fertiger Abschluss ist

Die Aussage „für jede Teilfrage gibt es ein Modell mit einem passenden Resultat“ hat die Form forall i exists M_i. Der benötigte Abschluss lautet exists M forall i. Der Tausch der Quantoren ist unzulässig. Man benötigt konkrete gemeinsame Algebra-, Zustands-, Feld-, Zeit- und Kopplungsabbildungen, die alle betreffenden Gleichungen gleichzeitig erhalten.

Der direkte Originalvergleich liefert drei besonders klare Brüche:

1. `EM.WARD.01` führt die physischen Ward-/Quillen-Ursprungsannahmen ausdrücklich separat von den Koeffizientenidentitäten. Bündelstabilität liefert deshalb nicht automatisch einen Yukawaoperator oder Polmassen aus demselben W[J].
2. `v1035_matter_coupling.py` behandelt eine bedingte Kopplung bis erster Ordnung und benennt Hindernisse des entfernten homogenen Gravitationsblocks. Die Behauptung, die angegebene homogene Lösung kompensiere den vierteiligen Gitterdefekt, enthält keine entsprechende gemeinsame Rechnung.
3. Der Full-Process-Contract `tfpt-full-process-20261002` verlangt dieselbe Modellwahl und denselben Zustand für alle Gates. Sein geprüfter Theta0-Zweig scheitert; Theta1 besitzt dort kein Gesamturteil. Dieser dokumentierte Stand darf nicht mit Resultaten aus anderen Parameterwelten als erledigt überschrieben werden.

Bei der modularen Route ist außerdem „strikt positiver Generator“ zu stark: Der gemeinsame invariante Referenzvektor liegt im Kern des Translationsgenerators. Die relevante Aussage ist P>=0. Nichttrivialität P!=0 muss zusätzlich sichergestellt sein; etwa N=M ergibt keine nichttriviale Translation. Mehrere Inklusionen ergeben ohne die passenden gemeinsamen Netz-, Kommutations-, Lokalisierungs- und Spektralbedingungen keine 3+1D-Poincaré-Theorie.

## Konsequenz für die laufende Herleitung

Die brauchbare Verbindung ist die **gemeinsame geordnete Quellenwirkung**, einschließlich der ursprünglichen komplexen Phasen und der nicht wegkomprimierbaren Korrelationen. Der Atlas enthält im g=0-Ladungs-/Skalarmodell bereits eine solche regulierte Geschichte samt Determinantengedächtnis. Seine Originalchecker bestanden erneut mit 39 und 25 Prüfgruppen. Diese positiven Ergebnisse dürfen nicht hinter einem leeren Atlas-Summary verschwinden.

Der parallel erarbeitete Beweis `PROOF.md` zeigt zugleich, warum man diese Quellenwirkung nicht aus einem endlichen Mittelkanal zurückgewinnen kann: Sogar stationäre, kovariante, ergodische Prozesse mit **allen** gleichen Phi^n und demselben Fixzustand können verschiedene Vier-Zeit-Antworten besitzen. Den fehlenden physikalischen Zusammenhang allein durch einen zusätzlich gewählten Puls, Reset oder eine Gibbspräparation zu ersetzen, würde die ursprüngliche Aufgabe nicht lösen.

Die eingefügte Zusammenstellung ist deshalb als Sammlung möglicher Bausteine brauchbar, nach Korrektur ihrer Formeln und Geltungsbereiche. Ihr Abschnitt 8 ist ein Konstruktionsprogramm. Er ist kein Nachweis, dass die gemeinsame ursprüngliche Quelle bereits konstruiert wurde.

## Entscheidende Originalquellen und Reproduktion

- `verification/v3_em_alpha.py`, Zeilen 17–42: beide Originalimplementierungen exp(-2 alpha).
- `verification/v1024_parent_internal_boundaries.py`, Zeilen 239–318: deklarierte rationale Modellkopplungen, tatsächlicher Hamiltonoperator und exakter heller Vakuum-/Paarblock.
- `verification/status_ledger.csv`, EM.WARD.01: Identität und physischer Ursprung getrennt.
- `experiments/theory-contracts/rr-continuous-clock-20260921/QUELLENPRUEFUNG.md`: echte RR/W-Kovarianz, aber keine automatische Gleichsetzung mit Hamiltonzeit.
- `experiments/theory-contracts/source-localized-spinor-20260919/PROOF.txt`: positive Randwortkonstruktion unter acht Kopien und expliziten Quellprämissen.
- `experiments/theory-contracts/tfpt-full-process-20261002/contract_index.json`: gemeinsames Modell als tatsächlicher Vertrag.
- `pasted_formula_check.py` und `pasted-formula-result.json`: reproduzierbare direkte Prüfungen dieser Sitzung.

Keine Änderung der ursprünglichen Verifikationsmodule, keine Ledger- oder Paper-Promotion. Keine vollständige TFPT-Lösung behauptet.


## Ergänzung zur Quellenpriorität, 14. Oktober 2026

Der Vergleich mit dem engeren Vertrag vom 21. September war korrekt, die Einordnung als neue Verschärfung im Projekt war unvollständig: `source-continuous-marked-time-search-20260922` beweist bereits den Ausschluss sämtlicher zeitlich homogener GKLS-Generatoren auf demselben M4 für alle positiven ganzzahligen Potenzen des angegebenen Kanals. Der hier geprüfte Choi-Rangsprung ist ein alternativer Beweis für Phi, keine erstmalige Schließung dieser Projektfrage. Sein mathematischer Geltungsbereich bleibt unverändert. Größere Umgebungen, GNS-unitäre Zeit und zeitabhängige Generatoren sind dadurch nicht ausgeschlossen.
