# Positive Quellenkerne automatisch konstruieren

## Adaptive Auswahl bereits zertifizierter Antworten

`adaptive_discriminator.py runs/source-clustering --out /tmp/meine-auswahl.json`
prüft zuerst die vorhandenen exakten Endpunkte und das ganze gespeicherte
Segment. Es ordnet anschließend die unterscheidenden Wortfragen nach der Zahl
neu getrennter Kandidatenpaare pro ausdrücklich gewählter Kostenangabe;
Vorgabe ist die Wortlänge. Auf den geprüften Kandidaten affine Duplikate
werden gemeinsam aufgeführt; der exakte affine Rang berücksichtigt auch
lineare Abhängigkeiten zwischen verschiedenen Gruppen. Für das ursprüngliche
Segment liefern zehn Fragen ein einziges Signal. Die Wertung ist eine lokale
Heuristik und keine Wahrscheinlichkeit oder Aussage über alle Prozesse.

Optional nimmt `--hypotheses FILE` eine JSON Liste mit den Feldern
`hypothesis`, `source_ids`, `premises`, `additional_assumptions` und `tests`
entgegen. Diese Vorschläge werden ausschließlich strukturell geprüft und als
unverifiziert ausgegeben. Sie verändern keine Suchbedingungen; Quellen und
Prämissen müssen separat gegen Atlas und Originale geprüft werden. Die
Identifikation dieser optionalen endlichen Klasse mit dem ursprünglichen P1
Quellenaufruf ist damit nicht bewiesen. Die vorhandene gemeinsame Originalkette
Randkern → Carrier → S/Z → OS/CAR → Gamma bleibt erhalten; diese Suchklasse ist
keine notwendige globale TFPT Ursprungsvoraussetzung. Siehe Bestandskorrektur
L20261102C und R101 im kanonischen Forschungsplan.

## Lineare Beziehungen zwischen Momenten

Der optionalen `--source-control` Datei kann `linear_relations` hinzugefügt
werden. Jede Gleichung beschreibt eine komplexe Summe von Wortmomenten mit
gaußschen rationalen Koeffizienten. `-1` steht für den Transfer T; andere
Zahlen bezeichnen die Majorana Generatoren. Die folgende Gleichung ist
**ausschließlich eine synthetische Softwarekontrolle**, keine TFPT Quellrelation:

```json
{
  "linear_relations": [{
    "terms": [
      {"word": [0, -1, 0], "re": "1"},
      {"word": [1, -1, 1], "re": "1"}
    ],
    "re": "1", "im": "0"
  }],
  "source_observables": {"label": "synthetische Kontrolle"}
}
```

Jeder Term hat `word`, rationales `re` und optional rationales `im`;
`re` und `im` am Ende der Gleichung bilden die rechte Seite. Ganze Zahlen
und Bruchstrings werden akzeptiert, Fließkommazahlen und ausführbare Ausdrücke
nicht. CAR Reduktion bewahrt das Vorzeichen und fasst gleiche Wörter zusammen.
Die normierten Gleichungen und die Herkunftsmetadaten erscheinen in der
Kampagnenkonfiguration sowie in der Momentdatei. Die SDP Suche rechnet weiterhin
numerisch. Beide endlichen Vervollständigungen prüfen die Gleichungen danach
zusätzlich exakt auf den rational rekonstruierten Momenten; ein fehlgeschlagener
Test bleibt unentschieden. Der Offlineprüfer prüft sie ebenfalls exakt.
Alte Eingaben ohne `linear_relations` bleiben zulässig.

Die Software liefert keine ursprüngliche Quellenrelation und keine vollständige
E8, Higgs oder Geometrieabbildung. Die jeweils benannte physikalische Bedeutung
und Zulässigkeit einer eingegebenen Beziehung müssen unabhängig begründet werden.

Aktualisierung 16. Oktober: Der Suchpfad verwendet jetzt standardmäßig reine
Zulässigkeit (`--objective feasibility`) mit Floor und Cluster-Gap 0. Die
nachstehenden Ergebnisse vom 15. Oktober sind historische Extremwertkontrollen;
`diversity` und min/max-Spur wählen keine physische Quelle aus. Die Transferklasse
ist ausdrücklich `vacuum_effect`, `thermal_gns` oder `bounded_positive`; ein
thermischer GNS-Transfer darf größer als eins sein. Die neue positive Kragenkontrolle
und die Mehrfeld-Herleitung stehen im benachbarten Contract
`seam-collar-multifield-20261016/PROOF.txt`. Historische Ausgabedateien bleiben erhalten.

Für unnormierte Quellenamplituden ist außerdem `T Omega = lambda Omega` mit
`--reference-eigenvalue` explizit möglich. Standard bleibt lambda=1. Zum Beispiel
akzeptiert `--transfer-kind bounded_positive --reference-eigenvalue 0.0625
--transfer-ceiling 0.0625` eine passende skalierte Quellenkontrolle. Lambda ist
hier eine Eingabe, kein gefundener Vakuumfaktor. `omega(1)=1` normiert den
bedingten Zustand und bedeutet nicht `Z[0]=1`. Eine Fock-Vervollständigung für
lambda ungleich1 wird als unentschieden ausgewiesen. Auch der direkte
Vervollständigungsaufruf prüft Klasse, Referenzamplitude und Deckel, bevor er ein
exaktes Zertifikat erzeugt. Cacheidentität und
Resultate behalten die neue Normierungsannahme.

Stand: 15. Oktober 2026. Forschungs-ID `UR.SEAM.POSITIVE_KERNEL_SEARCH.01`.

Die Suche konstruiert einen unbekannten gemeinsamen Wortkern: aus ihm werden
Zustand und Wirkung eines Transferoperators gemeinsam realisiert. Ein fertiges H,
eine Dichtematrix oder ein Zielkern werden dem Optimierer nicht vorgegeben.
Die endliche CAR-Algebra, deklarierte Referenzeigenwirkung und interne Vierteldrehung sind
deklarierte Suchannahmen. Sie wurden hier nicht aus P1/P2 hergeleitet.

## Was wirklich ausgeführt wurde

`runs/certified-positive/summary.json` enthält den vollständigen Fünf-Moden-Lauf.
Er sucht 2245 komplexe Wortmomente bis Grad 4 und vervollständigt sie auf zwei
unterschiedliche vollständige, positive Kerne in Dimension 32. Beide geben alle
2245 rational rekonstruierten Momente exakt wieder, darunter 508 mit mehreren
Transferfaktoren. Das sind keine bloß numerisch positiven Tabellen.

Die beiden Vervollständigungen stimmen für jedes Wort bis Länge 4 überein.
Die automatisch gelieferte kürzeste Unterscheidung hat Länge 5:

    omega(gamma_0 gamma_2 T gamma_2 gamma_0) = 1/4  beziehungsweise 1.

Dies ist die Transferantwort nach einer zusammengesetzten Zweifeld-Anregung.
Die Übereinstimmung aller kürzeren Antworten wird durch vollständigen Replay
geprüft. Zehn solche Zweifeld-Proben unterscheiden dieselben beiden Treffer.
Es geht um diese expliziten Kerne, nicht um eine Klassifikation aller Quellen.

Der Lauf verwendet die Suchschranke T >= (1/4) I, damit ein endlicher Logarithmus
existiert. Diese 1/4 ist eine frei gewählte numerische Schranke und hat **keinen
hergeleiteten Zusammenhang mit der ursprünglichen Viertelholonomie**. Der
separate Lauf ohne positive Untergrenze liefert ebenfalls vollständige Kerne,
aber mit Nullmoden; für diese wird kein endliches H=-log(T) behauptet.

Die Besetzungen werden aus den gefundenen Momenten gelesen. Die lineare
Suchzielfunktion bevorzugt in diesem Lauf einen reinen Besetzungszustand; diese
Auswahl ist keine physikalische Zustandsregel. Die beiden anschließenden
Optimierungsziele min/max tr(T) erkunden unbestimmte Richtungen. Keines
der Ziele wird als ursprüngliches Wirkungsprinzip ausgegeben.

## Nachprüfung einer tatsächlichen Anforderung der Originalquelle

Die operational seed category verlangt in der Originalquelle zusätzlich
Clustering. Beim ersten positiven Lauf besitzt der Maximalspur-Treffer einen
27-dimensionalen stationären Raum. Er erfüllt damit das endliche Gegenstück
dieser Anforderung nicht. Das Programm weist dies jetzt automatisch aus:

    finite_clustering <=> rank(I-T) = dimension-1.

Für einen positiven selbstadjungierten Transfer mit T Omega=Omega bedeutet dies
exakt T^n -> |Omega><Omega|. Der zusätzliche Lauf `runs/source-clustering`
verlangt dieses Verhalten konstruktiv über

    T <= (1-gap) I + gap |Omega><Omega|,  gap=1/4.

Die positive numerische Gap ist eine ausgewiesene Suchschranke, keine aus der
Quelle berechnete physikalische Skala. Sie liefert für alle endlichen A,B eine
Abschätzung der verbundenen Antwort durch
(1-gap)^n ||A*Omega|| ||B Omega||. Beide neuen Kerne besitzen nur noch Omega
als stationäre Richtung; sie unterscheiden sich in der Zweifeld-Probe durch
1/4 gegenüber 3/4. Diese Kontrolle berücksichtigt somit eine Originalanforderung,
ohne sie mit einer vollständigen lokalen OS-Netz-Realisierung gleichzusetzen.
Die freie Auswahl zwischen diesen beiden verbleibenden Kernen ist weiter offen.

## Die Suchbedingungen und die exakte Annahmegrenze

Das Alphabet besteht aus 2m selbstadjungierten Majoranas und einem zusätzlichen
selbstadjungierten T. Die Majoranas erfüllen CAR. Zwischen T und den Majoranas
wird keine Kommutatorrelation und für T keine Projektoreigenschaft eingesetzt.
Insbesondere bleibt T T ein Wort mit zwei Transferfaktoren.

Für eine unbekannte Momentfunktion omega wird M[u,v]=omega(u* v) gesucht.
Die Bedingungen sind Normierung, Sternverträglichkeit, M >= 0, T Omega=Omega
und Invarianz unter gamma_(2j)->gamma_(2j+1), gamma_(2j+1)->-gamma_(2j), T->T.
Stationarität entfernt ausschließlich äußere T in Erwartungswerten.

Für a=transfer_floor enthalten die endlichen Relaxationen getrennte Localizer
für T-aI, I-T und zusätzlich (T-aI)(I-T)=(1+a)T-T²-aI. Die letzte Bedingung ist
für einen echten beschränkten T notwendig und verschärft die endliche Stufe.
Ohne sie passierte ein Gegenfall die beiden ersten lokalen Schranken, obwohl
omega(gamma T² gamma)=13/16 > omega(gamma T gamma)=1/2 war. Das ist nun eine
Regression. Ein früherer Floor-Lauf mit falscher zweiter Transferantwort wird
ebenfalls verworfen; seine Daten bleiben als Entwicklungsgegenfall erhalten.

Die Fünf-Moden-SDP-Matrix ist nicht flach: Ränge 6 und 16 auf den beiden
Wortstufen. Deshalb wird aus ihr kein eindeutiger GNS-Abschluss behauptet.
`finite_completion.py` sucht stattdessen ausdrücklich in einer zusätzlichen
minimalen irreduziblen CAR-Darstellung der Dimension 2^m. Sie ist nur zulässig,
wenn die gelieferten CAR-Momente zum inferierten reinen Besetzungszustand passen.
Der Transfer bleibt eine unbekannte hermitesche Matrix. Nach der Optimierung
werden CAR, Positivität, Untergrenze, Stationarität, Viertelwirkung und sämtliche
gelieferten Momente exakt geprüft. Ein Fehlversuch bleibt unentschieden, kein
Ausschluss aller Darstellungen.

Ein separater Ein-Moden-Lauf demonstriert die andere, allgemeine Route:
`runs/m1-exact-certificate.json` rekonstruiert ohne vorgegebenen Fockabschluss
einen Rang-2-Kern durch exakte Gram-Faktorisierung auf drei aufeinanderfolgenden
Wortstufen. Alle 240 Momente werden exakt wiedergegeben; alle Generatoren
erhalten die Nullrelationen. Für Momente wie 1/256 beträgt die dokumentierte
Rationalisierungsgrenze 256. Eine zu kleine Grenze 64 wird nicht als physische
Unmöglichkeit interpretiert.

## Warum die erfolgreichen Vervollständigungen vollständige Kerne sind

Die gespeicherten exakten Matrizen und Omega definieren für beliebig lange Wörter

    K(u,v) = <A_u Omega, A_v Omega>.

Jede endliche Gram-Matrix ist deshalb positiv. Multiplikation und Adjungierung
werden durch tatsächliche Matrizen realisiert. Ein Nullvektor bleibt unter jeder
Generatorfortsetzung null. Die CAR-Wörter erzeugen im minimalen Abschluss den
gesamten Raum. Die Vollständigkeit dieses endlichen Kerns folgt daraus für alle
Wortlängen, nicht aus einem numerischen Rangtest.

Der GNS-Implementer der internen Vierteldrehung ist R=i^(N-N_Omega), mit R^4=I
und R Omega=Omega. Das ursprüngliche interne W=i^N wird separat ausgewiesen.
Die globale Referenzphase ändert die Feldkonjugation nicht. Diese interne
Wirkung ist noch keine Identifikation mit einer chiralen Nahtdomäne.

Wenn T strikt positiv ist, ist H=-log(T) ein endlicher selbstadjungierter Operator.
U_t=T^(it) und alpha_t(A)=U_t A U_t* liefern eine unitäre Feldzeit und eine
UCP-Automorphismusgruppe; Omega bleibt invariant. Die Zeiteinheit ist dabei
deklariert. Dies beweist weder einen euklidischen UCP-Kanal mit GNS-Implementer T
noch ein physisches Ereignisinstrument, eine native Zeitskala oder Raumzeit.

## Wie die Suche künftig statt im Kreis weiterarbeitet

`run_search.py` verbindet unbekannte Momente, exakte Vervollständigung und
Unterscheidung der Treffer in einem Aufruf. Kandidatenvorhersagen werden niemals
automatisch zu Quellenannahmen gemacht. Code, Konfiguration und Ausgabedateien
werden für die Wiederverwendung geprüft; geänderte Bedingungen brauchen einen
neuen Laufordner. Ungeänderte abgeschlossene Läufe werden wiederverwendet.

Die erzeugte nächste Frage lautet konkret: Welche **ursprüngliche** Nahtrelation
koppelt Transfer nach zwei Feldoperationen an die Einfeld-Antwort? Ein
Wick-/quasifreier Abschluss oder eine additive Zweitquantisierung wäre eine
solche zusätzliche Bindung, muss aber für dieselbe Quelle begründet sein.
Ein Treffer allein rechtfertigt nicht, diese Regeln einzusetzen.

Im Atlas existiert bereits die vollständige Produktkragen-Rekonstruktion:
der gesamte Zweiflächen-Graph bestimmt unter seinen Hypothesen B; die zusätzliche
CAR-/Wick-Identifikation verbindet den eingeschränkten Zustand mit derselben
Zeit. Diese positive Vorarbeit wird wiederverwendet. Die neue Methode sucht
unbekannte Kerne, statt nur den bekannten inversen Leser erneut auszuführen.

Die originalen Quelltypen bleiben getrennt: OS-Reflexionspositivität auf einem
positiven Zeitnetz ist mehr als die Gram-Positivität hier. P2 markiert einen
Einteilchenträger; eine 3+2-Besetzung von Moden wäre keine Herleitung seiner
Markierung. Eine nichttriviale Nahtindexklasse hat in dieser vollständig
endlichen Darstellung ebenfalls noch keine Realisierung. Somit ist die Methode
implementiert, der ursprüngliche gemeinsame TFPT-Quellenvertrag nicht geschlossen.

## Aufrufen und ohne Optimierer nachprüfen

Im Verzeichnis dieses Contracts:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_search.py --out runs/my-search
.venv/bin/python run_search.py --out runs/my-explicit-calibration --objective diversity --transfer-floor 0.25 --cluster-gap 0.25
.venv/bin/python proof_check.py runs/my-search
.venv/bin/python -m unittest discover -p 'test_*.py'
```

Die SDP arbeitet hier mit einer 78x78-Momentmatrix und die Vervollständigung mit
32x32-Matrizen. Der gemessene erste vollständige positive Lauf dauerte rund
6,5 Sekunden auf der lokalen CPU. Für diese Größe wurden keine H200 belegt.
Große parallele Kandidatenräume benötigen zuerst einen geeigneten numerischen
Engpass; GPU-Aktivität ist kein eigenständiges Forschungsresultat.

## Quellen und Beziehungen

- GNS-Original im Repo: `universalraum-primitive-source-audit-20260915/RESULTS.md`,
  Zeilen 145–194: vollständiger positiver Kern und Nullstabilität.
- `_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex`, Zeilen 480–520:
  ursprüngliche operational seed category mit OS-Netz, Reflexion und Nahtklasse.
- `_archive/tfpt-45/source_extracts/02_carrier_source.tex`, Zeilen 172–205:
  tatsächliche P2-Trägermarkierung und Hyperladung.
- `seam-full-collar-reconstruction-20261015/PROOF.txt`, Zeilen 63–143:
  vollständiger Graph sowie bedingte Zustands-/Zeitverbindung.
- Pironio, Navascués, Acín: [Convergent relaxations of polynomial optimization
  problems with non-commuting variables](https://arxiv.org/abs/0903.4368).
  Die Arbeit begründet das Verfahren hierarchischer positiver Momentrelaxationen
  und Rekonstruktionen unter passenden Voraussetzungen. Ihre Konvergenz wird
  hier nicht als TFPT-Kontinuumsbeweis beansprucht.

Verdict: **PARTIAL für den ursprünglichen Forschungsauftrag; exakte endliche
Konstruktionen und ein funktionierender automatischer Such-/Unterscheidungsweg.**
Kein neuer P1/P2-Beweis und keine Änderung des Beweisledgers.


## Withheld response search added on 2 November 2026

The existing campaign accepts `--withheld-response FILE`. The file has only
`terms`, for example `{"terms":[{"word":[1,-1,1],"re":"1"}]}`. A fitted
target is rejected. This example is a software control, not a physical P1
input. Source controls remain separate. The same relaxation is minimized and
maximized for this response. Numerical endpoints are estimates only. A
rational affine identity has an explicit replayable linear-combination
certificate, conditional on existence in the declared class. It can also
hold in a positivity-empty class; existence is not claimed by the identity.

```sh
python response_search.py --m 1 --degree 2 --response QUERY.json --source-control ORIGINAL_CONTROL.json --out NEW_RESULT.json
python response_search.py --verify NEW_RESULT.json
python run_search.py --out NEW_CAMPAIGN --source-control ORIGINAL_CONTROL.json --withheld-response QUERY.json
python certify_response_segment.py runs/source-clustering --out NEW_SEGMENT.json
```

Use the project's installed numerical dependencies. The current execution used
`/private/tmp/tfpt-kernel-relations-venv/bin/python`; that temporary environment
is not a permanent installation promise. `check_response_search.py` executes
complete exact 2D witness controls, tampering, contradictory inputs, complex
coefficients and positivity-empty cases. The CLI integration emitted its
conditional identity but retained an inconclusive finite rational recovery.

The stored five-mode source-clustering endpoints extend to an exact full
segment: all 2245 supplied moments are constant, while all ten existing
separating responses are `1/4+s/2`, `0<=s<=1`. These are attained values in
that constructed family, not global extrema. The original source must supply
an independent relation for this composite-field transfer or exclude the
class; no candidate answer is recycled as a physical premise. Full scope and
source reminders, including the existing conditional v534 selector, are in
[RESPONSE-SEARCH-20261102.txt](RESPONSE-SEARCH-20261102.txt).


## Source-indexed AI origin hypotheses, 2 November2026

`origin_hypothesis_search.py` uses the existing local Ollama backend to propose a typed hypothesis, validate its declared IDs and source premises, execute fixed curated checks and feed counterexamples into the next proposal. All proposals remain unverified hypotheses. Passing a check is not a fundamental physical solution. `GENESIS-SYNTHESIS-20261102.txt` records the common P1, double-cover, information, horizon and compiler context, the entire variable-mark conditional lattice class and the stronger primitive-projection selector.

The index includes all707 original registered Python sources with exact hashes, docstrings and check/assert text, the complete current ledger and reviewed Atlas records, and canonical research documents. This is static extraction, not a new run or proof review of707 modules. The model receives selected original excerpts; full source inventory is saved separately. Do not run indexing with `-OO`, which can strip docstrings.

```sh
python3 -B origin_hypothesis_search.py --catalog genesis_checks.json --rounds 4 --output runs/NEW_RUN.json
```

Run from this directory. A saved full inventory can be reused with `--reuse-inventory runs/OLD_RUN.inventory.json`; the driver checks the Atlas, registry, ledger and all original script hashes before reuse. Output names must be new. `--question` changes the mathematical proposal task. `--rounds` is an explicit positive number, not a scheduled service. The current curated tests cover four source/class questions; arbitrary new physical hypotheses still require an actual derivation or a new reviewed checker. No model text or source ID is silently inserted as a constraint, and no model-generated command is executed.

Completed run: `runs/genesis-20261102/qwen-structured-completion.json` contains two actual local qwen3.6:35b rounds, including a wrong primitive-projection prediction and its correction after the exact checker. Older failed/invalid proposals are retained. `curated-check-results.json` records all four computations. The raw broad hypothesis remains unproved.

```sh
python3 -B -m unittest test_origin_hypothesis_search test_adaptive_discriminator
python3 -B -OO genesis_checks.py --check variable-mark-completion
python3 -B -OO adaptive_discriminator.py runs/source-clustering --out runs/NEW_RANKING.json
```

The adaptive discriminator replays the existing full segment certificates and groups affinely equivalent response questions. It ranks separation per declared cost only on supplied candidates; it is not Shannon information or a search over all physical processes.
