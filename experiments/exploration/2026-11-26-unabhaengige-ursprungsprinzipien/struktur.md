# Strukturtest: vierpolige Netzantwort als unabhängiger Ursprung von 3 + 2

## Ergebnis

Die Spur liefert einen echten neuen Strukturtreffer, aber noch keinen Gesamtursprung.

Aus einer passiven reziproken Netzantwort mit `n` abstrakten Anschlusspolen folgt ein
Antwortvektor auf den ungeordneten Polpaaren. Seine Dimension ist `n(n−1)/2`. Die
gemeinsame Leitwertskala ist physikalisch nur eine Skalierung der ganzen Antwort. Gesucht
wurde deshalb nicht `5`, sondern der kleinste Polzahlwert, für den die skalenfreie
Antwortform unter allen Polvertauschungen in mehr als einen nichttrivialen irreduziblen
Sektor zerfällt.

Dieser erste Wert ist exakt `n = 4`:

\[
\mathbb R^{\binom{4}{2}}
  = \mathbf 1 \oplus \mathbf 3 \oplus \mathbf 2,
\qquad
\mathbb P(\text{Antwortskala entfernt})
  = \mathbf 3 \oplus \mathbf 2.
\]

Die `3 + 2` sind hier keine bloßen Dimensionsnamen. `strukturtest.py` konstruiert drei
orthogonale Projektoren auf den Skalensektor sowie die beiden S4 Module, prüft ihre
Idempotenz und Vertauschung mit allen 24 Polpermutationen und rekonstruiert die sechs
effektiven Paarleitwerte aus einer positiven Dirichlet zu Neumann Matrix.

## Ausgangsprinzip und kleinster entscheidender Test

Ausgangspunkt ist nur:

1. Das physisch Zugängliche ist die vollständige lineare Randantwort eines passiven
   reziproken Verbunds.
2. Unbenannte Anschlusspole sind vor einer weiteren geometrischen Wahl gleichwertig.
3. Eine gemeinsame Multiplikation aller Leitwerte ändert nur die Antwortskala.
4. Gesucht ist der kleinste Verbund, dessen skalenfreie Form mehr als einen
   nichttrivialen Symmetriesektor besitzt.

Der Erfolgstest war vor der Rechnung: Die minimale Antwort muss unter ihrer tatsächlichen
Permutationswirkung invariant in zwei irreduzible Sektoren der Dimensionen drei und zwei
zerfallen. Abbruchkriterium war: nur ein Zahlentreffer ohne äquivariante Projektoren oder
derselbe Treffer bereits bei `n = 3` oder `n = 5`.

Die Gegenfamilien sind scharf:

| Pole | Paarantwort | skalenfreie irreduzible Sektoren |
|---:|---:|---:|
| 2 | 1 | keine |
| 3 | 3 | 2 |
| **4** | **6** | **3 + 2** |
| 5 | 10 | 4 + 5 |

Allgemein zerfällt die Permutationsdarstellung auf Zweiermengen für `n ≥ 4` als

\[
[n]\oplus[n-1,1]\oplus[n-2,2]
\]

mit Dimensionen `1`, `n−1` und `n(n−3)/2`. Vier ist damit der erste mehrsektorige Fall,
nicht ein nachträglich gewählter Treffer.

## Was die Rechnung tatsächlich trägt

Eine Netzantwort ist ein konkreter Operator. Für positive Leitwerte ist ihre
Dirichlet zu Neumann Matrix symmetrisch, positiv semidefinit und hat genau die konstanten
Spannungen als Kern. Innere Knoten sind bereits durch das Schurkomplement eliminiert. Die
sechs Offdiagonalen bestimmen die effektiven Paarleitwerte. Die drei plus zwei Sektoren
sind damit verschiedene Formen derselben vollständigen Randantwort, keine separat
eingesetzten Teilräume.

Das ist gegenüber dem bisherigen TFPT Bestand neu: Die Atlas Suche fand Calderón Operatoren,
Netzattraktoren und frühere E6 plus A2 Rechnungen, aber keinen gespeicherten Versuch, der
vier Pole und die 3 plus 2 Teilung aus der Zweiermengen Darstellung einer reziproken
Antwort herleitet. Ein fehlender Treffer allein wäre kein Neuheitsbeweis. Entscheidend ist,
dass der konkrete Mechanismus in den gefundenen Versuchen nicht vorkommt.

## Ehrliche Gegenkontrollen

### Das natürliche Gitter ist falsch

Die ganzzahligen Antwortformen mit verschwindender Gesamtsumme bilden natürlich das
Gitter

\[
A_5=\{z\in\mathbb Z^6:\sum z_i=0\},\qquad \det A_5=6.
\]

Der Dreisektor ist ein um Faktor zwei skaliertes A3 Gitter mit Determinante 32. Der
Zweisektor ist ein skaliertes A2 Gitter mit Determinante 12. Ihre orthogonale Summe hat
Determinante 384 und liegt mit Index 8 in A5. Das TFPT Trägergitter D5 hat dagegen
Determinante 4. Die Netzantwort erzeugt also nicht heimlich D5. Damit sind auch E8 Glue,
die 16 Komponenten und die physische Standardmodellwirkung aus dieser Rechnung nicht
hergeleitet.

### Der genaue H6 Abgleich: V4 trifft, D4 braucht eine Signumlinie

Der vollständige Vergleich mit der tatsächlichen H6 Trägerdarstellung aus Z34 ist stärker
als der Gittervergleich. In der Klassenreihenfolge `(e, r², r, s, rs)` hat die skalenfreie
Netzantwort den D4 Charakter

\[
\chi_{\rm Netz}=(5,1,-1,1,1).
\]

Ihr Dreisektor ist exakt die Familienrepräsentation

\[
\chi_F=(3,-1,-1,1,-1),
\]

während ihr Zweisektor den Charakter `(2,2,0,0,2)` trägt. Der tatsächliche V4 Schnitt des
H6 Trägers hat dasselbe schwache Dublett, aber als Farbtripel nicht `F`, sondern
`Λ²F`:

\[
\chi_{\Lambda^2F}=(3,-1,1,-1,-1),\qquad
\chi_{\rm H6}=(5,1,1,-1,1).
\]

Das Skript löst die Intertwiner Gleichungen exakt. Unter dem geraden Klein Vierer gibt es
einen invertierbaren Intertwiner vom ganzen Netzraum auf den H6 Träger. Unter vollem D4
hat jeder unverdrehte Intertwiner höchstens Rang vier. Das ist ein echter Ausschluss, kein
fehlender Basiswechsel.

Die Abweichung ist jedoch vollständig typisiert. Sei

\[
\delta=\det F=\operatorname{sgn}(\text{Polpermutation}),
\]

also `δ(r)=δ(s)=−1` und `δ(r²)=δ(rs)=+1`. Dann gilt exakt

\[
\Lambda^2F\cong\delta\otimes F,
\qquad
V_{\rm H6}\cong\delta\otimes V_{\rm Netz}.
\]

Für die verdrehte Netzform findet das Skript einen D4 Intertwiner vollen Rangs mit
Determinante 256. Die Signumlinie ist eine reelle eindimensionale Determinantenlinie. Sie
ist weder die komplexe μ4 Phase noch ein Spin(10) Lift. Bildlich ist sie das Vorzeichen
eines alternierenden Vierpolvolumens. Die kanonische dreidimensionale Identität
`F ⊗ det(F) → Λ²F` ist genau die Hodge Abbildung ohne vorab gewählte Orientierung.

Damit ist die neue Prämisse scharf: Die Paarantwort müsste nicht skalar, sondern in der
alternierenden Determinantenlinie der vier Pole leben. Ein gewöhnliches passives
Widerstandsnetz liefert diese Verdrehung nicht. Auch die Kreisorientierung genügt nicht,
denn die Vierteldrehung `r` erhält die Orientierung der Scheibe, besitzt als Permutation
der vier Pole aber trotzdem Signum `−1`. Die Prämisse ist eine alternierende Polhebung,
keine bloße Wahl von Uhrzeigersinn. Unter V4 ist sie unsichtbar, weshalb dort der Treffer
vollständig war.

### Statische Antwort ist noch keine Zeit

Passivität gibt Positivität der Randenergie. Ein Widerstandsnetz liefert aber nur eine
statische Antwort. Zustand, kontinuierliche Zeit, 3+1 dimensionale Raumzeit, Dynamik,
Gravitation und die Normierung `c3 = 1/(8π)` sind nicht gewonnen. Eine Erweiterung zu
positiv reellen Impedanzfunktionen würde neue Speicherbauteile und ihre Werte einführen.
Das wäre ein neuer Ansatz und keine automatische Fortsetzung dieses Treffers.

## Urteil

**Kandidat Stufe 1.** Der Ansatz erzeugt aus einem unabhängigen Inversproblemprinzip die
Vier und eine tatsächliche 3 plus 2 Modulzerlegung derselben positiven Randantwort. Er liegt
über der Zahlenbasislinie, weil Projektoren, Gruppenwirkung, Positivität, Rekonstruktion und
Gegenfamilien gemeinsam geprüft sind.

Er erzeugt P2 noch nicht vollständig. Die fünf Antwortformen sind Formrichtungen, keine
fünf bereits identifizierten Trägerplätze. Das natürliche Gitter ist A5 statt D5. P1 fehlt
bis auf den allgemeinen positiven Randkern vollständig. Der stärkste nächste Test wäre
daher nicht eine weitere Netzvariante. Zuerst müsste die alternierende Determinantenlinie
aus demselben reziproken Antwortprinzip folgen. Danach braucht es eine kanonische Abbildung
der verdrehten Antwortform auf fünf physische Trägerplätze mit D5 Integralstruktur.
Existiert die erste Hebung nicht, endet der Weg bereits am vollen D4. Existiert sie, bleibt
das A5 gegen D5 Hindernis. Windung oder Dynamik vorher hinzuzufügen würde diese beiden
entscheidenden Engpässe nur verdecken.

## Neuheitsprüfung gegen den Bestand

Der Atlas wurde nach elektrischen Netzen, Dirichlet zu Neumann Antworten, vier Polen,
Paarleitwerten und der S4 Zweiermengen Darstellung durchsucht. Gefunden wurden vorhandene
Calderón Operatoren, K4 Verbindungen in einer Twistorzentrum Rechnung, allgemeine
Netzattraktoren und die H6 V4 Trägerprüfung. Der Inhaltsindex lieferte für den kombinierten
Vierterminal Netzantwort Mechanismus keinen Abschnittstreffer. Die Codegraph Suche in
OriginLab fand gerichtete unitäre Transportnetze und verschiedene Quellantworten, aber
keine passive Paarleitwertantwort oder ihre `1 + 3 + 2` Zerlegung.

`K4` kam also bereits als Graph oder Konnektivitätswort vor. Neu ist nicht der vollständige
Graph auf vier Knoten. Neu ist die Kette

\[
\text{reziproke Randantwort}
\longrightarrow \mathbb R^{\binom42}
\longrightarrow \mathbf1\oplus\mathbf3\oplus\mathbf2
\longrightarrow \delta\text{ Verdrehung}
\longrightarrow V_{\rm H6}.
\]

Die Atlasabwesenheit allein beweist keine weltweite Neuheit. Sie belegt hier nur, dass der
konkrete Mechanismus nicht in den geprüften aktuellen TFPT Versuchen gespeichert ist.

## Primärquellen

- Curtis, Ingerman und Morrow, [Circular planar graphs and resistor networks](https://doi.org/10.1016/S0024-3795(98)10087-3), 1998. Definiert die Netzantwort, gibt ihre Schurkomplementform und rekonstruiert kritische kreisplanare Netze aus der Antwort.
- Lam und Pylyavskyy, [Electrical networks and Lie theory](https://arxiv.org/abs/1103.3475), 2011. Ordnet positive kreisplanare Netzoperationen einer konkreten Gruppenwirkung und Stern Dreieck Relationen zu.

Die S4 Zerlegung, die exakten Projektoren, die Polzahl Gegenfamilien und die Gitterkontrolle
sind die neue Rechnung in `strukturtest.py`. Die Literatur wird nur für den Netzantwort und
Inverse Problem Rahmen verwendet.
