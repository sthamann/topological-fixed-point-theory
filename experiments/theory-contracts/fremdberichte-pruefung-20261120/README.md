# Prüfung zweier Fremdberichte, 20. November 2026

Anlass: Der Nutzer hat zwei Berichte eines anderen Assistenten eingefügt und um Prüfung gebeten. Die dort verlinkten
Dateien lagen hier nicht vor; geprüft ist nur, was sich aus dem Berichtstext und dem Repository unabhängig nachrechnen
lässt. Rechnung: `check_reports.py`, Ausgabe `check_reports.json`. Kein Ledgerwechsel, kein Atlaseintrag.

## Bericht A: gemeinsamer Quantenprozess

Bestätigt (11):

- Q = (1 − T01·T23·T45)/2 projiziert auf die drei Differenzen; U6 und M6 wirken dort als die angegebene Vierteldrehung
  und Dreieroperation.
- u und m erzeugen eine Gruppe der Ordnung 24, und es gibt genau einen unitären Intertwiner zum Paar aus `v117`.
- Unter der Kleinschen Vierergruppe: H5 = 2·1 + χ1 + χ2 + χ3, der Dreierraum trägt χ1 + χ2 + χ3.
- B hat die Eigenwerte 1, 2/3, 1/3. Der angegebene Generator hat die Raten log(2/3) und log(1/3).
- Vollständige Positivität erzwingt 2/9 ≤ η ≤ 4/9; eine stetige kovariante Entwicklung 0,3205 ≤ η ≤ 0,3467.
- N = P1 + 3/4·P5 + 11/36·P9 + 1/4·P10 ist ein Quantenkanal; sein kovarianter Logarithmus fällt im Positivitätstest
  mit −0,087696 durch. Raten, Verhältnis exp(4t/9) und die Zahlen bei t = 2 stimmen.

Nicht geprüft (5), weil die Operatoren des Berichts hier nicht vorlagen: Norm 10/9 und M ≥ 2/405; die Umweltisometrie
zur Dreikindrekursion; Rang 25 mit drei Einstellungen und d_R; die native Rückkehr 59/84, 11/36, 1/4; der gesamte
Koide Teil einschließlich des Intervalls für h.

Einordnung: Die geprüften Aussagen sind richtig, darunter zwei, die man nicht errät (−0,087696 und das enge
η Intervall). Das stützt die Sorgfalt des Berichts, beweist aber die ungeprüften Teile nicht. Der Bericht sagt selbst,
dass die verlangte Kopplung auf dem Dreierraum nicht folgt und dass Koide nicht getroffen wird.

## Bericht B: P1 als Verlust der Unterscheidbarkeit

Bestätigt (8):

- Drei Kopplungsmuster (2,1,1), (1,2,1), (1,1,2): tr(P1·P2·P3) = −1/8.
- Matrix aus `v395`: Sektoren 1, 2/3, 1/3 geben bei zyklischer Versetzung 0, +1/8, −1/8.
- Stern und Dreieck haben dieselbe statische Randantwort; mit speicherndem Mittelpunkt 3s/(9 + s) je Anschluss.
- Zwei Ketten mit 8 und 9 Gliedern unterscheiden sich erst in Ordnung 18.
- Von 24 umkehrbaren Regeln auf zwei Bits erfüllen genau fünf die Austauschbedingung; es sind die Identität und vier
  Tausche mit festen Bitumkehrungen. Die Quantenphase auf 11 ist für jede Phase konsistent; zweimal tauschen gleich
  Identität lässt +1 und −1.

Nicht gezeigt: dass die Theorie die zyklische Versetzung selbst erzeugt. Der Bericht nennt das als Annahme.

Einordnung: Alle Rechnungen stimmen. Die Vorzeichen ±1/8 folgen allein daraus, dass die eine Mode die Summe und die
andere die Differenz zweier Kanäle ist; sie sind keine Eigenschaft der Zahlen 2/3 und 1/3. Der Befund ist ein
möglicher Mechanismus, kein Beleg.

## Verbindung zum Neuansatz am singulären Punkt

Beide Ansätze setzen an derselben Stelle an. Der Punkt gleicher Kopplungen im Bericht B ist ein Entartungspunkt mit
Vorzeichenwechsel nach einem Umlauf; τ = i ist ein Entartungspunkt mit Deckabbildung als Quadrat der Vierteldrehung.
Ob das dieselbe Entartung ist, ist offen und prüfbar: Die drei Schritte vom Grad 2 aus dem Torus (die drei Marken)
wären die drei Paarungen des Berichts B.


## Dritter und vierter Bericht (Nachtrag, 20. November)

Aufrufe: `check_report3.py` (14 von 14), `check_report4.py` (8 von 8). Die verlinkten Rechenpakete liegen in einem
privaten Konto und wurden nicht geöffnet; geprüft ist nur, was sich am Repository und exakt nachrechnen lässt.

**Bestätigt (dritter Bericht).** Am unveränderten Generator aus `v999` (sechs Sprungoperatoren, Zeiteinheit eine
Hand): Spektrum 1, 8/81 zweimal, 64/729, 2√2/81 viermal, 1/729; A = diag(1,−1,0) ist Eigenoperator mit Rate
Delta = 6·log(3/2), also trägt die Auslese r = tr(A·rho)/2 den ganzen Koidefluss; Delta − log(81/8) = log(9/8);
Exponenten 0,95158… und 2,70951…. Fluss dq/dt = (Delta/3)(q−2)(q−5) und Koide als Q = 1/(3|<s|x>|²). In Q(zeta_30):
C^30 = 1, G^4 = 1, G·C·G^-1 = C^7; invariante Formen 4 (nur C) und 1 (beide), die verbleibende ist positiv;
J = 2C^5 − 1 mit J² = −3.

**Bestätigt (vierter Bericht).** Der Kanal 1/2·id + 1/3·W3 + 1/6·D3 ist vollständig positiv und spurtreu, hat das
Spektrum 1, 2/3 dreimal, 1/3 fünfmal und wirkt in der Clockbasis der Vierteldrehung auf die Populationen exakt als B.
In der reellen Basis gibt er dagegen (1/3)I + (2/9)J, nicht B. End_0(F) zerfällt unter V4 wie angegeben; die
Kommutatorregeln [5,5] in 3, [3,5] in 5, [3,3] in 3 stimmen.

**Nicht geprüft (Operatoren oder Dateien fehlen).** Das paarunsichtbare Zustandspaar (56,82 % gegen 87,73 %), der
Kommutator 4√2/75·W_A, die Ereignisuhr (10/9, 2/405), die Umweltisometrie, die 30 Effekte vom Rang 25, der kritische
Code (121/2187, 173/77), 7/144, 7/6, die Zerlegung von N, die Higgsvertices am 96 dimensionalen Operator, das
Intervall für h.

**Einordnung.** Beide Berichte sagen selbst, dass die volle Theorie nicht gelöst ist; das stimmt mit der Nachrechnung
überein. Die Gewichte 1/2, 1/3, 1/6 des Familienkanals sind gesetzt, nicht hergeleitet (der Bericht sagt das).
