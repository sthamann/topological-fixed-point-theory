# Quellengebundene Suche: eine oder zwei neutrale Grad-2-Richtungen

## Ergebnis

In der explizit benannten Klasse eines einzelnen neutralen Vertices
\(F=aR+bR_Y\), mit reellen \(a,b\), gibt es **keinen** Operator, der gleichzeitig Familienantwort 1, Normquadrat 48 und eine nichtverschwindende hypergeladene chirale Wardantwort hat. Das ist ein globaler Ausschluss für alle reellen Koeffizienten, nicht bloß ein Fehlschlag des endlichen Koeffizientengitters.

Die nächstliegende zulässige Quellenkonstruktion braucht zwei getrennte Richtungen derselben Originalquelle:

\[
S_Y=R_Y-\frac5{24}R,\quad
\operatorname{Gram}(R,S_Y)=\begin{pmatrix}48&0\\0&115/12\end{pmatrix},
\quad \langle R,J_YJ_Y\rangle=0,\quad
\langle S_Y,J_YJ_Y\rangle=115/36.
\]

Damit ist die erforderliche Rang-2-Struktur der hier benannten Klasse exakt belegt. Weder folgt daraus, dass TFPT insgesamt einen rang-zwei primitiven Port benötigt, noch ein 4D-Photon-Polarisator.

## Reproduktion und Originaldaten

`source_joint_vertex_search.py` ruft direkt `tfpt_explorer.source_charge_response.exact_charge_response()` auf. Diese Funktion baut Hyperladungsgewichte aus denselben 48 retained E8-Wurzeln, berechnet den ladungsgewichteten Momenttensors, dessen spurfreien Rest \(R_Y\), den ursprünglichen Familienoperator \(R\), Gramformen und \(J_Y\)-Wardantwort. Der neue Code bildet lediglich \(S_Y=R_Y-(5/24)R\) und prüft die behaupteten Identitäten symbolisch.

Für die Familienantwort auf die schweren neutralen Paare ruft das Skript separat `source_majorana_pairs._native_pairs()` und `_joint_neutral_pair(...)` auf. Die \(R_Y\)-Antwort wird aus den tatsächlichen Quellwurzelmomenten und Paar-Oszillatortermen der sechs von `exact_charge_response` erzeugten Paarantworten ausgewertet. Die Werte für \(S_Y\) folgen anschließend durch die exakte Quellenidentität; es werden keine gespeicherten Resultatmetadaten als Eingabe benutzt.

Reproduktion:

```bash
python3 source_joint_vertex_search.py --outputpath results/vertices.json
```

Die ausführliche maschinenlesbare Ausgabe steht in `results/vertices.json`; alle elf exakten Kontrollbedingungen sind `true`.

## Beweis des Ein-Vertex-Ausschlusses

Die direkte Originalrechnung liefert

\[
R_Y=\frac5{24}R+S_Y,\quad
\|R\|^2=48,\quad \langle R,S_Y\rangle=0,\quad
\|S_Y\|^2=\frac{115}{12}.
\]

Außerdem sind die chiralen Wardantworten

\[
W_R=\langle R,J_YJ_Y\rangle=0,\qquad
W_{R_Y}=\frac{115}{36},\qquad W_{S_Y}=\frac{115}{36}.
\]

Für \(F=aR+bR_Y\) ist die Familienkomponente \(f=a+5b/24\). Die Norm folgt aus der orthogonalen Zerlegung:

\[
\|F\|^2=48\left(a+\frac5{24}b\right)^2+\frac{115}{12}b^2.
\]

Setze \(f=1\). Die Normbedingung \(\|F\|^2=48\) ergibt \(115b^2/12=0\), also \(b=0\). Dann ist \(W_F=0\). Anders ausgedrückt: Für jeden realen \(b\) ist bei Familienkomponente 1 der Wardwert \(115b/36\), aber die Norm wächst strikt über 48, sobald \(b\ne0\). Ein komplexer Koeffizientenraum wurde hier nicht beansprucht.

Zusätzlich wurden alle 8.649 Paare \((a,b)\) aus dem ausdrücklich vorgegebenen rationalen Raster (Zähler −12 bis 12, Nenner 1,2,3,4,6,12,24, Duplikate entfernt) exakt durchlaufen: null erfüllende Kandidaten. Dieser endliche Lauf dient als reproduzierbare Kontrolle; der reelle Beweis ist stärker.

## Zweirichtungs-Kandidat und Paarantwort

Die Quellmatrix gibt die Gramform

\[
\operatorname{Gram}(R,R_Y)=\begin{pmatrix}48&10\\10&35/3\end{pmatrix},
\qquad
R_Y=\frac5{24}R+S_Y,
\]

und daher für \((R,S_Y)\) die positive diagonale Gramform mit Determinante 460. Die Quellenantwort trägt damit zwei unabhängige Richtungen, Familienkomponente und Hyperladungsantwort getrennt. Der Rekonstruktionskoeffizient \(5/24\) kommt aus derselben exakten 48-Wurzelrechnung und ist nicht gefittet.

Auf allen sechs retained schweren Paaren berechnet die ursprüngliche Paarantwort

\[
W_R=4,\qquad W_{R_Y}=-\frac53,
\qquad W_{S_Y}=W_{R_Y}-\frac5{24}W_R=-\frac52.
\]

Die Paare sind selbst hyperladungsneutral; diese Antwort ist eine Operator-/Wardantwort, keine Behauptung, dass ein einzelnes Paar elektrische Ladung trage oder kondensiere. Der hypergeladene chirale Testkanal hat Quellenwert \(115/36\) nur im chiral-korrelatorischen Sinn.

## Filterkorrektur im vorigen Bericht

Die vorhandene neutrale Korrelatorform hat führenden Term \(48/(z-w)^4\). Für den zugehörigen Gewicht-2-Modus folgt über die Konturkommutator-Normierung

\[
\kappa_n=\frac{48}{6}n(n^2-1)=8n(n^2-1),\qquad n\ge2,
\]

nicht \(4n(n^2-1)\). Zweitseitiges \(D^{-2}\)-Smearing gibt daher

\[
\frac{\kappa_n}{n^4}=8\left(\frac1n-\frac1{n^3}\right).
\]

Ein Modengewicht, das stattdessen genau \(1/n\) erzeugt, müsste \(|S_n|^2=1/[8n^2(n^2-1)]\), also \(S_n=1/[2\sqrt2\,n\sqrt{n^2-1}]\), erfüllen. Das bleibt eine nachträgliche Inversion des Zielkerns, solange keine Originalregel den Filter auswählt.

## Erste offene Herkunftsbedingung

Die Rechnung liefert zwei quelleneigene Operatorantworten; sie beweist nicht, dass der primitive P1/P2-Randkern diese beiden Richtungen erhält und zusammenführt. Für diese Kandidatenklasse wäre die nächste echte Auswahlbedingung ein gemeinsamer Originalzustands-/Randoperator, dessen geladene Mehrpunktfunktion gleichzeitig die Familien-Wardantwort, die Hyperladungs-Wardantwort und die positive Gramform oben erzeugt. Keine weitere Anpassung von \(a,b\) in der Ein-Vertex-Klasse kann diese Bedingung erfüllen.
