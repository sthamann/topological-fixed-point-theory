# Fortsetzung: Grunddoublett und die tatsächlich zustandsauflösende Algebra

## Ergebnis

Für das unveränderte

\[
K_\alpha(p,N)=N+\frac{|p|^2}{2}+\delta_\alpha\cdot p,
\qquad
(\delta_\alpha)_k=\frac{k+\alpha}{8}-\frac12,
\]

mit (\alpha=1/4,3/4) gibt es auf

\[
\mathbb Z^8\ \cup\ (\mathbb Z^8+s),
\qquad s=(1/2)^8,
\]

genau zwei (N=0)-Grundvektoren:

\[
|0\rangle,\qquad
|r_0\rangle,\qquad
r_0=(1/2,1/2,1/2,1/2,-1/2,-1/2,-1/2,-1/2).
\]

Beide liegen auch im markierten E8-Gitter

\[
L_s=D_8\cup(D_8+s).
\]

Die Zustandsantwort hängt danach entscheidend von der Algebra ab:

| Gegebene Daten | Exakter Schluss |
|---|---|
| (K_\alpha\ge0) | keine Zustandsauswahl |
| normaler Zustand mit (\omega(K_\alpha)=0) | Träger im zweidimensionalen Grundraum (D) |
| zusätzlich dual-(\mathbb Z_2)-invariant | Dichte (\operatorname{diag}(a,1-a)); nicht eindeutig |
| zusätzlich (R_c)-invariant | eindeutig (I_2/2), aber diese Invarianz ist eine weitere Annahme |
| ursprüngliche vollständige Cartan-Vakuumwerte | eindeutig (|0\rangle\langle0|) |
| Cartan-Vakuumwerte und (R_c)-Invarianz zugleich | unvereinbar |

Auf der bloßen Verschiebungsalgebra sehen alle dual-invarianten
Grundmischungen sogar wie derselbe Zustandsfunktional aus. Erst die erhaltenen
Cartan-Nullmoden unterscheiden die beiden Dichten. In der lokalen
Lattice-VOA ist (|0\rangle) durch die Vakuumstruktur ausgezeichnet;
(|r_0\rangle) ist ein geladener Gewicht-eins-Primärzustand, obwohl er für das
getwistete (K_\alpha) ebenfalls Energie null besitzt.

Das genaue Verdict lautet

```text
EXACT_GROUND_DOUBLET_DUAL_Z2_INSUFFICIENT_CARTAN_VACUUM_UNIQUE
```

## 1. Quellenstatus

Der ursprüngliche operative Seed enthält den Zustand (\omega) bereits als
Eingabedatum zusammen mit lokaler Algebra, Zeit, Reflexion, Nahtklasse und
Collaroperator ([`01_boundary_kernel_source.tex`, Zeilen 480–520](../../../../../_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex)). Der Collar-Satz rekonstruiert aus diesem Tupel weitere Daten; er wählt (\omega) nicht.

Die Zweifachentartung ist keine neue Entdeckung dieser Fortsetzung. Der
ältere Quellenvertrag
[`quotient-locality-e8-20261014/PROOF.txt`](../../../quotient-locality-e8-20261014/PROOF.txt)
beweist bereits durch Koordinatenkonvexität den eindeutigen zweiten
Nullvektor und vergleicht ihn mit der unveränderten mikroskopischen
QWZ/CAR-Quelle (dort Abschnitte 4–5, Zeilen 105–176). Die spätere
Nahtrechnung beweist ebenfalls die Zweifachentartung und die
bedingte Cartan-Auswahl in ihrem eigenen Vergleichsrahmen
([`seam-joint-ward-state-20261014/PROOF.txt`, Zeilen 55–82 und 113–171](../../../seam-joint-ward-state-20261014/PROOF.txt)). Die folgende Fortsetzung isoliert zusätzlich:

1. die vollständige Klassifikation zuerst auf
   (\mathbb Z^8\cup(\mathbb Z^8+s)),
2. die Rolle der dualen (\mathbb Z_2)-Invarianz,
3. die Antikommutation mit dem bereits konstruierten (R_c),
4. den Unterschied zwischen Verschiebungsalgebra, Cartan-Erweiterung und
   lokaler VOA.

Es wird kein Hamiltonoperator, Zustand oder Symmetriepostulat ergänzt.

## 2. Universeller Grundzustandsbeweis ohne Raster

Schreibe koordinatenweise

\[
f_k(x)=\frac{x^2}{2}+\delta_kx
      =\frac12(x+\delta_k)^2-\frac{\delta_k^2}{2}.
\]

Für beide Viertelwerte gilt

\[
-\frac12<\delta_k<\frac12,
\]

wobei (\delta_0,\ldots,\delta_3<0) und
(\delta_4,\ldots,\delta_7>0).

### Ganzgitter

Für (m\in\mathbb Z\setminus\{0\}) gilt

\[
f_k(m)-f_k(0)
=\frac{m^2}{2}+\delta_km
\ge |m|\left(\frac{|m|}{2}-|\delta_k|\right)>0.
\]

Also ist in jeder Koordinate (m=0) der eindeutige Minimierer. Damit ist
(p=0) der eindeutige Ganzgitter-Minimierer und hat Energie null.

### Halbgitter

Auf (\mathbb Z+1/2) ist wegen strikter Konvexität der eindeutige nächste
Punkt zu (-\delta_k)

\[
x_k^*=-\frac12\operatorname{sign}(\delta_k).
\]

Das ergibt genau (r_0). Der minimale Koordinatenwert ist

\[
f_k(x_k^*)=\frac18-\frac{|\delta_k|}{2}.
\]

Direkt aus den acht festen Koeffizienten folgt für jedes
(0<\alpha<1), insbesondere für beide Viertelwerte,

\[
\sum_{k=0}^7|\delta_k|
=\frac{10-4\alpha}{8}+\frac{6+4\alpha}{8}=2.
\]

Somit

\[
K_\alpha(r_0,0)=1-\frac12\sum_k|\delta_k|=0.
\]

Jede Vorzeichenänderung einer Halbkoordinate kostet (|\delta_k|>0), jede
Änderung des Betrags von (1/2) auf (3/2) mindestens
(1-|\delta_k|>|\delta_k|). Daher ist (r_0) auf dem gesamten unbeschränkten
Halbgitter eindeutig. Weil (N\ge0) ist, kann kein positiver Oszillatorgrad
einen weiteren Grundzustand liefern.

Schließlich

\[
r_0-s=(0,0,0,0,-1,-1,-1,-1)\in D_8.
\]

Deshalb bleiben genau dieselben zwei Minima nach der E8-Paritätsrestriktion
auf (L_s) erhalten. Für den gegenüberliegenden Spinorcoset beträgt die
kleinste Energie in beiden Viertelfällen (1/32); das ändert den hier
festgehaltenen (L_s)-Befund nicht.

## 3. Energiepositivität und duale (\mathbb Z_2)

(K_\alpha\ge0) ist eine Operatoreigenschaft und schränkt einen beliebigen
Zustand nicht ein. Erst die zusätzliche Grundzustandsbedingung
(\omega(K_\alpha)=0) erzwingt für normale Zustände den Träger

\[
D=\operatorname{span}\{|0\rangle,|r_0\rangle\}.
\]

In dieser Basis hat jede normale Grunddichte die Form

\[
\rho=\begin{pmatrix}a&c\\\bar c&b\end{pmatrix},
\qquad a,b\ge0,
\quad a+b=1,
\quad |c|^2\le ab.
\]

Die duale Sektorwirkung ist

\[
\Gamma=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\]

Aus (\Gamma\rho\Gamma=\rho) folgt nur (c=0). Somit bleibt

\[
\rho_a=\operatorname{diag}(a,1-a),\qquad 0\le a\le1.
\]

Energie-null plus duale Invarianz liefert also keine eindeutige Dichte.

## 4. Warum die kanonische Erwartung auf der kleinen Algebra trotzdem eindeutig aussieht

Sei zunächst nur die algebraische Verschiebungsalgebra gegeben,

\[
\mathcal A_{\rm shift}
=\operatorname{span}\{E_nX^c:n\in\mathbb Z^8, c=0,1\}.
\]

Die duale Wirkung ist (+1) auf (c=0) und (-1) auf (c=1). Ihre
Grad-null-Erwartung ist

\[
\mathcal E_0(E_nX^c)=
\begin{cases}
E_n,&c=0,\\
0,&c=1.
\end{cases}
\]

Mit der Ganzgitter-Vakuumrestriktion
(\omega_0(E_n)=\delta_{n,0}) ist

\[
\widetilde\omega=\omega_0\circ\mathcal E_0
\]

der eindeutige dual-invariante Fortsetzungsfunktional mit dieser
Restriktion. Das ist die bereits im Suchvertrag angegebene algebraische
Eindeutigkeit.

Sie ist aber keine eindeutige Dichtematrix auf (D). Jede (\rho_a) liefert

\[
\operatorname{Tr}(\rho_a E_nX^c)
=\delta_{n,0}\delta_{c,0}.
\]

Die Algebra kann ihre verschiedenen diagonalen Repräsentationsdichten nicht
unterscheiden. Ohne duale Invarianz würden Kohärenzen über die Verschiebung
zwischen (0) und (r_0) sichtbar. Die kanonische Erwartung beweist daher eine
eindeutige **kleine-Algebra-Fortsetzung**, nicht die Auswahl eines eindeutigen
physischen Grundvektors.

## 5. Die vorhandene (R_c)-Wirkung macht die Alternativen messerscharf

Der vorhandene affine Operator hebt

\[
\rho(p)=Pp+r_0
\]

an. Er erfüllt (\rho^2=1), erhält (K,Q) und vertauscht (0\leftrightarrow r_0)
([`seam-joint-ward-state-20261014/PROOF.txt`, Zeilen 820–832](../../../seam-joint-ward-state-20261014/PROOF.txt)). Sein tatsächlicher globaler Stromgruppenlift vertauscht die beiden Grundvektoren bis auf feste Phasen und ist ausdrücklich kein vakuumerhaltender VOA-Automorphismus (ebenda, Zeilen 934–993). Der angepasste (R_c)-Lift besitzt dieselbe Ladungswirkung und (R_c^2=1) (Zeilen 1070–1077).

Nach einer harmlosen Basisrephasierung gilt auf (D)

\[
R_c=\begin{pmatrix}0&1\\1&0\end{pmatrix},
\qquad
\Gamma R_c=-R_c\Gamma.
\]

Dual-Invarianz macht (\rho) diagonal; (R_c)-Invarianz vertauscht danach die
beiden Diagonaleinträge. Gemeinsam folgt eindeutig

\[
\rho_{\Gamma,R_c}=\frac12I_2.
\]

Das ist ein exakter bedingter Selektor, aber keine Quellenherleitung:
Physische (R_c)-Invarianz des Zustands ist nirgends gegeben. Sie darf nicht
aus ([K,R_c]=0) gefolgert werden.

## 6. Cartan-Erweiterung und lokale VOA

Sobald die ursprünglichen Cartan-Nullmoden in der Algebra erhalten werden,
ist die Blindheit vorbei. Auf (D) gilt

\[
H_{r_0,0}=\begin{pmatrix}0&0\\0&2\end{pmatrix}.
\]

Die ursprüngliche Ganzgitter-Vakuumrestriktion enthält die vollständigen
Cartan-Wardwerte (\langle H_{i,0}\rangle=0). Daher

\[
0=\operatorname{Tr}(\rho H_{r_0,0})=2b.
\]

Also ist (b=0); Positivität erzwingt dann (c=0), und

\[
\boxed{\rho=|0\rangle\langle0|}.
\]

Der zugleich dual- und (R_c)-invariante Zustand (I_2/2) hat dagegen
(\langle H_{r_0,0}\rangle=1). Es gibt somit keinen normalen
Grundzustand, der gleichzeitig die vollständige originale
Cartan-Vakuumrestriktion und (R_c)-Invarianz erfüllt. Diese Unvereinbarkeit
zeigt genau, welche zusätzliche Zustandsannahme jeweils gewählt wurde.

In der lokalen Lattice-VOA sind die Cartan- und Wurzelmoden Teil der Algebra.
Die VOA-Vakuumaxiome zeichnen (|0\rangle) aus. Der Zustand (|r_0\rangle)
hat gewöhnliches konformes Gewicht eins und ist ein geladener Primärzustand;
nur das fest getwistete (K_\alpha) setzt ihn zusätzlich auf Energie null.
Der vorhandene (R_c)-Operator ist deshalb folgerichtig nicht
vakuumerhaltend. Diese lokale VOA-Struktur ist eine Voraussetzung des
stipulierten E8-Abschlusses, keine aus dem bloßen Verschiebungsfunktional
gewonnene Zustandswahl.

Die neuere neutrale Clockanalyse bestätigt dieselbe Typgrenze: Die dort
erlaubte Heisenberg-/Jordan-Algebra erhält alle acht Cartanladungen und der
von (|0\rangle) erzeugte Modul bleibt (M(1)\otimes|0\rangle)
([`clock-source-twist-transport-20261015/ADDENDUM.md`, Zeilen 15–61](../../../clock-source-twist-transport-20261015/ADDENDUM.md)). Sie setzt diesen Vakuumvektor voraus; sie leitet ihn nicht aus einer dualen Symmetrie ab.

## 7. Konsequenz für den konstruktiven Suchbetrieb

Der ältere mikroskopische Quellenvergleich beweist auf der unveränderten
QWZ/CAR-Darstellung einen eindeutigen normalgeordneten Vielteilchengrundzustand
und eine strikt positive skalierte Anregungslücke. Eine gemeinsame,
Hilbertraum-, Zustands- und Zeit-erhaltende Erweiterung kann deshalb keinen
zusätzlichen orthogonalen Nullvektor liefern
([`quotient-locality-e8-20261014/PROOF.txt`, Zeilen 151–173](../../../quotient-locality-e8-20261014/PROOF.txt)).

Alle gegenwärtigen gekreuzten Defektkandidaten verwenden dasselbe
(L_s,K_\alpha). Die Phasen (\chi) und (\kappa) treten in der Energie nicht
auf. Daher haben **alle** aktuellen Kandidaten das hier klassifizierte
Grunddoublett. Die Auswahl einer Dichte durch Cartan-Wardwerte behebt den
Quellenwiderspruch nicht: Sie wählt einen Anfangszustand innerhalb (D),
entfernt aber den zweiten Nullvektor nicht und stellt die ursprüngliche
Lücke nicht wieder her.

Folglich gilt für die Suche:

```text
keiner der aktuellen Kandidaten = kompatible unveränderte QWZ/CAR-Quellenlösung
alle aktuellen Kandidaten       = bedingte algebraische Defekterweiterungen
```

Ein Reopening muss mindestens eine tragende Prämisse sichtbar ändern und
aus der Quelle begründen:

1. eine andere/zusätzliche physische Sektor- oder Quellrepräsentation, die
   den neutralen Partner tatsächlich enthält, oder
2. ein anderes physisches Zeit-/Stressgesetz beziehungsweise einen anderen
   Hilbertraumtransfer, der die ursprüngliche Gap-Aussage nicht erhalten
   soll.

Nur (R_c)-Invarianz oder nur die kanonische duale Erwartung genügt nicht.
Ein frei gesetzter Sektorenergieoffset wäre ebenfalls keine Reparatur, weil
er gerade das vorgegebene (K_\alpha) ändert. Der Wechsel zum anderen
Spinorcoset (L_c) stellt zwar die Eindeutigkeit her, verliert aber die bereits
festgehaltene zweimal-Wurzel-/Wardmarkierung (q/2\in L_s); auch das ist ein
expliziter Wechsel der Prämisse, keine Auswahl innerhalb der aktuellen
Kandidaten.

## 8. Reproduzierbare exakte Prüfung

```bash
python3 experiments/theory-contracts/seam-constructive-search-20261015/continuation/ground_state/check_ground_state.py
```

Der Rechner prüft mit rationaler Arithmetik:

- die Vorzeichen und strikten Abstände aller acht separablen Koordinaten,
- (\sum_k|\delta_k|=2) und (K_\alpha(0,0)=K_\alpha(r_0,0)=0),
- die E8-Parität von (r_0-s) sowie die gegenüberliegende Spinorlücke (1/32),
- die Antikommutation (\Gamma R_c=-R_c\Gamma),
- die dual-invariante Dichtesimplex, die bedingte Auswahl (I_2/2) durch
  zusätzliche (R_c)-Invarianz und deren Cartan-Wardwert eins,
- die identischen Verschiebungsalgebra-Funktionale verschiedener
  Grunddichten und ihre Trennung durch (H_{r_0,0}).

Der universelle Gitterbeweis ist die separable Konvexitätsabschätzung oben;
der Rechner führt keinen endlichen Ersatzscan aus.
