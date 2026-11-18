# Konstruktive Suchdomäne: markierte CAR–Defekt-Erweiterung

## Urteil vorab

Die erste ehrliche konstruktive Aufgabe ist **nicht**, einen neuen Hamiltonoperator, eine Gibbs-Dichte oder ein Messinstrument zu optimieren. Sie ist kleiner und näher an der offenen Nahtstelle:

> Erweitere die tatsächlich gegebene achtkanalige CAR-Ladungsalgebra mit ihrem geordneten Klein-Cocycle um den Halbsektor (s=(1/2)^8), so dass auch die ungeraden CAR-Felder als Defekt-Endpunkte wirken und zugleich die bereits feste Viertelphase, Zeitwirkung, Ladung, Parität und der Operator (W_4) erhalten bleiben.

Das Bild ist ein Reißverschluss mit zwei schon festen Seiten. Auf der einen Seite liegen alle ganzzahligen CAR-Ladungen, einschließlich der ungeraden Fermionfelder. Auf der anderen liegt der bereits konstruierte gerade (D_8/E_8)-Anschluss. Gesucht wird nur der Zahn, der beide Seiten als **Defektalgebra** zusammenhält. Er darf die ungeraden Felder nicht nachträglich zu gewöhnlich lokalen (E_8)-Feldern erklären.

Der bekannte Kandidat

\[
 X E_n X^{-1}=(-i)^{Q(n)}E_n,\qquad X^2=-E_q
\]

zeigt, dass die Domäne nicht leer ist. Der neue Gehalt liegt in der vollständigen Wortalgebra auf geraden **und** ungeraden Zuständen, ihrer Markierung durch (P, Z_a, U_s,R_c,W_4), der unveränderten Zeitwirkung sowie der noch zu zertifizierenden induzierten Zustandsauswertung. Die bekannte gerade Faktorisierung allein zählt nicht als neues Ergebnis.

## Warum diese Domäne die richtige erste Grenze trifft

Der Originalbeweis fixiert auf dem gesamten Roh-CAR-Raum

\[
 \epsilon_{\mathbb Z}(n,m)=(-1)^{\sum_{i>j}n_i m_j},\qquad
 W_4=Z_aE_aP_{\rm CAR}.
\]

Er konstruiert außerdem auf dem **geraden** (D_8)-Sektor einen Halbsektorgenerator (X=E_s), die beiden (P)-Liftzweige und (U_sR_c=W_4). Derselbe Beweis sagt ausdrücklich, was fehlt: Die ungeraden Vektorfelder brauchen einen Spin-/Defektmodultransport; eine gewöhnlich lokale CAR/(E_8)-Inklusion wurde nicht konstruiert. Genau dort setzt die Suche an. Sie wiederholt weder den bekannten geraden Beweis noch den bereits ausgeschlossenen ordinary-local Ansatz.

Ein Dimensionsvergleich (128\leftrightarrow16) spielt keine Rolle. Die Suchdomäne verwendet die tatsächliche Gruppenwirkung, den tatsächlichen Cocycle und die tatsächlichen Ladungen. Der frühere (V_{16})-Pullback bleibt eine bedingte inverse Spinwirkung und wird hier nicht als Rohquellen-Wörterbuch ausgegeben.

## Feste Algebra und einziger Suchgegenstand

Für (n,m\in\mathbb Z^8) sind die ganzzahligen Wörter bereits fest:

\[
 E_nE_m=\epsilon_{\mathbb Z}(n,m)E_{n+m}.
\]

Setze

\[
 q=(1,\ldots,1),\quad s=q/2,\quad
 a=(1,1,1,1,0,0,0,0),\quad b=Pa=(0,0,0,0,1,1,1,1),
\]

\[
 r_0=s-b=(1/2,1/2,1/2,1/2,-1/2,-1/2,-1/2,-1/2).
\]

Ein Kandidat enthält genau zehn Exponenten in \(\mathbb Z/4\mathbb Z\):

```json
{
  "schema": "tfpt.crossed_defect.v1",
  "chi_exp": [3, 3, 3, 3, 3, 3, 3, 3],
  "zeta_exp": 2,
  "kappa_exp": 2
}
```

Dabei bedeutet

\[
 X E_{e_j}X^{-1}=i^{\lambda_j}E_{e_j},\qquad
 X^2=i^zE_q,\qquad
 \widehat P X\widehat P^{-1}=i^kX.
\]

Der Kandidat enthält **keine** Zielmatrix, keinen Hamiltonoperator, keinen Zustand und kein frei gewähltes (W_4). Für beliebige ganzzahlige Ladungen ist

\[
 L(n)=\sum_j\lambda_jn_j\pmod4,\qquad
 X E_n=i^{L(n)}E_nX.
\]

Jedes Wort besitzt danach die Normalform (E_nX^c) mit (c\in\{0,1\}). Das ist eine gekreuzte, \(\mathbb Z_2\)-graduierte Defekt-Erweiterung. Sie behauptet keine gewöhnliche Lokalität.

## Die kleinsten entscheidenden Gleichungen

Die Suchprüfung arbeitet mit exakten Kongruenzen, nicht mit angenäherten Matrixnormen.

1. **Das Quadrat des Defekts muss die tatsächliche CAR-Innenwirkung liefern.** Aus dem geordneten Cocycle folgt

   \[
   \frac{\epsilon(q,n)}{\epsilon(n,q)}=(-1)^{Q(n)}.
   \]

   Deshalb muss für jeden Generator

   \[
   2\lambda_j=2\pmod4
   \]

   gelten. Ein lokal kommutierendes (X) mit \(\lambda_j=0\) scheitert bereits hier.

2. **Das Quadrat muss assoziativ mit (E_q) sein.** Da (X^2=i^zE_q), ist

   \[
   \sum_j\lambda_j=0\pmod4
   \]

   erforderlich. Dann kommutiert (X) mit (E_q).

3. **Die bekannte gerade Konstruktion wird als Quelle, nicht als Zielwert, benutzt.** Auf (D_8=\{d:Q(d)\text{ gerade}\}) gilt bereits

   \[
   X E_dX^{-1}=(-1)^{Q(d)/2}E_d.
   \]

   Auf den Wurzeln reduziert das auf

   \[
   \lambda_i-\lambda_j=0,\qquad \lambda_i+\lambda_j=2\pmod4.
   \]

   Damit bleiben nur die beiden Orientierungen \(\lambda_j=1\) und \(\lambda_j=3\). Die Suche prüft, ob eine **bereits feste** Quellidentität die Originalorientierung \(-i\) tatsächlich auswählt. Solange eine solche Brücke fehlt, sind beide Orientierungen zulässig und die Auswahl bleibt offen; der Checker darf \(\lambda_j=3\) nicht vorab als Erfolg fest eintragen.

4. **Der Permutationslift bleibt derselbe.** Für den Halbtausch (P) muss \(\lambda_j=\lambda_{j+4}\) gelten. Aus \(\widehat P^2=1\) folgt (2k=0\pmod4). Die Zweige (k=0,2) bleiben getrennt sichtbar, weil der gerade Quellbeweis sie nicht auswählt. (k=2) ist die dort verwendete Konvention, keine aus P1/P2 hergeleitete Auswahl.

5. **Die volle Operatorfaktorisierung wird als Wortgleichheit geprüft.** Definiere

   \[
   E_{r_0}:=E_{-b}X,\quad
   R_c:=Z_aE_{r_0}\widehat P,\quad
   U_s:=i^{-z}X,\quad
   W_4:=Z_aE_a\widehat P.
   \]

   Der Faktor (i^{-z}) in (U_s) ist durch die gewählte Sektion (X^2=i^zE_q) bestimmt. Er ist kein Suchparameter. Weil \(\epsilon(-b,q)=1\), gilt

   \[
   XE_{r_0}=i^{L(-b)+z}E_a.
   \]

   Daher ist

   \[
   U_sR_c=W_4\quad\Longleftrightarrow\quad L(-b)=0\pmod4.
   \]

   Diese Gleichung gilt als Algebraelement und damit auf beiden Graden. Ein bloßes (R_c=W_4) scheitert bereits an der Ladung: R_c ist neutral, W_4 verschiebt Q um vier. Die Parität allein unterscheidet diese beiden geraden Operatoren nicht.

6. **Alle Wörter, nicht nur zehn Generatorwerte, müssen assoziativ sein.** Lean reduziert beide Klammerungen von

   \[
   (E_nX^c)(E_mX^d)(E_\ell X^f),\qquad c,d,f\in\{0,1\},
   \]

   symbolisch auf Ladung, Grad und Phasenexponent. Die Aussage bleibt universell in (n,m,\ell\in\mathbb Z^8).

Die unbekannte odd-sector Wirkung wird dadurch nicht in eine Annahme geschrieben. Sie folgt, wenn sie folgt, aus einem einzigen Defektgenerator, dessen Quadrat, CAR-Wirkung und Markierungen bereits fest miteinander verträglich sein müssen.

### Analytische Vorklassifikation

Die Suche kennt ihr erwartetes Ende schon algebraisch. Aus \(\lambda_i-\lambda_j=0\) sind alle acht Exponenten gleich. Aus \(2\lambda_j=2\) sind sie ungerade. Daher bleiben exakt \(\lambda=1\) und \(\lambda=3\). Aus \(2k=0\) bleiben \(k=0\) und \(k=2\). Somit sind höchstens vier markierte algebraische Kandidatenklassen möglich.

Die Sektionsphase \(z\) ist Gauge. Schreibe \(t=\exp(\pi i g/4)\) und \(X'=tX\). Dann gilt

\[
 z'=z+g\pmod4,\qquad E_{r_0}'=tE_{r_0},\qquad R_c'=tR_c.
\]

Wegen der festen Normalisierung \(U_s=i^{-z}X\) transformiert gleichzeitig

\[
 U_s'=i^{-z'}X'=t^{-1}U_s.
\]

Damit bleibt \(U_s'R_c'=U_sR_c=W_4\) exakt invariant. Weder \(\chi\) noch \(\kappa\) ändern sich. Deshalb darf \(z\) auf null kanonisiert und nicht als neue Klasse gezählt werden.

Der GPU-Lauf ist deshalb eine vollständige ausführbare Gegenprüfung dieser Klassifikation und der implementierten Zeitidentitäten; die vollständige Zustandsprüfung bleibt getrennt. Er ist keine blinde Hoffnungssuche. Lean muss anschließend beweisen, dass diese Reduktion für alle ganzzahligen Ladungen und alle Wörter gilt.

## Statistik: absichtlich Defekt statt ordinary local

Für den kalibrierten Ast \(\lambda_j=3\) gilt

\[
 XE_{e_j}X^{-1}=-iE_{e_j}.
\]

Das ist gerade die benötigte gemischte Defektstatistik. Zweimaliges Überkreuzen ergibt (-1), also die CAR-Parität. Die ungeraden Felder bleiben Endpunkte eines Spin-/Defektsektors. Ein Kandidat mit gegenseitiger Phase (+1) wird nicht als „einfachere Lösung“ akzeptiert, sondern verletzt die feste Innenwirkung von (E_q).

## Viertelphase und Zeit werden nicht neu gewählt

Die Zeitfunktion bleibt die bereits vorgegebene quadratische Funktion

\[
 K_\alpha(p,N)=N+\frac{|p|^2}{2}+\delta_\alpha\!\cdot p,
 \qquad (\delta_\alpha)_j=\frac{j+\alpha}{8}-\frac12,
 \]

für \(\alpha\in\{1/4,3/4\}\). Neu ist nur, dass dieselbe Formel auf allen Ladungen (p=n+cs) ausgewertet wird. Das ist eine erklärte Erweiterungsvoraussetzung; es ist keine Herleitung dieser Zeit aus P1/P2.

Die Prüfung verlangt weiterhin

\[
 [Q,W_4]=4W_4,
 \qquad
 [K_\alpha,W_4]=W_4\left(\frac Q2+\frac34+\frac\alpha2\right).
\]

Die echte Viertelmarkierung lebt in \(\mu_8\), obwohl der kleine Cocycle-Genotyp nur \(\mu_4\) benötigt. Ihre festen Exponenten sind

\[
 7+4(Q\bmod2)\pmod8\quad(\alpha=1/4),
 \]

\[
 1+4(Q\bmod2)\pmod8\quad(\alpha=3/4).
\]

Damit kann weder ein Kandidat \(\alpha=0\) einsetzen noch die ganzzahlige Ladung auf vier zyklische Zustände reduzieren. Der Carry (T_1^8=T^2) und die ganzzahlige Ladung bleiben erhalten.

Diese Prüfung hält die Viertelphase fest, behauptet aber keine noch unbewiesene Gleichung zwischen ihr und dem Vorzeichen von \(\chi\). Überleben \(+i\) und \(-i\), ist das Ergebnis eine exakte Zweideutigkeit der Markierung und kein Grund, den kalibrierten \(-i\)-Ast automatisch „physikalisch“ zu nennen.

### Was die bekannten Originalmarkierungen tatsächlich unterscheiden

Die zitierten Originalidentitäten wählen weder \(\chi=+i\) gegen \(\chi=-i\) noch \(\kappa=+1\) gegen \(\kappa=-1\) aus:

- Beide \(\chi\)-Orientierungen haben dasselbe Quadrat \(\chi(n)^2=(-1)^{Q(n)}\), dieselbe Restriktion auf (D_8) und dieselbe Ladungs-/Paritätswirkung. Der feste (W_4)-Viertelcharakter enthält \((-1)^Q\), also genau das gemeinsame Quadrat, nicht die einzelne Phase \(\pm i\).
- (K_\alpha\) und seine Wardidentität hängen von Ladung und Modenverschiebung ab, nicht vom Vorzeichen der einfachen Defektüberkreuzung. Daher bleiben die geprüften Zeitwerte in beiden \(\chi\)-Ästen identisch.
- S-, T- und Arf-Daten sehen die Sektorklasse und die doppelte Überkreuzung; auch sie bleiben unter \(\chi\mapsto\overline\chi\) unverändert. Die kanonische Erwartung hat in beiden Fällen dieselbe Ganzgitterrestriktion und löscht den ungeraden Defektgrad.
- Der Originalbeweis sagt ausdrücklich, dass beide \(\kappa\)-Lifts zulässig sind. Auf der ganzzahligen CAR-Algebra wirken sie identisch, und die gegenphasige Normalisierung von \(U_s,R_c\) hält ihr Produkt (W_4) fest. \(\kappa\) erscheint erst beim Überkreuzen von \(\widehat P\) mit dem zusätzlichen (X); dafür enthält die vorhandene Quelle keinen fest ausgewerteten gemischten Marker.

Damit sind die vier Ausgaben algebraisch verschiedene markierte Kandidaten innerhalb der gewählten Präsentation. Eine physikalische Auswahl oder physikalische Inequivalenz folgt aus den bekannten Wardmarkierungen nicht.

## Noch zu zertifizierende Zustandsfortsetzung

Der Zustand wird nicht gesucht. Sei \(\omega_0\) der bereits feste positive CAR-Zustand auf dem Ganzgitterteil. Die Gradprojektion

\[
 \mathcal E_0(E_nX^c)=
 \begin{cases}
 E_n,&c=0,\\
 0,&c=1
 \end{cases}
\]

liefert die induzierte Auswertung

\[
 \widetilde\omega=\omega_0\circ\mathcal E_0.
\]

Für ein beliebiges endliches Wort in Normalform (y=a_0+a_1X) ist der Grad-null-Anteil

\[
 \mathcal E_0(y^*y)
 =a_0^*a_0+\phi^{-1}(a_1^*a_1),
 \qquad \phi(E_n)=i^{L(n)}E_n.
\]

Das ist eine exakte Summe positiver Quadrate im festen Ausgangssystem. Die Suchfunktion optimiert also keinen gewünschten Gram-Kern. Die noch ausstehende C11-Zertifizierung muss prüfen, ob der Kandidat diese positive Struktur transportiert; der aktuelle GPU-Filter prüft das nicht. Die Aussage betrifft endliche algebraische Wörter. Sie behauptet weder eine neue (C^*)-Vervollständigung noch eine Gibbs-Auswahl. Die Erweiterung ist nur in der dual-(\mathbb Z_2)-invarianten Klasse mit fester Restriktion kanonisch eindeutig.

Die vier Zustandssektoren werden dabei nicht auf zwei Grade verkürzt. In der Reihenfolge \([0,s,v,c]\) entsprechen sie

\[
 (c,Q(n)\bmod2)=(0,0),(1,0),(0,1),(1,1).
\]

Der Defektgenerator \(X\) vertauscht \(0\leftrightarrow s\) und \(v\leftrightarrow c\); ein ungerades \(E_{e_j}\) vertauscht \(0\leftrightarrow v\) und \(s\leftrightarrow c\). Die geplante vollständige C12-Prüfung soll außerdem die festen S-/T-Daten und die Arf-Rekombination

\[
 Z_F[\rho_x,\rho_t]=\frac12\sum_{a,b\in\mathbb F_2}
 (-1)^{ab+a\rho_t+b\rho_x}Z_B[a,b]
\]

unverändert erhalten. Bisher sind nur die Sektortoggles, T-Vorzeichen und die fest definierte Arf-Formel zertifiziert; das beweist weder die vollständige Zustandskompatibilität noch eine Auswahl von \(\omega_0\).

## GPU-Aufgabe und warum sie trotzdem klein bleiben darf

Das Eingabeformat erlaubt zehn rohe \(\mu_4\)-Exponenten. Für die ausführbare Suche wird die reine Sektionsphase auf \(z=0\) gesetzt; alle vier Werte von \(k\) werden zunächst mitgeführt. Damit werden \(4^9=262{,}144\) Zeilen ausgewertet. Die Gleichung \(2k=0\) verwirft die ungeraden Werte als eingebaute Negativfälle und lässt vor den übrigen Bedingungen \(4^8\cdot2=131{,}072\) quellgültige Zeilen. Nach allen Quellrelationen sollten höchstens vier algebraische Kandidatenklassen übrig bleiben: zwei \(\chi\)-Orientierungen mal zwei \(P\)-Liftzweige. Eine physikalische Inequivalenz dieser vier Klassen ist nicht bewiesen. Die GPU klassifiziert sie reproduzierbar auf zwei H200 und protokolliert Gerät, UUID, Datentyp, Batchgröße, synchronisierte Laufzeit und Spitzenspeicher.

Die Zeitfixtures verwenden in der tatsächlichen Implementierung die ganzzahlige Skalierung `energy64`: \(\delta\) wird mit Nenner 32 und die Ladung verdoppelt dargestellt, die vollständige Energie-/Wardgleichung ist mit 64 skaliert. Es gibt hier keine Nenner-16-Prüfung.

Das ist kein Geschwindigkeitsengpass. Sobald die vier Klassen analytisch klassifiziert und in Lean zertifiziert sind, endet dieser Lauf; weitere Millionen künstlicher Varianten wären keine zusätzliche Forschung.

Die Domäne wird nicht künstlich mit quadratischen Phasenkoeffizienten vergrößert. Sobald \(X\) durch Konjugation einen Algebraautomorphismus liefern soll, muss \(\chi(n+m)=\chi(n)\chi(m)\) gelten; scheinbar allgemeinere quadratische Terme kollabieren deshalb wieder auf den linearen Charakter oder scheitern exakt. Mehr GPU-Arbeit wäre hier kein größerer Erkenntniswert.

Die vier implementierten Zielfunktionen zählen Verletzungen der Phasen-/Assoziatorrelationen, der D8-/P-Verträglichkeit, der UsRc-Komposition und der festen Zeit-/Viertelproben. Die Wort-SOS-Identität ist eine anschließende Beweisaufgabe, keine bereits ausgewertete Fitnessfunktion.

Ein numerischer Nullrest ist kein Beweis. Ein Kandidat wird erst nach exakter Wiederholung und Lean-Zertifikat als konstruiert markiert.

**Implementierungsgrenze am 15.10.2026:** Der GPU-Filter prüft aktuell die Phasenrelationen, 200 Normalform-Assoziator-Koeffizienten und 36 feste Zeit-/Viertelproben in vier Zielfunktionsgruppen. Sein Flag bedeutet deshalb nur `provisional_feasible`. Der universelle Allwort-Beweis ist inzwischen im festen Lean-Zertifikat implementiert und für alle vier Klassen geprüft. Ebenfalls bewiesen sind der affine Ladungscarry und die fest vorgegebene rationale Wardwirkung. Für das volle Forschungs-Verdict fehlen die vollständige Paritätsoperatorwirkung, die Wort-SOS-Aussage und die vollständige Viersektoren-S-/T-/Arf-Prüfung. Diese Punkte sind Zertifikatspflichten, keine stillschweigenden Annahmen.

## Lean-Zertifikat

Der Checker akzeptiert ausschließlich das strikte Vier-Schlüssel-JSON. Kandidaten dürfen keine Annahmen oder Zielbehauptungen mitsenden. Die noch nicht vollständig implementierte Zielaussage hat die Form

```text
check g = true ->
  CrossedWordAssociative g /\
  D8Restriction g /\
  DefectOddBraiding g /\
  UsRcEqW4 g /\
  ChargeParity g /\
  QuarterTimeWard g /\
  CanonicalExpectationSOS g /\
  AllSectorArfCompatibility g
```

Die Universalaussagen über Ladungen werden auf ganzzahlige Exponentenidentitäten zurückgeführt; nur die drei Gradbits werden endlich ausgeschöpft. Es werden kein Fließkomma-Premissensatz, kein frei formulierter Kandidatensatz und kein rational erfundenes (Q^{-1/2}) oder \(-\log Q\) zugelassen.

## Pflicht-Gegenproben

Die folgenden Fehler müssen sichtbar und aus dem richtigen Grund scheitern:

- `ordinary_local_X`: \(\lambda_j=0\) verletzt (X^2={\rm Ad}(E_q)), die gerade Restriktion und die Defektstatistik.
- `drop_half_transport`: (R_c\) allein besitzt nicht die Ladung und Wirkung von (W_4).
- `wrong_U_normalization`: Eine Änderung von (X^2=i^zE_q) bei fälschlich festgehaltenem (U_s=-X) verletzt die Operatorwortgleichheit. Mit der vorgeschriebenen Regel (U_s=i^{-z}X) ist dieselbe Änderung nur ein Sektions-/Gaugewechsel und keine neue Lösung.
- `erase_quarter` oder zyklisches Reduzieren der Ladung löscht die feste \(\mu_8\)-Markierung beziehungsweise den Carry und verletzt die Zeit-/Ladungsidentität.
- `flip_twisted_Arf_sign` verändert den festen twisted-sector-Beitrag; diese Negativkontrolle gehört zur noch ausstehenden vollständigen C12-Prüfung.

## Aussagekraft eines Ergebnisses

Ein Erfolg trägt genau das Verdict `CONDITIONAL_MARKED_DEFECT_LIFT`: Unter der schon gewählten achtkanaligen CAR-Quelle, ihrem Cocycle, dem Halbsektor und derselben Zeit-/Zustandsklasse existiert eine kohärente volle Defektalgebra, in der (U_sR_c=W_4) auch die ungeraden Felder korrekt mitführt.

Das wäre ein echter Operatorfortschritt gegenüber der bisherigen even-(D_8)-Faktorisierung. Es würde weiterhin nicht erklären, warum P1/P2 gerade diese Quelle, diesen Halbsektor, diese Zeit oder diesen Zustand auswählt. Auch eine vollständige Klassifikation kann den (k=0/2)-Ast nur dann physikalisch auswählen, wenn eine zusätzliche **bereits aus der Quelle stammende** Markierungsidentität einen Ast ausschließt.

Ein negativer Ausgang schließt nur diese gekreuzte Ein-Generator-Domäne aus. Er ist kein allgemeines No-Go gegen Spin-/Defektmodule.

## Verbindliche Originalstellen

- `seam-joint-ward-state-20261014/PROOF.txt`, Zeilen 752–861: (T_1,W_4), Ladung, Parität, (K_\alpha)-Wardidentität und Viertelcharakter.
- dieselbe Datei, Zeilen 1003–1087: tatsächlicher geordneter CAR-Cocycle, (P_{\rm CAR}), (Z_aE_aP_{\rm CAR}), gerade (D_8/E_8)-Erweiterung, (\kappa)-Zweige und die explizite odd-sector-Grenze.
- `charged-source-time-audit-20260920/PROOF.txt`, Zeilen 154–290: geladene Zeitdetektoren, Halbspinorfeld und die Grenze neutraler Stromtests.
- `source-sector-reconstruction-20260920/PROOF.txt`, Zeilen 1–293: vier Sektoren, Arf-Rekombination, odd fields als Defekt-Endpunkte und gemeinsame Zeitvoraussetzung.
- `v16-operator-pullback-20261014/PROOF.md`: keine dimensionsbasierte Rohquellen-Abbildung.
- `native-source-reconstruction-20261014/PROOF.md`: Zustand und Zeit bleiben bedingte Quellvoraussetzungen; ihre Auswahl ist nicht durch diese Konstruktion gelöst.

Die zugehörigen SHA-256-Werte und maschinenlesbaren Grenzen stehen in `domain_spec.json`.
