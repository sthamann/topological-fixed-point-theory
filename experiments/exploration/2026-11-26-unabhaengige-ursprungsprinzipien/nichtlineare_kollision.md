# Nichtlineare Wellenkollision: Der Paaransatz besteht den Strukturtest nicht

## Ergebnis

Die exakte Khan Penrose Kollision erhält \(J=0\) nicht. Zwei einzelne ebene Vakuumimpulse sind jeweils vom Petrov Typ N. Für sie gilt \(I=J=0\). Auch eine bloß lineare Addition der beiden transversalen Weyl Anteile lässt \(J=0\), weil dabei keine longitudinale Komponente entsteht. In der exakten nichtlinearen Wechselwirkungszone erzeugen die Einstein Gleichungen dagegen \(\Psi_2\). Dadurch wird \(J\) bereits an einem regulären, algebraisch exakt ausgerechneten Punkt ungleich null.

Das ist ein scharfer Negativbefund für die erhoffte Übersetzung aus Z174. Die Kollision zeigt zwar eine echte Antwort, die nur bei zwei vorhandenen Eingangswellen entsteht. Sie leitet aber nicht her, warum es überhaupt genau zwei physische Partner geben soll. Die Lösung setzt zwei beliebige charakteristische Randdaten voraus und verarbeitet sie. Sie wählt weder den zweiten Puls noch sein Amplitudenverhältnis, seine Polarisation oder einen TFPT Partnerzustand aus.

## Exakte Lösung und zwei getrennte Schalter

Halil gibt den Khan Penrose Grenzfall in den Gleichungen (51) bis (54) explizit an. In der glatten Wechselwirkungszone \(u>0,v>0\) schreiben wir

\[
p=a u,\qquad q=b v,\qquad
r=\sqrt{1-p^2},\quad w=\sqrt{1-q^2},\quad
t=\sqrt{1-p^2-q^2}.
\]

Die Konstanten \(a\) und \(b\) sind die unabhängig schaltbaren Impulsstärken. Die Metrik lautet

\[
ds^2=2e^{-M}du\,dv-e^{-U}\left(e^Vdx^2+e^{-V}dy^2\right),
\]

mit

\[
e^{-U}=t^2,\qquad
e^{-M}=\frac{t^3}{r w(pq+rw)^2},\qquad
e^V=\frac{r+q}{r-q}\frac{w+p}{w-p}.
\]

Das Einsetzen von \(p=au,q=bv\) ist kein zusätzlicher Lösungsansatz. Aus der normierten Lösung folgt es durch \(P=au,Q=bv,X=\sqrt{ab}x,Y=\sqrt{ab}y\) und eine anschließende konstante Skalierung der gesamten Metrik mit \(1/(ab)\). Eine konstante Metrikskalierung erhält die Vakuumgleichung. Direktes Nachrechnen des Ricci Tensors im Skript ergibt am Testpunkt zusätzlich exakt \(R_{ab}=0\).

Dabei ist eine wichtige physikalische Grenze sichtbar: Ein longitudinaler Boost skaliert \(a\) und \(b\) entgegengesetzt. Ihr Verhältnis ist daher kein invariantes Messdatum. Das Produkt \(ab\) bleibt invariant. Die beiden Grenzfälle \(a=0\) und \(b=0\) sind trotzdem getrennte, wohldefinierte Ausschaltkontrollen.

## Der entscheidende Test

Für die Rosen Szekeres Tetrade gilt in der linear polarisierten Lösung \(\Psi_1=\Psi_3=0\). Die drei verbleibenden Komponenten werden direkt aus \(U,V,M\) differenziert:

\[
\begin{aligned}
\Psi_0&=-\frac{e^M}{2}\left(V_{vv}-U_vV_v+M_vV_v\right),\\
\Psi_2&=\frac{e^M}{2}M_{uv},\\
\Psi_4&=-\frac{e^M}{2}\left(V_{uu}-U_uV_u+M_uV_u\right).
\end{aligned}
\]

Damit sind die tetradenunabhängigen Weyl Polynome

\[
I=\Psi_0\Psi_4+3\Psi_2^2,
\qquad
J=\Psi_2\left(\Psi_0\Psi_4-\Psi_2^2\right).
\]

Wir wählen \(a=2,b=3,u=1/4,v=1/6\). Dann ist \(p=q=1/2\) und \(t^2=1/2>0\), also liegt der Punkt regulär vor der Fokussierungsfläche. Die Rechnung ergibt exakt

\[
\Psi_0=-54\sqrt6,\qquad
\Psi_2=39\sqrt2,\qquad
\Psi_4=-24\sqrt6,
\]

und somit

\[
I=16902,\qquad
J=184626\sqrt2\ne0.
\]

Zusätzlich ist

\[
I^3-27J^2=2987837829504\ne0.
\]

Der Testpunkt ist daher sogar algebraisch allgemein, Petrov Typ I. Das Resultat beruht nicht auf einer Benennung der Lösung, sondern auf den aus der publizierten Metrik neu berechneten Invarianten.

## Gegenfälle

Erstens wird \(b=0\) oder \(a=0\) gesetzt. Abseits der distributionswertigen Wellenfront verschwindet die entsprechende Einzelwellenkrümmung in der hier getesteten Zone exakt. Auf der Front ist der Impuls vom Typ N, also bleiben seine skalaren Weyl Polynome ebenfalls null. Die Produkte distributionswertiger Größen werden dabei ausdrücklich nicht als gewöhnliche Funktionen behandelt.

Zweitens wird als lineare Negativkontrolle \(\Psi_0=9,\Psi_4=4,\Psi_2=0\) gesetzt. Dann sind beide gegenläufigen Typ N Anteile vorhanden und \(I=36\), aber \(J=0\). Genau diesen Befund kippt die nichtlineare Einstein Wechselwirkung durch das neu erzeugte \(\Psi_2\).

Die Amplitudenstruktur zeigt zugleich die Grenze des Ergebnisses. Wegen \(\partial_u=a\partial_p\) und \(\partial_v=b\partial_q\) skalieren

\[
\Psi_0\sim b^2,\qquad \Psi_2\sim ab,\qquad \Psi_4\sim a^2,
\qquad J\sim a^3b^3
\]

bei festem \(p,q\). Die Antwort hat also eine saubere UND Struktur. Sie verschwindet, wenn ein Puls fehlt. Aus diesem Produktgesetz folgt aber keine Auswahl eines bestimmten Paars.

## Abgrenzung zu Z174

Der Atlas beschreibt Z174 als endliches Standardformmodell aus Naht und Tomita Partner. Er hält ausdrücklich fest, dass \(P_{\mathrm{std}}\) zwei Prämissen zusammenlegt, der Oszillatoranteil nur als Phasenarithmetik geprüft ist und \(P_{\mathrm{paar}}\), also warum genau zwei Nähte summiert werden, offen bleibt.

Die Khan Penrose Raumzeit schließt diese Lücke nicht. Ihre zwei Wellen sind charakteristische Anfangsdaten. Die Feldgleichungen liefern zu diesen Daten eine eindeutige nichtlineare Wechselwirkungszone. Sie erzeugen den zweiten Eingang nicht und unterscheiden keinen kanonischen Partner unter mehreren möglichen Eingängen. Außerdem ist das hier geprüfte \(J\) der kubische Weyl Invariant der vierdimensionalen Raumzeit. Eine Gleichsetzung mit einem gleich bezeichneten internen oder modularen TFPT Objekt wäre eine neue, bislang unbegründete Identifikation.

Damit ist der stärkste Einwand zweigeteilt. Der gewünschte Nullinvariant überlebt die echte nichtlineare Physik nicht. Und selbst die robuste bilineare Paarantwort ist eine Antwort auf vorgegebene zwei Randdaten, keine physische Herleitung der Paarzahl.

## Reproduktion

```bash
python3 experiments/exploration/2026-11-26-unabhaengige-ursprungsprinzipien/nichtlineare_kollision.py
```

Das Skript prüft die exakten Weyl Werte, \(I\), \(J\), die Petrov Diskriminante, beide Ausschaltkontrollen, den linearen Gegenfall und den Ricci Tensor am selben Punkt.

## Primärquellen

1. K. A. Khan und R. Penrose, „Scattering of Two Impulsive Gravitational Plane Waves“, *Nature* 229 (1971), 185 bis 186, DOI [10.1038/229185a0](https://doi.org/10.1038/229185a0). Originalarbeit zur exakten Zweipulskollision.
2. M. Halil, „Colliding plane gravitational waves“, *Journal of Mathematical Physics* 20 (1979), 120 bis 125, DOI [10.1063/1.523951](https://doi.org/10.1063/1.523951), [offenes PDF](https://i-rep.emu.edu.tr/xmlui/bitstream/handle/11129/2067/Journal%20of%20Mathematical%20Physics%20Volume%2020%20issue%201%201979%20%5Bdoi%2010.1063%252F1.523951%5D%20Halil%2C%20Mustafa%20--%20Colliding%20plane%20gravitational%20waves.pdf?sequence=1). Gleichungen (51) bis (54) liefern die hier eingesetzte Khan Penrose Metrik, Gleichung (42) die verwendete Nulltetrade, der Anhang bestätigt \(\Psi_1=\Psi_3=0\).
3. E. Newman und R. Penrose, „An Approach to Gravitational Radiation by a Method of Spin Coefficients“, *Journal of Mathematical Physics* 3 (1962), 566 bis 578, DOI [10.1063/1.1724257](https://doi.org/10.1063/1.1724257). Primärquelle des eingesetzten Spin Koeffizienten Formalismus.
4. P. Szekeres, „Colliding Plane Gravitational Waves“, *Journal of Mathematical Physics* 13 (1972), 286 bis 294, DOI [10.1063/1.1665972](https://doi.org/10.1063/1.1665972). Allgemeine exakte Lösung für kollidierende linear polarisierte ebene Vakuumwellen.

Die Webquellen wurden nur für die externe Physik verwendet. Der Atlas wurde gelesen, aber nicht verändert.
