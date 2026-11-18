# TFPT als markierte Reflexionsgeometrie
## Vierfach punktierte Sphäre, Viertelclock und zentrale Liftobstruktion

**Datum:** 23. Oktober 2026  
**Art:** Neue Definitionen, explizite Ableitungen und vollständige endliche Klassifikationen in einer angegebenen Klasse.  
**Originalbasis:** `sthamann/tfpt`, Commit `7e7fefa42f0e26cc4229c12b35543fe34ad56f18`.

## 0. Ergebnis und genaue Reichweite

Wir untersuchen die ursprüngliche 60er Reflexionsklasse in ihrer tatsächlichen vierdimensionalen Darstellung. Für vier geordnete Punkturen und Produktrelation eins ergeben sich:

* 14.760 geordnete Reflexionsvierer;
* 17 Klassen unter simultaner Konjugation durch die volle Matrixgruppe;
* fünf Bahnen unter den zusätzlich zugelassenen Hurwitzoperationen;
* 62 Klassen bei zusätzlich gewählter, strikt äquivarianter Faserclock **exakter Ordnung vier**; ohne Forderung exakter Ordnung sind es 80;
* ein nichttrivialer, genau zweitorsionärer zentraler Erweiterungskokykel;
* 64 zentrale Lifts jeder festen ungerahmten Basismonodromie, wenn die vier peripheren zentralen Phasen frei sind;
* kein Lift, der zugleich die festgelegte zyklische Viertelclockstruktur streng erhält;
* zwei nichtkonjugierte Einbettungstypen von C2 × C2, die durch integrale lokale Kohomologie unterschieden werden: nur einer besitzt Z/2 Torsion in H¹.

Die endlichen Klassifikationen sind neu in dieser Ausarbeitung berechnet. Eine Prioritätsbehauptung gegenüber der gesamten Literatur wird nicht erhoben. Die bekannte Konstruktion der Gruppe G31 und die allgemeine Theorie zentraler Erweiterungen sind Ausgangswerkzeuge, keine neuen Sätze dieser Ausarbeitung.

Die Ergebnisse gelten **nicht** für beliebige Konjugationsklassen in W(E8), beliebige Randwindungen oder sämtliche denkbaren Quadratikmarkierungen. Der erste abgeschlossene Modellfall verwendet ausdrücklich vier ursprüngliche komplexe Reflexionen in G31. Eine veränderte Klasse ist ein anderes, klar anschließbares Klassifikationsproblem.

## 1. Quellen und neue Definitionen auseinanderhalten

### 1.1 Aus den Originalen übernommene Daten

[S1] gibt 240 Gaussian E8 Wurzeln in verdoppelten reellen Koordinaten und ihre 60 komplexen Linien an. Auf V = C⁴ lauten die Reflexionen

\[
r_z=I-\frac12zz^*,\qquad z^*z=4.
\]

Dies ist ein hermitescher Operator mit Eigenwerten −1,1,1,1. Die erzeugte Gruppe heißt im Original G31. Ihre Ordnung 46.080, ihr skalares Zentrum μ4 und die projektive Ordnung 11.520 werden hier durch vollständige Enumeration erneut bestätigt. Die Darstellung von G31 über die zweite äußere Potenz ist auch in [L1] untersucht.

[S2] gibt den fünfdimensionalen Quartikraum Q ⊂ Sym⁴(V), seine invariante Metrik und die Darstellung R(g) an. Es gilt

\[
\det R(g)=\det g=(-1)^{\epsilon(g)}.
\]

Für die dort festgelegte Sektion L und das dortige zentrale Element z von Ordnung vier gilt

\[
L(g)L(h)=z^{-\epsilon(g)\epsilon(h)}L(gh).
\]

Zentralität ist hier relativ zum angegebenen Liftunterverband zu verstehen; z ist nicht als Zentrum der gesamten E8 Gruppe zu lesen.

### 1.2 Neuer präziser Ausgangsgegenstand

Um die entfernten Punkte nicht zugleich als Punkte des bereits punktierten Raums zu behandeln, verwenden wir

\[
(\bar X,\Sigma),\qquad U=\bar X\setminus\Sigma.
\]

Im ersten Modell ist \(\bar X=\mathbb{CP}^1\) und \(\Sigma=(1,i,-1,-i)\) geordnet. Ein Objekt der Kategorie markierter flacher Reflexionsgeometrien besteht aus:

1. einer orientierten kompakten Fläche mit geordneten Punkturen und gegebenenfalls markiertem Basispunkt;
2. einem lokalen System \(\mathcal L\) positiver gerader Gitter mit festem Fasertyp \(\Lambda\) und einer Reduktion seiner Strukturgruppe auf eine angegebene Untergruppe \(\Gamma\leq O(\Lambda)\);
3. vorgegebenen peripheren Konjugationsklassen \(\mathcal C_j\subseteq\Gamma\);
4. gegebenenfalls einer zusätzlichen Quadratikmarkierung und einer angehobenen Transportstruktur.

Ein lokales System ist äquivalent zu einem Transportfunktor aus dem Fundamentalgruppoid. Nach Wahl eines Faserrahmens ist seine Monodromie \(\rho:\pi_1(U)\to\Gamma\). **Lokales System, Monodromie und flacher Transport sind daher nicht drei unabhängige Eingaben.**

Für das ursprüngliche gerade Gitter kann man intrinsisch

\[
q_{\rm lat}(\bar v)=\frac{\langle v,v\rangle}{2}\pmod2,
\qquad \bar v\in\Lambda/2\Lambda
\]

verwenden. Die Wohlbestimmtheit folgt aus Geradheit und Integrität. Diese Form ist automatisch unter O(Λ) invariant. Eine der früher betrachteten zusätzlichen quadratischen Verfeinerungen auf einem anderen endlichen Quotienten wäre eine **separate** Markierung. Soll sie global parallel sein, muss die Monodromie ihren Stabilisator treffen. Ohne diese Festlegung dürfen beide Quadratikbegriffe nicht gleichgesetzt werden.

Ein Morphismus besteht aus einer Abbildung der markierten Basen, welche die festgelegten peripheren Klassen erhält, und einer isometrischen Isomorphie zum zurückgezogenen lokalen System. Zusätzliche Quadratik- und Liftdaten müssen mittransportiert werden. Bei festem Basismodell werden für die untenstehenden Modulräume nur Morphismen über der Identität verwendet. Das ist genau simultane Konjugation der Monodromien.

Homotopien werden relativ zu den markierten Daten betrachtet. Die induzierten Transportfunktoren sind dann natürlich isomorph. Eine zusätzliche Quotientierung nach Abbildungsklassen oder Hurwitzbewegungen wird ausdrücklich gesondert angegeben.

### 1.3 Produkte und Nahtverklebung

Disjunkte Vereinigung definiert eine monoidale Operation. Das äußere Produkt zweier Gitterlokalsysteme kann über einem kartesischen Produkt mit Gitter \(\Lambda_1\oplus\Lambda_2\) gebildet werden. Es verlässt im Allgemeinen die Kategorie der Flächen mit festem E8 Fasertyp: E8 ⊕ E8 ist nicht E8. Es wird daher **nicht** als kategoriales Produkt dieser engen Klasse ausgegeben.

Für die Verklebung ersetzt man Punkturen durch kleine Randkreise. Eine orientierungsumkehrende Randidentifikation muss die beiden lokalen Systeme durch eine Isometrie identifizieren. In Rahmenform lautet die Bedingung

\[
M_+=uM_-^{-1}u^{-1}.
\]

Für angehobene Systeme muss dieselbe Gleichung einschließlich der zentralen Faktoren gelten. Die universelle Transportzusammensetzung folgt aus dem van Kampen Pushout der Fundamentalgruppoide. Ohne Randidentifikation ist die Verklebung noch nicht spezifiziert.

## 2. Das endliche Monodromieproblem

Wähle positiv orientierte Meridiane \(m_1,\ldots,m_4\) mit

\[
\pi_1(U)=\langle m_1,m_2,m_3,m_4\mid m_1m_2m_3m_4=1\rangle\cong F_3.
\]

Es sei \(\Gamma=G_{31}\) in den ursprünglichen Matrizen und \(\mathcal R\) die 60er Reflexionsklasse. Das erste Problem lautet

\[
\mathcal M_{\mathcal R}
=
\{(r_1,r_2,r_3,r_4)\in\mathcal R^4:r_1r_2r_3r_4=I\}/\Gamma.
\]

Die vier Punkturen bleiben geordnet. Eine Surjektivität der Monodromie wird nicht vorausgesetzt.

### Satz 1: Jede solche Monodromie besitzt mindestens zwei gemeinsame Fixrichtungen

Für vier unitäre Reflexionen vom Typ (−1,1,1,1) mit Produkt eins gilt

\[
\dim\bigcap_{j=1}^4\ker(r_j-I)\geq2.
\]

**Beweis.** Sind r1=r2, so folgt r3=r4. Die vier Normalrichtungen spannen dann höchstens zwei Dimensionen. Sind r1 und r2 verschieden, haben ihre Normalvektoren u1,u2 unabhängige Richtungen. Aus r1 r2 x=x folgt r2 x=r1 x. Die Differenz der beiden Rang eins Beiträge verschwindet nur, wenn x auf beiden Normalen senkrecht steht. Somit

\[
\operatorname{Fix}(r_1r_2)
=(\operatorname{span}\{u_1,u_2\})^\perp.
\]

Wegen r1 r2 = r4 r3 besitzt auch das zweite Paar genau diesen bewegten Zweierraum. Alle vier Normalen liegen darin. ∎

Insbesondere kann kein solcher Vierer die irreduzible ursprüngliche Viererdarstellung der vollen Gruppe erzeugen. Das ist eine Einschränkung dieses gewählten Modellfalls, kein allgemeiner Satz über vier beliebige Monodromien.

### Satz 2: Vollständige endliche Klassifikation

Es gibt genau 14.760 geordnete Vierer und genau 17 simultane Konjugationsklassen. Die Bildgruppen besitzen die Ordnungen 2,4,6 oder 8.

**Exakter endlicher Beweis.** Alle ursprünglichen Reflexionen werden als Matrizen \(r=A/2\) mit Gaussian ganzzahligem A gespeichert. Jeder Multiplikationsschritt ist ganzzahlig, und jede notwendige Division durch zwei oder vier wird auf exakte Teilbarkeit geprüft. Die sechs im Zertifikat angegebenen ursprünglichen Reflexionen erzeugen durch vollständigen Abschluss 46.080 verschiedene Matrizen. Alle 60 Ausgangsreflexionen liegen darin. Ihre Konjugationswirkung auf der 60er Menge hat Ordnung 11.520. Ihr Kern ist genau das nachgeprüfte skalare μ4.

Bilde für alle 3.600 Paare die Produktklassen. Die Multiplikationen haben folgende Kollisionsgrößen:

| Größe eines Paarproduktblocks | Anzahl solcher Blöcke |
|---:|---:|
| 1 | 0 |
| 2 | 360 |
| 3 | 640 |
| 4 | 180 |
| 6 | 30 |
| 60 | 1 |

Da r1r2=(r3r4)⁻¹=r4r3, ist jedes zulässige Vierertupel genau ein geordnetes Paar aus demselben Produktblock. Deshalb

\[
60^2+30\cdot6^2+180\cdot4^2+640\cdot3^2+360\cdot2^2=14760.
\]

Der Abschluss unter den sechs tatsächlichen Konjugationsgeneratoren zerlegt diese Menge in die 17 in Anhang A aufgeführten disjunkten Orbits. Jeder Orbit wird vollständig erzeugt, kein weiterer Vertreter bleibt übrig. ∎

Die Isotropiegruppe einer Darstellung ist \(C_\Gamma(\operatorname{im}\rho)\). Die angegebenen Stabilisatorordnungen folgen aus der exakten Orbitformel.

### Geordnete Moduli und Hurwitzorbits sind verschieden

Die zusätzlichen Hurwitzoperationen lauten

\[
(r_i,r_{i+1})\mapsto(r_ir_{i+1}r_i^{-1},r_i).
\]

Sie erhalten das Produkt und erzeugen auf den 17 Klassen genau fünf Orbits:

| Bildtyp | Anzahl geordneter Klassen | Größe jedes Konjugationsorbits | Stabilisator in G31 |
|---|---:|---:|---:|
| C2 | 1 | 60 | 768 |
| Diedergruppe der Ordnung 8 | 6 | 720 | 64 |
| S3 | 4 | 1920 | 24 |
| C2 × C2, Typ A | 3 | 180 | 256 |
| C2 × C2, Typ B | 3 | 720 | 64 |

Die beiden Vierergruppen sind abstrakt isomorph, aber verschieden in die ursprüngliche Gruppe eingebettet. Die fünf Hurwitzorbits sind **keine fünf Zusammenhangskomponenten des festen geordneten Modulraums**: dieser ist eine diskrete Menge aus 17 Punkten. Sie beschreiben einen zusätzlichen Quotienten nach Punkturbewegungen.

Als reduzierte komplexe Nullraumgeometrie ist die endliche grobe Modulimenge glatt und nulldimensional. Ihre nichttrivialen Automorphismen müssen durch den Quotientgruppoid erhalten bleiben. Sie sind nicht mit gewöhnlichen positiven dimensionalen Singularitäten zu verwechseln.

## 3. Die Viertelclock wird als zusätzliche äquivariante Struktur klassifiziert

Die Basisrotation c(z)=iz hat Ordnung vier und permutiert die vier Punkturen. Mit einem festen Basispunkt 0 kann ihr äquivarianter Monodromievertrag durch

\[
C^4=I,\qquad r_{j+1}=Cr_jC^{-1}
\]

geschrieben werden. Dabei ist C ein Element der **ursprünglichen Gruppe** Γ. Ein Element geometrischer Ordnung vier darf auf der Faser auch Ordnung eins oder zwei haben; für eine strikte Faserquarterclock verlangen wir zusätzlich C²≠I.

Nach Eliminieren von r2,r3,r4 lautet die zusätzliche Relation

\[
(r_1C)^4=I.
\]

Der zugehörige Orbifoldgruppoid besitzt damit die fundamentale Gruppenpräsentation

\[
\langle C,D\mid C^4=D^4=1\rangle\cong C_4*C_4,
\qquad D=r_1C.
\]

### Satz 3: Genau 62 strikte Faserclockklassen

Von den 17 Viererklassen sind nur die Nummern 1,13,16 unter zyklischer Verschiebung der Punkturen stabil. Sie besitzen die Muster

\[
(r,r,r,r),\qquad(r,s,r,s),
\]

wobei die alternierenden Muster zu den beiden verschiedenen Vierergruppeneinbettungen gehören.

| Grundvierer | Tatsächliche Clocks C der exakten Ordnung vier über dem festen Vertreter | Klassen unter seinem Zentralisator |
|---|---:|---:|
| Klasse 1 | 256 | 44 |
| Klasse 13 | 96 | 15 |
| Klasse 16 | 24 | 3 |
| **Summe der Klassen** | | **62** |

Unter der gesamten simultanen Konjugation gibt es 49.920 geordnete Paare aus Vierer und strikter Faserclock. Lässt man auch Faserordnungen eins und zwei zu, erhält man insgesamt 80 Klassen.

**Beweis.** Die vollständige Gruppe mit 46.080 Matrizen wird für jeden der drei stabilen Vertreter nach den exakten Gleichungen Crj=rj+1C und C⁴=I durchsucht. Die Clocks werden unter dem ebenfalls vollständig erzeugten Viererzentralisator quotientiert. Die Liste aller 62 Vertreter einschließlich Matrizen, Ordnung, Stabilisator und charakteristischem Polynom steht in `verified_results.json`. ∎

Die **zentrale ursprüngliche Clock** C=iI ist ein spezieller Fall. Sie konjugiert trivial. Wird gerade dieser Operator vorgeschrieben, müssen alle vier rj gleich sein; es bleibt eine einzige Grundviererklasse. Die Zahl 62 darf nicht als Klassifikation nur dieses einen festgehaltenen zentralen Operators ausgegeben werden.

## 4. Der vollständige zentrale Kokykel und seine Klasse

Bezeichne den zentralen Generator hier mit z, damit er nicht mit der Reflexionsgruppe Γ verwechselt wird. Aus [S2] stammt

\[
c(g,h)=-\epsilon(g)\epsilon(h)\pmod4.
\]

Die Erweiterung ist

\[
\widehat\Gamma=(\mathbb Z/4)\times_c\Gamma,
\]

\[
(a,g)(b,h)=(a+b-\epsilon(g)\epsilon(h),gh).
\]

Die drei Paritätsvariablen genügen, um die volle Kokykelgleichung zu verifizieren; sie folgt auch direkt aus \(\epsilon(gh)=p+q-2pq\).

### Satz 4: Die Klasse ist nichttrivial und hat genau Ordnung zwei

Für jede ursprüngliche Reflexion r gilt

\[
(0,r)^2=(-1,e)=z^{-1}.
\]

Jeder andere Lift besitzt die Form zᵏ(0,r), dessen Quadrat z²ᵏ⁻¹ ist. Dieses ist niemals eins. Über der Untergruppe ⟨r⟩≅C2 ist die Erweiterung daher C8 mit Kern C4 und nicht gespalten. Die globale Kokykelklasse ist somit nicht null.

Andererseits gilt für b(g)=−ε(g) modulo vier

\[
\delta b(g,h)=-2\epsilon(g)\epsilon(h)=2c(g,h).
\]

Also

\[
\boxed{[c]\neq0,\qquad 2[c]=0.}
\]

Diese Klasse wird hier vollständig bestimmt; eine Berechnung der gesamten Gruppe H²(G31,Z/4) wird nicht behauptet.

### Genaues Spaltkriterium auf Untergruppen

Es gibt eine konkrete alternative Beschreibung:

\[
\widehat\Gamma\cong
\{(g,t)\in\Gamma\times\mathbb Z/8:t\bmod2=\epsilon(g)\},
\qquad t=2a-\epsilon(g).
\]

Daraus folgt: Die Einschränkung auf eine Untergruppe H spaltet genau dann, wenn der Charakter ε|H zu einem Homomorphismus H→Z/8 angehoben werden kann.

Insbesondere:

* Auf ker ε verschwindet die Klasse.
* Auf jeder Untergruppe mit einer ursprünglichen Reflexion verschwindet sie nicht.
* Auf einer zyklischen Gruppe ⟨g⟩ mit ε(g)=1 spaltet sie genau dann, wenn die Ordnung von g durch acht teilbar ist.

### Die Koeffizientengruppe gehört zur Aussage

Nach Erweiterung der zentralen Phasen von μ4 auf μ8 wird c durch die Einbettung Z/4→Z/8 als 2c geschrieben. Dann ist bereits der Kozyklus selbst ein Korand:

\[
2c=\delta(-\epsilon)\quad\text{in }\mathbb Z/8.
\]

Die Obstruktion ist also nicht durch μ4 zulässige Umphasungen entfernbar, wohl aber nach einer **Änderung** der Koeffizientengruppe. Diese Änderung ist nicht bloß eine neue Sektion desselben ursprünglichen μ4 Objekts.

### Abstrakte Erweiterung und effektives E8 Bild

Die abstrakte Erweiterung über Γ hat Ordnung 184.320. Die in [S2] angegebene Abbildung (a,g)↦zᵃL(g) hat den viergliedrigen diagonalen Kern

\[
\{(-k,(iI_4)^k):k\in\mathbb Z/4\}.
\]

Ihr effektives Bild besitzt Ordnung 46.080 und ist eine Erweiterung der projektiven Gruppe. Es darf nicht ohne weiteres als dieselbe Erweiterung über dem ursprünglichen G31 behandelt werden. Die nachstehenden Liftzählungen beziehen sich ausdrücklich auf \(\widehat\Gamma\to\Gamma\).

## 5. Der Hauptsatz: Freie Lifts existieren, äquivariante Lifts nicht

### 5.1 Auf der punktierten Sphäre besteht keine absolute H² Obstruktion

U zieht sich auf einen Strauß aus drei Kreisen zurück. Daher gilt für konstante Z/4 Koeffizienten

\[
H^2(U;\mathbb Z/4)=0.
\]

Man kann dies ohne Kohomologiesatz direkt benutzen: Wähle beliebige Lifts der ersten drei freien Meridiane und definiere den vierten als Inverses ihres Produkts. Jede Γ Monodromie auf U ist somit liftbar.

Die Nichttrivialität von [c] auf Γ widerspricht dem nicht. Entscheidend ist der Rückzug auf die Fundamentalgruppe; die allgemeine Form dieses Liftkriteriums ist in [L2] behandelt.

### 5.2 Exakt 64 Lifts bei festen Basismonodromien

Für einen ursprünglichen Vierer r1r2r3r4=I liefern die kanonischen Lifts

\[
L(r_1)L(r_2)L(r_3)L(r_4)=z^{-2}=z^2.
\]

Mit \(\widehat r_j=z^{k_j}L(r_j)\) lautet die notwendige und hinreichende Produktbedingung

\[
\boxed{k_1+k_2+k_3+k_4=2\pmod4.}
\]

Es gibt 4³=64 Lösungen. Beispielsweise funktioniert (2,0,0,0). Die natürliche unveränderte Sektion mit (0,0,0,0) funktioniert dagegen nicht.

Weil der Kokykel nur vom Konjugationsinvarianten ε abhängt und symmetrisch in den beiden Paritäten ist, erhält Konjugation in \(\widehat\Gamma\) diese zentralen Labels. Deshalb entstehen über den 17 geordneten Basisklassen insgesamt 17·64=1.088 Konjugationsklassen, wenn beliebige zentrale Randlifts erlaubt sind. Diese Zahl ist nicht die Zahl der Klassen im effektiven E8 Bild nach dessen weiterem Quotienten.

### 5.3 Neue äquivariante Obstruktion

**Satz 5.** Eine C4 äquivariante ursprüngliche Vierreflexionsmonodromie kann nicht zu einer strikt C4 äquivarianten Monodromie in \(\widehat\Gamma\) gehoben werden, solange die ursprüngliche zentrale μ4 Struktur beibehalten wird.

**Erster Beweis.** Für jeden angehobenen Clockoperator ist

\[
\widehat C\,L(r_j)\,\widehat C^{-1}=L(r_{j+1}).
\]

Zentrale Faktoren ändern sich bei Konjugation nicht. Eine strikt zyklische Äquivarianz der vier angehobenen Meridiane erzwingt daher k1=k2=k3=k4. Deren Summe ist null modulo vier, die Produktbedingung verlangt aber zwei. Widerspruch. ∎

**Unabhängiger Beweis über die Orbifoldrelation.** Die beiden Generatoren C und D=r1C haben vierte Potenz eins. Ihre Determinantenparitäten sind entgegengesetzt. Für einen beliebigen Lift eines Elements g mit g⁴=I gilt

\[
(z^aL(g))^4=z^{-2\epsilon(g)}.
\]

Mindestens einer der beiden Generatoren hat ε=1 und somit vierte Potenz z²≠1. Also kann die Darstellung von C4*C4 nicht geliftet werden. ∎

Für den konstanten Koeffizientenvertrag gilt

\[
H^2(C_4*C_4;\mathbb Z/4)
\cong\mathbb Z/4\oplus\mathbb Z/4.
\]

Die rückgezogene Klasse hat je nach Clockparität die Werte (0,2) oder (2,0). Sie ist in beiden Fällen nichttrivial von Ordnung zwei.

**Bedeutung:** Derselbe Grundraum erlaubt 64 Lifts, aber seine zusätzlich verlangte Viertelclockäquivarianz erzeugt eine echte Obstruktion. Absolute, relative und äquivariante Liftfragen sind deshalb getrennte mathematische Probleme.

Eine Clock, die zusätzlich eine andere Wirkung auf dem Zentrum hätte, oder eine veränderte Koeffizientengruppe definiert eine andere Kategorie und wird von diesem Satz nicht ausgeschlossen.

## 6. Vierer, Quartikfünfer und äußerer Sechser als gemeinsame Bündeldaten

Aus jeder ursprünglichen Monodromie entstehen die assoziierten Bündel

\[
\mathcal V=\widetilde U\times_\rho\mathbb C^4,
\quad
\mathcal W=\Lambda^2\mathcal V,
\quad
\mathcal Q\subseteq\operatorname{Sym}^4\mathcal V.
\]

Die Ränge sind vier, sechs und fünf. Auf dem unnormierten ursprünglichen Quartikrahmen lautet die Metrik

\[
D=\operatorname{diag}(4,12,12,12,24).
\]

Alle sechs verwendeten Gruppengeneratoren wurden exakt auf die fünf Quartikbasisvektoren eingeschränkt; sie erfüllen R*DR=D, Spur 3 und Determinante −1. Damit werden keine verschiedenartigen Fünferräume gleichgesetzt.

Die determinantkorrigierte Darstellung ist

\[
R_0(g)=\det(g)^{-1}R(g).
\]

Sie beschreibt den im Original benannten parallelen Quartiktransport und liegt im gegebenen reellen Chart in SO(5).

### Kein universeller Intertwiner vom Quartikfünfer in den äußeren Sechser

Der zentrale ursprüngliche Operator iI wirkt auf Sym⁴(V), also auf Q, als +I. Auf Λ²(V) wirkt er als −I. Für einen Γ Intertwiner T gilt damit

\[
-T=T,
\qquad T=0.
\]

Also

\[
\boxed{
\operatorname{Hom}_\Gamma(Q,\Lambda^2V)=0,
\qquad
\operatorname{Hom}_\Gamma(\det^{-1}Q,\Lambda^2V)=0.
}
\]

Die zweite Aussage benutzt det(iI)=1. Ein zusätzlicher exakter Kommutantentest ergibt

\[
\dim\operatorname{End}_\Gamma(\Lambda^2V)=1,
\]

also Irreduzibilität des äußeren Sechsers über C. Die allgemeine äußere Konstruktion wird auch in [L1] behandelt.

Für eine konkrete kleine Monodromiebildgruppe können zusätzliche Intertwiner entstehen. Der universelle Γ Ausschluss wird nicht als Ausschluss aller Untergruppenabbildungen ausgegeben.

## 7. Kohomologie und eine integrale Unterscheidungsinvariante

Für eine Darstellung V der freien Gruppe F3 lautet der zelluläre Kokettenkomplex

\[
0\longrightarrow V
\xrightarrow{d^0}
V^3\longrightarrow0,
\qquad
v\mapsto((M_1-I)v,(M_2-I)v,(M_3-I)v).
\]

Daher

\[
H^0(U;V)=V^{\operatorname{im}\rho},
\qquad
\dim H^1(U;V)=2\dim V+\dim H^0(U;V),
\qquad
H^{k\ge2}(U;V)=0.
\]

Für V=C⁴ ergibt das H¹ Dimension 11 im zyklischen Typ und Dimension 10 in allen anderen Viererklassen. Für Q und Λ²V stehen sämtliche Werte in Anhang A und im JSON Zertifikat.

Diese H¹ Räume sind Kohomologie mit den angegebenen **linearen Koeffizienten**. Ihre Dimensionen sind nicht die Dimensionen des endlichen G31 Modulraums.

### Das ursprüngliche integrale Gitter

Die 240 ganzzahligen Wurzeln spannen ein Rang acht Gitter mit Bilinearform x·y/2. Eine durch Hermitenormalform gewonnene Basis B hat det B=16. Ihre Grammatrix

\[
G_\Lambda=B^TB/2
\]

ist gerade, integral und hat Determinante eins. Alle 60 Reflexionen wirken in dieser Basis ganzzahlig und erhalten GΛ. Basis und Grammatrix stehen im Ergebniszertifikat.

### Satz 6: Die beiden Vierergruppentypen sind integral topologisch verschieden

Für Typ A gilt

\[
\boxed{H^1(U;\mathcal L_A)\cong\mathbb Z^{20}\oplus\mathbb Z/2.}
\]

Für Typ B gilt

\[
\boxed{H^1(U;\mathcal L_B)\cong\mathbb Z^{20}.}
\]

**Beweis.** Der vollständige ganzzahlige Komplex wird in der genannten Gitterbasis gebildet. Seine nichtverschwindenden Smithinvarianten sind für Typ A (1,1,1,2), für Typ B (1,1,1,1). Der Kokettenzielraum hat Rang 24. Damit folgen freie Ränge 20 und der jeweils angegebene Torsionsanteil. Die drei geordneten Klassen jedes Typs liefern dieselben Smithinvarianten. ∎

Für den zyklischen Typ ergibt sich Z²² ohne Torsion, für alle anderen Typen Z²⁰ ohne Torsion.

Bereits die komplexe Quartikdarstellung unterscheidet die zwei Vierergruppeneinbettungen: dim Q^H beträgt bei Typ A vier, bei Typ B drei. Die integrale Torsion ist eine zusätzliche, unmittelbar topologische Unterscheidung.

Da U auf einen Graphen zurückgezogen werden kann, verschwinden gewöhnliche positive Chernklassen dieser komplexen Bündel. Das macht ihre flachen Monodromien nicht trivial. Nichttriviale Daten liegen hier in Transport, Isotropie, lokalen Koeffizienten und der äquivarianten zentralen Obstruktion.

## 8. Was abgeschlossen ist und was gesondert bleibt

Abgeschlossen sind der definierte vierreflexive Modellfall, seine Clockverfeinerung, der konkrete μ4 Kokykel, die angegebenen Liftobstruktionen und die lokalen Kohomologien auf den ursprünglichen Darstellungen.

Nicht als erledigt ausgegeben werden:

* die gesamte Untergruppen- und zentrale Erweiterungsklassifikation von G31;
* andere periphere Konjugationsklassen oder eine variable Anzahl von Punkturen;
* eine Quotientierung nach einem anderen, zusätzlich festgehaltenen Quadratikstabilisator;
* eine freie Gleichsetzung des Quartikfünfers und eines äußeren Fünfers;
* eine bereits konstruierte allgemeine Fusions- oder Bordismustheorie.

Die nächste präzise mathematische Erweiterung wäre die Vergrößerung des peripheren Klassenvertrags oder die Verwendung anderer Punkturzahlen, jeweils mit demselben vollständigen Gruppen- und Liftwörterbuch. Der jetzige Modellfall bietet dafür einen fertigen Basistest statt einer nur programmatischen Definition.

## 9. Reproduzierbarkeit

Das Paket führt fünf Prüfer aus:

1. vollständige Gaussian ganzzahlige Gruppenenumeration und Monodromieorbits;
2. vollständige Clockklassifikation und Zentralisatorquotienten;
3. exakte Quartik- und Außendarstellungen sowie komplexe Kohomologie;
4. integrale Gitterbasis und Smithnormalformen;
5. modulare Kokykel-, Lift- und äquivariante Obstruktionszertifikate.

Die analytischen Sätze über beliebige Untergruppen, Homotopie oder die universelle Rangschranke beruhen auf den oben ausgeschriebenen Beweisen. Die Programme prüfen die endliche konkrete Instanz und ihre vollständige Enumeration, nicht sämtliche Aussagen der Originaltheorie.

Ausführung:

```bash
python run_checks.py --optimized
```

Es werden keine Repositoryzugänge benötigt und keine externen Dienste angesprochen. Die fünf Skripte laufen in einem frischen temporären Verzeichnis. Sämtliche entscheidenden Gruppenvergleiche verwenden Gaussian Ganzzahlarithmetik, die Bündelrechnungen rationale komplexe Arithmetik, die Gitterrechnung ganzzahlige Normalformen. Pythonassertions werden nicht als Prüfmechanismus verwendet.

## Quellen

**[S1]** `experiments/theory-contracts/compiler-quartic-carry-20260919/core.py`, gleicher gepinnter Commit. Übernommen werden die Wurzel- und Quartikbasisdefinitionen, nicht die dortigen Laufresultate. Git Blob SHA: `67dbda074e1564b9bb98c237f0b11fe921224b12`.

**[S2]** `experiments/theory-contracts/compiler-quartic-holonomy-20260919/PROOF.txt`, gleicher gepinnter Commit. Verwendet werden die tatsächliche Sektion L, ihre Kompositionsrelation, die Determinantencharaktere und der Quartiktransport. Git Blob SHA: `89b58a7fa555aaf42a5e4a21a833fadcb150ed01`.

**[L1]** C. Bonnafé, *Une (nouvelle ?) construction du groupe de réflexion complexe G31*, arXiv:math/0302042. Vergleich für Gruppe und äußere Darstellung, nicht Quelle der hier neu enumerierten 17 und 62 Klassen.

**[L2]** R. Joshi, S. Spallone, *Central Extensions and Cohomology*, arXiv:2307.14658. Allgemeiner Rahmen für die Kohomologieklasse einer zentralen Erweiterung und das Rückzugskriterium für Lifts. Die hier verwendeten Kokykel- und Nichtspaltungsbeweise sind vollständig oben ausgeführt.


## Anhang A: Sämtliche 17 geordneten Klassen

Die Reflexion s_v ist I−2vv*/(v*v). Die sieben benötigten Linienvertreter sind

```text
a  = (1,0,0,0)
b  = (1,1,0,0)       b′ = (1,−1,0,0)
c  = (1,1,1,1)       c′ = (1,−1,−1,−1)
d  = (0,1,0,0)
e  = (0,1,1,0)
```

Sie sind exakt proportional zu den im Originalcode erzeugten Linien mit Indizes 0,1,12,13,44,45,46. Die Skalierung des Linienvektors verändert seine Reflexion nicht.

| Nr. | Monodromievertreter | Orbitgröße | Bildordnung | Stabilisator | dim H⁰(V₄) | dim H⁰(Q₅) | dim H⁰(Λ²V₄) | Torsion H¹(Λ) |
|---:|---|---:|---:|---:|---:|---:|---:|---|
| 1 | (a, a, a, a) | 60 | 2 | 768 | 3 | 4 | 3 | 0 |
| 2 | (a, a, b, b) | 720 | 8 | 64 | 2 | 3 | 1 | 0 |
| 3 | (a, a, c, c) | 1920 | 6 | 24 | 2 | 3 | 1 | 0 |
| 4 | (a, a, d, d) | 180 | 4 | 256 | 2 | 4 | 1 | Z/2 |
| 5 | (a, a, e, e) | 720 | 4 | 64 | 2 | 3 | 1 | 0 |
| 6 | (a, b, a, b′) | 720 | 8 | 64 | 2 | 3 | 1 | 0 |
| 7 | (a, b, b, a) | 720 | 8 | 64 | 2 | 3 | 1 | 0 |
| 8 | (a, b, b′, d) | 720 | 8 | 64 | 2 | 3 | 1 | 0 |
| 9 | (a, b, d, b) | 720 | 8 | 64 | 2 | 3 | 1 | 0 |
| 10 | (a, c, a, c′) | 1920 | 6 | 24 | 2 | 3 | 1 | 0 |
| 11 | (a, c, c, a) | 1920 | 6 | 24 | 2 | 3 | 1 | 0 |
| 12 | (a, c, c′, c) | 1920 | 6 | 24 | 2 | 3 | 1 | 0 |
| 13 | (a, d, a, d) | 180 | 4 | 256 | 2 | 4 | 1 | Z/2 |
| 14 | (a, d, b, b′) | 720 | 8 | 64 | 2 | 3 | 1 | 0 |
| 15 | (a, d, d, a) | 180 | 4 | 256 | 2 | 4 | 1 | Z/2 |
| 16 | (a, e, a, e) | 720 | 4 | 64 | 2 | 3 | 1 | 0 |
| 17 | (a, e, e, a) | 720 | 4 | 64 | 2 | 3 | 1 | 0 |

Die vollständigen Vertreter der Clockklassen stehen mit exakten Matrizen (jeweils als zweimal Real- und Imaginärteil) im Ergebniszertifikat.

## Ergänzende unabhängige Orbitzählung

Die 17 Konjugationsorbits wurden zusätzlich über Burnside bestätigt. Für jede der 11.520 projektiven Konjugationswirkungen wurde aus ihren fixierten Reflexionslabels die Zahl fixer Vierertupel bestimmt. Die Summe beträgt 195.840, also 195.840 / 11.520 = 17. Diese Zählung verwendet nicht die zuvor erzeugten Orbitlisten als Klassifikation. Außerdem wurde die Zugehörigkeit aller 60 ursprünglichen Reflexionen und aller vier zentralen Skalare zur vollständig enumerierten Matrixgruppe gesondert geprüft.
