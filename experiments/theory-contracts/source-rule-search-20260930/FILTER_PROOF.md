# Neutrale Rand-zu-Feld-Suche: Standardnorm und inverse Casimirregel

## Ergebnis

Aus dem ursprünglichen neutralen Feld mit \(\langle R(z)R(w)\rangle_0=48/(z-w)^4\) folgt in der üblichen Laurent-Konvention exakt

\[
\kappa_n=\langle0|R_nR_{-n}|0\rangle=8n(n^2-1),\qquad n\ge2.
\]

Damit gibt es im erweiterten Suchraum eine kurze inverse Casimirregel, die auf dem gesamten aktiven Modenturm die Zielantwort erzeugt:

\[
s_n=\frac{1}{\sqrt8}\,n^{-1}(n-1)^{-1/2}(n+1)^{-1/2}
=\frac{1}{\sqrt{8n^2(n^2-1)}} ,
\qquad \kappa_ns_n^2=\frac1n.
\]

Die Faktorisierung erzwingt diese Regel eindeutig in der Klasse \(c n^a(n-1)^b(n+1)^d\), \(a,b,d\) halbzahlig in \([-3,3]\), \(c>0\), \(c^2=2^j\), \(-8\le j\le8\): \(a=-1\), \(b=d=-1/2\), \(j=-3\). Die reine D-only-Klasse \(c n^a\) hat keinen Treffer.

## Modennorm und Primarität

Die Quelle definiert \(R=4T_{A_2}-8T_{U(1)}\) als neutralen Gewicht-2-Zustand im Grad-2-Heisenberg-Vakuumsektor. Dort gilt \(\operatorname{Tr}(A_R)=0\), also verschwindet \(L_2|R\rangle\); \(L_1\) und höhere positive Virasoro-Moden verschwinden unter den angegebenen Vakuumannahmen. Damit ist \(R\) in diesem Sektor Virasoro-primär. Die Quelle warnt zugleich, dass Neutralität und Primarität allein \(R\) nicht eindeutig auswählen; der 48-Zweig ist durch die orthogonale Projektion festgelegt.

Mit \(R(z)=\sum_mR_mz^{-m-2}\) gilt für \(|z|>|w|\):

\[
\frac{48}{(z-w)^4}=48z^{-4}\sum_{k\ge0}\binom{k+3}{3}(w/z)^k.
\]

Der Koeffizient mit \(m=k+2=n\) liefert

\[
\langle0|R_nR_{-n}|0\rangle
=48\binom{n+1}{3}=8n(n^2-1).
\]

Das prüft insbesondere \(\kappa_2=48=\|R_{-2}|0\rangle\|^2\). Eine Halb-Normierung wird im Hauptbefund nicht verwendet. Bei \(n=1\) verschwindet der Faktor \(n^2-1\); der aktive Turm beginnt bei \(n=2\).

## Suchklassen

Für den erweiterten Ansatz

\[
s_n=c n^a(n-1)^b(n+1)^d
\]

ist die Zielgleichung nach Quadrieren

\[
8c^2 n^{1+2a}(n-1)^{1+2b}(n+1)^{1+2d}=n^{-1}.
\]

Die eindeutige Faktorisierung der rationalen Funktionen erzwingt \(a=-1\), \(b=d=-1/2\), \(8c^2=1\). Mit \(c^2=2^j\) ergibt sich \(j=-3\), also der oben angegebene positive Gram-Wurzelfaktor \(c=2^{-3/2}=1/\sqrt8\). Das ist ein echter Treffer in der verlangten erweiterten Koeffizientenklasse und gilt identisch für alle ganzzahligen \(n\ge2\), nicht nur für eine endliche Anpassung.

Für die reine Potenzregel \(s_n=cn^a\) enthält die quadrierte linke Seite zusätzlich die nichtkürzbaren Faktoren \((n-1)(n+1)\), während die Zielseite nur \(n^{-1}\) enthält. Kein Exponent \(a\) und kein positiver Vorfaktor kann diese Faktoren beseitigen. Somit gibt es keinen D-only-Potenzfilter, auch außerhalb des endlichen Exponentenfensters.

Ein endlicher lokaler Differentialfilter hat im Modenraum einen Polynomialmultiplikator \(P(n)\). Für \(\deg P=m\ge0\) wächst \(8n(n^2-1)P(n)^2\) wie eine von null verschiedene Konstante mal \(n^{3+2m}\), während \(1/n\) abfällt. Kein endlicher nichtnull Differentialoperator erreicht daher die Identität auf allen aktiven Moden.

## D-only \(D^{-2}\) mit führender Gram-Amplitude

Setze \(s_n=8^{-1/2}n^{-2}\). Dann

\[
\kappa_ns_n^2=\frac1n-\frac1{n^3},
\qquad
\frac{\kappa_ns_n^2-1/n}{1/n}=-\frac1{n^2}.
\]

Der signierte relative Fehler bei \(n=2,4,8\) ist exakt \(-1/4,-1/16,-1/64\); sein Betrag ist \(1/4,1/16,1/64\).

## Produktverträglichkeit des ganzen Modenturms

Ein skalarer Modenfilter, der dieselbe ursprüngliche Produktabbildung unverändert lässt, müsste \(s(p+q)=s(p)s(q)\) erfüllen. Für die gefundene Regel reicht der exakte Gegencheck \(p=q=2\), um diese Forderung zu widerlegen:

\[
s_4^2=\frac1{1920},\qquad (s_2s_2)^2=\frac1{9216}.
\]

Beide Multiplikatoren sind positiv, aber ihre Quadrate unterscheiden sich. Also ist dieser Filter keine Multiplikative Charakterabbildung auf dem Quellfeld-Modenturm. Der Geltungsbereich ist begrenzt: Eine Rand-zu-Feld-Portabbildung muss kein Algebra-Homomorphismus sein. Der Test verbietet nur die ungeänderte gemeinsame Produktkarte; er widerlegt keine separat definierte Portabbildung mit eigener Kopplung, Smearing- oder Rekonstruktionsregel.

## Auswahlgrenze der ursprünglichen Mark-Green-Karte

Die betrachtete Mark-Green-Abbildung lebt auf vier Marklabels. Ihre orbit-gemittelte Matrix ist \(-c_3\log(2)C\); die gerade Kompression ergibt \(4P_0\) und den Determinanten \(1-4\epsilon c_3\log(2)\). Diese Daten spezifizieren keine Modenabbildung \(n\ge2\), keinen Casimir \(n(n^2-1)\) und keine Auswahl des neutralen R-Vertex. Die Originaldatei `source_flavor_correlator.py` benennt die Rand-zu-R-Vertexabbildung samt Referenzkontraktion ausdrücklich als fehlend. Ein rückwärts angepasster Filter ist daher ein exakter bedingter Kandidat, aber keine Auswahl aus P1 oder aus der logarithmischen Mark-Green-Karte.

## Reproduzierbarkeit

Das Skript [source_filter_search.py](source_filter_search.py) verwendet nur die Python-Standardbibliothek und exakte `Fraction`-Arithmetik. Es erzeugt die Resultate in [results/filter.json](results/filter.json):

```sh
python3 source_filter_search.py --outputpath work/results/filter.json
```

**Verdict:** Die inverse Casimirregel ist der eindeutige exakte bedingte Monomfilter in der angegebenen erweiterten Klasse für die Standardnorm (K=8). D-only-Potenzen und endliche lokale Differentialfilter sind ausgeschlossen. Die Regel erhält die ursprüngliche gemeinsame Produktstruktur nicht automatisch, und die Mark-Green-Quelle wählt weder die Modenabbildung noch den Vertex. Keine physische Ursprungs- oder TOE-Aussage.
