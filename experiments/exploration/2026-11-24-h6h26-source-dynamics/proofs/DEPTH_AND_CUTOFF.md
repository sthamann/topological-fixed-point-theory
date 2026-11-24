# Armtiefe, lokale Kopplung und die Grenze des effektiven Modells

## Ergebnis

Die fünf Schichten aus H6h27 sind für das angegebene Größenordnungsbeispiel
nachvollziehbar. Sie sind weder eine hergeleitete physische Endtiefe noch
eine allgemeine UV Grenze. Mit den tatsächlichen Kantenmassen und relativen
Planckmassen bestehen beim Beispiel
\(M_c=2{,}435\times10^{18}\,\mathrm{GeV}\), \(\mu=1\,\mathrm{TeV}\)
in beiden Modellen fünf Armlagen den lokalen Test \(m_{\rm edge}<\Lambda_3\).
Die sechste Lage besteht ihn nicht.

Die fünfte Lage besitzt dabei nur noch
\(\Lambda_3/m_{\rm edge}\simeq3{,}02\) im flachen und \(2{,}31\) im
gewichteten Modell. Verlangt man einen Faktor zehn Abstand zwischen
Massenschwelle und lokaler Kopplungsskala, bleiben vier Lagen; bei einem
Faktor hundert bleiben drei.

Der Faktor \(535{,}49\) pro Schicht bezeichnet das Wachstum von \(1/M_d\).
Ein quartischer Kontakt wächst stattdessen mit \(q^{-2}\simeq286751\);
im gewichteten Modell kommt der Quotient der Tiefengewichte hinzu.
Der kleine räumliche Nullmodenschwanz ist kein Nachweis, dass alle
Quantenkorrekturen aus tiefen Schichten klein sind.

Die algebraischen Formeln und Tiefenzertifikate stehen in
[process_depth_control.py](../../../../research/originlab/process_depth_control.py).
Die zehn gezielten Tests bestanden in 0,292 Sekunden. Entscheidungen über
die native Tiefe benutzen rationale Schranken für \(e^{-2\pi}\), keine
Gleitkomma Proben.

## 1. Welche Tiefe gezählt wird

In Z51 trägt eine Armmetrik den Index \(n=0,1,\ldots\) und die Planckmasse
\(M_{j,n}=M_cq^{n+1}/2\). Hier bezeichnet

\[
d=n+1\ge1,\qquad M_d=\frac{M_c}{2}q^d
\]

die Anzahl der Schritte vom Zentrum zu dieser Armmetrik. Fünf vollständige
Armlagen bedeuten \(d=1,\ldots,5\), also \(n=0,\ldots,4\). Ein so
abgeschnittener Stern hat 21 Metriken und 20 Kanten. Die fünf Lagen sind
nicht die Gesamtzahl seiner massiven Moden.

Die normierte Kopplung einer Kante hängt von beiden anliegenden Metriken
ab. Setze

\[
\rho_d=\frac{M_d}{M_{\rm parent}}=
\begin{cases}
q/2,&d=1,\\
q,&d\ge2.
\end{cases}
\]

Die erste Kante ist wegen der vier Arme besonders zu behandeln.

## 2. Exakte lokale Kantenformeln

Für den flachen Stern ist \(w_d=1\), für den neuen gewichteten Stern
\(w_d=d\). Die gewählte Steifigkeit lautet

\[
k_d=\mu^2w_dM_d^2.
\]

Mit der Konvention \(M_v^2R[g_v]/2\) und
\(g_v=\eta+2h_v/M_v\) ergibt das isolierte Paar:

\[
M_{{\rm rel},d}
=\left(M_{\rm parent}^{-2}+M_d^{-2}\right)^{-1/2}
=\frac{M_d}{\sqrt{1+\rho_d^2}},
\]

\[
m_{{\rm edge},d}^2
=\frac{k_d}{M_{{\rm rel},d}^2}
=\mu^2w_d(1+\rho_d^2).
\]

Dies ist der einzige nichtverschwindende Fierz Pauli Eigenwert der
isolierten Zweimetrikkante. Im vollständigen Graphen sind die Eigenmoden
Linearkombinationen der Vertexfelder; eine Kante wird damit nicht zu einem
eigenständigen globalen Teilchen erklärt.

Die übliche lokale Potenzzählung verwendet die Skala

\[
\Lambda_{3,d}
=\left(m_{{\rm edge},d}^2M_{{\rm rel},d}\right)^{1/3}
=\left[\mu^2w_dM_d\sqrt{1+\rho_d^2}\right]^{1/3}.
\]

Der Zusammenhang von massiven Gravitonen und einer solchen
\(\Lambda_3\) Skala ist in der Primärliteratur ausgearbeitet. Die
multimetrische Kopplungsskala kann darüber hinaus von der Graphstruktur
und der Mischung longitudinaler Moden abhängen. Insbesondere sind
perturbative Kopplungsskala und physische Abschneideskala nicht automatisch
identisch:
[Scargill und Noller, Abschnitte 1 und 2](https://arxiv.org/abs/1511.02877);
[Bonifacio und Hinterbichler](https://arxiv.org/abs/1806.10607).
Die hier folgende Prüfung ist daher eine benannte lokale Abschätzung,
keine vollständige Streurechnung des Sterns.

## 3. Das tatsächlich geprüfte Kriterium

Für eine vorgegebene Trennung \(\gamma>0\) ist

\[
\frac{\Lambda_{3,d}}{m_{{\rm edge},d}}>\gamma
\quad\Longleftrightarrow\quad
M_d>\gamma^3\mu\sqrt{w_d}(1+\rho_d^2).
\]

Die schwächste hier verwendete Schwelle ist \(\gamma=1\).
Der in Z51 verwendete Knotentest \(M_d\ge\mu\) ist beim flachen Modell
für sehr kleines \(q\) eine gute Näherung. Für das gewichtete Modell
fehlt ihm der Faktor \(\sqrt d\); für allgemeines \(q\) fehlt beiden
Lesarten der Faktor \(1+\rho_d^2\). Gleichheit an der Schwelle ergibt
noch kein offenes Intervall perturbativer Energien oberhalb der Masse;
der Helfer verwendet deshalb ein striktes Größerzeichen.

Für innere Kanten \(d\ge2\), mit
\(a=\log(1/q)\) und \(R=M_c/(2\mu)\), ist die kontinuierliche Grenztiefe
bei \(\gamma=1\)

\[
d_*^{\rm flach}
=\frac{\log[R/(1+q^2)]}{a},
\]

\[
d_*^{\rm gewichtet}
=\frac{W\!\left(2a[R/(1+q^2)]^2\right)}{2a}.
\]

Hier ist \(W\) die positive Lambert Funktion. Zulässig sind die ganzen
Tiefen \(d<d_*\); die erste Kante wird mit ihrer eigenen Normierung
getrennt geprüft. Eine höhere Trennung ersetzt \(R\) durch \(R/\gamma^3\).

Der lokale Quotient nimmt mit jeder Lage strikt ab. Im Inneren gilt

\[
\frac{(M_{\rm rel}/m_{\rm edge})_{d+1}^2}
     {(M_{\rm rel}/m_{\rm edge})_d^2}
=q^2\frac{w_d}{w_{d+1}}<1.
\]

Beim ersten Schritt kommt der zusätzliche Faktor
\([(1+q^2/4)/(1+q^2)]^2<1\) hinzu. Deshalb reicht eine zertifizierte
bestandene Lage und die erste zertifiziert nicht bestandene Lage,
um die gesamte lokale Grenze festzuhalten.

## 4. Zahlen und ihre Voraussetzungen

Der reduzierte Planckwert folgt aus
\(\overline M_{\rm P}=(8\pi G_N)^{-1/2}\). Die aktuelle PDG Tabelle gibt
\(G_N=6{,}70883(15)\times10^{-39}\,\mathrm{GeV}^{-2}\) in natürlichen
Einheiten an:
[PDG, Physical Constants 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-phys-constants.pdf).
Für die folgenden Zahlen wird ausdrücklich
\(M_c=2{,}435\times10^{18}\,\mathrm{GeV}\) angenommen.
Im unendlichen Stern ist die gemessene effektive Planckmasse
\(M_{\rm eff}=M_c/\sqrt{1-q^2}\); eine exakte Anpassung an Newtons Konstante
setzt daher \(M_c=\overline M_{\rm P}\sqrt{1-q^2}\).
Beim nativen \(q\) ändert diese Korrektur keine der angegebenen Tiefen.

Für \(q=e^{-2\pi}\) und \(\mu=1\,\mathrm{TeV}\):

| Armlage \(d\) | \(M_d\) in GeV | Flache Kantenmasse in GeV | Flache \(\Lambda_3\) in GeV | Gewichtete Kantenmasse in GeV | Gewichtete \(\Lambda_3\) in GeV |
| --- | ---: | ---: | ---: | ---: | ---: |
| 3 | \(7{,}929\times10^9\) | \(1000{,}002\) | \(199406\) | \(1732{,}054\) | \(287593\) |
| 4 | \(1{,}481\times10^7\) | \(1000{,}002\) | \(24556\) | \(2000{,}003\) | \(38980\) |
| 5 | \(27650{,}66\) | \(1000{,}002\) | \(3023{,}91\) | \(2236{,}072\) | \(5170{,}81\) |
| 6 | \(51{,}636\) | \(1000{,}002\) | \(372{,}38\) | \(2449{,}494\) | \(676{,}66\) |

Die kontinuierlichen Schwellen liegen bei \(5{,}5283\) beziehungsweise
\(5{,}3942\). Die ganzzahlige Zahl bestandener Lagen ist also jeweils fünf.

| Geforderter Abstand \(\Lambda_3/m_{\rm edge}\) | Flacher Stern | Gewichteter Stern |
| --- | ---: | ---: |
| \(>1\) | 5 Lagen | 5 Lagen |
| \(>10\) | 4 Lagen | 4 Lagen |
| \(>100\) | 3 Lagen | 3 Lagen |

Z51 verwendet \(\mu/M_c=10^{-16}\) als Größenordnungsbeispiel.
Das ist nicht genau dieselbe Planckkonvention; bei dem oben verwendeten
reduzierten Wert beträgt das Verhältnis etwa \(4{,}11\times10^{-16}\).
Beide Beispiele ergeben fünf Lagen.

Die Unterschiede der genaueren Normierung werden bei \(q=1/2\) sichtbar:

| \(\mu/M_c\) | Flach, genaue lokale Prüfung | Gewichtet, genaue lokale Prüfung |
| --- | ---: | ---: |
| \(10^{-3}\) | 8 Lagen | 7 Lagen |
| \(10^{-10}\) | 31 Lagen | 29 Lagen |
| \(10^{-16}\) | 51 Lagen | 49 Lagen |

Der angenäherte Z51 Knotentest zählt in den beiden letzten flachen Fällen
je eine Lage mehr. Bei \(q=e^{-2\pi}\) sind die Unterschiede zu klein,
um die Beispiele \(0,3,5\) bestandener Lagen zu ändern.

## 5. Welche Kopplung um welchen Faktor wächst

Im Inneren eines Arms gilt exakt

\[
\frac{M_{d+1}^{-1}}{M_d^{-1}}=q^{-1}
=e^{2\pi}\simeq535{,}4917.
\]

Das ist beispielsweise der Faktor eines kanonischen kubischen
Einstein Vertexkoeffizienten. Ein quartischer Koeffizient mit
\(M_d^{-2}\) wächst um

\[
q^{-2}=e^{4\pi}\simeq286751{,}313.
\]

Für den lokalen Lorentzschen Quartikrest, mit kanonischem Relativfeld
\(H_{{\rm rel},d}\), gilt

\[
\Delta\mathcal L_{4,d}
=\frac{\tau_d k_d}{2M_{{\rm rel},d}^4}
 \det H_{{\rm rel},d}
=\frac{\tau_d\mu^2w_d(1+\rho_d^2)^2}{2M_d^2}
 \det H_{{\rm rel},d}.
\]

In der Transferzeile \((Bh)_d\) lautet dieselbe Wechselwirkung

\[
\Delta\mathcal L_{4,d}
=\frac{\tau_d\mu^2w_d}{2M_d^2}\det(Bh)_d.
\]

Für konstantes \(\tau_d\) ist der Tiefenfaktor deshalb
\(w_dq^{-2d}\). Mit \(n=d-1\) ist er im gewichteten Modell

\[
(n+1)q^{-2(n+1)}
=q^{-2}(n+1)q^{-2n}.
\]

Bis auf einen tiefenunabhängigen Vorfaktor ist die gesuchte Abhängigkeit
also genau \((n+1)q^{-2n}\). Ein Ansatz \(q^{-4n}\) übersieht den
Faktor \(M_d^2\) in der Kantensteifigkeit.

Die lokale \(\Lambda_3\) selbst sinkt pro innerer Lage dagegen um
\(q^{1/3}\) im flachen beziehungsweise
\([q(d+1)/d]^{1/3}\) im gewichteten Modell.
Ein pauschaler Faktor 535 für alle Kopplungen und Skalen wäre falsch.

## 6. Der kleine Nullmodenschwanz und seine Grenze

Die normierte Nullmode des unendlichen Sterns besitzt über alle vier
Knoten einer Lage zusammen das Gewicht

\[
\sum_{j=0}^3|\eta_{j,d}|^2=(1-q^2)q^{2d}.
\]

Jenseits von \(D\) vollständigen Armlagen ergibt sich exakt

\[
\sum_{d=D+1}^{\infty}\sum_j|\eta_{j,d}|^2
=q^{2(D+1)}.
\]

Für \(D=5\) ist das
\(q^{12}=1{,}7987363\times10^{-33}\).
Z51 berechnet diese Zahl richtig. Seine verbale Formel
\(q^{2(n_{\max}+1)}\) ist beim verwendeten Index \(n_{\max}=4\)
jedoch um eine Lage verschoben: korrekt ist
\(q^{2(n_{\max}+2)}\).

Der Ausdruck misst die räumliche Norm des Nullvektors. Er ist keine
allgemeine Abschätzung von Schleifen, Streuamplituden oder unbekannter
starker Dynamik. Insbesondere gilt an jedem Vertex

\[
\frac{\eta_v}{M_v}=\frac1{M_{\rm eff}}.
\]

Die kleine kanonische Wellenfunktion wird in der physikalischen
Quellenkopplung durch denselben kleinen lokalen Planckfaktor kompensiert.
Ebenso bleibt für orthonormale massive Moden

\[
\sum_v\frac{\eta_v}{M_v}\overline{u_a(v)}u_b(v)
=\frac{\delta_{ab}}{M_{\rm eff}}.
\]

Eine kleine räumliche Nullmodennorm allein macht massive Schleifenmoden
daher nicht automatisch unsichtbar. Ein spezieller Randprozess kann
zusätzliche Unterdrückungen besitzen; diese benötigen seinen tatsächlichen
Propagator, seine Quellen und seine Vertexsummen.

## 7. Flaches Band und gewichtete Leiter sind verschiedene Modelle

Der flache unendliche Stern hat ein beschränktes massives Spektralband.
Darin liegen unendlich viele massive Moden. Der gewichtete Stern hat
dagegen im symmetrischen Sektor die diskrete Leiter
\(\mu^2(1-q^2)n\) und in den drei weiteren Sektoren je eine strikt
interlaciernde Leiter. Unter jeder festen endlichen Grenze der Masse zum
Quadrat liegen nur endlich viele Eigenwerte.

Die alte Aussage über unendlich viele Moden in einem beschränkten Band
darf deshalb nicht auf den gewichteten Stern übertragen werden.
Die fallenden lokalen Planckmassen bleiben in beiden Konstruktionen
vorhanden, sodass der hier benannte lokale Test bei endlicher Tiefe endet.
Aus diesem Test folgt noch keine physische Terminierung des Graphen.
Auch ein endlicher Graph, eine neue Beschreibung ab einer bestimmten
Tiefe oder eine kontrollierte starke Fortsetzung wären zusätzliche
dynamische Entscheidungen.

## 8. Was das neue Zertifikat tatsächlich beweist

Die Identitäten für \(M_{\rm rel}\), Kantenmasse und Quartikgewicht sind
symbolisch geprüft. Für das native \(q\) benutzt der Helfer die streng
bewiesenen rationalen Schranken

\[
\frac{1867442}{10^9}<e^{-2\pi}<\frac{1867443}{10^9}.
\]

Der Beweis verwendet die Machin Identität
\(\pi=16\arctan(1/5)-4\arctan(1/239)\), alternierende Restschranken der
Arctan Reihen und eine positive Exponentialreihe mit geometrisch
beschränktem Rest. Anschließend wird die äquivalente rationale
Ungleichung

\[
\frac{q^{2d}}
 {4(\mu/M_c)^2w_d(1+\rho_d^2)^2}>\gamma^6
\]

durch obere und untere Schranken entschieden. Eine unentschiedene
Intervallschranke erhält ausdrücklich keinen Wahrheitswert aus einer
numerischen Probe.

Damit ist die lokale Tiefenzählung reproduzierbar und ihre Bedeutung
klar begrenzt. Der tatsächliche Abschneider der gesamten Theorie, ein
kollektiver Test der longitudinalen Moden und eine dynamisch gewählte
Randbedingung bleiben offen.

## 9. Ein positiver endlicher Elternoperator mit dem richtigen Rand

Für jede gewählte Tiefe \(D\ge1\) gibt es auch eine exakte endliche
Präparationskonstruktion. Definiere den abgeschnittenen Oszillator

\[
a_D=\sum_{n=1}^D\sqrt n\,|n-1\rangle\langle n|
\]

auf \(D+1\) Zuständen sowie einen zweiten unabhängigen Modus \(b_D\).
Setze

\[
\mathcal H_D=\frac{\lambda}{2}
\left[
(a_D-qb_D^\dagger)^\dagger(a_D-qb_D^\dagger)
 +(b_D-qa_D^\dagger)^\dagger(b_D-qa_D^\dagger)
\right],\qquad\lambda>0.
\]

Dieser endliche Operator ist positiv. Auf dem Paarsektor
\(\{|n,n\rangle:0\le n\le D\}\) gilt im Inneren

\[
(\mathcal H_D)_{nn}
=\frac{\lambda}{2}[2n+2q^2(n+1)]
=\lambda[n+q^2(n+1)],
\]

\[
(\mathcal H_D)_{n,n+1}=-\lambda q(n+1).
\]

Am Endpunkt verschwinden dagegen beide Anhebenormen:
\(a_Da_D^\dagger|D\rangle=b_Db_D^\dagger|D\rangle=0\).
Deshalb wird der letzte Diagonaleintrag exakt

\[
(\mathcal H_D)_{DD}=\frac{\lambda}{2}(D+D)=\lambda D.
\]

Die Paarsektorrestriktion ist somit genau
\(\lambda A_D^\dagger\operatorname{diag}(1,\ldots,D)A_D\),
also der physische rechteckige Graphenabschluss und nicht die
Hauptkompression des unendlichen Jacobi Operators.

Beide Summanden vernichten den normierten Vektor

\[
|\eta_{q,D}\rangle
=\sqrt{\frac{1-q^2}{1-q^{2(D+1)}}}
  \sum_{n=0}^D q^n|n,n\rangle.
\]

Sein Grundzustand ist im gesamten endlichen Zweimodenraum eindeutig.
Zum Beweis schreibe einen gemeinsamen Nullvektor als vektorisierte Matrix
\(X=q^NY\). Die beiden Nullgleichungen sind
\(a_DX=qXa_D\) und \(Xa_D^\dagger=qa_D^\dagger X\);
nach Einsetzen werden sie
\([a_D,Y]=[a_D^\dagger,Y]=0\).
Da \(a_D,a_D^\dagger\) die volle Matrixalgebra erzeugen, ist \(Y\)
ein Vielfaches der Identität.

Der reduzierte Zustand ist exakt der endliche Gibbs Zustand

\[
\rho_D=\frac{1-q^2}{1-q^{2(D+1)}}
       \sum_{n=0}^D q^{2n}|n\rangle\langle n|.
\]

Der zusätzliche Preis dieser Konstruktion ist ebenso exakt:

\[
[a_D,a_D^\dagger]=I-(D+1)|D\rangle\langle D|.
\]

Das endliche Modell ist daher keine reguläre Darstellung der
unendlichen CCR. Der frühere Eindeutigkeitssatz für einen gaußschen
Elternoperator auf zwei regulären CCR Moden lässt sich nicht unverändert
auf diese Randmodifikation übertragen.

Der Helfer enthält einen unabhängigen \(D=2\) Zeugen aus dem gesamten
neundimensionalen Zweimodenraum, seine eindimensionale Nullmode, den
thermischen reduzierten Zustand und die exakte Kommutatorabweichung. Die allgemeinen
Matrixkoeffizienten folgen aus den oben ausgewerteten Anhebenormen.
Diese Konstruktion beseitigt die Unsicherheit über den endlichen
terminalen Eintrag. Sie liefert für jede vorgegebene Tiefe einen positiven
Elternoperator und wählt die Tiefe selbst nicht.
