# Wellenpaar unter dem Z174 Standardform und Z4 Klebevertrag

## Ergebnis

Der neueste Quellvertrag trennt zwei verschiedene Paarungen, die bisher leicht
gleich aussehen konnten.

Der Standardvektor paart eine Naht mit ihrem Tomita Partner tatsächlich aus dem
Zustand. Für eine Mode und die in Z174 verwendete KMS Normierung
`beta = 4 pi`, Energie `1/4`, lautet er in der Reihenfolge
`00, 10, 01, 11`

\[
 \Omega={|00\rangle+e^{-\pi/2}|11\rangle\over\sqrt{1+e^{-\pi}}}.
\]

Die Gewichte von `10` und `01` sind exakt null. Außerdem besitzt der Zustand
eine echte Kohärenz zwischen `00` und `11`. Für die nativen kommutierenden
Standardform Operationen `X_A=X tensor 1` und `X_B=1 tensor X^t` ist

\[
 {\langle\Omega,X_AX_B\Omega\rangle\over
  \sqrt{\langle X_A^2\rangle\langle X_B^2\rangle}}
 =\operatorname{sech}(\pi/2)=0.398536815338\ldots .
\]

Dieser Koeffizient ist zunächst nur eine **graduierte Kontrolle**. Ein
Einmoden Pauli `X` ändert Ladung und Fermionparität. Es ist als ungerade
Feldeinfügung in einem insgesamt geraden Korrelator zulässig, aber nicht als
bosonische Präparation einer Gravitationswelle. Der frühere Schluss von diesem
Koeffizienten auf einen physisch admissiblen Wellenquellenkern war daher zu
stark.

Sie wählt jedoch **kein Paar ganzer Nähte**. Für zwei unabhängige Kopien von
`Omega` haben die vier Besetzungswörter vor jeder Klebeprojektion die Gewichte

| Wort | Gewicht |
|---|---:|
| `00` | 0.918868269539 |
| `10` | 0.0397078982950 |
| `01` | 0.0397078982950 |
| `11` | 0.00171593387135 |

Alle vier sind positiv. Der in Z174 modulverträgliche Sektor `X = -1` mit
normaler Statistik `G=X^2=+1` entsteht erst durch die diagonale
Z4 Sektorprojektion. Im reduzierten Besetzungszweig bleibt dann `11`. Das ist
Konditionierung, keine Auswahl durch den unprojizierten Standardzustand.

Die Zahl 0.0017159 und ihr Renormierungsfaktor 582.773 gehören nur zu diesem
endlichen Besetzungszweig. Z174 summiert vollständige Modulbahnen mit Spin und
Oszillatoren. Die Zahl ist daher kein physisches Z174 Sektorgewicht. Die
strukturelle Aussage ist lediglich, an welcher Rechenoperation die Auswahl
entsteht. Insbesondere ist diese Vierworttafel **keine Nachrechnung der
ursprünglichen Z4 Summe**. Deren Aussage folgt aus den vollständigen
Modulbahnen, Charakteren und Positivitätsprüfungen B2 bis B6 von Z174. Die
Tafel ist nur eine Besetzungskontrolle im kleinsten Standardformzweig.

## Ungerade Einmodenkontrolle und Reihenfolge

Für die beiden schwachen nativen Quellen

\[
 U_A(\epsilon)=e^{-i\epsilon X_A},\qquad
 U_B(\epsilon)=e^{-i\epsilon X_B}
\]

wurde das Vakuumpersistenzfunktional für `00`, `10`, `01`, `11` direkt
ausgewertet. Bei `epsilon = 0.173` ergibt sich

| Wort | Amplitude |
|---|---:|
| `00` | 1 |
| `10` | 0.985072785496 |
| `01` | 0.985072785496 |
| `11` | 0.958559106327 |

Die verbundene gemischte Differenz ist

\[
 Z_{11}-Z_{10}Z_{01}
 =-\sin^2(\epsilon)\operatorname{sech}(\pi/2)
 =-0.0118092863971\ldots .
\]

Die vertauschte Reihenfolge liefert exakt denselben Wert. Das folgt nicht aus
einer Näherung, sondern aus der Standardform: linke und rechte Algebra
kommutieren. Auch der Tausch der beiden vollständigen Nahtkopien ändert weder
Zustand noch Sektorprojektor.

Verwendet man stattdessen die rohen ungeraden CAR Erzeuger `a^+` und `b^+`,
gilt

\[
 a^+b^+\Omega=-b^+a^+\Omega .
\]

Die Norm beider Wörter ist gleich. Das Vorzeichen ist die CAR Ordnung und keine
Paarselektion. Die physisch geraden Paaroperationen und die kommutierenden
Standardform Operationen dürfen damit nicht vermischt werden.

Der Wellenpaarversuch aus `wellenpaar.md` braucht neben der symmetrischen
gemischten Antwort auch eine geordnete Tensorantwort. Eine solche
Reihenfolgeabhängigkeit erzeugt der Z174 Standardzustand nicht. Das passt zum
Negativbefund von Z174: Die Standardform trägt den Zustand und seine Uhr, aber
`P_ord`, die Kette aus mehr als zwei Nähten, bleibt offen.

## Physisch gerade Zweimodenpräparation

Der tragende Admissibilitätstest verwendet zwei Fermionmoden auf jeder Seite
und den nativen neutralen Bilinearoperator

\[
 H=c_1^\dagger c_2+c_2^\dagger c_1.
\]

Er ist selbstadjungiert, strikt fermiongerade und erhält die Seitenladung:

\[
 [H,(-1)^Q]=0,\qquad [H,Q]=0.
\]

Damit ist `exp(-i epsilon H)` eine zulässige neutrale Seitenpräparation im
endlichen CAR Modell. Sie ist noch kein Spin 2 Operator und deshalb nicht
automatisch eine Gravitationswelle.

Auf der Z174 Zweimoden Standardform seien

\[
 H_L=H\otimes1,\qquad H_R=JH_LJ=1\otimes\bar H.
\]

Der Zustand hat im einfach besetzten Sektor das Gewicht

\[
 p_{Q=1}={2e^{-\pi}\over(1+e^{-\pi})^2}
         =0.0794157965900\ldots .
\]

Direkt berechnet gilt

\[
 \langle H_L^2\rangle=\langle H_R^2\rangle
 =\langle H_LH_R\rangle=p_{Q=1}.
\]

Der normierte Original Partner Koeffizient ist daher exakt eins. Das ist eine
echte, gerade und ladungserhaltende Standardformpaarung. Sie ist zugleich eine
allgemeine Tomita Identität für einen selbstadjungierten Operator, der mit
`rho=f(Q)` kommutiert. Der Zustand wählt also den J Partner von `H`, aber er
wählt `H` nicht gegenüber allen anderen neutralen Bilinearen aus.

### Vier normierte Schaltungen

Alle fünf Zustandsnormen `00`, `10`, `01`, `11` und `11 vertauscht` sind in
beiden Versuchen exakt eins. Für `epsilon=0.173` ergibt die gleiche
Parameterorientierung auf beiden Seiten

| Schaltung | Amplitude |
|---|---:|
| `00` | 1 |
| `10` | 0.998814543369 |
| `01` | 0.998814543369 |
| `11` | 0.995293564608 |
| `11 vertauscht` | 0.995293564608 |

Die verbundene Differenz ist `-0.00233692743808`.

Für die wirkliche J konjugierte **Schaltung** muss das Vorzeichen des
Parameters umkehren, weil J antiunitär ist:

\[
 J e^{-i\epsilon H_L}J=e^{+i\epsilon H_R}.
\]

Dann gilt

| Schaltung | Amplitude |
|---|---:|
| `00` | 1 |
| `10` | 0.998814543369 |
| `01` | 0.998814543369 |
| `11` | 1 |
| `11 vertauscht` | 1 |

Die verbundene Differenz ist `+0.00236950795406`. Die gemeinsame Schaltung
führt den Standardvektor exakt auf sich selbst zurück. Beide Reihenfolgen sind
identisch, weil die beiden Standardformseiten kommutieren. Das ist eine
Purifikationsidentität und kein Wechselwirkungs oder Zeitordnungssignal.

Als gerade, aber nicht ladungserhaltende Kontrolle wurde zusätzlich

\[
 D=c_1^\dagger c_2^\dagger+c_2c_1
\]

geprüft. Ihr normierter J Paarwert ist
`sech(pi)=0.0862667383341`. Gerade Parität allein reicht somit nicht für den
perfekten neutralen Treffer. Die Ladungserhaltung und das Kommutieren mit dem
Zustandsgewicht sind tragend.

## Was der Zustand auswählt

Die Rechnung trennt drei Stufen.

1. **Aus dem Zustand:** Innerhalb einer Naht trägt `Omega` nur gleiche
   Original und Partnerbesetzung. Zulässige gerade Seitenoperatoren werden mit
   ihrem J Partner gepaart. Für das neutrale `H` ist der normierte Paarwert
   exakt eins.
2. **Nicht aus den Diagonalgewichten allein:** Der unkorrelierte gemischte
   Zustand mit denselben KMS Randgewichten hat
   `Tr(rho X_A X_B)=0`. Auch die Konditionierung dieses Zustands auf den geraden
   Raum `span{00,11}` erzeugt keine Kohärenz. Die Antwort stammt aus dem
   Standardvektor, nicht bloß aus seiner Besetzungsstatistik. Diese Kontrolle
   benutzt ungerade Felder und belegt nur die Herkunft der Kohärenz.
3. **Aus Projektor und Modultest:** Welche zwei vollständigen Nähte gemeinsam
   die Z4 Statistik tragen, folgt nicht aus `Omega tensor Omega`. Z174 findet
   innerhalb der geprüften Kopienzahl nur bei zwei Paaren einen positiven,
   modulinvarianten Sektor mit normaler Statistik. Dafür werden die diagonale
   Klebesumme und der Sektor `X=-1` benutzt. Die neue Prämisse `P_paar` ist
   dadurch präzise lokalisiert und nicht beseitigt.

## Weitere native Zustandshypothese

Z174 prüft außer `K=pi Q` auch das native Ladungsgewicht

\[
 K=2\pi Q^2-\pi Q.
\]

Für zwei Moden bleibt die Original Partner Paarmatrix diagonal. Die
Paarselektion gleicher Moden ist also gegen diese zweite in Z174 bereits
verwendete Zustandshypothese stabil. Ihre Normierung ist es nicht:

| Gewicht | ungerade Diagonalkontrolle | gerades H roh | gerades H normiert | Besetzungskovarianz |
|---|---:|---:|---:|---:|
| `K=pi Q` | 0.398536815338 | 0.0794157965900 | 1 | 0 |
| `K=2pi Q^2-pi Q` | 0.382715434742 | 0.0795523026046 | 1 | -0.00158213669495 |

Das zweite Gewicht korreliert die beiden Besetzungen und verschiebt die rohe
Antwort. Die normierte J Paarung des neutralen Bilinears bleibt eins. Der
Standardformvertrag allein wählt somit noch nicht zwischen den absoluten
Normierungen. Ein Treffer der projektiven Wellengeometrie wäre von dieser
offenen Wahl nicht automatisch absolut normiert.

## Uhrenzeit

Der native Moduloperator ist

\[
 \Delta=\rho\mathbin{\otimes}\rho^{-1},\qquad
 \Delta^{-i/2}=i^{Q_A-Q_B}.
\]

Er lässt `Omega` fest. Für die ungerade Einmodenkontrolle rotiert die Antwort
zwischen `X` und `Y`:

| modulare Zeit | X Antwort | Y Antwort | Quadraturbetrag |
|---:|---:|---:|---:|
| 0 | 0.398536815338 | 0 | 0.398536815338 |
| -1/4 | 0.281808084678 | 0.281808084678 | 0.398536815338 |
| -1/2 | 0 | 0.398536815338 | 0.398536815338 |
| -1 | -0.398536815338 | 0 | 0.398536815338 |

Bei der in Z174 ausgezeichneten halben Modularzeit verschwindet also nur die
gleichphasige X Auslesung. Die Antwort ist vollständig in die konjugierte
Quadratur gewandert. Eine Messung eines einzigen reellen Kanals würde an dieser
Stelle fälschlich „kein Paar“ melden. Uhrkovariant ist der Betrag beider
Quadraturen. Diese Aussage betrifft den ungeraden Feldtest.

Die physisch bessere neutrale Präparation `H` kommutiert mit `rho=f(Q)`. Ihre
Antwort ist für alle modularen Zeiten konstant `0.0794157965900`. Sie liefert
also **kein Uhrenzeitsignal**. Die gerade, ladungsändernde Kontrolle `D` rotiert
dagegen zwischen ihren zwei Paarquadraturen. Bei `t=-1/4` ist ihre X Antwort
null und ihre Y Antwort `0.0794157965900`; der Betrag bleibt gleich. Im
kleinsten Modell sind Ladungserhaltung und eine nichttriviale Modularuhr damit
nicht zugleich durch denselben Seitenoperator realisiert.

Diese Zeit ist modulare Zeit des gesetzten Standardzustands. Die Rechnung
identifiziert sie nicht zusätzlich mit OS Zeit, Lorentzzeit oder der zeitlichen
Reihenfolge zweier Gravitationspulse.

## Folge für den Wellenpaarweg

Der aktuelle Vertrag verbessert den Wellenpaaransatz an einer engeren Stelle:
Es gibt eine physisch gerade, ladungserhaltende und exakt normierte
Original Partner Paarung. Sie ist eine Tomita Paarung im Standardzustand. Sie
ist weder ein Wechselwirkungsnachweis noch eine Auswahl eines bestimmten
Wellenoperators. Die nichttriviale Uhrrotation der früheren Rechnung gehörte
zum ungeraden Feldtest und darf diesem neutralen Quellenoperator nicht
zugeschrieben werden.

Er liefert aber nur die skalare Original Partner Kovarianz. Er wählt weder die
zwei primitiven Weyl Hessians `Q0` und `Q1` noch die projektive Tensorantwort
`[X^2 : tr(TX)X]`. Auch die Zahl von zwei ganzen Nähten kommt erst aus
Modularität, Positivität und Sektorprojektion. Eine geordnete verbundene Antwort
fehlt. Damit bleibt der zentrale Quellenvertrag des Wellenpaars offen:

\[
 \text{Standardform Paarung}
 \quad\not\Rightarrow\quad
 \text{ausgewählte Weyl Paarantwort}.
\]

Der nächste entscheidende Test ist deshalb nicht ein zusätzlicher frei
gewählter Hamiltonoperator. Er muss die beiden Weyl Eingriffe als benannte
Operatoren derselben Nahtalgebra abbilden und zeigen, dass ihr uhrkovarianter
verbundener Kern gerade die beiden primitiven zweiten Antworten und die
Tensorinformation trägt. Scheitert diese Abbildung oder erscheint ein dritter
gleichrangiger Kanal, bleibt die geometrische Auswahl in `wellenpaar.md` eine
Konditionierung.

## Quellen und Ausführung

Verwendet wurden der registrierte aktuelle Z174 Originalquelltext
[`z174_standardform_statistik_komposition.py`](../2026-11-23-phasenraum-f5-holocode/z174_standardform_statistik_komposition.py),
seine [Ausgabe](../2026-11-23-phasenraum-f5-holocode/ausgabe_z174.txt), der
aktuelle [TFPT Stand](../../../docs/research/TFPT-STAND.md) und der zu prüfende
[Wellenpaarvertrag](wellenpaar.md). Es wurde kein neuer Hamiltonoperator als
Quellbefund eingeführt.

Der Test läuft mit

```text
python3 -B experiments/exploration/2026-11-26-unabhaengige-ursprungsprinzipien/standardform_paar.py
```

und endet mit `STANDARDFORM PAIR TEST: PASS`.
