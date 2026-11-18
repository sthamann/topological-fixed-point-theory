# Welche ursprüngliche Dynamik wählt Markierung und Viertel-Clock?

2026-10-16. Fortsetzung von `UR.NATIVE.JOINT_OPERATOR.01`, Verdict PARTIAL.
Der gemeinsame native Vierer-/Quartik-/Erasure-Kern bleibt erhalten.
Untersucht werden die tatsächliche Markierung und zwei schon vorhandene
Quellkonstruktionen; keine neue passende Kopplung wird eingesetzt.

## 1. Ergebnis

Die richtige Viertelrelation, der markierte Clockfixpunkt und die native
Symmetrieverträglichkeit bestimmen den positiven Polartransport nicht.
Die ausdrücklich geprüfte Klasse besitzt eine vollständige CP1-Faser von
Markierungen, bei reellen Projektoren RP1. Die originale Polarregel wählt
darin einen konkreten Punkt. Der Parameter wird bei festgehaltenem nativen
Clock und Quellpunkt nicht durch einen nativen Framewechsel entfernt.

Die Polarregel besitzt eine eindeutige geometrische Charakterisierung:
Sie ist die nächstgelegene Isometrie zur ursprünglichen Differenzmatrix im
ungewichteten Frobeniusabstand. Dies erklärt ihre mathematische Kanonizität.
Eine aus P1 stammende physische Minimierung gerade dieses Abstands folgt
daraus nicht. Der ursprüngliche Moment-Hamiltonoperator kann diese Auswahl
nicht leisten: Er ist auf dem gesamten Fünfercode null. Die vorhandene
geometrische Verbindung lokaler Produkt-Codewege liefert ebenfalls keinen
markierten Viertel-Clock; ihre Holonomie bleibt projektiv in S6.

Das sind eng bestimmte Aussagen über vorhandene Operatorbedingungen und
Quellkandidaten. Die alternativen Markierungen werden dadurch weder als
vollständige P1-Quellenmodelle noch als physisch erlaubte Eingriffe ausgegeben.

## 2. Vollständige Markierungsfaser im tatsächlichen Originalrahmen

Benutze unverändert die nativen Simplexspalten n_j des vorigen Nachtrags,
mit n_i* n_j=delta_ij-1/6, sigma=(234), tau=(01), q=5 und
F=(I+sigma+sigma²)/3. Setze m=n2+n3+n4 und

    a=(n0-n1)/sqrt(2),
    e0=sqrt(6/5)n5,
    e1=sqrt(5/6)[m+3n5/5].

Diese drei Vektoren sind eine orthonormale Basis von Fix(sigma). Tau wirkt
darauf als diag(-1,+1,+1). Die negative Linie ist also fest; die positive
Ebene E+=span(e0,e1) ist zweidimensional.

Klassifiziere hermitesche orthogonale Rang-2-Projektoren P mit Ran(P) in
Fix(sigma), [P,tau]=0 und je einem tau-Eigenwert +1 und -1 auf Ran(P).
Die negative Linie Ca muss enthalten sein. In E+ muss genau eine Linie
verbleiben. Damit sind ALLE Projektoren dieser Klasse exakt

    P=F-cc*,     [c] in CP1(E+).

Für reelle orthogonale Projektoren ist die Klasse RP1, weiterhin nicht ein
einziger Punkt. Ein vollständiger komplexer Chart ist
c proportional e0+(x+iy)e1, x,y reell, ergänzt um c=e1 bei Unendlich.
Es wird kein Raum außerhalb der tatsächlichen nativen Darstellung eingeführt.

Für jedes dieser P hat der bestehende Ausdruck

    R(P)=(I-P)(I-2F)+i P tau

dieselben orthogonalen Wirkungen: +I auf Ran(I-F), -1 auf Cc, -i auf Ca
und +i auf E+ orthogonal zu c. Also gelten für die ganze Klasse

    spec R(P)={1,1,-1,i,-i},       R(P)²=I-2P,       R(P)^4=I.

Der ursprüngliche positive Polartransport ist genau der Punkt

    c_pol=-sqrt(3/5)e0+sqrt(2/5)e1,
    x_pol=-sqrt(2/3), y_pol=0.

Die reelle Gegenmark c=e0 erfüllt dieselben genannten Clockbedingungen.
Sie ist von P_pol operatorisch verschieden:

    ||P_pol-P_e0||²_HS=4/5,
    e0* P_pol e0=2/5,            e0* P_e0 e0=0.

Diese Größen sind mathematische Antworten auf dem bereits vorhandenen
markierten Strahl; es wird keine tatsächliche physische Projektormessung
oder Präparation als aus P1 erlaubt behauptet.

Ersetze q=5 durch q in {0,1,5} in e0(q),e1(q). Bei festem CP1-Parameter
ist die ganze Familie unter dem tatsächlichen nativen Clockzentralisator-
Quotienten S3_fix äquivariant. Der native Stabilisator eines festgehaltenen
q und sigma ist auf E+ trivial; selbst dessen größerer logischer
Normalisator S3_color×S2_otherfixed wirkt auf e0,e1 als Identität. Globale
Phasen wirken skalar. Sie beseitigen den Faserparameter nicht. Ein beliebig
erweiterter kontinuierlicher logischer Zentralisator wäre eine andere
Äquivalenzregel und wird hier nicht als ursprünglich gegeben eingesetzt.

## 3. Was die positive Polarregel zusätzlich festlegt

Die ursprüngliche Differenzmatrix D_q=(n_j-n_q)_(j≠q) erfüllt

    D_q* D_q=I5+J5,
    S=(D_q* D_q)^(1/2)=I5+(sqrt(6)-1)J5/5,
    U_pol=D_q S^(-1).

Für jede unitäre Matrix U ist mit Z=U_pol* U

    ||D_q-U||²_F=15-2 Re tr(S Z) >= 7-2sqrt(6).

Denn S hat positive Eigenwerte 1,1,1,1,sqrt(6), und Re Z_ii<=1 in seiner
Eigenbasis. Gleichheit verlangt alle Z_ii=1 und damit Z=I. Also ist U_pol
der eindeutige Minimierer. Das ist eine elementare Anwendung der üblichen
Polar-/Procrustes-Eigenschaft; siehe [Higham 1986, Theorem 2.2 und Corollary 2.3](https://eprints.maths.manchester.ac.uk/694/1/high86p.pdf#page=4).

Diese Aussage ersetzt keine physische Herkunftsbedingung: Ohne Nachweis,
dass die Originalwirkung gerade diese Abstandsminimierung auf denselben
Slot-/Codeabbildungen erzeugt, bleibt dies die Charakterisierung einer
vorhandenen geometrischen Vorschrift. Reflexionspositivität eines Randkerns
ist eine Aussage über einen anderen Operator-/Funktionaltyp; sie ist nicht
bereits die Bedingung U*D_q>0 eines endlichen Polartransports.

## 4. Der tatsächliche Moment-Hamiltonoperator trifft diese Auswahl nicht

Die Originaldokumentation gibt 40 M4=Pi5+Psym4 und

    H_nat=I-(Pi5+Psym4)/2=I-20M4.

Die tatsächlichen Codevektoren liegen im Bild beider Projektoren. Daher

    H_nat V=0,        exp(-it H_nat)V=V für alle reellen t.

Alle obigen Markierungen liegen in diesem selben fünfdimensionalen
Grundraum und haben dieselbe Energie null. Das ist keine neue Familie
ausgewählter Hamiltonians, sondern ein Vergleich am vorhandenen Parent.
Auf E+ wirken zusätzlich sigma und tau jeweils als I2. Deshalb wirkt die
von H_nat|Code, sigma und tau erzeugte Algebra auf E+ nur skalar; sie kann
c_pol und seine orthogonale Richtung nicht mit den benötigten verschiedenen
Viertelphasen versehen. Auch die zusätzliche Forderung minimaler Energie
dieses Parents hebt den Faserfreiheitsgrad nicht auf.

Dies schließt keine andere hergeleitete Wechselwirkung, zeitabhängige
kombinierte Antriebe mit zeitweiligem Codeverlust oder allgemeinen kollektiven
Prozess aus. Der physische Einsatz des Moment-Parents selbst wäre ebenfalls
eine zusätzliche Identifikation, keine Folgerung allein aus seiner Existenz.

## 5. Vorhandenen geometrischen Transport wiederverwenden

Der ursprüngliche Vertrag `compiler-quartic-holonomy-20260919/PROOF.txt`
beweist bereits V* dGamma4(X)V=tr(X)I5 und für gemeinsame lokale Wege
U(t)^⊗4 die Verbindung A=tr(U* dU)I5. Diese Rechnung wird nicht neu als
Entdeckung ausgegeben. Der inzwischen vollständig bestimmte lokale
Stabilisator K=U(1)G31 schließt nun die globale Holonomieklasse ab.

Für allgemeine Produktwege T(t)=U1(t)⊗...⊗U4(t) liefert die tatsächliche
Einregister-Erasure-Korrektur an jeder Stelle j sogar

    V* X_j V=(tr X)/4 I5,
    A=(1/4) sum_j tr(U_j* dU_j) I5.

Die Verbindung ist projektiv flach. Endet ein solcher Weg wieder beim
ursprünglichen Code, ist T(1) ein Produkt-Codeerhalter und wirkt nach dem
vorigen vollständigen Stabilisatorsatz logisch als Phase×rho5(g).
Paralleltransport fügt nur eine skalare Phase hinzu. Die Holonomie liegt
deshalb projektiv in S6; der tatsächliche R5 liegt außerhalb dieser Klasse.
Für gemeinsame U(t) liefert der alte Vertrag sogar die schärfere Form
det(g)^(-1)rho5(g).

Dies klassifiziert geometrischen Paralleltransport genau dieser lokalen
Produkt-Codewege mit der vorhandenen Verbindung. Es ist kein allgemeiner
Ausschluss kollektiver Codebewegungen, zusätzlicher Dynamik im Unterraum,
Ancillas oder physischer adiabatischer Protokolle mit anderen Hamiltonians.

## 6. Entscheidende nächste Originalquellenfrage

Die Auswahl kann nicht allein aus dem richtigen Viertelquadrat, dem
vorhandenen Moment-Parent oder dieser lokalen geometrischen Holonomie
kommen. Benötigt wird die tatsächliche Originalantwort, die im selben
E+-Zweiraum eine Orientierung bzw. kohärente Wirkung unterscheidet.
Eine konkret zu prüfende Herkunftsbehauptung wäre, dass die ursprüngliche
Wirkungsvariation die oben charakterisierte Polar-Abstandsregel liefert.
Diese Behauptung ist hier nicht eingesetzt oder bewiesen. Eine passende
neue Energie, ein Steuerfeld oder eine frei erlaubte Gruppenalgebra würden
die gesuchte Auswahl voraussetzen.

Quellen: Originaldokumentation vom 27.09., Zeilen 8528–8557 und 8671–8701;
voriger Nachtrag `MARKED-CLOCK-LIFT.md`; der vollständige alte geometrische
Vertrag, insbesondere Abschnitte 1,2 und5. Der neue Checker
`check_polar_fibre.py` prüft den ganzen rationalen CP1-Chart symbolisch,
seinen fehlenden Unendlichkeitspunkt, die Originalpolarstelle, die native
Fixpunktäquivarianz und sämtliche 64 Einregister-Matrixeinheiten exakt.
Die Vollständigkeit der Faser und der Holonomieschluss sind analytische
Sätze, keine Folgerungen aus einer endlichen Parameterstichprobe.
