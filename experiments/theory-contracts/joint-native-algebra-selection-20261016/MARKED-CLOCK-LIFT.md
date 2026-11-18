# Gemeinsamer Ursprung: tatsächliche Markierung und Clocklift

2026-10-16. Fortsetzung von `UR.NATIVE.JOINT_OPERATOR.01`, Verdict **PARTIAL**.
Dieser Nachtrag ergänzt `PROOF.md`, ohne dessen belegte Aussagen umzuschreiben.
Die physische P1-Quellenauswahl, Felder und Kopplungsnormierungen bleiben offen.

## 1. Ergebnis und Herkunft

Die vorhandene 3+2-Markierung lässt sich ohne Standardmodell-Ladungen aus dem
ursprünglichen Simplex, seinem Dreier-Clock, einem markierten Fixpunkt und der
vorhandenen positiven Polarregel rekonstruieren. Die drei Fixpunktwahlen sind
unter dem tatsächlichen nativen Clockzentralisator äquivalent. Die Polarregel
selbst wird dabei vorausgesetzt; sie folgt nicht aus dieser Rechnung aus P1.

Der vorhandene Vierteloperator hat dagegen keinen codeerhaltenden Lift als
Produkt beliebiger vier lokaler Unitaries. Auch klassische Mischungen nativer
Gruppenwortkanäle erzeugen seinen vollständigen logischen Kanal nicht.
Die ursprüngliche Erasure-Rückabbildung liefert einen kollektiven logischen
Operator. Das Einsetzen der tatsächlichen nativen Clockwirkung in dasselbe
Viertelpolynom ergibt zwar eine unitäre Operation mit richtigem Quadrat,
verletzt aber diese Erasure-Wirkung: Ein exakter Eingang hat Spurabstand 1/2.

Originalstellen: `_newest2/TFPT_Gesamtdokumentation2_20260927.md`, Zeilen
8628–8697, insbesondere der positive Polartransport und die tatsächliche
Dreiregister-Rückabbildung am Ende; die zugehörigen Konstruktoren sind
`process._code_basis`, `marked_source_process._source_operators`,
`marker_selection.build_marker_data` und `charge_selection._process_frame`.
Der frühere markierte R-Abschluss ist eine kohärente Polynomkonstruktion,
keine bereits aus P1 bestimmte primitive Intervention.

## 2. Derselbe Rahmen auf allen Trägern

Sei C der vorhandene 6×5-Prozessframe, G=diag(4,12,12,12,24) die originale
Rohmetrik und W6 der dortige 5×6-Simplex. Die nativen Spalten sind
N=G^(1/2)W6/sqrt(48). Exakt gilt

    N=O C^T,       O=-C^T P_reverse C,       det(O)=1.

Der unveränderte originale Vierer-Clock sigma4 komprimiert im N-Rahmen zu
(2 3 4), im C-Rahmen zu (1 3 2). Ebenso ist (0 1)_N gleich (4 5)_C und der
markierte Punkt 5_N gleich 0_C. Alle Labels sind hier nullbasiert.
Nur Y mit O zu transportieren und die nominalen N-Permutationen im C-Rahmen
beizubehalten, behält diesen unveränderten sigma4 nicht bei.

Die quellenkonsistente Variante ist R_source=O R_previous_C O^T. Weil
O in U(1)rho5(G31) liegt, bleiben die früheren abstrakten Nichtwort-,
Stabilisator- und Kommutantenresultate erhalten. Für die neue gemeinsame
Clockprüfung werden durchgehend N-Labels und der unveränderte sigma4 benutzt.
Dies ist eine explizite Rahmentransformation, keine neue physische Zuordnung.

## 3. Geschlossene Form der tatsächlichen 3+2-Markierung

Seien n_j die Spalten von N. Ihre Gram-Matrix ist I6-J6/6. Setze

    sigma=sigma5,       F=(I+sigma+sigma²)/3,
    m=n_2+n_3+n_4,      q in {0,1,5}.

F hat Rang 3, sein Komplement Q=I-F Rang 2. Die Originalregel nimmt
D_q=(n_j-n_q)_(j≠q) und U_q=D_q[I+(1/sqrt(6)-1)J5/5]. Aus der angegebenen
Gram-Matrix folgt U_q^* U_q=I5. Die drei bewegten Slots bilden unter U_q
den Farbraum. Ihr konstanter Einheitsvektor wird zu

    c_q=[m-3(sqrt(6)-1)n_q/5]/sqrt(3).

Beweis: Für den konstanten bewegten Slotvektor v gilt
D_q v=(m-3n_q)/sqrt(3), D_q 1=-6n_q und 1^T v=sqrt(3).
Damit ergibt das ursprüngliche Polarpolynom genau die obige Formel.
c_q ist normiert und liegt in Fix(sigma). Die anderen beiden Farbrichtungen
sind gerade Q. Deshalb ist die originale Markierung exakt

    P2(q)=F-c_q c_q^*,          P3(q)=I-P2(q).

Es werden keine Ladungswerte verwendet. Für q=5 stimmt dies mit der
vorhandenen Markierung (6Y5+2I)/5 überein. Alle drei P2 sind Rang-2-Projektoren.
Für jede Transposition der drei Fixpunkte existiert in ihrem viergliedrigen
nativen Reflexionslift genau ein mit dem tatsächlichen sigma4 kommutierender
Lift. Seine Fünferwirkung vertauscht die P2(q) wie die Fixpunkte. Somit ist
diese Familie **eine native Clockzentralisator-Bahn**. Die Wahl eines q erzeugt
hier keine drei inequivalenten Markierungen. Daraus folgt keine Invarianz
einer beliebigen anisotropen Quelle oder eines physischen Zustands.

Die Orientierung von P2 relativ zu Simplex und Clock ist damit berechnet.
Die positive Polarvorschrift, der symmetrische Code und die physische
Interpretation seiner beiden Teilräume werden dadurch nicht aus P1 ausgewählt.
Insbesondere bleibt P2 verschieden von allen 15 nativen Matching-Projektoren.

## 4. Vollständiger Ausschluss von Produktunitär-Lifts

Für die tatsächliche Isometrie V:C5→(C4)^⊗4 gilt für jedes Registerpaar

    Tr_others(VV^*)=(I16+Swap)/4.

Erhält U1⊗U2⊗U3⊗U4 den Code, muss jedes Ui⊗Uj diese Paarmatrix erhalten.
Aus der Swap-Identität folgt (Ui Uj^*)⊗(Uj Ui^*)=I16. Somit sind Ui und Uj
bis auf Phase gleich. Jeder solche Produktoperator ist also Phase×U^⊗4.

Der tatsächliche identische lokale Stabilisator ist K=U(1)G31. Beweis:
Für die 15 spurfreien hermiteschen Paulis ist der exakte Leckage-Gram 32 I15;
die Lie-Algebra von K ist daher allein iR I4. Die Determinantenscheibe
F_K={U∈K:det(U)=±1} ist endlich und enthält die tatsächliche komplexe
Cliffordgruppe X2=<G31,zeta8 I>. Nach dem endlichen Übergruppensatz von
Nebe–Rains–Sloane ist sie nur um skalare Wurzeln größer. Durch skalare
Determinantenkorrektur folgt K=U(1)G31. Der Satz wird auf F_K angewandt,
nicht auf die unendliche Gruppe K. Die tatsächliche Clifford-Identifikation
ist eine ausdrückliche vorhandene Beweisvoraussetzung.

Primärquelle: [Theorem 6.5, S.18; X2/G31, S.19](https://arxiv.org/pdf/math/0001038v2#page=18).
Die vollständige Beweisprüfung steht in `transversal_audit.md`.

Nun ist R5=P3(I-2F)+iP2 tau, mit tau=(0 1)_N. Es hat Spur 1 und
Spektrum {1,1,-1,i,-i}. Die nativen S6-Matrizen sind reell mit ganzzahliger
Spur. Eine Gleichheit R5=lambda rho5(g), |lambda|=1, würde durch die Spur
lambda=±1 erzwingen. Das widerspricht ||Im R5||²_HS=2. Kein codeerhaltender
Produktunitär-Lift existiert. Auch eine konvexe Mischung nativer Wortkanäle
ist ausgeschlossen: Der Zielkanal Ad_R5 hat Choi-Rang 1; jeder positiv
gewichtete Summand müsste derselbe Kanal und damit ein ausgeschlossenes Wort
bis auf Phase sein. Dies gilt für den vollständigen Kanal, nicht nur einzelne
ausgewählte Eingangszustände.

Kollektive Operationen, Ancillas, allgemeines LOCC, Postselektion und
Protokolle mit zeitweiligem Verlassen des Codes werden hier nicht ausgeschlossen.

## 5. Originaldecoder und tatsächliches gemeinsames Clockpolynom

Die originale Dreiregister-Rückabbildung ist, in den bereits bewiesenen
Erasure-Koordinaten W:C5⊗bar(C4)→(C4)^⊗3, genau

    D(A)=W(A⊗I4)W^*,       D(I5)=WW^*.

Damit wird Y5 zu Y5⊗I4 und der logische Viertel-Clock zu R5⊗I4.
Die unitäre Dreiregister-Erweiterung ist W(R5⊗I4)W^*+I64-WW^*.
Diese Wirkung folgt aus dem **vorhandenen Decoder**; ihre physische Erlaubnis
als primitive Intervention folgt daraus nicht.

Die tatsächliche native Gruppenwirkung auf dem Zwanziger lautet dagegen

    sigma20=sigma5⊗bar(sigma4),   tau20=tau5⊗bar(tau4).

Verwendet man die durch denselben Decoder festgelegte Markierung P220=P2⊗I4,
F20=(I+sigma20+sigma20²)/3 und dasselbe vorhandene R-Polynom, erhält man

    R20_native=(I-P220)(I-2F20)+iP220 tau20.

Für alle vier ursprünglichen Reflexionslifts von tau5 gilt exakt:
R20_native ist unitär und R20_native²=I-2P220. Trotzdem gilt nicht
Ad_R20_native E=E Ad_R5 für E(rho)=rho⊗I4/4.

Ein entscheidender Eingang ist rho=Q(I+s23)/2 mit s23=(2 3)_N und Q=I-F.
Er ist ein reiner Zustand im bewegten Standardzweiraum, P2 rho=0 und
R5 rho R5^*=rho. Die vier Gauge-Clocklinien von bar(sigma4) haben Eigenwerte
1,1,omega,omega². Zwei lassen rho unverändert; die beiden anderen spiegeln
ihn auf denselben orthogonalen Zustand. Die Ausgangsdifferenz hat Eigenwerte
{+1/4,+1/4,-1/4,-1/4,0×16}, Rang 4 und Spurabstand 1/2.
Der tau-Term wirkt auf diesem Eingang nicht; deshalb gilt das Ergebnis für
alle vier tatsächlichen tau4-Lifts. Ein richtiges Quadrat ist hier folglich
keine ausreichende Übereinstimmung der tatsächlichen Operationswirkung.

Der Test betrifft genau die angegebene originale Decoder-Markierung und das
gleiche native Clockpolynom. Er schließt weder jeden kohärenten gemeinsamen
Lift noch beliebige andere Definitionen einer Zwanzigermarkierung aus.

## 6. Konsequenz für die laufende Herleitung

Die gemeinsame native Struktur liefert jetzt zusätzlich die tatsächliche
Markierungsorientierung relativ zum Compiler-Clock und ihre äquivalente
Fixpunktfamilie. Der Viertel-Clock benötigt für diesen Code eine kollektive
Wirkung; lokale Produktoperationen und bloßes klassisches Mischen reichen
nachweislich nicht. Der vorhandene Decoder liefert ihre logische Form,
aber das direkte native Polynom liefert eine andere Operation.

Der konkrete verbleibende Quellenschritt ist deshalb die ursprüngliche
Erlaubnis und Auswahl dieser kollektiven Decoderwirkung einschließlich ihrer
Zeitverkettung und Normierung. Ein freier kollektiver Hamiltonoperator oder
eine frei erklärte kohärente Gruppenalgebra schließt diese Herkunftslücke
nicht. Aus diesem Nachtrag folgen keine neue Kopplung, kein ausgewählter
physischer Zustand und keine neue alpha-/Flavor-Normierung.

Die exakten Checker heißen `marked_lift_checker.py` und
`transversal_checker.py`; ihre Ergebnisdateien sind gleichnamig mit Suffix
`_results.json` bzw. `transversal_checker.json`. Numerische Toleranzen werden
für die angegebenen Identitäten nicht benutzt. Der Übergruppensatz bleibt
eine externe Beweisvoraussetzung, kein Ergebnis eines endlichen Checkers.
