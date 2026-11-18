# Originaldecoder: Produkt, Zeit und alle unitären Erasure-Lifts

2026-10-16. Nachtrag zu `UR.NATIVE.JOINT_OPERATOR.01`, weiterhin PARTIAL.
Vorausgesetzt sind nur der tatsächliche vorhandene Erasure-Intertwiner
W:C5⊗bar(C4)→C64 mit W*W=I20, der originale Decoder und eine ausdrücklich
gegebene logische Operation. Keine neue Feldidentifikation oder Kopplung.

## 1. Die Komposition muss nicht zusätzlich eingesetzt werden

Setze Pi=WW* und D(A)=W(A⊗I4)W*. Dies ist dieselbe Rückabbildung wie in der
Originalquelle `_newest2/TFPT_Gesamtdokumentation2_20260927.md:8694–8697`;
ihre Erasure-Koordinaten wurden im bestehenden Vertrag bestimmt. Dann gilt
für beliebige A,B in M5 unmittelbar aus W*W=I20:

    D(A)D(B)=D(AB),     D(A*)=D(A)*,     D(I5)=Pi.

Außerdem ist D injektiv, weil W*D(A)W=A⊗I4. D ist also eine treue
*-Darstellung in der durch Pi bestimmten Unteralgebra mit Einheit Pi.
Es ist keine unital auf ganz M64 abgebildete Matrixalgebra: Pi ist nicht I64.
Diese Unterscheidung verhindert eine falsche globale Identifikation.

Für unitäres U ist die Erweiterung mit festgehaltenem Komplement

    Uhat=D(U)+I64-Pi

unitär. Weil D(A)(I-Pi)=(I-Pi)D(A)=0, gilt exakt

    Uhat Vhat = (UV)hat.

Jede endliche Folge vorhandener logischer Operationen besitzt damit dieselbe
dekodierte Produktwirkung. Zusätzliche Zwischenprojektionen oder Messregeln
werden für diese Identität nicht eingeführt. Das Festhalten des Komplements
ist eine angegebene Erweiterung; die physische Quelle könnte dort andere
Operationen erlauben. Der Decoder allein wählt diese Physik nicht aus.

Für einen gegebenen hermiteschen H5 gilt entsprechend

    exp(-it D(H5))=D(exp(-it H5))+I64-Pi.

Auf dem kodierten Träger ist die Entwicklung genau exp(-it H5)⊗I4. Für den
gegebenen diskreten Vierteloperator gilt Rhat^n=D(R5^n)+I64-Pi. Keine dieser
Identitäten wählt H5, einen Logarithmuszweig von R5 oder eine physische
Zeitskala aus. Sie liefern die Zeitverkettung relativ zu bereits gegebenen
Operationen, keine Herkunft der erlaubten Eingriffe.

## 2. Sämtliche unitären Lifts mit voller branchweiser Erasure-Verträglichkeit

Sei R eine gegebene logische Unitäre und T eine beliebige Unitäre auf
C5⊗C4. Fordere nur die vollständige branchweise Erasure-Gleichung

    T(rho⊗I4/4)T*=(R rho R*)⊗I4/4   für jeden Zustand rho.

Diese Bedingung erzwingt bereits die Faktorform; sie wird nicht als Ansatz
eingesetzt. Beweis: Setze S=(R*⊗I4)T. Die Gleichung besagt
S(rho⊗I4)S*=rho⊗I4. Zustände spannen M5 komplex linear auf, also kommutiert
S mit M5⊗I4. Dessen Kommutant ist I5⊗M4. Folglich gilt genau

    T=R⊗V,                 V in U4.

Die Umkehrung ist unmittelbar, da V I4 V*=I4. Dies klassifiziert alle
unitären Lifts auf dem Erasure-Träger, nicht sämtliche CPTP-Kanäle oder
physisch erlaubten Instrumente. Der Originaldecoder gibt V=I4; Erasure
allein lässt V frei. Der erste Nachtrag widerlegt das direkte native
Viertelpolynom bereits auf einem Zustand. Weil jede V-Variante denselben
Erasure-Ausgang hat, kann eine andere Wahl von V diesen Fehler nicht beheben.

Auf dem vollständigen Dreiregisterraum erzwingt bereits der Eingang I5/5
die Erhaltung von Pi, da sein Erasure-Zustand auf genau diesem Träger treu ist.
Alle unitären Erasure-Lifts auf C64 haben somit die Blockform
W(R⊗V)W* plus eine beliebige Unitäre auf dem 44-dimensionalen Komplement.

## 3. Der genaue Unterschied zur nativen Gruppenoperation

Aus dem ursprünglichen Intertwiner gilt

    W* g^⊗3 W = rho5(g)⊗bar(g),       g in G31.

Der Originaldecoder liefert stattdessen rho5(g)⊗I4. Die Differenz ist exakt
die zusätzliche Gauge-Kompensation I5⊗bar(g)* auf diesem Träger. Sie gehört
als Matrix zum vorhandenen gemeinsamen Algebraabschluss. Ihre physische
Erlaubnis folgt nicht aus der bloßen Matrixzugehörigkeit. Insbesondere ist
ein physischer Zugriff auf diesen Faktor nicht automatisch ein ursprünglicher
Gruppenwort-Eingriff auf den drei Registern.

Die native Operation rho5(g)⊗bar(g) erfüllt selbst bereits die vollständige
Erasure-Gleichung für R=rho5(g). Die Gauge-Kompensation wird nur benötigt,
wenn genau der Decoder-Repräsentant V=I4 umgesetzt werden soll; sie ist keine
zusätzliche notwendige Bedingung der Erasure-Verträglichkeit.

Damit ist die nächste Herkunftsfrage kleiner geworden: Die mathematische
Komposition und die unitäre Erasure-Liftklasse sind bestimmt. Zu bestimmen
bleibt, welche kollektiven Decoder-/Gaugeoperationen die Originalquelle
tatsächlich erlaubt, wie sie deren Generator und Zeitskala auswählt und ob
sie die positive Polarregel begründet. Das Ergebnis ersetzt keine dieser
Auswahlen durch einen frei erklärten Hamiltonoperator.

## Beweisprüfung

Die Produkt- und Zeitidentitäten folgen algebraisch allein aus W*W=I20;
der Wert und die tatsächliche Quellengleichheit von W wurden im vorhandenen
exakten `marked_lift_checker.py` geprüft. Der Kommutantenschluss ist ein
analytischer Satz für alle Eingaben, keine endliche Zustandsstichprobe.
Der Anwendungsumfang ist der bereits gegebene symmetrische Code und
Originaldecoder. Kein P1-Feld, keine physische Kopplung und keine unabhängige
empirische Vorhersage werden hier behauptet.
