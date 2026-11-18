# Originaler R064 Kanal: Spektrum und kontinuierlicher Parent

Ausführen im TFPT Repository:

```
python3 experiments/theory-contracts/compiler-correlated-event-clock-20260919/synthesis_20261102/check_source_parent.py
```

Der Check lädt unverändert den vorhandenen `source_channel.py` Builder, erzeugt die 60 originalen Gaussian Code Strahlen und die zugehörigen unitären Involutionen U auf Sym²(C⁴) und verwendet deren originalen Basiswechsel. Er betrifft den unraffinierten Kanal Φ(X)=Σ U X U*/60. Dieser ist nicht der in R088 untersuchte fünfgliedrige Hodge Kanal auf Λ²(C⁵); gleiche Trägerdimension ersetzt kein Wörterbuch.

Für die 100×60 Matrix M der vektorisierten U ist G=M* M reell ganzzahlig. Die int64 Berechnung prüft mit expliziter Überlaufschranke exakt `G(G−8)(G−16)(G−96)=0` und die ersten drei Spuren. Wegen einfacher Polynomwurzeln bestimmen diese die Multiplizitäten 5,45,9,1. Das normierte Choiobjekt MM*/600 hat folglich 0⁴⁵, (1/75)⁴⁵, (2/75)⁹ und (4/25)¹. Der normierte Identitätsvektor ist wegen ΣU=24I Eigenvektor zum größten Eigenwert. Für jeden unitären Kanal ergibt die Projektion auf dessen reinen Choivektor den halben Diamantabstand mindestens 1−4/25; die Zerlegung Φ=(4/25)Id+(21/25)Ψ liefert die passende obere Schranke. Der kleinste Abstand ist somit exakt 21/25.

Für die 3.600 Produkte UiUj ist der Rang mindestens 100 über C: Nach Abbildung i→10 in F_101 ergibt eine explizite modulare Zeilenreduktion Rang 100. Zusammen mit dem obigen Rang 55 des ursprünglichen Choioperators bestätigt dies den Rangwechsel bei zwei Schritten. Die Konsequenz für zeitlich homogene Lindbladhalbgruppen benötigt zusätzlich den entsprechenden Halbgruppensatz; der Test selbst beweist den Rangwechsel.

Ein gesonderter expliziter Parent nimmt einen zusätzlichen 60 dimensionalen Recordraum und |+>=Σz|z>/√60. Für W=Σz Uz⊗|z><z| gilt W*=W und W²=I; H=π(I−W)/(2τ) erzeugt zum Zeitpunkt τ genau W. Der Check verifiziert die Recordoperatoren, ihre Spurtreue, den Zwischenzeitkanal und das Originalinstrument bei τ. Dieselben Records erneut gekoppelt ergeben den Ausgangszustand zurück, frische Records dagegen Φ². Für die im Check benutzte reine Eingabe beträgt der Frobeniusabstand 0,843274042712. Andere Eingaben können andere Werte liefern.

Die positiven Hamiltonianformeln sind algebraisch; der Quellencode verifiziert die angegebenen Endpunkt und Zwischenzeitformeln numerisch. Recordpräparation, Zellvorrat, Ereigniszeitskala und Zugehörigkeit zur ursprünglichen TFPT Quelle werden nicht hergeleitet. Ein kontinuierlicher Parent zeigt Existenz, keine physische Auswahl. Die Bibliotheksverweise der Nutzervorlage werden hier nicht als lokal ausgeführtes Notebook bezeichnet.
