# Quellengebundene Regelsuche vom 30.09.2026

**Verdict: PARTIAL. Die vollständige gemeinsame Ursprungsregel wurde nicht gefunden.**

Der Nutzerauftrag war, das hinter den gemeinsamen TFPT-Antworten liegende Gesetz mit automatischer Suche zu suchen. Dieser Durchgang untersucht ausdrücklich benannte Operator-/Filterklassen an der vorhandenen Quelle; er ersetzt das Ursprungsproblem nicht durch ein neues Hilfsmodell. Das finite Q/P-Wörterbuch und dessen Auswahl von eta=-3 werden als bestehende Randbedingungen behandelt.

Die exakte Suche liefert einen eindeutigen inversen Modenfilter in einer deklarierten Monomklasse. Die gemeinsame Zeit-/Markprüfung widerlegt dessen wörtliche Ein-Port-Auswertung als vollständigen Alpha-Anschluss. Getrennte Auslesungen derselben Quelle bleiben konsistent. Die zweite Suche schließt den benannten einzelnen neutralen Vertex mit unveränderter Familienkomponente, Norm48 und nichtverschwindender Hyperladungs-Wardantwort aus. Die bereits vorhandenen quelleneigenen Richtungen R und SY tragen getrennte Antworten; ihre erneute Berechnung ist kein neuer Herkunftsbeweis.

## Ausführen

Python 3 mit den vorhandenen TFPT-Projektabhängigkeiten (insbesondere NumPy und SymPy); für die Originalmodule muss das TFPT-Repository verfügbar sein. Keine neuen Pakete wurden installiert.

```sh
TFPT_REPO=/Users/stefanhamann/Projekte/tfpt-theoryv4 python3 checker.py
```

Der Checker berechnet die tatsächlichen Quellwurzelmomente und schweren Paarantworten neu, bestätigt die R-Zustandsnorm gegen die Laurent-Modennorm und schreibt die drei Teilergebnisse sowie `results/summary.json`. Ein PASS der ausführbaren Kontrollen bestätigt die mathematischen Befunde innerhalb ihres Scopes, keinen Abschluss der TFPT-Herleitung.

## Erste zusätzliche Wahlen

1. Die Modenfilter-Grammatik und das Ziel 1/n werden für die Suche vorgegeben.
2. Für den bedingten Markvergleich werden vier Kreispositionen, ein gemeinsamer Diagonal-Referenzabzug und die positive Sewing-Zeit alpha gewählt.
3. Für die Ein-Vertex-Suche werden die Bedingungen Familienkomponente1 und Norm48 gesetzt. Ihre gleichzeitige physikalische Notwendigkeit wurde hier nicht abgeleitet.
4. Für den positiven Mehr-Port-Kontrollfall werden die Grad2-Projektion und eine separate Rand-Auslesung eingesetzt. Deren gemeinsame Auswahl aus P1/P2 fehlt.

## Quellen und Beweisgrenzen

Originalmodule: `tfpt_explorer/neutral_source_response.py`, `source_charge_response.py`, `source_majorana_pairs.py`, `source_flavor_correlator.py`; Quellhashes in der Zusammenfassung. Kontext: vorhandene Contracts `source-flavor-joint-dictionary-20260921`, `source-native-clock-flavor-bridge-20260922` und `source-variation-origin-gate-20260922`.

Die vollständigen Rechnungen und Geltungsbereiche stehen in `FILTER_PROOF.md`, `PORT_PROOF.md`, `VERTEX_PROOF.md`. Keine globale TFPT-Unmöglichkeit, keine Auswahl des geladenen Collar-Operators aus P1/P2, kein physikalischer 4D-Photon-Kern und keine Auswahl der Massenkopplungen werden beansprucht. Der nächste gemeinsame Suchgegenstand bleibt die source-owned Zustands-/Kopplungs-/Randregel, welche diese Antworten ohne nachträgliche Anpassung auswählt.

**Firewall: experiments; keine Promotion, keine empirische Scorecard, keine Änderung von Ledger, Verifikationssuite oder Papers.**
