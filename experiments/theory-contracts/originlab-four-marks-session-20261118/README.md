# OriginLab · Sitzung zu den vier Marken · 17./18. November 2026

Prüfskripte der Sitzung, verschoben aus
`research/originlab/runs/session-scenarios-20261117/` (dort von Git
ausgeschlossen). Nichts hier ist ein Ledger- oder Atlas-Eintrag.

- Erläuterung aller Runden: `research/originlab/THEORY-BRIEF.md`
- Vorschlag für einen Atlas-Eintrag:
  `research/originlab/ATLAS-PROPOSAL-20261118.md`
- Prüfsummen und Ergebnis jedes Laufs: `MANIFEST.txt`
- Ausgaben: `output/<Skriptname>.txt`

## Ausführen

```
sh run_all.sh <python> <Wurzel des Repositorys>
```

Benötigt werden `sympy`, `numpy`, `mpmath` und `scipy`. Vier Skripte
lesen Originale aus dem Repository und erwarten dessen Wurzel als
Argument (`check_monodromy_e8.py`, `check_pasted_and_hypotheses.py`,
`check_native_consumer_control.py`, `check_scalar_marker_control.py`);
die übrigen ignorieren das Argument. `cover_lattice.py` ist eine
Bibliothek für mehrere Skripte.

## Geltung

Die Rechnungen sind exakt oder mit hoher Stellenzahl geprüft. Sie
betreffen aufgezeichnete Objekte und die freie Quelle. Keine davon
leitet P1 her; F1 bleibt offen. Welche Aussagen belegt, bedingt oder
ausgeschlossen sind, steht im Vorschlag.
