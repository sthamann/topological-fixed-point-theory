# Naht-Erreichbarkeit und absolute Quellen-Normierung

**PARTIAL.** Der Test verbindet die vorhandene No-Spectator-Prüfung mit der
Quellen-Normierungsfrage. Am Originalmodell bleiben Erreichbarkeit und
Kommutant unter einer skalaren Energieverschiebung gleich. Eine absolute
Normierung folgt aus diesen Kriterien allein nicht.

Der ergänzende Basiswechsel-Test findet eine Grenze der ursprünglichen
Kommutanten-Hilfsfunktion: Sie ist für die verwendeten reellen Matrizen korrekt,
aber nicht für allgemeine komplexe Darstellungen. Die korrigierte komplexe
Nullraumrechnung liegt im Experiment-Checker.

Ausführen vom Repository-Verzeichnis:

```sh
python3 experiments/theory-contracts/spectator-normalization-separation-20261013/checker.py
```

Benötigt NumPy und die unveränderten Originaldateien aus `verification/`.
`results.json` enthält Parameter, Quellhash und tatsächliche Messwerte;
`PROOF.txt` die algebraische Begründung, Vorarbeiten und Grenzen.

**Experiments-Firewall:** Diese Studie erzeugt keine neue physikalische Quelle,
ändert keinen Ledgerstatus und liefert keine T1–T8-Gesamtlösung. Sie beweist
insbesondere nicht, dass die untersuchte Energieverschiebung innerhalb des
ursprünglichen vollständigen TFPT-Quellenfunktionals zulässig ist. Für die
festgehaltene Fermi-Energie wird keine Gleichheit der Besetzungsprojektoren
behauptet. Die ursprünglichen reellen v1020-Prüfwerte bleiben gültig.
