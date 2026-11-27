"""Knowledge map of the ledger: what is on record for every part of the theory, with a normalised class.

Why this exists (finding of 27 November 2026): the entry document of exploration sessions,
docs/research/TFPT-STAND.md, named 12 of 784 active ledger claims; the Atlas reported 547 of 789 claims as
"mixed" because the ledger status column is free text (724 distinct strings).  Whole branches of established work
(twistor bridge, 4D quantum field theory policy, gravity, predictions) were therefore invisible to a session that
followed the reading order of AGENTS.md, and an exploration could call something new that the ledger had decided.

The map is generated, never edited by hand.  It classifies nothing scientifically: the class is read off the
status text of the ledger, and the ledger with its original sources stays authoritative.

    ../../.venv/bin/python ledger_map.py --write          regenerate docs/research/LEDGER-KARTE.md and the Atlas summary
    ../../.venv/bin/python ledger_map.py --check          fail when either file is out of date
    ../../.venv/bin/python ledger_map.py --find 'twistor|birefring'
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / "verification" / "status_ledger.csv"
TARGET = ROOT / "docs" / "research" / "LEDGER-KARTE.md"
# compact per-area summary that the Atlas attaches to theory_overview, read live by the server
ATLAS_SUMMARY = ROOT / "semantic_explorer" / "ledger_areas.json"

PARTIAL = "exakt, Teil bedingt oder offen"
CLASSES = ["exakt", PARTIAL, "bedingt", "offen", "Kill-Test", "gesetzt", "Regel oder Synthese"]
EXACT = ("identity", "lattice", "numerical", "formal", "exact", "proven", "theorem", "float-track", "sheet generator",
         "bridge [e]", "abstract master theorem")
RULE = ("firewall", "decision", "branch policy", "guard", "analysis", "roadmap", "synthesis", "capstone", "audit",
        "discipline", "inventory", "method", "review", "correction", "registry", "policy")
COND = ("physical", "conditional", "candidate", "conjecture", "undecided", "success per", "data confrontation",
        "contract [c]", "bridge", "reading")

# reading order of the theory; a block that is not listed lands in the last area
AREAS = [
    ("Eingaben, Anker und Selbstkonsistenz", "AX ANCHOR ARCH BOOT P1 P2 TFPT OVERDET SYNTH FIRSTPRINCIPLES FOUR META "
                                              "CAR CARRIER K5 PG32 DOILY SEL CENTER"),
    ("Naht, Rand und Transfer", "SEAM BOUNDARY SHEET MARKLOCAL CMNORM DHR GNET QFT FTR FTRANSFER TRANSFER RATE RES "
                                "RESIDUAL DYN DYNAMICS DECOUPLING"),
    ("E8, Gitter, Uhren und Codes", "E8 COX COXETER MCKAY GAUSSIAN GOLD QEC ARF CURVE GALOIS ARITH H2 CAS CFIN FORM "
                                    "NORMALFORM DUAL"),
    ("Standardmodell, Eichstruktur und 4D", "SM PS GAUGE HIGGS QCD QFT4D CHIRAL4D TFPT4D PARENT DIMENSION GATE"),
    ("Feinstruktur und Elektromagnetismus", "ALPHA EM"),
    ("Flavor, Massen und CP", "FLAV FLAVOR CP LEP THETA13 SEED DIAMOND"),
    ("Twistor und Himmelskugel", "WOIT CELEST"),
    ("Gravitation, Horizont und Geometrie", "GRAV GRAVITY HOR QG QGEO QGAMB FLATAWAY CAUSAL BULK VGEO HYP SCALE"),
    ("Kosmologie, Vorhersagen und Daten", "COSMO LAMBDA DM CCC PRED DATA KILL AS OBS OBSERVATORY METROLOGY FR "
                                          "UVBRANCH"),
    ("Methodik, Pruefung und Sonstiges", ""),
]


def classify(status):
    s = status.strip().lower()
    head, span = s[:70], s[:160]
    if head.startswith("open"):
        return "offen"
    if head.startswith("kill"):
        return "Kill-Test"
    if head.startswith("axiom"):
        return "gesetzt"
    if head.startswith(RULE):
        return "Regel oder Synthese"
    if head.startswith(COND) and not head.startswith(EXACT):
        return "bedingt"
    exact = head.startswith(EXACT) or any(t in span for t in ("[e]", "[i]", "[l]", "[n]", "[f]", "[e-", "numerical",
                                                              "identity", "exact"))
    part = any(t in span for t in ("[c]", "[p]", "[o]", "conditional", "candidate", "conjecture", " open"))
    if exact:
        return PARTIAL if part else "exakt"
    if "[o]" in span:
        return "offen"
    return "bedingt" if part else "Regel oder Synthese"


def load():
    with LEDGER.open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r["active"].strip().lower() != "false"]
    for r in rows:
        r["block"] = r["claim_id"].split(".")[0]
        r["class"] = classify(r["status"])
    return rows


def area_of(block):
    for index, (_, blocks) in enumerate(AREAS[:-1]):
        if block in blocks.split():
            return index
    return len(AREAS) - 1


def short(text, limit=190):
    text = re.sub(r"\s+", " ", text).strip()
    if len(text) <= limit:
        return text
    return text[:limit].rsplit(" ", 1)[0] + " …"


def find(pattern, rows=None, limit=400):
    rx = re.compile(pattern, re.I)
    return [(r["claim_id"], r["class"], short(r["claim"], limit), r["script"])
            for r in (rows or load()) if rx.search(r["claim_id"] + " " + r["claim"])]


def counts(rows=None):
    rows = rows or load()
    out = {c: 0 for c in CLASSES}
    for r in rows:
        out[r["class"]] += 1
    return out


def render(rows=None):
    rows = rows or load()
    total = counts(rows)
    lines = ["# Ledger-Karte: was zu welchem Teil der Theorie bereits feststeht", "",
             "Erzeugt von `research/originlab/ledger_map.py --write` aus `verification/status_ledger.csv`. "
             "Nicht von Hand ändern. Die Klasse ist aus dem Statustext des Ledgers abgelesen und keine neue "
             "Bewertung; maßgeblich bleiben das Ledger und die Originalquellen.", "",
             f"**{len(rows)} aktive Behauptungen:** " + ", ".join(f"{total[c]} {c}" for c in CLASSES if total[c]) + ".",
             "",
             "Zweck: Vor jeder Aussage „neu“ oder „offen“ hier nachsehen. Suche im Volltext mit "
             "`ledger_map.py --find 'Muster'`.", "",
             "## Übersicht", "", "| Gebiet | Blöcke | " + " | ".join(CLASSES) + " |",
             "|---|---|" + "---|" * len(CLASSES)]
    grouped = {}
    for r in rows:
        grouped.setdefault(area_of(r["block"]), {}).setdefault(r["block"], []).append(r)
    for index, (title, _) in enumerate(AREAS):
        blocks = grouped.get(index, {})
        c = {k: 0 for k in CLASSES}
        for members in blocks.values():
            for r in members:
                c[r["class"]] += 1
        names = ", ".join(f"{b} {len(m)}" for b, m in sorted(blocks.items(), key=lambda kv: -len(kv[1])))
        lines.append(f"| {title} | {names} | " + " | ".join(str(c[k]) for k in CLASSES) + " |")
    for index, (title, _) in enumerate(AREAS):
        lines += ["", f"## {title}"]
        for block, members in sorted(grouped.get(index, {}).items(), key=lambda kv: (-len(kv[1]), kv[0])):
            lines += ["", f"### {block} ({len(members)})", ""]
            for r in members:
                script = f" · `{r['script'].split(';')[0].strip()}`" if r["script"].strip() else ""
                lines.append(f"- **{r['claim_id']}** [{r['class']}] {short(r['claim'])}{script}")
    return "\n".join(lines) + "\n"


def atlas_summary(rows=None):
    """Per-area counts for the Atlas overview, so a broad status shows the whole ledger and not a selection."""
    rows = rows or load()
    grouped = {}
    for r in rows:
        grouped.setdefault(area_of(r["block"]), {}).setdefault(r["block"], []).append(r)
    areas = []
    for index, (title, _) in enumerate(AREAS):
        blocks = grouped.get(index, {})
        classes = {c: 0 for c in CLASSES}
        for members in blocks.values():
            for r in members:
                classes[r["class"]] += 1
        areas.append({"title": title, "claims": sum(classes.values()),
                      "classes": {c: n for c, n in classes.items() if n},
                      "blocks": ", ".join(f"{b} {len(m)}" for b, m in
                                          sorted(blocks.items(), key=lambda kv: (-len(kv[1]), kv[0])))})
    total = counts(rows)
    return {"claims": len(rows), "classes": {c: n for c, n in total.items() if n}, "areas": areas,
            "map": "docs/research/LEDGER-KARTE.md",
            "lookup": "research/originlab/ledger_map.py --find 'pattern'; search_knowledge with type claim",
            "ledger_path": "verification/status_ledger.csv",
            "ledger_sha256": hashlib.sha256(LEDGER.read_bytes()).hexdigest(),
            "notice": "Counts of active ledger claims per area. The class is read off the ledger status text and "
                      "is no new assessment; 'exakt, Teil bedingt oder offen' means an exact core with a "
                      "conditional or open part. Look a subject up here before calling it new or open."}


def summary_text(rows=None):
    return json.dumps(atlas_summary(rows), ensure_ascii=False, indent=1) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--find")
    args = parser.parse_args(argv)
    rows = load()
    if args.find:
        for cid, cls, text, script in find(args.find, rows):
            print(f"{cid} [{cls}] {text} ({script.split(';')[0].strip()})")
        return 0
    text, summary = render(rows), summary_text(rows)
    if args.check:
        ok = all(path.exists() and path.read_text() == wanted
                 for path, wanted in ((TARGET, text), (ATLAS_SUMMARY, summary)))
        print("current" if ok else "out of date")
        return 0 if ok else 1
    if args.write:
        TARGET.write_text(text)
        ATLAS_SUMMARY.write_text(summary)
        print(f"{TARGET.relative_to(ROOT)}: {len(rows)} claims, {len(text.splitlines())} lines; "
              f"{ATLAS_SUMMARY.relative_to(ROOT)}: {len(AREAS)} areas")
        return 0
    for name, value in counts(rows).items():
        print(f"{value:4d}  {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
